'''
From Stern et al. (2019), power-law model with cooling.
Using Jonathan's cooling-flow package for reasonable estimates of 
the cooling rate Lambda
'''
import numpy as np

import cooling_flow.WiersmaCooling as wcool

import fire_an.mstar_mhalo.analytical as an
import fire_an.ionrad.ion_utils as iu
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu
import fire_an.utils.opts_locs as ol

# copied from m12q AGN-CR
cosmopars_base_fire = {'h': 0.702, 
                       'omegab': 0.0455,
                       'omegam': 0.272,
                       'omegalambda': 0.728}

class PLmodel:
    def __init__(self, mvir_msun, redshift, mdot_msunperyr, z_sol, pli_vc):
        '''
        Parameters:
        -----------
        mvir_msun: float
           halo mass (BN98 definition), units: solar mass.
        redshift: float
           redshift (affects Rvir given Mvir)
        mdot_msunperyr: float
           inflow rate (solar masses / yr)
        z_sol: float
           metallicity in solar units
        pli_vc: float
           power law index for the circular velocity profile
           0. : isothermal profile
           0.5: point mass
        '''
        self.mvir_cgs = mvir_msun * c.solar_mass
        self.cosmopars = cosmopars_base_fire.copy()
        self.cosmopars['z'] = redshift
        self.cosmopars['a'] = 1. / (1. + redshift)
        self.mdot_cgs = mdot_msunperyr * c.solar_mass / c.sec_per_year
        self.z_sol = z_sol 
        self.pli_vc = pli_vc
        self.pli_entropy = 1. + (4. / 3.) * self.pli_vc
        self._setup_model()

    def _setup_model(self):
        self.mu = 0.59 # about right for ionised (hot) gas, primordial
        self.hmassfrac = 0.752 # roughly primordial
        self.rvir_cgs = cu.rvir_from_mvir(self.mvir_cgs, self.cosmopars,
                                          meandef='BN98')
        self.vc_rvir_cgs = np.sqrt(c.gravity * self.mvir_cgs 
                                   / self.rvir_cgs)
        self.tvir_cgs = self.mu * c.u * c.atomw_H * self.vc_rvir_cgs**2 \
                        / (2. * c.boltzmann)
        self._A = 0.9 * (1. - 2 * self.pli_vc)
        self._estimate_Lambda()
    
    def _estimate_Lambda(self):
        # estimate a representative Lambda: 
        # use a consistent metallicity, redshift, high nH (CIE)
        # temperature ~ middle in the CGM
        self.cooling = wcool.Wiersma_Cooling(self.z_sol, self.cosmopars['z'])
        self._reprT = self.t_K(0.5 * self.rvir_cgs) * wcool.un.K
        self._reprnH = 1e-3 * wcool.un.cm**-3
        self._Lambda = self.cooling.LAMBDA(self._reprT, self._reprnH)
        self._Lambda_cgs = self._Lambda.to('erg * cm**3 / s').value 

    def vc_cmps(self, r3d_cm):
        # eq 19
        return self.vc_rvir_cgs * (r3d_cm / self.rvir_cgs)**self.pli_vc
    def t_K(self, r3d_cm):
        # eq 21, 24
        self._t = 6. / (5. * self._A) \
                  * (r3d_cm / self.rvir_cgs)**(2. * self.pli_vc) \
                  * self.tvir_cgs
        return self._t
    def nH_cm3(self, r3d_cm):
        self._nH_temp = np.sqrt(9. * self.pli_entropy * self.mdot_cgs
                                / (40. * np.pi * self._A * r3d_cm**3 
                                   * self._Lambda_cgs)) \
                        * self.vc_cmps(r3d_cm)
        return self._nH_temp
    
    def coldensprof(self, ion, impactpars_pkpc, loslen_rvir=4.):
        # check reasonability of nH, T, abundance, sightline length,
        # sum of dT probabilities
        # ion fractions seem to be lower than expected, even with scatter
        self._ion = ion
        self._impactpars_cm = impactpars_pkpc * (c.cm_per_mpc * 1e-3)
        self._los0_cm = -0.5 * loslen_rvir * self.rvir_cgs
        self._los1_cm = 0.5 * loslen_rvir * self.rvir_cgs
        self._losbins_cm = np.linspace(self._los0_cm, self._los1_cm, 
                                       500)
        self._r3d_cm = np.sqrt((self._losbins_cm**2)[np.newaxis, :] +
                               (self._impactpars_cm**2)[:, np.newaxis])
        self._tab = iu.Linetable_PS20(self._ion, self.cosmopars['z'], 
                                      emission=False, lintable=True)
        #self._dTsample = np.linspace(-3. * sigma_logT, 3. * sigma_logT, 21)
        #self._normob = spstat.norm()
        #self._dTweights = self._normob.pdf(self._dTsample / sigma_logT)
        #self._dTweights *= 1. / np.sum(self._dTweights)
        print(self.t_K(self._r3d_cm))
        self._logT_K = np.log10(self.t_K(self._r3d_cm))
        #self._logT_K = self._logT_K[..., np.newaxis] \
        #               + self._dTsample[np.newaxis, np.newaxis, :]
        self._logZ = np.log10(self.z_sol * self._tab.solarZ) \
                     * np.ones(np.prod(self._logT_K.shape))
        print(self.nH_cm3(self._r3d_cm))
        self._lognH_cm3 = np.log10(self.nH_cm3(self._r3d_cm))
        #self.__lognH_cm3 = np.tile(self._lognH_cm3, 
        #                           (1, 1, len(self._dTsample)))
        self._interpdct = {'logT': self._logT_K.flatten(), 
                           'lognH': self._lognH_cm3.flatten(),
                           'logZ': self._logZ}
        self._ionfrac = self._tab.find_ionbal(self._interpdct, log=False)
        self._ionfrac = self._ionfrac.reshape(self._logT_K.shape)
        #self._ionfrac = np.sum(self._ionfrac 
        #                       * self._dTweights[np.newaxis, np.newaxis, :],
        #                       axis=2)
        self._interpdct = {'logZ': np.array([self._logZ[0]])}
        self._eltabund = self._tab.find_assumedabundance(self._interpdct, 
                                                         log=False)
        self._eltabund = self._eltabund[0]
        self._iondens = self._ionfrac * self._eltabund * 10**self._lognH_cm3
        self._dl = np.average(np.diff(self._losbins_cm)) # uniform bins     
        self.coldens = np.sum(self._iondens * self._dl, axis=1)
        
        #del self._ion, self._impactpars_cm, self._los0_cm, self._los1_cm
        #del self._losbins_cm, self._r3d_cm, self._tab, self._logZ
        #del self._logT_K, self._lognH_cm3, self._interpdct
        #del self._ionfrac, self._eltabund, self._iondens, self._dl
        return self.coldens
         