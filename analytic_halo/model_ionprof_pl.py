'''
Use power law models from Stern et al.
'''

import numpy as np
import scipy.stats as spstat

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
    def __init__(self, mvir_msun, redshift, fcgm, z_sol, pli):
        '''
        Parameters:
        -----------
        mvir_msun: float
           halo mass (BN98 definition), units: solar mass.
        redshift: float
           redshift (affects Rvir given Mvir)
        fcgm: float between 0 and 1
           fraction of halo mass in the warm/hot phase of the CGM,
           relative to Omega_baryon / Omega_matter
        z_sol: float
           metallicity in solar units
        pli: float
           power law index (for the circular velocity profile)
           0. : isothermal profile
           0.5: point mass
        '''
        self.mvir_cgs = mvir_msun * c.solar_mass
        self.cosmopars = cosmopars_base_fire.copy()
        self.cosmopars['z'] = redshift
        self.cosmopars['a'] = 1. / (1. + redshift)
        self.fcgm = fcgm
        self.z_sol = z_sol 
        self.pli = pli

        self._setup_model()

    def _setup_model(self):
        self.mu = 0.59 # about right for ionised (hot) gas, primordial
        self.hmassfrac = 0.752 # roughly primordial
        self.rvir_cgs = cu.rvir_from_mvir(self.mvir_cgs, self.cosmopars,
                                          meandef='BN98')
        self.cgmmass_cgs = self.mvir_cgs * self.fcgm \
                           * self.cosmopars['omegab'] \
                           / self.cosmopars['omegam']
        
        self.vc_rvir_cgs = np.sqrt(c.gravity * self.mvir_cgs 
                                   / self.rvir_cgs)
        self.tvir_cgs = self.mu * c.u * self.vc_rvir_cgs**2 \
                        / (2. * c.boltzmann)
        # eq 20
        self.nH_rvir_cgs = (1.5 + self.pli) * self.hmassfrac * self.mvir_cgs \
                           / (4. * np.pi * c.atomw_H * c.u * self.rvir_cgs**3)
        
    def vc_cmps(self, r3d_cm):
        # eq 19
        return self.vc_rvir_cgs * (r3d_cm / self.rvir_cgs)**self.pli
    def t_K(self, r3d_cm):
        # eq 21, 24
        self._a = 0.9 * (1. - 2. * self.pli)
        self._t = 6. / (5. * self._a) \
                  * (r3d_cm / self.rvir_cgs)**(2. * self.pli) \
                  * self.tvir_cgs
        return self._t
    def nH_cm3(self, r3d_cm):
        # eq 20
        return self.nH_rvir_cgs * (r3d_cm / self.rvir_cgs)**(-1.5 + self.pli)
    
    def coldensprof(self, ion, impactpars_pkpc, loslen_rvir=4.,
                    sigma_logT=0.3):
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
        self._dTsample = np.linspace(-3. * sigma_logT, 3. * sigma_logT, 21)
        self._normob = spstat.norm()
        self._dTweights = self._normob.pdf(self._dTsample / sigma_logT)
        self._dTweights *= 1. / np.sum(self._dTweights)
        self._logT_K = np.log10(self.t_K(self._r3d_cm))
        self._logT_K = self._logT_K[..., np.newaxis] \
                       + self._dTsample[np.newaxis, np.newaxis, :]
        self._logZ = np.log10(self.z_sol * self._tab.solarZ) \
                     * np.ones(np.prod(self._logT_K.shape))
        self._lognH_cm3 = np.log(self.nH_cm3(self._r3d_cm))
        self.__lognH_cm3 = np.tile(self._lognH_cm3, 
                                   (1, 1, len(self._dTsample)))
        self._interpdct = {'logT': self._logT_K.flatten(), 
                           'lognH': self.__lognH_cm3.flatten(),
                           'logZ': self._logZ}
        self._ionfrac = self._tab.find_ionbal(self._interpdct, log=False)
        self._ionfrac = self._ionfrac.reshape(self._logT_K.shape)
        self._ionfrac = np.sum(self._ionfrac 
                               * self._dTweights[np.newaxis, np.newaxis, :],
                               axis=2)
        self._interpdct = {'logZ': np.array([self._logZ[0]])}
        self._eltabund = self._tab.find_assumedabundance(self._interpdct, 
                                                         log=False)
        self._eltabund = self._eltabund[0]
        self._iondens = self._ionfrac * self._eltabund * 10**self._lognH_cm3
        self._dl = np.average(np.diff(self._losbins_cm)) # uniform bins     
        self.coldens = np.sum(self._iondens * self._dl, axis=1)
        
        del self._ion, self._impactpars_cm, self._los0_cm, self._los1_cm
        del self._losbins_cm, self._r3d_cm, self._tab, self._logZ
        del self._logT_K, self._lognH_cm3, self._interpdct
        del self._ionfrac, self._eltabund, self._iondens, self._dl
        return self.coldens
         


