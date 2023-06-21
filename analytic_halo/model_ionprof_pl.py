'''
Use power law models from Stern et al.
'''

import numpy as np

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
        self.tvir_cgs = self.mu * c.u * self.vc_rvir_cgs * self.rvir_cgs**2 \
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
    
    def coldensprof(self, ion, impactpars_pkpc, loslen_rvir=4.):
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
        self._logZ = np.log10(self.z_sol * self._tab.solarZ) \
                     * np.ones(np.prod(self._r3d_cm.shape))
        self._logT_K = np.log10(self.t_K(self._r3d_cm))
        self._lognH_cm3 = np.log(self.nH_cm3(self._r3d_cm))
        self._interpdct = {'logT': self._logT_K.flatten(), 
                           'lognH': self._lognH_cm3.flatten(),
                           'logZ': self._logZ}
        self._ionfrac = self._tab.find_ionbal(self._interpdct, log=False)
        self._ionfrac = self._ionfrac.reshape( self._r3d_cm.shape)
        self._eltabund = self._tab.find_assumedabundance(self._interpdct, 
                                                         log=False)
        self._eltabund = self._eltabund.reshape( self._r3d_cm.shape)
        self._iondens = self._ionfrac * self._eltabund * 10**self._lognH_cm3
        self._dl = np.diff(self._losbins_cm)
        self.coldens = np.sum(self._iondens * self._dl, axis=1)

        del self._ion, self._impactpars_cm, self._los0_cm, self._los1_cm
        del self._losbins_cm, self._r3d_cm, self._tab, self._logZ
        del self._logT_K, self._lognH_cm3, self._interpdct
        del self._ionfrac, self._eltabund, self._iondens, self._dl
        return self.coldens
         


