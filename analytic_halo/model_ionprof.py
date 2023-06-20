'''
uses Jonathan Stern's cooling flow package to set up an analytical
halo model,
then calculates ion profiles from the Ploeckinger & Schaye (2020) 
tables
'''

import h5py
import numpy as np
import sys
import scipy.interpolate as spi

import fire_an.ionrad.ion_utils as iu
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu
import fire_an.utils.opts_locs as ol

sys.path.insert(0, ol.path_jscoolingflow)
import cooling_flow.cooling_flow as cf
import cooling_flow.WiersmaCooling as wcool
import cooling_flow.HaloPotential as halopot

outdir_profiles = ('/Users/nastasha/ciera/projects_lead/fire3_ionabs/'
                   'analytical/')
# Planck 2015 from astropy, used in Jonathan Stern's model
cosmo_base = {'h': 0.6774, 'omegam': 0.3075, 'omegab': 0.0486}
cosmo_base['omegalambda'] = 1. - cosmo_base['omegam']

def solutionset(logmvirs_msun, redshift, mdot=1. * cf.un.Msun / cf.un.yr,
                zsol=0.1, plind=-0.1):
    '''
    get density/temperature solutions for a given halo mass set (BN98
    halo masses), mass flow rate (mdot), metallicity (zsol), 
    v_circ slope (plind)

    plind: float
        power law index for the v_c profile
    zsol: float
        metallicity [solar units]
    '''
    mvirs_msun = 10**logmvirs_msun
    cosmopars = cosmo_base.copy()
    cosmopars['z'] = redshift
    cosmopars['a'] = 1. / (1. + redshift)
    rvirs_cm = cu.rvir_from_mvir(mvirs_msun * c.solar_mass,
                                  cosmopars, meandef='BN98')
    rvirs_kpc = rvirs_cm / (c.cm_per_mpc * 1e-3) * cf.un.kpc
    vcs_cmps = np.sqrt(c.gravity * mvirs_msun * c.solar_mass / rvirs_cm)
    print(vcs_cmps / 1e5)
    potentials = [halopot.PowerLaw(plind, vc_cmps * cf.un.cm / cf.un.s, 
                                   rvir_kpc)
                  for vc_cmps, rvir_kpc in zip(vcs_cmps, rvirs_kpc)]
    cooling = wcool.Wiersma_Cooling(zsol, redshift)
    # want to shoot from Rcirc, since the halos considered here will
    # generally be fully virialized. (The simulated ones sure seem to
    # be.)
    solutions = {}
    for lmv, potential in zip(logmvirs_msun, potentials):
        print(lmv)
        rcirc = 0.02 * potential.Rvir # inner stalled radius
        rmax = 2. * potential.Rvir
        #print(rcirc, rmax)
        #print(potential.vc(rcirc))
        #print(potential.vc(rmax))
        solution = cf.shoot_from_R_circ(potential, cooling, 
                                        rcirc, mdot, rmax, 
                                        v0=1. * cf.un.km/cf.un.s,
                                        max_step=0.1, 
                                        T_low=1e4 * cf.un.K, 
                                        T_high=1e6 * cf.un.K,
                                        tol=1e-6, epsilon=0.1,
                                        terminalUnbound=True,
                                        pr=True, 
                                        return_all_results=False)
        solutions[lmv] = solution
    return solutions

def calc_cgmfrac(solution, mvir_msun):
    '''
    get the CGM gas mass fraction. (Inflow rates as an input parameter
    seem uncertain, but gas fractions I can check.)
    built-in Mgas includes all radii, I want 0.1 -- 1 Rvir
    assumes there are at least two solution points before and after
    the radial range edges
    '''
    _rmin_rvir = 0.1
    _rmax_rvir = 1.0
    rmin = _rmin_rvir * solution.potential.Rvir
    rmax = _rmax_rvir * solution.potential.Rvir
    rs = solution.Rs()
    rimin = np.searchsorted(rs, rmin, side='right')
    rimax = np.searchsorted(rs, rmax, side='left')
    frac_rimin_m1 = (rmin - rs[rimin - 1]) / (rs[rimin] - rs[rimin - 1])
    frac_rimax = (rs[rimax] - rmax) / (rs[rimax] - rs[rimax - 1])
    print(frac_rimin_m1, frac_rimax)

    rsel = slice(rimin - 1, rimax + 1, None)
    rsel_diff = slice(rimin - 2, rimax + 2, None)
    drs = 0.5 * (rs[rsel_diff][2:] - rs[rsel_diff][:-2])
    print(drs)
    drs[0] = drs[0] * frac_rimin_m1
    drs[-1] = drs[-1] * frac_rimax
    print(rs[rsel])
    dms = 4. * np.pi * rs[rsel]**2 * drs * solution.rhos()[rsel]
    mtot = np.sum(dms)
    return mtot.to('Msun') / (mvir_msun * cf.un.Msun)

def calcionprof(solution, ion, redshift, zsol,
                impactpars_cm, lossample_cm, 
                truncate_inner_cm):
    r_cm = solution.Rs().to('cm').value 
    logT_K = np.log10(solution.Ts.to('K').value)
    lognH_cm3 = np.log10(solution.nHs.to('cm**-3').value)
    tab = iu.Linetable_PS20(ion, redshift, emission=False,
                            lintable=True)
    logZ = np.log10(zsol * tab.solarZ) * np.ones(len(r_cm))
    interpdct = {'logT': logT_K, 'lognH': lognH_cm3, 'logZ': logZ}
    ionfrac = tab.find_ionbal(interpdct, log=False)
    eltabund = tab.find_assumedabundance(interpdct, log=False)
    iondens = ionfrac * eltabund * 10**lognH_cm3
    rimin = np.where(r_cm > truncate_inner_cm)[0][0]
    _iondens_prof = spi.interp1d(np.log10(r_cm), np.log10(iondens), 
                                 kind='linear',
                                 fill_value=0.)
    id0 = _iondens_prof(np.log10(truncate_inner_cm))
    _r_cm = np.append(truncate_inner_cm, r_cm[rimin:])
    _iondens = np.append(id0, iondens[rimin:])
    iondens_prof = spi.interp1d(np.log10(_r_cm), np.log10(_iondens), 
                                 kind='linear',
                                 fill_value=0.)

    logr3d_cm = 0.5 * np.log10(impactpars_cm[:, np.newaxis]**2,
                               lossample_cm[np.newaxis, :]**2)
    dl = 0.5 * (lossample_cm[2:] - lossample_cm[:-2])
    dl = np.append(dl[0], dl)
    dl = np.append(dl, dl[-1])
    coldens = np.sum(iondens_prof(logr3d_cm) * dl[np.newaxis, :], axis=1)
    return coldens
    
def runsolsgrid(outfile='test.hdf5'):
    logmvirs_msun = np.arange(11.5, 13.65, 0.1)
    #redshifts = np.arange(0.5, 1.05, 0.1)
    redshifts = np.array([0.75])
    #zs_sol = np.array([0.1, 0.2, 0.5, 1.])
    zs_sol = np.array([0.1, 1.])
    #plinds = np.array([0., -0.1, -0.2, -0.5])
    plinds = np.array([-0.1])
    #mdots = np.array([0.01, 0.03, 0.1, 0.3, 1., 3., 10., 30., 100.])
    mdots = np.array([0.01, 0.1, 1., 10., 100.])
    
    with h5py.File(outdir_profiles + outfile, 'a'):
        for redshift in redshifts:
            for zsol in zs_sol:
                for plind in plinds:
                    for mdot in mdots:
                        sols = solutionset(logmvirs_msun, redshift, 
                                           mdot=mdot * cf.un.Msun / cf.un.yr,
                                           zsol=zsol, plind=plind)


