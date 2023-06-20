# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 14:43:18 2017

@author: wijers

General cosmological utility functions; initially copied from make_maps to be
loaded without the entire read_Eagle machinery
"""

import numpy as np
import numbers as num # for instance checking

import fire_an.utils.opts_locs as ol # needed for some ion data
import fire_an.utils.constants_and_units as c

def comoving_distance_cm(cosmopars=None): # assumes Omega_k = 0
    '''
    input:
    ------
    z: redshift (cosmopars redshift is used if cosmopars is None)
    cosmopars: dictionary of cosmological parameter containing
               'omegam': (baryon + dark) matter density / rho_critical at z=0
               'omegam': baryon density / rho_critical at z=0
               'omegalambda': dark energy density / rho_critical at z=0
               'z': redshift
               'a': expansion factor 1 / (1 + z)
               'h': hubble factor z=0 in units of km/s/Mpc
               'boxsize': size of the simulation box in cMpc/h (not used here)
               
    returns:
    --------
    comoving distance in cm
    '''
    if cosmopars is None:
        raise ValueError('Must specify cosmopars for comoving_distance_cm')
    z = cosmopars['z']
    if z < 1e-8:
        print('Using 0 comoving distance from z. \n')
        return 0.
    hpar = cosmopars['h']
    omega0 = cosmopars['omegam']
    omegalambda = cosmopars['omegalambda']
    def integrand(zi):
        return (omega0 * (1. + zi)**3 + omegalambda)**0.5
    zi_arr = np.arange(0., z + z / 512., z / 512.)
    com = np.trapz(1. / integrand(zi_arr), x=zi_arr)
    return com * c.c / (c.hubble * hpar)

def ang_diam_distance_cm(cosmopars=None):
    z = cosmopars['z']
    return comoving_distance_cm(cosmopars=cosmopars) / (1. + z)

def lum_distance_cm(cosmopars=None):
    z = cosmopars['z']
    return comoving_distance_cm(cosmopars=cosmopars) * (1. + z)

def Hubble(cosmopars=None): 
    hpar = cosmopars['h']
    z = cosmopars['z']
    omega0 = cosmopars['omegam']
    omegalambda = cosmopars['omegalambda']
        
    return (c.hubble * hpar) * (omega0 * (1. + z)**3 + omegalambda)**0.5

def rhocrit(cosmopars=None):
    '''
    critical density at z; units: g / cm^-3
    cosmopars z overrides input z if cosmopars are given
    '''
    rhoc = 3. / (8. * np.pi * c.gravity) * (Hubble(cosmopars=cosmopars))**2
    return rhoc

def rhom(cosmopars=None):
    '''
    mean matter (DM + baryons) density at z; units: g / cm^-3
    cosmopars z overrides input z if cosmopars are given 
    '''
    cp = cosmopars.copy()
    cp['z'] = 0.
    cp['a'] = 1.
    rhoc0 = rhocrit(cosmopars=cp)
    omegam0 = cosmopars['omegam']
    _z = cosmopars['z']
    rhom = rhoc0 * omegam0 * (1. + _z)**3
    return rhom

def expfactor_t(time, cosmopars=None):
    '''
    calculate the expansion factor at time t, for a 
    matter + cosmological constant cosmology

    Parameters
    ----------
    time : float
        cosmic time (> 0, since aexp = 0).
    cosmopars : dict or object with cosmopars dict keys as attributes
        cosmological parameters. 
        (The redshift and expansion factor of the cosmopars or object
        object are ignored.)
        
    Raises
    ------
    NotImplementedError
        the cosmological parameters are not for a flat LambdaCDM
        cosmology.
        
    Returns
    -------
    the expasion factor a (float)
    '''
    if isinstance(cosmopars, dict):
        hpar = cosmopars['h']
        omega0 = cosmopars['omegam']
        omegalambda = cosmopars['omegalambda']
    else:
        hpar = cosmopars.h 
        omega0 = cosmopars.omegam
        omegalambda = cosmopars.omegalambda
    if not np.isclose(omegalambda + omega0, 1.):
        msg = ('expfactor_t only works for flat universes;'
               ' received Omega_0 = {o0}, Omega_Lambda = {ol}')
        raise NotImplementedError(msg.format(o0=omega0, ol=omegalambda))
    hzero = hpar * c.hubble
    ft = np.tanh(1.5 * np.sqrt(omegalambda) * time * hzero)**2
    aexp = (omega0 / omegalambda * ft / (1. - ft))**(1./3.)
    return aexp

def t_expfactor(aexp, cosmopars=None):
    '''
    calculate the cosmic time (since the big bang) at a given expansion
    factor (cosmological parameter a)

    Parameters
    ----------
    aexp : float
        cosmological expansion factor a.
    cosmopars :  dict or object with cosmopars dict keys as attributes
        cosmological parameters. 
        (The redshift and expansion factor of the cosmopars or object
        object are ignored.)

    Raises
    ------
    NotImplementedError
        the cosmological parameters are not for a flat LambdaCDM cosmology.

    Returns
    -------
    time : float
        the cosmological time (since the big bang) in seconds.

    '''
    if cosmopars is None:
        hpar = c.hubbleparam 
        omega0 = c.omega0
        omegalambda = c.omegalambda
    elif isinstance(cosmopars, dict):
        hpar = cosmopars['h']
        omega0 = cosmopars['omegam']
        omegalambda = cosmopars['omegalambda']
    else:
        hpar = cosmopars.h 
        omega0 = cosmopars.omegam
        omegalambda = cosmopars.omegalambda
    if not np.isclose(omegalambda + omega0, 1.):
        msg = ('expfactor_t only works for flat universes;'
               ' received Omega_0 = {o0}, Omega_Lambda = {ol}')
        raise NotImplementedError(msg.format(o0=omega0, ol=omegalambda))
    hzero = hpar * c.hubble
    oa = omegalambda * aexp**3
    fa = np.arctanh(np.sqrt(oa / (omega0 + oa)))
    time = 1. / (1.5 * np.sqrt(omegalambda) * hzero) * fa
    return time
    
def conc_mass_MS15(Mh, cosmopars=None):
    '''
    Schaller et al. 2015 Eagle concentration-mass relation: 
    DM in full hydro Eagle fit
    Note: fit is for z=0, so avoid at high z!!
    '''
    hpar = cosmopars['h']
    
    return 5.699 * (Mh / (1e14 * hpar * c.solar_mass))**-0.074

def rho_NFW(r, Mh, delta=200, ref='rhocrit', z=0., 
            cosmopars=None, c='Schaller15'):
    '''
    returns: density (g /cm^3)
    Mh: halo mass (g)
    r: cm, physical
    delta: overdensity threshold
    c: concentration - number or 'Schaller15' for that relation (z=0 fit)
    ref: reference density for delta ('rhocrit' or 'rhom')
    '''
    if c == 'Schaller15':
        c = conc_mass_MS15(Mh, cosmopars=cosmopars)
    elif not isinstance(c, num.Number):
        raise ValueError('Value %s for c is not a valid option'%(c))
    if ref == 'rhocrit':
        rho_ref = rhocrit(cosmopars=cosmopars)
    elif ref == 'rhom':
        rho_ref = rhom(cosmopars=cosmopars)
    else:
        raise ValueError('Value %s for ref is not a valid option'%(ref))
    Redge = (Mh * 3. / (4. * np.pi * delta * rho_ref)) ** (1. / 3.)
    rnorm = c * r / Redge
    rhoval = (delta * rho_ref * c**3 ) / \
             (3. * (np.log(1. + c) - c / (1. + c)) * rnorm * (1. + rnorm)**2)
    return rhoval

def Rhalo(Mh, delta=200, ref='rhocrit', z=0., cosmopars=None):
    if ref == 'rhocrit':
        rho_ref = rhocrit(cosmopars=cosmopars)
    elif ref == 'rhom':
        rho_ref = rhom(cosmopars=cosmopars)
    else:
        raise ValueError(f'Value {ref} for ref is not a valid option')
    Redge = (Mh * 3. / (4. * np.pi * delta * rho_ref)) ** (1. / 3.)
    return Redge

def Tvir_hot(Mh, delta=200, ref='rhocrit', z=0., cosmopars=None):
    '''
    Mh in g
    '''
    mu = 0.59 # about right for ionised (hot) gas, primordial
    Rh = Rhalo(Mh, delta=delta, ref=ref, z=z, cosmopars=cosmopars)
    return (mu * c.protonmass) / (3. * c.boltzmann) * c.gravity * Mh / Rh

# 200c, 200m tested against cosmology calculator at z=0, 2.8
# BN98 values lay between the two at those redshifts
def getmeandensity(meandef, cosmopars):
    if meandef == 'BN98':
        # Bryan & Norman (1998)
        # for Omega_r = 0: Delta_c = 18*np.pi**2 + 82 x - 39x^2
        # x = 1 - Omega(z) = 1 - Omega_0 * (1 + z)^3 / E(z)^2
        # E(z) = H(z) / H(z=0)
        # 
        _Ez = Hubble(cosmopars=cosmopars) \
            / (cosmopars['h'] * c.hubble)
        _x = cosmopars['omegam'] * (1. + cosmopars['z'])**3 / _Ez**2 - 1.
        _Deltac = 18. * np.pi**2 + 82. * _x - 39. * _x**2
        meandens = _Deltac * rhocrit(cosmopars=cosmopars)
    elif meandef.endswith('c'):
        overdens = float(meandef[:-1])
        meandens = overdens * rhocrit(cosmopars=cosmopars)
    elif meandef.endswith('m'):
        overdens = float(meandef[:-1])
        _csm = cosmopars.copy()
        _csm['z'] = 0.
        cosmo_meandens = rhocrit(cosmopars=_csm) \
                         * cosmopars['omegam'] * (1. + cosmopars['z'])**3
        meandens = cosmo_meandens * overdens
    return meandens

def Tvir_hot_meandef(mh, cosmopars, meandef='BN98'):
    '''
    mh in g, Tvir in K
    '''
    meandens = getmeandensity(meandef, cosmopars)
    rh = (mh * 3. / (4. * np.pi * meandens)) ** (1. / 3.)
    mu = 0.59 # about right for ionised (hot) gas, primordial
    return (mu * c.protonmass) / (3. * c.boltzmann) * c.gravity * mh / rh

def mvir_from_rvir(rh, cosmopars, meandef='BN98'):
    '''
    rh in cm, mvir in g
    '''
    meandens = getmeandensity(meandef, cosmopars)
    return (4. * np.pi / 3.) * meandens * rh**3

def rvir_from_mvir(mvir, cosmopars, meandef='BN98'):
    '''
    mvir in g, rvir in cm
    '''
    meandens = getmeandensity(meandef, cosmopars)
    return (mvir / ((4. * np.pi / 3.) * meandens))**(1./3.)

def solidangle(alpha, beta): 
    # alpha = 0.5 * pix_length_1/D_A, 
    # beta = 0.5 * pix_length_2/D_A
    #from www.mpia.de/~mathar/public/mathar20051002.pdf
    # citing  A. Khadjavi, J. Opt. Soc. Am. 58, 1417 (1968).
    # stored in home/papers
    # using the exact formula, with alpha = beta, 
    # the python exact formula gives zero 
    # for alpha = beta < 10^-3--10^-4
    # assuming the pixel sizes are not so odd that the exact formula
    # is needed in one direction and gives zero for the other,
    # use the Taylor expansion to order 4 
    # testing the difference between the Taylor and python exact
    # calculations shows that for square pixels, 10**-2.5 is a
    # reasonable cut-off
    # for rectangular pixels, the cut-off seems to be needed in both
    # values
    if alpha < 10**-2.5 or beta <10**-2.5:
        return 4 * alpha * beta - 2 * alpha * beta * (alpha**2 + beta**2)
    else: 
        return 4 * np.arccos(((1. + alpha**2 + beta**2) 
                              /((1. + alpha**2) * (1. + beta**2)))**0.5)

def Tvir(m200c, cosmopars=None, mu=0.59):
    '''
    input: 
    m200c in solar masses
    z directly or via cosmopars (cosmopars 'wins')
    output:
    Tvir in K       
    '''
    # formula from LSS + Gal. Form. notes, page 35
    # mu = means molecular weight / hydrogen mass. 0.59 is for primordial, fully ionised, gas
    h = cosmopars['h']
    z = cosmopars['z']
    return 4.1e5 * (mu / 0.59) * (m200c / (1e12 / h))**(2. / 3.) * (1. + z)


def R200c_pkpc(M200c, cosmopars):
    '''
    M200c: solar masses
    '''
    M200c = np.copy(M200c)
    M200c *= c.solar_mass # to cgs
    rhoc = (3. / (8. * np.pi * c.gravity) * Hubble(cosmopars=cosmopars)**2)
    R200c = (3 * M200c / (4. * np.pi* 200. * rhoc))**(1./3.)
    return R200c / c.cm_per_mpc * 1e3 

def Teos_eagle(rho):
    '''
    rho = density (g / cm^-3)
    '''
    X = 0.752 # (primodial also used in code, but for T <-> entropy)
    nstar = 0.1 # cm^-3
    Tstar = 8.0e3 # K
    gamma_eos = 4. / 3.
    rhostar = nstar * c.atomw_H * c.u / X
    Teos = Tstar * (rho / rhostar)**(gamma_eos - 1.)
    return Teos
        
def hasSFR_eagle(rho, T, Z, cosmopars):
    '''
    rho = density (g / cm^-3)
    Z = metallcitity (mass fraction)
    T = temperature (K) !! temperature including EOS, not imposed 10^4 K !!
    
    Eagle paper (Schaye et al. 2015), section 4.3
    
    seems to mostly match recorded SFR values when using particle Z values
    '''
    X = 0.752
    nH = X * rho / (c.atomw_H * c.u)

    nHmin = 0.1 * (Z / 0.002)**-0.64
    nHmin = np.minimum(10., nHmin)
    nHmin = np.maximum(nHmin, 57.7 * rhocrit(cosmopars=cosmopars)
                              * cosmopars['omegab']) 
    
    Tmax = 10**0.5 * Teos_eagle(rho) # 0.5 dex above EOS
    
    hassfr = np.logical_and(nH >= nHmin, T <= Tmax)
    return hassfr 

def getdX(L_z, cosmopars=None):
    # assuming L_z is smaller than the distance over which H varies
    # significantly; assumed in single-snapshot projection anyway 
    redshift = cosmopars['z']
    hpar = cosmopars['h']
    dz = Hubble(cosmopars=cosmopars) / c.c * L_z * c.cm_per_mpc
    dX = dz * (1. + redshift)**2 * c.hubble * hpar \
         / Hubble(cosmopars=cosmopars)
    return dX

def getdz(L_z, cosmopars=None):
    # assuming L_z is smaller than the distance over which H varies
    # significantly; assumed in single-snapshot projection anyway 
    dz = Hubble(cosmopars=cosmopars) / c.c * L_z * c.cm_per_mpc
    return dz

def luminosity_to_Sb(Ls, Axis1, Axis2, Axis3, npix_x, npix_y, ion,
                     cosmopars, ps20tables=True):
    '''
    converts cgs luminosity (erg/s) to cgs surface brightness
    (photons/s/cm2/steradian)
    ion needed because conversion depends on the line energy
    
    Parameters
    ----------
    vardict: Vardict instance 
        used to get cosmological parameters
    Ls: array of floats            
        the dimensions of the projected box: Ls[0] is the full extent along 
        the x axis, Ls[1] along y, Ls[2] is along z (diameter, not radius)
    Axis1, Axis2: int  
        axes perpendicular to the line of sight (0=x, 1=y, 2=z
    Axis3: int 
        axis along the line of sight (int)
    npix_x, npix_y: int
        number of pixels along Axis1 and Axis2, respectively 
    ion: str           
        name for the line to get the conversion for, as used in 
        make_maps_opts_locs
    ps20tables: bool    
        if True, the ion refers to a PS20 table line (get energy from the 
        line name, not make_maps_opnts_locs)
    In Ls, the the indices match the simulation axes every time: indices 0, 1, 
    and 2 always correspond to the X, Y, and Z axes respectively.
    Axis1, Axis2, and Axis3 and as used as indices for Ls. Axis1 and Axis2, 
    are used with Ls and npix_<x/y> to determine the dimensions of each pixel, 
    and Axis3 is used to impose a minimum comoving distance of half the extent
    along the line of sight.

    returns:
    --------
    a number (float) which can be multiplied by the emission line luminosity 
    in a pixel in erg/s to get the surface brightness in 
    photons/cm^2/s/steradian            
        '''
    zcalc = cosmopars['z']
    comdist = comoving_distance_cm(cosmopars)
    longlen = max(Ls) * 0.5 * c.cm_per_mpc
    # even at larger values, 
    # the projection along z-axis = projection along sightline 
    # approximation will break down
    if comdist > longlen: 
        ldist = comdist * (1. + zcalc) # luminosity distance
        adist = comdist / (1. + zcalc) # angular size distance
    else:
        ldist = longlen * (1. + zcalc)
        adist = longlen / (1. + zcalc)

    # x, y are axis placeholders and may actually represent different
    # axes in the simulation, as with numpix_x, and numpix_y
    halfangle_x = 0.5 * Ls[Axis1] / (1. + zcalc) / npix_x \
                  * c.cm_per_mpc / adist
    halfangle_y = 0.5 * Ls[Axis2] / (1. + zcalc) / npix_y \
                  * c.cm_per_mpc / adist

    #solidangle = 2*np.pi*(1-np.cos(2.*halfangle_x))
    #print("solid angle per pixel: %f" %solidangle)
    # the (1+z) is not so much a correction to the line energy 
    # as to the luminosity distance:
    # the 1/4 pi dL^2 luminosity 
    # -> flux conversion is for a broad spectrum and includes energy
    # flux decrease to to redshifting
    # multiplying by (1+z) compensates for this: 
    # the number of photons does not change from redshifting
    if ps20tables:
        # units: Å (A), cm (c), and μm (m)
        _wl = (ion[4:]).strip()
        _unit = _wl[-1]
        wl = float(_wl[:-1])
        unit = 1. if _unit == 'c' else 1e-8 if _unit == 'A' \
               else 1e-4 if _unit == 'm' else np.NaN
        wl_cm = wl * unit
        eng_erg = c.planck * c.c / wl_cm
    else:
        raise ValueError('only PS20tables ions can be used here')
    return 1. / (4 * np.pi * ldist**2) * (1. + zcalc) / eng_erg *\
           1. / solidangle(halfangle_x, halfangle_y)