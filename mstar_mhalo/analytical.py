'''
various analytical Mstar-Mhalo relations
'''

import numpy as np
import scipy.special as sps
import scipy.interpolate as spi

import fire_an.utils.math_utils as mu
import fire_an.utils.cosmo_utils as cu
import fire_an.utils.constants_and_units as c

def mstar_moster_etal_2013(mvir_msun, redshift):
    # uses M200c for the Mvir definition
    # fits as a function of redshift, with uncertainty in M*
    fp_m10 = 11.590
    fp_m11 = 1.195
    fp_n10 = 0.0351
    fp_n11 = -0.0247
    fp_beta10 = 1.376
    fp_beta11 = -0.826
    fp_gamma10 = 0.608
    fp_gamma11 = 0.329
    fitpar_m1_msun = 10**(fp_m10 + fp_m11 * redshift / (1. + redshift))
    fitpar_n = fp_n10 + fp_n11 * redshift / (1. + redshift)
    fitpar_beta = fp_beta10 + fp_beta11  * redshift / (1. + redshift)
    fitpar_gamma = fp_gamma10 + fp_gamma11 * redshift / (1. + redshift)
    out = mvir_msun * 2. * fitpar_n \
          / ((mvir_msun / fitpar_m1_msun)**-fitpar_beta \
             + (mvir_msun / fitpar_m1_msun)**fitpar_gamma)
    return out


# numerically inverted to Mstar -> Mhalo in the paper
def mstar_burchett_etal_2019(mvir_msun, redshift):
    # modified Moster et al. (2013)
    fp_m10 = 11.590
    fp_m11 = 1.195
    fp_n10 = 0.0351
    fp_n11 = -0.0247
    fp_beta10 = 1.376
    fp_beta11 = -0.826
    fitpar_m1_msun = 10**(fp_m10 + fp_m11 * redshift / (1. + redshift))
    fitpar_n = fp_n10 + fp_n11 * redshift / (1. + redshift)
    fitpar_beta = fp_beta10 + fp_beta11  * redshift / (1. + redshift)
    out = mvir_msun * 2. * fitpar_n \
          / ((mvir_msun / fitpar_m1_msun)**-fitpar_beta + 1.)
    return out

def cumulgauss(xsig):
    return 0.5 * (1. + sps.erf(xsig / np.sqrt(2.)))

def calcdist(monorel_x, monorel_y, xbins, xprob,
             filter_monorels=False, loind_cut=20,
             hiind_cut=-20, cutoff_xbins=False):
    '''
    for a scatter-free relation sampled by monorel_x, monorel_y,
    calculate the y probability density functions
    corresponding to x probabilties of xprob in xbins.

    Parameters:
    -----------
    monorel_x: array of floats
        x values defining the relation. Must be monotonically
        increasing or decreasing.
    monorel_y: array of floats
        y values corresponding to the x values. Must be monotonically
        increasing or decreasing.
    xbins: array of floats. must be monotonically increasing.xs
        bins in which the x probabilities are calculated. All bin edges
        must fall in the range of monorel_x points.
    xprob: array of floats
        probability value (NOT probability density) in each x bin.
    filter_monorels: bool
        If x_monorel and y_monorel might not be monotonic at the edges
        of the range (e.g., distribution medians), cut off the end of
        those distributions. Note that the xbins then need to fall 
        within the *filtered* x_monorel range.
        The filtering is done only if this parameter is True. The 
        default is False. Note the midpoint_x and midpoint_y parameters
        if filter_monorels is True. Typically, only one of the arrays
        will actually need filtering.
    midpoint_x: float
        if filter_monorels is True, midpoint_x sets the x point below
        which the maximum non-monotonic bin is cut off, and above which
        the minimum non-monotonic bin is cut off.
    midpoint_y: float
        if filter_monorels is True, midpoint_y sets the y point below
        which the maximum non-monotonic bin is cut off, and above which
        the minimum non-monotonic bin is cut off.
    loind_cut: int
        index below which non-monotonic bins are cut
    hiind cut: int
        index above which non-monotonic bins are cut. 
        (Values < 0 are allowed.)
    cutoff_xbins: bool
        if True, cut off the xbins outside the xbins range. A warning
        is printed if the cut-off bins have non-zero probability.

    
    Returns:
    --------
    ybins: array of floats
        bins in y. Will generally be of different sizes, even if the x
        bins were evenly speced.
    yprob: array of floats
        probability DENSITY in the y bins.
    '''
    if filter_monorels:
        # x filtering
        _dx = np.diff(monorel_x)
        _x_sgn = np.sign(np.median(_dx[np.isfinite(_dx)]))
        badbins = np.where(np.logical_or(_x_sgn * _dx <= 0.,
                                         np.logical_not(np.isfinite(_dx))))[0]
        _xi_min = 0
        _xi_max = len(monorel_x)
        if hiind_cut < 0:
            hiind_cut = len(monorel_x) + hiind_cut
        if len(badbins) >= 1:
            if np.any(badbins < loind_cut):
                _xi_min = np.max(badbins[badbins < loind_cut])
            if np.any(badbins > hiind_cut):
                _xi_max = np.min(badbins[badbins > loind_cut])
            if np.any(np.logical_and(badbins >= loind_cut,
                                     badbins <= hiind_cut)):
                msg = (f'loind_cut {loind_cut}, hiindcut {hiind_cut}'
                       f' miss non-monotonic steps in {monorel_x}')
                raise ValueError(msg)
        print(_xi_min, _xi_max)
        # y filtering; should be increasing of decreasing
        _dy = np.diff(monorel_y)
        _y_sgn = np.sign(np.median(_dy[np.isfinite(_dy)]))
        badbins = np.where(np.logical_or(_y_sgn * _dy <= 0.,
                                         np.logical_not(np.isfinite(_dy))))[0]
        _yi_min = 0
        _yi_max = len(monorel_y)
        if len(badbins) >= 1:
            if np.any(badbins < loind_cut):
                _yi_min = np.max(badbins[badbins < loind_cut])
            if np.any(badbins > hiind_cut):
                _yi_max = np.min(badbins[badbins > loind_cut])
            if np.any(np.logical_and(badbins >= loind_cut,
                                     badbins <= hiind_cut)):
                msg = (f'loind_cut {loind_cut}, hiindcut {hiind_cut}'
                       f' miss non-monotonic steps in {monorel_y}')
                raise ValueError(msg)
        print(_yi_min, _yi_max)
        
        selmin = max(_xi_min, _yi_min)
        selmax = min(_xi_max, _yi_max)
        monorel_x = monorel_x[selmin : selmax]
        monorel_y = monorel_y[selmin : selmax]
        print(monorel_x)
        print(monorel_y)
    if cutoff_xbins:
        xb_min = np.where(xbins < np.min(monorel_x))[0]
        if len(xb_min) >= 1:
            xb_min = xb_min[-1] + 1
        else:
            xb_min = 0
        xb_max = np.where(xbins > np.max(monorel_x))[0]
        if len(xb_max) >= 1:
            xb_max = xb_max[0]
        else:
            xb_max = len(xbins)
        xbins = xbins[xb_min : xb_max]
        xprob_excl = np.sum(xprob[:xb_min]) + np.sum(xprob[xb_max - 1:])
        if xprob_excl > 0.:
            msg = (f'x bins cutoffs at {xbins[0]:.2f}, {xbins[-1]:.2f}'
                   f' exclude a region with probability {xprob_excl:.2e}')
            print(msg)
        xprob = xprob[xb_min : xb_max - 1]
    
    print(monorel_x)
    print(monorel_y)
    ybins = np.array([mu.linterpsolve(monorel_x, monorel_y, xp)
                      for xp in xbins])
    yprob = xprob / np.diff(ybins)
    return ybins, yprob

def concentration_mass_m200c_dm14(hm_msun, cosmopars):
    '''
    Dutton and Maccio (2014) relation
    from their Table 3; listed overdensities there are relative 
    to critical. concentrations are r_{-2} / R_{Delta crit};
    classic R_s = r_{-2} (by my math at least)
    '''
    zvals = [0., 0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
    avals_z = [0.905, 0.814, 0.728, 0.612, 0.557, 0.528, 0.539]
    bvals_z = [-0.101, -0.086, -0.073, -0.050, -0.021, 0.000, 0.027]
    afunc = spi.interp1d(zvals, avals_z, kind='linear')
    bfunc = spi.interp1d(zvals, bvals_z, kind='linear')
    aval = afunc(cosmopars['z'])
    bval = bfunc(cosmopars['z'])
    c200 = 10**(aval + bval * np.log10(hm_msun * cosmopars['h']/ 1e12))
    return c200


# convert M200c to Mvir (Bryan and Norman):
def m200c_to_mvir(m200c_msun, cosmopars, c200):
    meandens_bn98 = cu.getmeandensity('BN98', cosmopars)
    meandens_200c = cu.getmeandensity('200c', cosmopars)
    m200c_g = m200c_msun * c.solar_mass
    r200c_cgs = (m200c_g / meandens_200c * 3. / (4. * np.pi))**(1./3.)
    
    rbins = np.linspace(0., 5. * r200c_cgs, 5000)
    rcens = 0.5 * (rbins[:-1] + rbins[1:])
    dr = np.average(np.diff(rbins))
    rhovals = cu.rho_NFW(rcens, m200c_g, delta=200, ref='rhocrit', 
                         z=cosmopars['z'],
                         cosmopars=cosmopars, c=c200)
    menc = np.cumsum(4. * np.pi * rhovals * rcens**2 * dr)
    rhomean = menc / (4. * np.pi / 3. * rbins[1:]**3)
    rbn98_cgs = mu.linterpsolve(rhomean, rbins[1:], meandens_bn98)
    mbn98_g = 4. * np.pi / 3. * meandens_bn98 * rbn98_cgs**3
    mbn98_msun = mbn98_g / c.solar_mass
    return mbn98_msun 