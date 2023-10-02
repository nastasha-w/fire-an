import astropy.io.fits as apfits
import numpy as np
import pandas as pd
import scipy.special as sps

import fire_an.mstar_mhalo.analytical as msmhan
import fire_an.mstar_mhalo.loader_smdpl_sfr as ldsmdpl
import fire_an.utils.math_utils as mu

## From Zhijie:
#Two notes to clarify the table:
# 1) For each measurements, there are three columns representing 
#    the value, 1sigma lower uncertainty, and upper uncertainty. 
#    In particular, a lower uncertainty of -1 means it is an upper 
#    limit in observation, while an upper uncertainty of -1 is for 
#    lower limits. These limits are all 2 sigma. If both 
#    uncertainties are -1, this value is unconstrained. 
# 2) The velocity of each ion in the fits file is relative to the 
#    column of “zfit”, instead of relative to the redshift of the 
#    associated galaxy.
## follow-up from Zhijie:
# (1) I think it okay to account for the calibration scatter using a 
#     quadrature summation. Now the reported uncertainty is just the 
#     propagated uncertainty from photometric and color measurements.
# (2) Yes, stellar mass and NeVIII column density are in log10 with
#     units of Msun and cm^-2.
# (3) Ms_rmin is the stellar mass of the galaxy with the smallest
#     dproj/rvir, while Ms_tot is the total stellar mass of all 
#     selected nearby galaxies in idv_gals.fits.
# (4) Actually, it is tricky to do the galaxy association with NeVIII
#     absorption. In Figure 3, we show different association scenarios
#     for OVI. If we consider the most massive galaxy is more or less
#     tracing the group properties (e.g., the group center), these
#     massive galaxies don't show a good correlation with observed OVI
#     properties. Therefore, we suggest that detected OVI is more likely
#     to be associated with nearby galaxies rather than the group medium.
#     Not sure whether this conclusion can be applied to NeVIII, but
#     considering these two ions are not so different from each other,
#     it is likely that NeVIII behaves similarly. Is it a question can
#     be answered by simulations?
# Checked in table: no LL/UL stellar masses Ms_rmin or Ms_tot

#cubsdf = ('/Users/nastasha/ciera/projects_lead/fire3_ionabs/'
cubsdf = ('/projects/b1026/nastasha/extdata/'
          'CUBSVII_qu_etal_2022_draft_table2.fits')
oddir = '/projects/b1026/nastasha/extdata/'
ofilen = oddir + 'data_burchett_etal_2019_table1.txt'

def readin_cubsdata():
    with apfits.open(cubsdf, 'readonly') as hdul:
        data = hdul[1].data
        # TODO: check MS_tot
        mstar = data['Ms_rmin'] # log10 Msun
        mstar_2s_loerr = data['Ms0_rmin'] 
        mstar_2s_hierr = data['Ms1_rmin'] 
        impactpar = data['dproj_rmin'] # kpc
        ne8col = data['NeVIII'] # log10 cm**-2
        ne8col_2s_loerr = data['NeVIII0']
        ne8col_2s_hierr = data['NeVIII1']
        z_gal = data['z_rmin']

        ne8col[ne8col == -1.] = np.NaN # missing values
        isul_ne8 = ne8col_2s_loerr == -1.
        isll_ne8 = ne8col_2s_hierr == -1.
        # I don't think there are any, but check.
        if np.any(isll_ne8[np.logical_not(np.isnan(ne8col))]):
            print('Warning: lower limits on Ne VIII col. in data.')
    out = {'mstar_log10Msun': mstar,
           'mstar_2s_loerr_dex': mstar_2s_loerr,
           'mstar_2s_hierr_dex': mstar_2s_hierr,
           'impactpar_kpc': impactpar,
           'z_gal': z_gal,
           'ne8col_logcm2': ne8col,
           'ne8col_2s_loerr_dex': ne8col_2s_loerr,
           'ne8col_2s_hierr_dex': ne8col_2s_hierr,
           'isul_ne8': isul_ne8,
           }
    print('CUBS data read')
    return out

def cumulgauss_lohalf(xsig):
    if 0. in xsig or np.max(xsig) < 0. or np.min(xsig) > 0.:
        ins = False
        _xsig = xsig
    else:
        ins = True
        insi = np.searchsorted(xsig, 0.)
        _xsig = np.array(list(xsig[:insi]) + [0.] + list(xsig[insi:]))
    pcumul = 0.5 * (1. + sps.erf(_xsig / np.sqrt(2.)))
    pcumul[_xsig >= 0.] = 0.5
    if ins:
        pcumul = np.array(list(pcumul[:insi]) + list(pcumul[insi + 1:]))
    return pcumul

def cumulgauss_pastehalfs(xpoints, mu, siglo, sighi, nsig_err=1.):
    xsig_lo = (xpoints - mu) / (nsig_err * siglo)
    xsig_hi = (xpoints - mu) / (nsig_err * sighi)

    pcumul_lo = cumulgauss_lohalf(xsig_lo)
    pcumul_hi = cumulgauss_lohalf(-1. * xsig_hi)
    pcumul_all = pcumul_lo + 0.5  - pcumul_hi
    return pcumul_all

def calchalomassdist_cubs(cubsdatadict):
    sigmas_tar = (1, 2)
    sigmas_tar = np.asarray(sigmas_tar)
    sig2t = msmhan.cumulgauss(sigmas_tar) - msmhan.cumulgauss(-sigmas_tar)
    cumulP_lo = 0.5 - 0.5 * sig2t
    cumulP_hi = 0.5 + 0.5 * sig2t

    redshifts = cubsdatadict['z_gal']  
    histobj = ldsmdpl.SMHMhists(np.array(redshifts), binsize=0.1)

    mss = cubsdatadict['mstar_log10Msun']
    mss_hierr_1s = 0.5 * cubsdatadict['mstar_2s_hierr_dex']
    mss_loerr_1s = 0.5 * cubsdatadict['mstar_2s_loerr_dex']
    calscatter = 0.2 # scatter in SED-fit to few-band M* calibration
    mss_hierr_1s = np.sqrt(mss_hierr_1s**2 + calscatter**2)
    mss_loerr_1s = np.sqrt(mss_loerr_1s**2 + calscatter**2)

    logmvir_msun_bestest = []
    logmvir_msun_lo = []
    logmvir_msun_hi = []
    logmvir_msun_loer = []
    logmvir_msun_hier = []

    for redshift, ms, ms_hierr_1s, ms_loerr_1s in zip(
        redshifts, mss, mss_hierr_1s, mss_loerr_1s):
        
        _msbins_hist = histobj.getbins(redshift, mode='ms')
        _Pcbin_ms = cumulgauss_pastehalfs(_msbins_hist, ms,
                                          ms_loerr_1s,
                                          ms_hierr_1s,
                                          nsig_err=1.)
        _Pbin_ms = np.diff(_Pcbin_ms)
        _mhP, _mhbins = histobj.matrixconv(_Pbin_ms, redshift, 
                                           mode='mstomh')    
        _mhpd = _mhP / np.diff(_mhbins)
        _mhcens = 0.5 * (_mhbins[:-1] + _mhbins[1:])
        logmvir_msun_bestest.append(_mhcens[np.argmax(_mhpd)])
        mlo = mu.linterpsolve(np.cumsum(_mhP), _mhbins[1:], cumulP_lo[0])
        logmvir_msun_lo.append(mlo)
        mhi = mu.linterpsolve(np.cumsum(_mhP), _mhbins[1:], cumulP_hi[0])
        logmvir_msun_hi.append(mhi)
        if len(sigmas_tar) > 1:
            mloer = mu.linterpsolve(np.cumsum(_mhP), _mhbins[1:],
                                    cumulP_lo[1])
            logmvir_msun_loer.append(mloer)
            mhier = mu.linterpsolve(np.cumsum(_mhP), _mhbins[1:],
                                    cumulP_hi[1])
            logmvir_msun_hier.append(mhier)
    out = {'logmvir_msun_bestest': np.array(logmvir_msun_bestest),
           'logmvir_msun_lo': np.array(logmvir_msun_lo),
           'logmvir_msun_hi': np.array(logmvir_msun_hi),
           'logmvir_msun_loer': np.array(logmvir_msun_loer),
           'logmvir_msun_hier': np.array(logmvir_msun_hier),
           }
    print('Halo masses and errors calculated')
    return out

def getplotdata_cubs():
    data = readin_cubsdata()
    _data = calchalomassdist_cubs(data)
    data.update(_data)
    return data

def calcmhalodist_casbah(logmstar_msun, sigmalogmstar, redshift):
    histobj = ldsmdpl.SMHMhists(np.array([redshift]), binsize=0.1)
    msbins_hist = histobj.getbins(redshift, mode='ms')
    _Pcbin_ms = msmhan.cumulgauss((msbins_hist - logmstar_msun) 
                                  / sigmalogmstar)
    _Pbin_ms = np.diff(_Pcbin_ms)
    mhp, _mhbins = histobj.matrixconv(_Pbin_ms, redshift, 
                                      mode='mstomh')    
    return _mhbins, mhp

def readdata_b19(nsigmas=(1, 2)):
    if not hasattr(nsigmas, '__len__'):
        nsigmas = np.array([nsigmas])
    else:
        nsigmas = np.array(nsigmas)
    data_bur = pd.read_csv(ofilen, comment='#', sep='\t')
    #cosmopars_bur = {'h': 0.677, 'omegam': 0.31, 'omegalambda': 0.69}
    sig2t = msmhan.cumulgauss(nsigmas) - msmhan.cumulgauss(-nsigmas)
    cumulP_lo = 0.5 - 0.5 * sig2t
    cumulP_hi = 0.5 + 0.5 * sig2t

    logmstars_msun = data_bur['log_Mstar_Msun']
    sigmalogmstars = data_bur['log_Mstar_Msun_err']
    #isul = data_bur['log_N_Ne8_isUL']
    reshifts = data_bur['zgal']

    logmvir_msun_bestest = []
    logmvir_msun_lo = []
    logmvir_msun_hi = []
    logmvir_msun_loer = []
    logmvir_msun_hier = []
    for lgmstar_msun, slgmstar, _z in zip(logmstars_msun, sigmalogmstars,
                                          reshifts):
        mhbins, mhP = calcmhalodist_casbah(lgmstar_msun, slgmstar, _z)
        mhpd = mhP / np.diff(mhbins)
        mhcens = 0.5 * (mhbins[:-1] + mhbins[1:])
        logmvir_msun_bestest.append(mhcens[np.argmax(mhpd)])
        mlo = mu.linterpsolve(np.cumsum(mhP), mhbins[1:], cumulP_lo[0])
        logmvir_msun_lo.append(mlo)
        mhi = mu.linterpsolve(np.cumsum(mhP), mhbins[1:], cumulP_hi[0])
        logmvir_msun_hi.append(mhi)
        if len(nsigmas) > 1:
            mloer = mu.linterpsolve(np.cumsum(mhP), mhbins[1:], cumulP_lo[1])
            logmvir_msun_loer.append(mloer)
            mhier = mu.linterpsolve(np.cumsum(mhP), mhbins[1:], cumulP_hi[1])
            logmvir_msun_hier.append(mhier)
    logmvir_msun_bestest = np.array(logmvir_msun_bestest)
    logmvir_msun_lo = np.array(logmvir_msun_lo)
    logmvir_msun_hi = np.array(logmvir_msun_hi)
    if len(nsigmas) > 1:
        logmvir_msun_loer = np.array(logmvir_msun_loer)
        logmvir_msun_hier = np.array(logmvir_msun_hier)

    data_bur = data_bur.assign(logmvir_msun_bestest=logmvir_msun_bestest,
                               logmvir_msun_lo=logmvir_msun_lo,
                               logmvir_msun_hi=logmvir_msun_hi)
    if len(nsigmas) > 1:
        data_bur = data_bur.assign(logmvir_msun_loer=logmvir_msun_loer,
                                   logmvir_msun_hier=logmvir_msun_hier)
    return data_bur






