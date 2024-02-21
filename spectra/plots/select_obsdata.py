# stub for obs. selection based on halo mass and redshift
import numpy as np
import scipy.special as sps

import fire_an.mstar_mhalo.loader_smdpl_sfr as ldsmdpl
import fire_an.utils.math_utils as mu

def getranges_Mh_z(logMh_sim: np.ndarray[float],
                   z_sim: np.ndarray[float],
                   Mhmargin_dex: float = 0.2,
                   zmargin: float = 0.05):
    minlogMh = np.min(logMh_sim) - Mhmargin_dex
    maxlogMh = np.max(logMh_sim) + Mhmargin_dex
    minz = np.min(z_sim) - zmargin
    maxz = np.max(z_sim) + zmargin
    return minlogMh, maxlogMh, minz, maxz

def cumulgauss(xsig):
    return 0.5 * (1. + sps.erf(xsig / np.sqrt(2.)))

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

def getMhest(logMstar_meas: float,
             err_dex_meas: float | tuple[float | float],
             z_meas,
             sigmas_tar: tuple[float] = (1., 2.)):
    if not hasattr(err_dex_meas, '__len__'):
        err_dex_meas = (err_dex_meas, err_dex_meas)
    
    sigmas_tar = np.asarray(sigmas_tar)
    sig2t = cumulgauss(sigmas_tar) - cumulgauss(-sigmas_tar)
    cumulP_lo = 0.5 - 0.5 * sig2t
    cumulP_hi = 0.5 + 0.5 * sig2t

    histobj = ldsmdpl.SMHMhists(np.array([z_meas]), binsize=0.1)

    _msbins_hist = histobj.getbins([z_meas], mode='ms')
    _Pcbin_ms = cumulgauss_pastehalfs(_msbins_hist, 
                                        logMstar_meas,
                                        err_dex_meas[0],
                                        err_dex_meas[1],
                                        nsig_err=1.)
    _Pbin_ms = np.diff(_Pcbin_ms)
    _mhP, _mhbins = histobj.matrixconv(_Pbin_ms, z_meas, 
                                        mode='mstomh')    
    _mhpd = _mhP / np.diff(_mhbins)
    _mhcens = 0.5 * (_mhbins[:-1] + _mhbins[1:])
    bestest = _mhcens[np.argmax(_mhpd)]
    mlo = [mu.linterpsolve(np.cumsum(_mhP), _mhbins[1:], cp)
            for cp in cumulP_lo]
    mhi = [mu.linterpsolve(np.cumsum(_mhP), _mhbins[1:], cp)
            for cp in cumulP_hi]
    return bestest, mlo, mhi

def getMhest_cubs(cubsdatadict):
    mss = cubsdatadict['mstar_log10Msun']
    mss_hierr_1s = 0.5 * cubsdatadict['mstar_2s_hierr_dex']
    mss_loerr_1s = 0.5 * cubsdatadict['mstar_2s_loerr_dex']
    calscatter = 0.2 # scatter in SED-fit to few-band M* calibration
    mss_hierr_1s = np.sqrt(mss_hierr_1s**2 + calscatter**2)
    mss_loerr_1s = np.sqrt(mss_loerr_1s**2 + calscatter**2)