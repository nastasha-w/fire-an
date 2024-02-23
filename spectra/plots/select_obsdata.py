# stub for obs. selection based on halo mass and redshift
import numpy as np
import pandas as pd
import scipy.special as sps

import fire_an.mstar_mhalo.loader_smdpl_sfr as ldsmdpl
import fire_an.utils.math_utils as mu

obsdatadir = '/projects/b1026/nastasha/extdata/highuv_velocities/'
filen_cubs_components = obsdatadir + 'ne8_detected_components_cubs.dat'
filen_cubs_systems = obsdatadir + 'ne8_detected_systems_cubs.dat'

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

def getMhest_cubs(cubssys: pd.DataFrame,
                  err_sigma: tuple[float] = (1, 2)):
    '''
    Returns: 
    --------
    pandas Dataframe
        input dataframe with errors at requested sigma levels added.
        columns are 
        'Mvir_bestest [log10 Msun]'
        'err_Mvir_<err_sigma>_plus [dex]'
        'err_Mvir_<err_sigma>_minus [dex]'
        for all requested err_sigma levels
    '''
    # assuming 1-sigma error bars in preprint Table B.1
    mss = cubssys['Mstar [log10Msun]']
    mss_hierr_1s = cubssys['err_Mstar_plus [dex]']
    mss_loerr_1s = cubssys['err_Mstar_minus [dex]']
    zsgal = cubssys['zgal']
    #mss_hierr_1s = 0.5 * cubsdatadict['mstar_2s_hierr_dex']
    #mss_loerr_1s = 0.5 * cubsdatadict['mstar_2s_loerr_dex']
    calscatter = 0.2 # scatter in SED-fit to few-band M* calibration
    mss_hierr_1s = np.sqrt(mss_hierr_1s**2 + calscatter**2)
    mss_loerr_1s = np.sqrt(mss_loerr_1s**2 + calscatter**2)
    
    mhs = []
    mh_loerrs = []
    mh_hierrs = []
    for ms, loerr_1s, hierr_1s, zgal in zip(
            mss, mss_hierr_1s, mss_loerr_1s, zsgal):
        mh, mh_hierr, mh_loerr = getMhest(
            ms, (loerr_1s, hierr_1s), zgal, sigmas_tar=err_sigma)
        mhs.append(mh)
        mh_loerrs.append(mh_loerr)
        mh_hierrs.append(mh_hierr)
    mhs = np.array(mhs)
    mh_hierrs = np.array(mh_hierrs).T
    mh_loerrs = np.array(mh_loerrs).T
    cubssys['Mvir_bestest [log10 Msun]'] = mhs
    for si, s in enumerate(err_sigma):
        cubssys[f'err_Mvir_{s:.1f}_plus'] = mh_hierrs[si] 
        cubssys[f'err_Mvir_{s:.1f}_minus'] = mh_loerrs[si]
    return cubssys

def selectsystems(logMh_sim: np.ndarray[float],
                  z_sim: np.ndarray[float],
                  dataset: {'cubs', 'casbah'},
                  Mhmargin_dex: float = 0.2,
                  zmargin: float = 0.05,
                  method: str = 'bestest_inrange'):
    '''
    Parameters:
    -----------
    method: str
        'inrange_bestest': zgal and mode of the Mvir distribution
                           lie within the margins for both
        'inrange_sigma_<n>': zgal lies in z range (with margin),
                             Mvir n-sigma range overlaps the 
                             simulation Mvir range (with margin)
    '''
    minlogMh, maxlogMh, minz, maxz =  getranges_Mh_z(
            logMh_sim, z_sim, Mhmargin_dex=Mhmargin_dex, zmargin=zmargin)
    if method == 'bestest_inrange':
        nsig = 1. # just use some default
    elif method.startswith('inrange_sigma_'):
        nsig = float(method.split('_')[-1])

    if dataset == 'cubs':
        df = pd.read_csv(filen_cubs_systems, sep='\t', comment='#',
                         converters={'sigmav [km/s]': np.float64,
                                     'QSO': str,
                                     'system': str,
                                     'impact_parameter [kpc]': np.float64,
                                     })
        getMhest_cubs(df, err_sigma=(nsig,))
    elif dataset == 'casbah':
        raise NotImplementedError('Casbah read-in: TODO')
    
    if method == 'bestest_inrange':
        sel = np.ones(len(df['Mvir_bestest [log10 Msun]']), dtype=bool)
        sel &= df['zgal'] <= maxz
        sel &= df['zgal'] >= minz
        sel &= df['Mvir_bestest [log10 Msun]'] <= maxlogMh
        sel &= df['Mvir_bestest [log10 Msun]'] >= minlogMh
    elif method.startswith('inrange_sigma_'):
        sel = np.ones(len(df['Mvir_bestest [log10 Msun]']), dtype=bool)
        sel &= df['zgal'] <= maxz
        sel &= df['zgal'] >= minz
        sel &= df[f'err_Mvir_{nsig:.1f}_minus'] <= maxlogMh
        sel &= df[f'err_Mvir_{nsig:.1f}_plus'] >= minlogMh
    return df[sel]


        