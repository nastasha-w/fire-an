import numpy as np
import scipy.interpolate as scpi

from fire_an.ionrad.ion_utils import Linetable_PS20
import fire_an.utils.math_utils as mu

# from Lide D.R., ed. 2003, CRC Handbook of Chemistry and Physics,
#  84 edn. CRC Press LLC, Boca Raton:
# O6, O7, O8, Ne8, Ne9, Fe17
#_ionization_energies_ev = {
#    'O6': 138.12,
#    'Ne8': 239.10,
#    'O7': 739.29,
#    'O8': 871.41,
#    'Ne9': 1195.83,
#    'Fe17': 1266.,
#}
#
#ionization_energies_erg = {key: val * c.ev_to_erg 
#                           for key, val in _ionization_energies_ev.items()}

ionclasses = {'CIE': 0,
              'PIE': 1,
              'C+PIE': 2}

def getCPIEtrans(ionfrac_nH, lognH_cm3, minjump_CPIEtrans=2.):
    '''
    ion fractions should be the actual values, not log values.
    '''
    icie = -6 # some weird stuff at the highest tabulated densities
    frac_cielim = ionfrac_nH[icie]
    cienH = lognH_cm3[icie]
    maxfrac_nhi = np.argmax(ionfrac_nH)
    maxfrac = ionfrac_nH[maxfrac_nhi]
    if np.isclose(maxfrac, 0.):
        msg = ('Input curve only contains ion fractions close to 0 '
               f'(<={maxfrac:e}).')
        print(msg)
        #raise ValueError(msg)
        # assume this is the high-temperature limit
        return -np.inf
    # classifications by Strawn et al. (2021)
    # CIE limit: ion fraction consistently declines with density
    #            return -np.inf: CIE at any density
    # PIE limit: ion fraction decreases as density increases away
    #            from the max ion fraction
    #            return np.inf: PIE at any density
    # transition: ~constant ion frac with density at high nH, but
    #            then increases as density decreases at some point
    #            return nH where 
    #            ion frac = minjump_CPIEtrans * cie fraction

    # will depend a bit on the table range, and could be an issue for 
    # low ions, but seems like it should be ok
    cielim_compnH = cienH - 2.
    cielim_altimar = 2
    deltamax_dex = 0.01
    
    ciecheck_minmaxfrac = (10**(-deltamax_dex) * frac_cielim,
                           10**(deltamax_dex) * frac_cielim)
    cielim_compnHi = np.argmin(np.abs(cielim_compnH - lognH_cm3))
    if cielim_compnHi == len(lognH_cm3) - 1:
        cielim_compnHi =  len(lognH_cm3) - 1 - cielim_altimar
    frac_ciecomp = ionfrac_nH[cielim_compnHi]
    hascielim = frac_ciecomp >= ciecheck_minmaxfrac[0] \
                and frac_ciecomp <= ciecheck_minmaxfrac[1] \
                and not np.isclose(frac_cielim, 0., atol=1e-10)
    if not hascielim: # always PIE
        return np.inf
    elif maxfrac < minjump_CPIEtrans * frac_cielim: # always CIE
        return -np.inf
    else:
        target = minjump_CPIEtrans * frac_cielim
        intercepts = mu.find_intercepts(ionfrac_nH, lognH_cm3, target, 
                                        xydct=None)  
        return  intercepts[-1]

def get_cie_pie_nHT_twolines(ion, redshift, useZ_log10sol=0.):
    '''
    single cut on nH/T
    untested
    '''
    # Clayton Strawn's minimum T criteria:
    # https://ui.adsabs.harvard.edu/abs/2023MNRAS.519....1S
    # https://ui.adsabs.harvard.edu/abs/2021MNRAS.501.4948S

    minjump_CPIEtrans = 2.
    minfactor_CPIEmix = 2.

    todoc = {}
    todoc['method_Tcut'] = ('Strawn et al. (2021) Tmin for CIE:'
                            ' largest tabulated temperature where '
                            ' ion fractions are >= minjump_CPIEtrans'
                            ' times the CIE ion fraction at some '
                            'density.')
    todoc['method_nHcut'] = ('at CIE max. temperature, the minimum '
                             'density '
                             'where the ion fraction differs from the '
                             'CIE limit by at least a factor '
                             'minfactor_CPIEmix')
    todoc['ion'] = ion
    todoc['redshift'] = redshift
    todoc['useZ_log10sol'] = useZ_log10sol
    todoc['minjump_CPIEtrans'] = minjump_CPIEtrans
    todoc['minfactor_CPIEmix'] = minfactor_CPIEmix

    iontab = Linetable_PS20(ion, redshift, emission=False, 
                            vol=True, lintable=True)
    iontab.findiontable()
    logTK = iontab.logTK
    lognH = iontab.lognHcm3
    logZsol = iontab.logZsol
    refZi = np.where(np.isclose(useZ_log10sol, logZsol))[0][0]    
    tab_T_nH = iontab.iontable_T_Z_nH[:, refZi, :]

    # find min. T where CIE limit exists
    icie = -6
    ciecurve = tab_T_nH[:, icie]
    ciemax_logTi = np.argmax(ciecurve)
    ciemax_logT = logTK[ciemax_logTi]
    guess_minlogT = ciemax_logT - 1.
    guess_Ti = np.argmin(np.abs(guess_minlogT - logTK)) 
    # we pretty much have to assume the value is in our tabulated range
    minval_minlogT = logTK[0]
    minval_minlogTi = 0
    # at CIE max, we can assume there is a CIE nH range
    maxval_maxlogT = ciemax_logT
    maxval_maxlogTi = ciemax_logTi
    
    # bisection search, I think
    # not the most efficient probably, but it's fast enough
    while minval_minlogTi < maxval_maxlogTi - 1:
        # integer division means guesses could miss values
        if guess_Ti == minval_minlogTi:
            guess_Ti = guess_Ti + 1
        elif guess_Ti == maxval_maxlogTi:
            guess_Ti = guess_Ti - 1
        nHtrans = getCPIEtrans(tab_T_nH[guess_Ti, :], 
                               lognH, 
                               minjump_CPIEtrans=minjump_CPIEtrans)
        if nHtrans < np.inf: # has a transition CIE -> PIE, or all CIE
            maxval_maxlogTi = guess_Ti
        elif nHtrans == np.inf: # only PIE
            minval_minlogTi = guess_Ti
        guess_Ti = (minval_minlogTi + maxval_maxlogTi) // 2
    # minimum T with a CIE transition is somewhere between
    # minval (PIE only) and maxval (CIE + PIE or just CIE)
    # for a 'round' number, just take the lower values, i.e.,
    # the highest tabulated temperature without a transition
    logTcut_K = logTK[minval_minlogTi]
    
    # nHcut: factor difference with CIE at CIE max T
    icie = -6
    maxfrac = ciecurve[ciemax_logTi]
    ionfrac_nH = tab_T_nH[ciemax_logTi]
    fracrange_cie = (maxfrac / minfactor_CPIEmix, 
                     maxfrac * minfactor_CPIEmix)
    nHi_lim = np.where(np.logical_or(ionfrac_nH < fracrange_cie[0],
                                     ionfrac_nH > fracrange_cie[1]))[0][-1]
    lognHcut_cm3 = lognH[nHi_lim + 1]

    todoc['ionbalfile'] = iontab.ionbalfile
    todoc['ionbaltable_vol'] = True
    todoc['ionbaltable_lintable'] = True
    todoc['lognHcut_cm3'] = lognHcut_cm3
    todoc['logTcut_K'] = logTcut_K 
    return lognHcut_cm3, logTcut_K, todoc

def get_cie_pie_nHT_strawn21(ion, redshift, useZ_log10sol=0.):
    '''
    get the CIE/PIE transition density from Strawn et al. (2021)
    at each temperature
    '''
        # Clayton Strawn's minimum T criteria:
    # https://ui.adsabs.harvard.edu/abs/2023MNRAS.519....1S
    # https://ui.adsabs.harvard.edu/abs/2021MNRAS.501.4948S

    minjump_CPIEtrans = 2.
    minfactor_CPIEmix = 2.

    todoc = {}
    todoc['method'] = ('Strawn et al. (2021) transition log10 nH '
                       '[cm**-3]'
                       ' for each tabulated temperature log10 T [K]. '
                       'All CIE at a given T: threshold is -np.inf; '
                       'all PIE at a given T: threshold in np.inf.'
                       ' Ion fractions are >= minjump_CPIEtrans'
                       ' times the CIE ion fraction at some '
                       ' at densities < the transition density.')
    todoc['ion'] = ion
    todoc['redshift'] = redshift
    todoc['useZ_log10sol'] = useZ_log10sol
    todoc['minjump_CPIEtrans'] = minjump_CPIEtrans
    todoc['minfactor_CPIEmix'] = minfactor_CPIEmix

    iontab = Linetable_PS20(ion, redshift, emission=False, 
                            vol=True, lintable=True)
    iontab.findiontable()
    logTK = iontab.logTK
    lognH = iontab.lognHcm3
    logZsol = iontab.logZsol
    refZi = np.where(np.isclose(useZ_log10sol, logZsol))[0][0]    
    tab_T_nH = iontab.iontable_T_Z_nH[:, refZi, :]

    lognHtrans_cm3 = np.array([getCPIEtrans(tab_T_nH[iT, :], 
                               lognH, 
                               minjump_CPIEtrans=minjump_CPIEtrans)
                               for iT in range(len(logTK))])
    
    # nHcut at high T (no PIE peak): factor difference with CIE at CIE max T
    icie = -6
    ciecurve = tab_T_nH[:, icie]
    ciemax_logTi = np.argmax(ciecurve)
    #ciemax_logT = logTK[ciemax_logTi]
    maxfrac = ciecurve[ciemax_logTi]
    ionfrac_nH = tab_T_nH[ciemax_logTi]
    fracrange_cie = (maxfrac / minfactor_CPIEmix, 
                     maxfrac * minfactor_CPIEmix)
    nHi_lim = np.where(np.logical_or(ionfrac_nH < fracrange_cie[0],
                                     ionfrac_nH > fracrange_cie[1]))[0][-1]
    lognHcut_hiT_cm3 = lognH[nHi_lim + 1]

    todoc['ionbalfile'] = iontab.ionbalfile
    todoc['ionbaltable_vol'] = True
    todoc['ionbaltable_lintable'] = True
    todoc['lognHcut_hiT_cm3'] = lognHcut_hiT_cm3
    todoc['logT_K'] = logTK
    todoc['lognHtrans_cm3'] = lognHtrans_cm3 
    return todoc


def get_ionclass_twolines(dct_lognH_logT, ion, redshift):
    todoc = {}
    lognHcut_cm3, logTcut_K, _todoc = get_cie_pie_nHT_twolines(ion, redshift)
    todoc.update(_todoc)
    todoc['ionclasses'] = ionclasses
    lognH = dct_lognH_logT['lognH_cm3']
    logT = dct_lognH_logT['logT_K']
    out = -1 * np.ones(lognH.shape, dtype=np.int8)
    
    hiT = logT > logTcut_K
    hinH = lognH > lognHcut_cm3
    out[np.logical_and(hiT, hinH)] = ionclasses['CIE']
    out[np.logical_and(hiT, np.logical_not(hinH))] = ionclasses['C+PIE']
    out[np.logical_not(hiT)] = ionclasses['PIE']
    #out[np.logical_and(np.logical_not(hiT), 
    #                   np.logical_not(hinH))] = ionclasses['lo']
    return out, 1, todoc
    
def get_ionclass_strawn21(dct_lognH_logT, ion, redshift):
    '''
    Strawn et al. (2021)-based division into CIE, PIE, and C+PIE gas
    only tested for Ne8 (last option for nH cuts relative to T_CIE_max
    transitions)
    '''
    todoc = {}
    _todoc = get_cie_pie_nHT_strawn21(ion, redshift)
    todoc.update(_todoc)
    todoc['ionclasses'] = ionclasses
    lognH_tocheck = dct_lognH_logT['lognH_cm3']
    logT_tocheck = dct_lognH_logT['logT_K']
    out = -1 * np.ones(lognH_tocheck.shape, dtype=np.int8)

    lognHcut_hiT_cm3 = _todoc['lognHcut_hiT_cm3']
    logTK = np.copy(_todoc['logT_K'])
    lognHtrans_cm3 = np.copy(_todoc['lognHtrans_cm3'])

    logTmaxi_PIE_K = np.where(lognHtrans_cm3 < np.inf)[0][0]
    logTmax_PIE_K = logTK[logTmaxi_PIE_K]
    todoc['logTmax_PIE_K'] = logTmax_PIE_K
    
    logTK = logTK[logTmaxi_PIE_K:]
    lognHtrans_cm3 = lognHtrans_cm3[logTmaxi_PIE_K:]
    logTmini_allCIE = np.where(lognHtrans_cm3 == -np.inf)[0][0]
    if lognHcut_hiT_cm3 > lognHtrans_cm3[logTmini_allCIE - 1]:
        # if the highest-T nH cut is to the right of the T_CIEmax
        # line in the phase diagram, just connect the line to all 
        # CIE temperatures
        logTKlast = logTK[logTmini_allCIE]
        lognHtranslast_cm3 = lognHcut_hiT_cm3
        logTK_interp = np.copy(logTK[:logTmini_allCIE + 1])
        lognHcm3_interp = np.copy(lognHtrans_cm3[:logTmini_allCIE + 1])
        logTK_interp[-1] = logTKlast
        lognHcm3_interp[-1] = lognHtranslast_cm3
    elif lognHcut_hiT_cm3 > np.max(lognHtrans_cm3):
        # the highest-T nH cut is to the right of all transition
        # densities -> effective just use the two-line case
        logTK_interp = np.array([logTK[0], logTK[-1]])
        lognHcm3_interp = np.array([lognHcut_hiT_cm3] * 2)
    else:
        # at each temperature, use the max CIE T transition as the
        # minimum tranisition density
        logTKlast = logTK[logTmini_allCIE]
        lognHtranslast_cm3 = lognHcut_hiT_cm3
        lognHtranslast_cm3 = np.maximum(lognHtranslast_cm3, lognHcut_hiT_cm3)
        logTK_interp = np.copy(logTK[:logTmini_allCIE + 1])
        lognHcm3_interp = np.copy(lognHtrans_cm3[:logTmini_allCIE + 1])
        logTK_interp[-1] = logTKlast
        lognHcm3_interp[-1] = lognHtranslast_cm3
    interpf = scpi.interp1d(logTK_interp, lognHcm3_interp, 
                            kind='linear', copy=True, bounds_error=False,
                            fill_value=(np.inf, lognHcut_hiT_cm3),
                            assume_sorted=False)
    
    hiT = logT_tocheck >= logTmax_PIE_K
    hinH = lognH_tocheck > interpf(logT_tocheck)
    out[np.logical_not(hiT)] = ionclasses['PIE']
    out[np.logical_and(hiT, hinH)] = ionclasses['CIE']
    out[np.logical_and(hiT, np.logical_not(hinH))] = ionclasses['C+PIE']
    #out[np.logical_and(np.logical_not(hiT), 
    #                   np.logical_not(hinH))] = ionclasses['lo']
    return out, 1, todoc


    


