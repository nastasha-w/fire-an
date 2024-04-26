import numpy as np

import fire_an.ionrad.ion_utils as iu
import fire_an.utils.math_utils as mu

def calc_cieranges(threshval, ions, zs, nHval_logpcm3=1.,
                   metval_logsolar=0.):
    cieranges = dict()
    ciemaxs = dict()
    for ion in ions:
        refrange = None
        refmax = None
        for z in zs:
            tab = iu.Linetable_PS20(ion, z, lintable=True)
            tab.findiontable()
            tab_logTK = tab.logTK
            tab_lognHpcm3 = tab.lognHcm3
            tab_logZsol = tab.logZsol 
            tab_T_Z_nH = tab.iontable_T_Z_nH
            nHi = np.argmin(np.abs(tab_lognHpcm3 - nHval_logpcm3))
            lZi = np.argmin(np.abs(tab_logZsol - metval_logsolar))
            tab_T = tab_T_Z_nH[:, lZi, nHi]
            mti = np.argmax(tab_T)
            maxlogTK = tab_logTK[mti]
            maxv = tab_T[mti]
            print('max ion fraction: ', maxv)
            minr = maxv * threshval
            logTKrange = mu.find_intercepts(tab_T, tab_logTK, minr)
            if refmax is None:
                refmax = maxlogTK
            elif not np.isclose(refmax, maxlogTK, atol=0.05, rtol=5e-2):
                msg = (f'for ion {ion}, redshifts {zs}, max. T difference'
                       f'was too large: {refmax}, {maxlogTK}')
                raise RuntimeError(msg)
            if refrange is None:
                refrange = logTKrange
            elif not np.allclose(refrange, logTKrange, atol=0.05, rtol=5e-2):
                msg = (f'for ion {ion}, redshifts {zs}, T range difference'
                       f' was too large: {refrange}, {logTKrange}')
                raise RuntimeError(msg)
        cieranges[ion] = refrange
        ciemaxs[ion] = refmax
    return ciemaxs, cieranges

# calc_cieranges(0.1, ['O6', 'O7', 'O8', 'Ne8', 'Ne9', 'Ne10', 'Fe17'], 
#                [0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 
#                nHval_logpcm3=1., metval_logsolar=0.)
ciemaxs1 = {'O6': 5.5, 'O7': 5.9, 'O8': 6.4, 
            'Ne8': 5.8, 'Ne9': 6.2, 'Ne10': 6.7, 
            'Fe17': 6.6}
cieranges1 = {'O6': np.array([5.30320019, 5.77078197]), 
              'O7': np.array([5.40875786, 6.48462019]), 
              'O8': np.array([6.09073196, 6.79483222]), 
              'Ne8': np.array([5.61993395, 6.15170053]),
              'Ne9': np.array([5.71938275, 6.78576501]), 
              'Ne10': np.array([6.33208108, 7.16666652]), 
              'Fe17': np.array([6.23186434, 6.99623786])}