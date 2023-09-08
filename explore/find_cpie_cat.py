import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import string

from fire_an.ionrad.ion_utils import Linetable_PS20
import fire_an.utils.constants_and_units as c
import fire_an.utils.math_utils as mu

mdir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/ionmech/tests/'

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

def test_transfinder(ion, redshift):
    outname = f'test_CPIE_transition_nH_finder_{ion}_z{redshift:.2f}'
    outname = outname.replace('.', 'p')
    outname = mdir + outname + '.pdf'

    useZ_log10sol = 0.
    iontab = Linetable_PS20(ion, redshift, emission=False, 
                            vol=True, lintable=True)
    iontab.findiontable()
    logTK = iontab.logTK
    lognH = iontab.lognHcm3
    logZsol = iontab.logZsol
    refZi = np.where(np.isclose(useZ_log10sol, logZsol))[0][0]    
    tab_T_nH = iontab.iontable_T_Z_nH[:, refZi, :]

    tnHs = []
    for ti in range(len(logTK)):
        tnH = getCPIEtrans(tab_T_nH[ti, :], lognH, minjump_CPIEtrans=2.)
        tnHs.append(tnH)
    tnHs = np.array(tnHs)
    print(tnHs)
    transis = np.where(np.isfinite(tnHs))[0]
    hiT = logTK[transis[0]]
    loT = logTK[transis[-1]]
    nHplot = tnHs.copy()
    Tplot = logTK.copy()
    nHplot[:transis[0]] = lognH[-1]
    nHplot[transis[-1] + 1:] = lognH[0]
    Tplot[:transis[0]] = hiT
    Tplot[transis[-1] + 1:] = loT
    
    deltaT = np.average(np.diff(logTK))
    deltanH = np.average(np.diff(lognH))
    extent = (lognH[0] - 0.5 * deltanH, lognH[-1] + 0.5 * deltanH,
              logTK[0] - 0.5 * deltaT, logTK[-1] + 0.5 * deltaT)
    xlabel = ('$\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
              '\\; [\\mathrm{cm}^{-3}]$')
    ylabel = '$\\log_{10} \\, \\mathrm{T} \\; [\\mathrm{K}]$'
    
    fig = plt.figure(figsize=(5.5, 5.))
    grid = gsp.GridSpec(ncols=2, nrows=1, width_ratios=(10., 1.))
    ax = fig.add_subplot(grid[0, 0])
    cax = fig.add_subplot(grid[0, 1])
    fontsize = 12

    clevels = [0.001, 0.01, 0.1]
    linewidths = [2., 1., 0.5]
    ctargs = {'colors': 'gray', 'linestyles': 'solid'}
    ax.tick_params(which='both', direction='in', labelsize=fontsize - 1,
                       top=True, right=True, labelbottom=True,
                       labelleft=True)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    img = ax.imshow(np.log10(tab_T_nH), extent=extent, cmap='viridis', 
                    interpolation='nearest', origin='lower',
                    aspect='auto')
    ax.contour(lognH, logTK, tab_T_nH, origin='lower',
                levels=clevels, linewidths=linewidths, **ctargs)
    ax.plot(nHplot, Tplot, color='red', linewidth=1.5, linestyle='dashed')
        
    plt.colorbar(img, cax=cax, orientation='vertical', aspect=15.)
    cax.set_ylabel(f'$\\log_{{10}} \\, f({ion})$', fontsize=fontsize)
    
    handles = [mlines.Line2D((), (), linewidth=lw, label=f'{cl:.1e}',
                             color=ctargs['colors'], 
                             linestyle=ctargs['linestyles'])
               for cl, lw in zip(clevels, linewidths)]
    legtitle = (f'$f(\mathrm{{{ion}}})$')
    ax.legend(handles=handles, loc='lower right',
              fontsize=fontsize, 
              title=legtitle,
              title_fontsize=fontsize,
              handlelength=1.5)
    plt.savefig(outname, bbox_inches='tight')


    
def get_cie_pie_nHT(ion, redshift, useZ_log10sol=0.):
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

def get_ionclass(dct_lognH_logT, ion, redshift):

    todoc = {}
    lognHcut_cm3, logTcut_K, _todoc = get_cie_pie_nHT(ion, redshift)
    todoc.update(_todoc)
    todoc['ionclasses'] = ionclasses
    lognH = dct_lognH_logT['lognH_cm3']
    logT = dct_lognH_logT['logT_K']
    out = -1 * np.ones(lognH.shape, dtype=np.int8)
    
    hiT = logT > logTcut_K
    hinH = lognH > lognHcut_cm3
    out[np.logical_and(hiT, hinH)] = ionclasses['C+PIE']
    out[np.logical_and(hiT, np.logical_not(hinH))] = ionclasses['CIE']
    out[np.logical_not(hiT)] = ionclasses['PIE']
    #out[np.logical_and(np.logical_not(hiT), 
    #                   np.logical_not(hinH))] = ionclasses['lo']
    return out, 1, todoc


def plot_iontab_Zcomp(ion, redshift):
    '''
    result:
    for Ne8 at z=0.5 -- 1.0, PS20 tables, comparing different Z:
    relative differences are largest where ion fractions are low,
    and below 0.1 dex, except for a few single points with ion 
    fractions < 1e-3, which still have differences < 0.2 dex.
    '''
    outname = f'tablecomp_Zvalues_{ion}_z{redshift:.2f}'
    outname = outname.replace('.', 'p')
    outname = mdir + outname + '.pdf'

    iontab = Linetable_PS20(ion, redshift, emission=False, 
                            vol=True, lintable=True)    
    iontab.findiontable()
    logTK = iontab.logTK
    lognH = iontab.lognHcm3
    logZsol = iontab.logZsol
    tab_T_Z_nH = iontab.iontable_T_Z_nH 

    deltaT = np.average(np.diff(logTK))
    deltanH = np.average(np.diff(lognH))
    extent = (lognH[0] - 0.5 * deltanH, lognH[-1] + 0.5 * deltanH,
              logTK[0] - 0.5 * deltaT, logTK[-1] + 0.5 * deltaT)
    
    numpanels = len(logZsol)
    ncols = 4
    nrows = (numpanels - 1) // ncols + 1
    panelsize = 2.5
    cbarwidth = 0.2
    width_ratios = [panelsize] * ncols + [cbarwidth]
    height_ratios = [panelsize] * nrows
    figsize = (sum(width_ratios), sum(height_ratios))

    fig = plt.figure(figsize=figsize)
    grid = gsp.GridSpec(ncols=ncols + 1, nrows=nrows,
                        width_ratios=width_ratios, 
                        height_ratios=height_ratios,
                        hspace=0., wspace=0.)
    cax = fig.add_subplot(grid[:, -1])
    
    clabel = (f'$\\log_{{10}} \\, f(\\mathrm{{{ion}}}) \\, /\\,'
              f'f(\\mathrm{{{ion}}}, '
              '{\\mathrm{Z} = \\mathrm{Z}_{\\odot}}) $')
    xlabel = ('$\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
              '\\; [\\mathrm{cm}^{-3}]$')
    ylabel = '$\\log_{10} \\, \\mathrm{T} \\; [\\mathrm{K}]$'
    cmap = 'RdBu'
    refZ = 0.
    fontsize = 12
    
    refZi = np.where(np.isclose(refZ, logZsol))[0][0]
    refZlin = 10**refZ
    ref = np.copy(tab_T_Z_nH[:, refZi, :])
    clevels = [0.001, 0.01, 0.1]
    linewidths = [2., 1., 0.5]
    ctargs = {'colors': 'gray', 'linestyles': 'solid'}

    _rat = np.log10(tab_T_Z_nH / ref[:, np.newaxis, :])
    vmax = np.max(np.abs(_rat[np.isfinite(_rat)]))
    vmin = -1. * vmax
    print(vmin, vmax)
    
    axes = []
    for zi in range(numpanels):
        xi = zi % ncols
        yi = zi // ncols
        doleft = xi == 0
        dobottom = zi >= numpanels - ncols

        ax = fig.add_subplot(grid[yi, xi])
        axes.append(ax)
        ax.tick_params(which='both', direction='in', labelsize=fontsize - 1,
                       top=True, right=True, labelbottom=dobottom,
                       labelleft=doleft)
        if dobottom:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        if doleft:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        toplot = np.log10(tab_T_Z_nH[:, zi, :] / ref)
        img = ax.imshow(toplot, extent=extent, cmap=cmap, vmin=vmin,
                        vmax=vmax, interpolation='nearest', origin='lower',
                        aspect='auto')
        ax.contour(lognH, logTK, ref, origin='lower',
                   levels=clevels, linewidths=linewidths, **ctargs)
        
        axtitle = (f'$\\mathrm{{Z}} = 10^{{{logZsol[zi]:.1f}}} \\,'
                    '\\mathrm{Z}_{\\odot}$')
        ax.text(0.05, 0.95, axtitle, color='black',
                fontsize=fontsize, horizontalalignment='left',
                verticalalignment='top', transform=ax.transAxes)
        
    plt.colorbar(img, cax=cax, orientation='vertical', aspect=15.)
    cax.set_ylabel(clabel, fontsize=fontsize)
    
    handles = [mlines.Line2D((), (), linewidth=lw, label=f'{cl:.1e}',
                             color=ctargs['colors'], 
                             linestyle=ctargs['linestyles'])
               for cl, lw in zip(clevels, linewidths)]
    if np.isclose(refZlin, 1.):
        legtitle = (f'$f(\mathrm{{{ion}}}, \\mathrm{{Z}} ='
                    f'\\mathrm{{Z}}_{{\\odot}})$')
    else:
        legtitle = (f'$f(\mathrm{{{ion}}}, \\mathrm{{Z}} ='
                    f'{refZlin:.1f} \\mathrm{{Z}}_{{\\odot}})$')
    axes[0].legend(handles=handles, loc='lower right',
                   fontsize=fontsize - 2, 
                   title=legtitle,
                   title_fontsize=fontsize - 2,
                   handlelength=1.5)
    plt.savefig(outname, bbox_inches='tight')

def checkZdep(ion, redshifts):
    '''
    result:
    exact ion balance values for Ne VIII (z=0.5 -- 1) do depend a
    bit on Z, but not very much
    '''
    for redshift in redshifts:
        print(redshift, ':')
        iontab = Linetable_PS20(ion, redshift, emission=False, 
                                vol=True, lintable=True)    
        iontab.findiontable()
        logTK = iontab.logTK
        lognH = iontab.lognHcm3
        logZsol = iontab.logZsol
        tab_T_Z_nH = iontab.iontable_T_Z_nH 

        ref = tab_T_Z_nH[:, -2, :]
        difflist = []
        smalldifflist = []
        for zi in range(len(logZsol)):
            #if zi == 0:
            #    print(f'\tmin. Z ({logZsol[zi]}) table range:')
            #    print('\t',
            #          np.min(tab_T_Z_nH[:, zi, :]),
            #          np.max(tab_T_Z_nH[:, zi, :]))
            if np.all(ref == tab_T_Z_nH[:, zi, :]):
                continue
            elif np.allclose(ref, tab_T_Z_nH[:, zi, :]):
                smalldifflist.append(logZsol[zi])
            else:
                maxreldiff = np.max(tab_T_Z_nH[:, zi, :][ref != 0.]
                                    / ref[ref != 0.])
                minreldiff = np.min(tab_T_Z_nH[:, zi, :][ref != 0.]
                                    / ref[ref != 0.])
                maxabsdiff = np.max(np.abs(tab_T_Z_nH[:, zi, :] - ref))
                difflist.append((logZsol[zi], 
                                 maxreldiff, minreldiff, maxabsdiff))
        print(f'\tsmall differences: {smalldifflist}')
        print(f'\tdifferences: {difflist}')
            

def check_cutvals(ion, redshift):
    '''
    result:
    for Ne8 at z=0.5 -- 1.0, PS20 tables, comparing different Z:
    looks generally ok, but there is a regime at the lowest rho/T of
    the CIE quadrant which deviates from CIE meaningfully, and would
    be classified as PIE but Clayton Strawn's definition. 
    Check phase diagrams to see if this actually matters for the haloes
    in our set.
    '''
    outname = f'check_cutvals_{ion}_z{redshift:.2f}'
    outname = outname.replace('.', 'p')
    outname = mdir + outname + '.pdf'

    iontab = Linetable_PS20(ion, redshift, emission=False, 
                            vol=True, lintable=True)    
    iontab.findiontable()
    logTK = iontab.logTK
    lognH = iontab.lognHcm3
    logZsol = iontab.logZsol
    tab_T_Z_nH = iontab.iontable_T_Z_nH 

    deltaT = np.average(np.diff(logTK))
    deltanH = np.average(np.diff(lognH))
    extent = (lognH[0] - 0.5 * deltanH, lognH[-1] + 0.5 * deltanH,
              logTK[0] - 0.5 * deltaT, logTK[-1] + 0.5 * deltaT)
    
    numpanels = len(logZsol)
    ncols = 5
    nrows = (numpanels - 1) // ncols + 1
    panelsize = 2.5
    cbarwidth = 0.2
    width_ratios = [panelsize] * ncols + [cbarwidth]
    height_ratios = [panelsize] * nrows * 2
    figsize = (sum(width_ratios), sum(height_ratios))

    fig = plt.figure(figsize=figsize)
    grid = gsp.GridSpec(ncols=ncols + 1, nrows=nrows * 2,
                        width_ratios=width_ratios, 
                        height_ratios=height_ratios,
                        hspace=0., wspace=0.)
    fcax = fig.add_subplot(grid[0, -1])
    rcax = fig.add_subplot(grid[1, -1])
    
    rclabel = (f'$\\log_{{10}} \\, f(\\mathrm{{{ion}}}) \\, /\\,'
               f'f(\\mathrm{{{ion}}}, '
               '\\mathrm{CIE})$')
    fclabel = (f'$\\log_{{10}} \\, f(\\mathrm{{{ion}}}) $')
    xlabel = ('$\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
              '\\; [\\mathrm{cm}^{-3}]$')
    ylabel = '$\\log_{10} \\, \\mathrm{T} \\; [\\mathrm{K}]$'
    diffcmap = 'RdBu'
    fraccmap = 'viridis'
    refZ = 0.
    fontsize = 12
    
    #refZi = np.where(np.isclose(refZ, logZsol))[0][0]
    #refZlin = 10**refZ
    ionfraclevels = [0.001, 0.01, 0.1]
    difflevels = [0.1, 0.3, 0.5]
    linewidths = [2., 1., 0.5]
    ctargs = {'colors': 'gray', 'linestyles': 'solid'}
    icie = - 6

    cierat = np.log10(tab_T_Z_nH / tab_T_Z_nH[:, :, icie][:, :, np.newaxis])
    vmaxcierat = np.max(np.abs(cierat[np.isfinite(cierat)]))
    vmaxcierat = min(vmaxcierat, 1.)
    vmincierat = -1. * vmaxcierat

    frac_T_Z_nH = np.log10(tab_T_Z_nH)
    vmaxfrac = np.max(frac_T_Z_nH[np.isfinite(frac_T_Z_nH )])
    vminfrac = np.min(frac_T_Z_nH[np.isfinite(frac_T_Z_nH )])
    vminfrac = max(vminfrac, vmaxfrac - 6.)

    lognHcut, logTcut, _ = get_cie_pie_nHT(ion, redshift, 
                                           useZ_log10sol=refZ)

    axes = []
    for zi in range(numpanels):
        xi = zi % ncols
        yi = zi // ncols
        doleft = xi == 0
        dobottom = zi >= numpanels - ncols

        fax = fig.add_subplot(grid[2 * yi, xi])
        rax = fig.add_subplot(grid[2 * yi + 1, xi])
        axes.append(fax)
        axes.append(rax)
        rax.tick_params(which='both', direction='in', labelsize=fontsize - 1,
                        top=True, right=True, labelbottom=dobottom,
                        labelleft=doleft)
        fax.tick_params(which='both', direction='in', labelsize=fontsize - 1,
                        top=True, right=True, labelbottom=False,
                        labelleft=doleft)
        if dobottom:
            rax.set_xlabel(xlabel, fontsize=fontsize)
        if doleft:
            rax.set_ylabel(ylabel, fontsize=fontsize)
            fax.set_ylabel(ylabel, fontsize=fontsize)

        fimg = fax.imshow(frac_T_Z_nH[:, zi, :], extent=extent,
                          cmap=fraccmap, 
                          vmin=vminfrac, vmax=vmaxfrac, 
                          interpolation='nearest', origin='lower',
                          aspect='auto')
        fax.contour(lognH, logTK, tab_T_Z_nH[:, zi, :], origin='lower',
                    levels=ionfraclevels, linewidths=linewidths, **ctargs)
        
        rimg = rax.imshow(cierat[:, zi, :], extent=extent,
                          cmap=diffcmap, 
                          vmin=vmincierat, vmax=vmaxcierat, 
                          interpolation='nearest', origin='lower',
                          aspect='auto')
        rax.contour(lognH, logTK, np.abs(cierat[:, zi, :]), origin='lower',
                    levels=difflevels, linewidths=linewidths, **ctargs)
        axtitle = (f'$\\mathrm{{Z}} = 10^{{{logZsol[zi]:.1f}}} \\,'
                    '\\mathrm{Z}_{\\odot}$')
        fax.text(0.05, 0.95, axtitle, color='red',
                 fontsize=fontsize, horizontalalignment='left',
                 verticalalignment='top', transform=fax.transAxes)
        
    plt.colorbar(fimg, cax=fcax, orientation='vertical', aspect=15.,
                 extend='min')
    fcax.set_ylabel(fclabel, fontsize=fontsize)
    plt.colorbar(rimg, cax=rcax, orientation='vertical', aspect=15.,
                 extend='both')
    rcax.set_ylabel(rclabel, fontsize=fontsize)

    for axi, ax in enumerate(axes):
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.plot(xlim, (logTcut,) * 2, linestyle='dashed', color='black')
        ax.plot((lognHcut,) * 2, (logTcut, ylim[1]), 
                linestyle='dashed', color='black')
        if axi == 0:
            ax.text(xlim[0], logTcut, 'PIE',
                    color='red', fontsize=fontsize,
                    horizontalalignment='left',
                    verticalalignment='top')
            ax.text(xlim[0], logTcut, 'C+PIE',
                    color='red', fontsize=fontsize,
                    horizontalalignment='left',
                    verticalalignment='bottom')
            ax.text(xlim[1], logTcut, 'CIE',
                    color='red', fontsize=fontsize,
                    horizontalalignment='right',
                    verticalalignment='bottom')
    
    handles = [mlines.Line2D((), (), linewidth=lw, label=f'{cl:.1e}',
                             color=ctargs['colors'], 
                             linestyle=ctargs['linestyles'])
               for cl, lw in zip(ionfraclevels, linewidths)]
    legtitle = (f'$f(\mathrm{{{ion}}})$')
    axes[0].legend(handles=handles, loc='lower right',
                   fontsize=fontsize - 2, 
                   title=legtitle,
                   title_fontsize=fontsize - 2,
                   handlelength=1.5)
    
    handles = [mlines.Line2D((), (), linewidth=lw, label=f'{cl:.1f}',
                             color=ctargs['colors'], 
                             linestyle=ctargs['linestyles'])
               for cl, lw in zip(difflevels, linewidths)]
    legtitle = (f'$|\\Delta \\log_{{10}} \\, f(\mathrm{{{ion}}})|$')
    axes[1].legend(handles=handles, loc='lower right',
                   fontsize=fontsize - 2, 
                   title=legtitle,
                   title_fontsize=fontsize - 2,
                   handlelength=1.5)
    
    plt.savefig(outname, bbox_inches='tight')

