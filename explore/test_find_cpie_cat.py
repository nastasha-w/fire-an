
import matplotlib.gridspec as gsp
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import numpy as np
import string

from fire_an.ionrad.ion_utils import Linetable_PS20
import fire_an.explore.find_cpie_cat as fcp
import fire_an.makeplots.tol_colors as tc
import fire_an.utils.constants_and_units as c
import fire_an.utils.math_utils as mu



mdir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/ionmech/tests/'


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
        tnH = fcp.getCPIEtrans(tab_T_nH[ti, :], lognH, minjump_CPIEtrans=2.)
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

# Can we just use the table values for solar Z?
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

# how much do ion balance values depend on Z (more quantitative)
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

def check_cutvals_twolines(ion, redshift):
    '''
    result:
    for Ne8 at z=0.5 -- 1.0, PS20 tables, comparing different Z:
    looks generally ok, but there is a regime at the lowest rho/T of
    the CIE quadrant which deviates from CIE meaningfully, and would
    be classified as PIE but Clayton Strawn's definition. 
    Check phase diagrams to see if this actually matters for the haloes
    in our set. 
    -> checked, and yes, for the m12 CGM, there seems to be quite a of
    gas in the CIE/PIE 'corner' region, including the bottom left corner
    of the 'two straight lines' CIE box.
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

    lognHcut, logTcut, _ = fcp.get_cie_pie_nHT_twolines(ion, redshift, 
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

def check_ionclass_strawn21_grid(ion, redshift):
    '''
    Looks like it's following the tables as intended
    (for 'Ne8', redshift 0.5)
    '''
    valsnH = np.arange(-8., 4., 0.025)
    valsT = np.arange(0.5, 9.5, 0.025)

    dct_lognH_logT = {'lognH_cm3': np.array([valsnH] * len(valsT)).flatten(),
                      'logT_K': np.array([[tv] * len(valsnH)
                                          for tv in valsT]).flatten()}
    ionclasses, _, todoc = fcp.get_ionclass_strawn21(dct_lognH_logT, 
                                                     ion, redshift)
    print(todoc)
    toplot = ionclasses.reshape((len(valsT), len(valsnH)))

    fig = plt.figure(figsize=(5.5, 5.))
    ax = fig.add_subplot(1, 1, 1)
    fontsize = 12
    ax.tick_params(which='both', top=True, right=True,
                   labelsize=fontsize - 1., direction='in')
    ax.set_xlabel('$\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
                  '\\; [\\mathrm{cm}^{-3}]$', fontsize=fontsize)
    ax.set_ylabel('$\\log_{10} \\, \\mathrm{T}'
                  '\\; [\\mathrm{K}]$', fontsize=fontsize)
    
    dnH = np.average(np.diff(valsnH))
    dT = np.average(np.diff(valsT))
    extent = (valsnH[0] - 0.5 * dnH, valsnH[-1] + 0.5 * dnH,
              valsT[0] - 0.5 * dT, valsT[-1] + 0.5 * dT)
    colors = tc.tol_cset('bright')
    labeltoval = fcp.ionclasses
    valtolabel = {val: key for key, val in labeltoval.items()}
    valtocolor = {key: col for key, col in zip(valtolabel.keys(), colors)}
    
    handles = []
    for catval in valtolabel:
        cval = np.array(mcolors.to_rgba(valtocolor[catval]))
        _toplot = np.zeros(toplot.shape + (4,), dtype=cval.dtype)
        valmatch = np.where(toplot == catval)
        _toplot[valmatch[0], valmatch[1], :] = cval
        plt.imshow(_toplot, interpolation='nearest', origin='lower',
                   extent=extent)
        handles.append(mpatch.Patch(color=cval, label=valtolabel[catval]))
    ax.legend(handles=handles, loc='lower right', fontsize=fontsize)
    
    # add Strawn et al. transition line
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
        tnH = fcp.getCPIEtrans(tab_T_nH[ti, :], lognH, minjump_CPIEtrans=2.)
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

    ax.plot(nHplot, Tplot, color='black', linewidth=1.5, linestyle='dashed')
    

def check_cutvals_strawn21(ion, redshift):
    '''
    result:
    for Ne8 at z=0.5 -- 1.0, PS20 tables, comparing different Z:
    '''
    outname = f'check_cutcurve_strawn21_{ion}_z{redshift:.2f}'
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

    lognHcut, logTcut, _ = fcp.get_cie_pie_nHT_twolines(ion, redshift, 
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

