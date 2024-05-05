import matplotlib as mpl
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import fire_an.makeplots.plot_utils as pu
from fire_an.makeplots.ionmech.plotpds import readpd
import fire_an.makeplots.tol_colors as tc
from fire_an.ionrad.ion_utils import Linetable_PS20
import fire_an.simlists as sl 

def plotpds_talk(physlabel='FIRE-2', ic='m12f', redshift=1.0,
                 rrange_rvir=(0.1, 1.)):
    '''
    plot phase diagram for an example halo: gas, Ne8
    '''
    ddir = '/projects/b1026/nastasha/hists/phasediagrams_all2/'
    fbase = ('hist_rcen_temperature_density_by_{wtstr}_'
             '{simname}_snap{snapnum}'
             '_bins1_v1.hdf5')
    filetemp = ddir + fbase
    
    mdir = '/projects/b1026/nastasha/imgs/pcie/test_minmaxTnHcuts/'
    outname = (f'phasediag_Ne8_mass_{physlabel}_{ic}'
               f'z{redshift:.1f}_{rrange_rvir[0]:.2f}_to_'
               f'{rrange_rvir[1]:.2f}_Rvir')
    outname = mdir + outname.replace('.', 'p') + '.pdf'

    zopts = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
    snapi = np.where(np.isclose(zopts, redshift))[0][0]
    snap_sr = sl.snaps_sr[snapi]
    snap_hr = sl.snaps_hr[snapi]
    snap_f2 = sl.snaps_f2md[snapi]
    sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    sims_f2 = sl.m12_f2md

    allsims = sl.m12_sr_all2 + sl.m12_hr_all2 + sl.m12_f2md + \
              sl.m13_sr_all2 + sl.m13_hr_all2
    simnames = [sn for sn in allsims
                if (sn not in sl.buglist2) 
                    and (sl.physlabel_from_simname(sn)) == physlabel
                    and (sl.ic_from_simname(sn)) == ic]
    snaps = [snap_f2 if sn in sims_f2
             else snap_hr if sn in sims_hr
             else snap_sr if sn in sims_sr
             else None
             for sn in simnames]
    simname = simnames[0]
    snap = snaps[0]
    
    ion = 'Ne8'
    otherions = [] #'O6']
    oilabels ={'O6': 'O VI / O'}
    wts = ['gasmass', 'Ne8'] #  'Neon', 
    wtlabels = {'gasmass': 'gas mass',
                'Neon': 'Ne',
                'Ne8': 'Ne VIII',
                'gasvol': 'Vol.'}
    colors = tc.tol_cset('vibrant')
    colors = [colors[0], colors[1], colors[3]]
    wtcolors = {wt: col for wt, col in zip(wts, colors)}
    wtlinestyles = {wt: 'solid' for wt in wts}
    ibpluscolors = [colors[len(wts) + i] for i in range(len(otherions))] 
    ibcolor = colors[len(wts) + len(otherions)]
    iblinestyle = 'dashed'
    encllevels = [0.99, 0.9, 0.5]
    iblevels = [0.001, 0.01, 0.1]
    linewidths = [0.5, 1.0, 1.5,]
    refZ = 0. # log solar mass fraction units

    iontab = Linetable_PS20(ion, redshift, emission=False, 
                            vol=True, lintable=True)    
    iontab.findiontable()
    logTK = iontab.logTK
    lognH = iontab.lognHcm3
    logZsol = iontab.logZsol
    tab_T_Z_nH = iontab.iontable_T_Z_nH 
    iZ = np.where(np.isclose(refZ, logZsol))[0][0]
    tab_T_nH = tab_T_Z_nH[:, iZ, :] 
    othertabs = {}
    for oi in otherions:
        othertabs[oi] = {}
        _iontab = Linetable_PS20(oi, redshift, emission=False, 
                                 vol=True, lintable=True)    
        _iontab.findiontable()
        _logTK = _iontab.logTK
        _lognH = _iontab.lognHcm3
        _logZsol = _iontab.logZsol
        _tab_T_Z_nH = _iontab.iontable_T_Z_nH 
        _iZ = np.where(np.isclose(refZ, _logZsol))[0][0]
        _tab_T_nH = _tab_T_Z_nH[:, _iZ, :]
        othertabs[oi]['logTK'] = _logTK
        othertabs[oi]['lognH'] = _lognH
        othertabs[oi]['tab_T_nH'] = _tab_T_nH

    npanels = len(simnames)
    ncols = 1
    nrows = 1
    #ncols = 4
    #nrows = (npanels - 1) // ncols + 1
    panelsize = 5.
    cbarwidth = 0.4
    width_ratios = [panelsize] * ncols + [cbarwidth]
    height_ratios = [panelsize] * nrows
    figsize = (sum(width_ratios), sum(height_ratios))

    fig = plt.figure(figsize=figsize)
    grid = gsp.GridSpec(ncols=ncols + 1, nrows=nrows,
                        width_ratios=width_ratios, 
                        height_ratios=height_ratios,
                        hspace=0., wspace=0.)
    cax = fig.add_subplot(grid[:, -1])
    
    #clabel = ('$\\log_{10} \\, \\mathcal{N}(\\mathrm{Ne\\,VIII\\,ions})$')
    clabel = ('$\\log_{10} \\, \\mathrm{gas\\,mass} \\; [\\mathrm{g}]$')
    xlabel = ('$\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
              '\\; [\\mathrm{cm}^{-3}]$')
    ylabel = '$\\log_{10} \\, \\mathrm{T} \\; [\\mathrm{K}]$'
    _cmap = mpl.cm.get_cmap('gist_yarg')
    cmap = pu.truncate_colormap(_cmap, minval=0., maxval=0.5)
    #print(cmap(0.), cmap(0.5), cmap(1.))
    fontsize = 12

    axes = []
    for axi, (simname, snap) in enumerate(zip(simnames, snaps)):
        xi = axi % ncols
        yi = axi // ncols
        doleft = xi == 0
        dobottom = axi >= npanels - ncols

        ax = fig.add_subplot(grid[yi, xi])
        axes.append(ax)
        ax.tick_params(which='both', direction='in', labelsize=fontsize - 1,
                       top=True, right=True, labelbottom=dobottom,
                       labelleft=doleft)
        if dobottom:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        if doleft:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        
        for wi, wt in enumerate(wts):
            dat = readpd(filetemp.format(simname=simname, snapnum=snap,
                                        wtstr=wt), 
                         rrange_rvir=rrange_rvir)
            if wi == 0:
                toplot = np.log10(dat['linpd'])
                vmax = np.max(toplot)
                vmin = max(vmax - 6., np.min(toplot))
                extent = (dat['nHbins'][0], dat['nHbins'][-1],
                          dat['Tbins'][0],  dat['Tbins'][-1])
                img = ax.imshow(toplot, extent=extent,
                          cmap=cmap, vmin=vmin, vmax=vmax, 
                          interpolation='nearest', origin='lower',
                          aspect='auto')
            pu.add_2dhist_contours(ax, dat['linpd'].T, 
                                   [dat['nHbins'], dat['Tbins']],
                                   (0, 1),
                                   histlegend=False, fraclevels=True, 
                                   levels=encllevels, colors=wtcolors[wt],
                                   linestyles=wtlinestyles[wt], 
                                   linewidths=linewidths)
           
        #if axi == 0:
        #    posy = 0.05
        #    va = 'bottom'
        #else:
        #    posy = 0.95
        #    va = 'top'
        #ax.text(0.05, posy, sl.ic_from_simname(simname), color='black',
        #        fontsize=fontsize, horizontalalignment='left',
        #        verticalalignment=va, transform=ax.transAxes)
        
    plt.colorbar(img, cax=cax, orientation='vertical', aspect=15.,
                 extend='min')
    cax.set_ylabel(clabel, fontsize=fontsize)

    xlims = [ax.get_xlim() for ax in axes]
    xlim = [min([l[0] for l in xlims]), max([l[1] for l in xlims])]
    xlim[1] = min([xlim[1], 1.])
    ylims = [ax.get_ylim() for ax in axes]
    ylim = [min([l[0] for l in ylims]), max([l[1] for l in ylims])]
    ylim[0] = max([ylim[0], 3.7])
    for axi, ax in enumerate(axes):
        #xlim = ax.get_xlim()
        #ylim = ax.get_ylim()
        for _color, oi in zip(ibpluscolors, otherions):
            ax.contour(othertabs[oi]['lognH'],
                       othertabs[oi]['logTK'],
                       othertabs[oi]['tab_T_nH'], 
                       origin='lower',
                       levels=iblevels[1:], linewidths=linewidths[1:], 
                       colors=_color,
                       linestyles=iblinestyle)
        ax.contour(lognH, logTK, tab_T_nH, origin='lower',
                   levels=iblevels, linewidths=linewidths, colors=ibcolor,
                   linestyles=iblinestyle)
        
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    
    phandles = [mlines.Line2D((), (), linewidth=lw, 
                              label=f'{100. * cl:.0f}%',
                              color='gray', 
                              linestyle='solid')
                for cl, lw in zip(encllevels, linewidths)]
    legtitle = 'encl. frac.'
    l0 = axes[0].legend(handles=phandles, loc='upper center',
                        fontsize=fontsize - 1, 
                        title=legtitle,
                        title_fontsize=fontsize - 1,
                        handlelength=1.0,
                        handletextpad=0.4)
    legtitle = 'weight'
    whandles = [mlines.Line2D((), (), linewidth=2., 
                             label=wtlabels[wt],
                             color=wtcolors[wt], 
                             linestyle='solid',)
               for wt in wts]
    l1 = axes[0].legend(handles=whandles, loc='upper left',
                        fontsize=fontsize - 1, 
                        title=legtitle,
                        title_fontsize=fontsize - 1,
                        handlelength=1.,
                        ncol=1, handletextpad=0.3,
                        columnspacing=0.7)
    
    legtitle = (f'Ne VIII / Ne')
    ihandles = [mlines.Line2D((), (), linewidth=lw, 
                              label=f'{ibl:.0e}'.replace('e-0', 'e-'),
                              color=ibcolor, 
                              linestyle=iblinestyle)
                for ibl, lw in zip(iblevels, linewidths)]
    ihandles = ihandles +\
               [mlines.Line2D((), (), linewidth=linewidths[1], 
                              label=oilabels[oi],
                              color=ibc,
                              linestyle='solid')
                for ibc, oi in zip(ibpluscolors, otherions)]
    axes[0].legend(handles=ihandles, loc='upper right',
                   fontsize=fontsize - 1, 
                   title=legtitle,
                   title_fontsize=fontsize - 1,
                   handlelength=1.5, handletextpad=0.4)
    axes[0].add_artist(l0)
    axes[0].add_artist(l1)

    plt.savefig(outname, bbox_inches='tight')
