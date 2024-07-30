import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
import fire_an.spectra.plots.select_obsdata as so


def plot_effect_varpars(filens, filelabels, title=None,
                        outname=None):
    cs = tc.tol_cset('bright')
    fontsize = 12
    #binsize = 0.15
    fmtsim = {'color': 'gray',
              'marker': 'o',
              's': 10. + 4. * len(filens),
              'alpha': 1.}

    fig = plt.figure(figsize=(9., 3.))
    grid = gsp.GridSpec(nrows=1, ncols=3, wspace=0.5)
    ax1 = fig.add_subplot(grid[0, 0])
    ax2 = fig.add_subplot(grid[0, 1])
    ax3 = fig.add_subplot(grid[0, 2])
    ax1.set_xlabel("$\\log_{10} \\, \\mathrm{N}_{\\mathrm{true}}"
                   "\\; [\\mathrm{cm}^{-2}]$", fontsize=fontsize)
    ax1.set_ylabel("det. fraction of N", fontsize=fontsize)
    ax2.set_xlabel("$\\log_{10} \\, \\mathrm{N}_{\\mathrm{det}}"
                   "\\; [\\mathrm{cm}^{-2}]$", fontsize=fontsize)
    ax2.set_ylabel("$|\\langle v - v_{\\mathrm{gal}} \\rangle_{\\tau}|"
                   "\\; [\\mathrm{km}\\,\\mathrm{s}^{-1}]$",
                   fontsize=fontsize)
    ax3.set_xlabel("$\\log_{10} \\, \\mathrm{N}_{\\mathrm{det}}"
                   "\\; [\\mathrm{cm}^{-2}]$", fontsize=fontsize)
    ax3.set_ylabel("$\\sigma(v) \\; [\\mathrm{km}\\,\\mathrm{s}^{-1}]$",
                   fontsize=fontsize)
    tickpars = {'which': 'both', 'direction': 'in',
                'labelsize': fontsize - 1, 'right': True,
                'top': True}
    ax1.tick_params(**tickpars)
    ax2.tick_params(**tickpars)
    ax3.tick_params(**tickpars)
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)

    for i, (label, filen) in enumerate(zip(filelabels, filens)):
        color = cs[i]
        fmtsim.update(color=color)
        s = fmtsim['s']
        fmtsim.update(s=s - 4)
        simdata = pd.read_csv(filen, sep='\t')
        ndet = simdata['Ndet_logcm2']
        ntot = simdata['Ntot_logcm2']
        vgal = np.abs(simdata['vgal_kmps'])
        sigmav = simdata['sigmav_kmps']

        xps = [ntot, ndet, ndet]
        yps = [ntot / ndet, vgal, sigmav]
        for xp, yp, ax in zip(xps, yps, [ax1, ax2, ax3]):
            ax.scatter(xp, yp,
                       label=label, rasterized=True, **fmtsim)
            # lolim = np.min(xp)
            # lolim = max(lolim, 13.0)
            # xlo = np.floor(lolim / binsize) * binsize
            # xhi = np.ceil(np.max(xp) / binsize) * binsize
            # xbins = np.arange(xlo, xhi + 0.5 * binsize, binsize)
            # xcen = 0.5 * (xbins[:-1] + xbins[1:])
            # yperc, _, xsel = pu.get_perc_and_points(xp, yp, xbins,
            #                     percentiles=(10., 50., 90.),
            #                     mincount_x=15,
            #                     getoutliers_y=False, 
            #                     getmincounts_x=False,
            #                     x_extremes_only=False)
            # ax.plot(xcen[xsel], (yperc[1:-1, 0])[xsel],
            #         color=color, linestyle='dashed',
            #         linewidth=1.2, label='10, 90%')
            # ax.plot(xcen[xsel], (yperc[1:-1, 2])[xsel], 
            #         color=color, linestyle='dashed',
            #         linewidth=1.2, label=None)
            # ax.plot(xcen[xsel], (yperc[1:-1, 1])[xsel], 
            #         color=color, linestyle='solid',
            #         linewidth=1.5, label='median')
    for ax in [ax1, ax2, ax3]:
        xl = ax.get_xlim()
        ax.set_xlim(max(13.25, xl[0]), xl[1])
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles=handles, fontsize=fontsize - 1)
    #ax2.legend(handles=handles[1:3],
    #           fontsize=fontsize - 1)

    if outname is not None:
        fig.savefig(outname, bbox_inches='tight')

def plotset_effect_varpars():
    datadir = '/projects/b1026/nastasha/imgs/spectra/test4/'
    outdir = datadir

    # minN effect
    _filens = ['sigmav_testset4_minN13.50_lsf10.0_maxb100_cutf0.05.dat',
               'sigmav_testset4_minN13.75_lsf10.0_maxb100_cutf0.05.dat',
               'sigmav_testset4_minN14.00_lsf10.0_maxb100_cutf0.05.dat']
    filens = [datadir + filen for filen in _filens]
    filelabels = ['N > 13.5', 'N > 13.75', 'N > 14.0']
    title = 'Varying component det. lim.'
    outname = outdir + ('roughcompdet_effect_minN_testset4'
                        '_lsf10.0_maxb100_cutf0.05.pdf')
    plot_effect_varpars(filens, filelabels, title=title,
                        outname=outname)
    
    _filens = ['sigmav_testset4_minN13.50_lsf10.0_maxb100_cutf0.01.dat',
               'sigmav_testset4_minN13.50_lsf10.0_maxb100_cutf0.05.dat',
               'sigmav_testset4_minN13.50_lsf10.0_maxb100_cutf0.25.dat']
    filens = [datadir + filen for filen in _filens]
    filelabels = ['cutfrac 0.01', 'cutfrac 0.05', 'cutfrac 0.25']
    title = 'Varying component edge cutoff, at det lim. log N 13.5'
    outname = outdir + ('roughcompdet_effect_cutf_testset4'
                        'minN13.50_lsf10.0_maxb100.pdf')
    plot_effect_varpars(filens, filelabels, title=title,
                        outname=outname)
    
    _filens = ['sigmav_testset4_minN14.00_lsf10.0_maxb100_cutf0.01.dat',
               'sigmav_testset4_minN14.00_lsf10.0_maxb100_cutf0.05.dat',
               'sigmav_testset4_minN14.00_lsf10.0_maxb100_cutf0.25.dat']
    filens = [datadir + filen for filen in _filens]
    filelabels = ['cutfrac 0.01', 'cutfrac 0.05', 'cutfrac 0.25']
    title = 'Varying component edge cutoff, at det lim. log N 14'
    outname = outdir + ('roughcompdet_effect_cutf_testset4'
                        'minN14.00_lsf10.0_maxb100.pdf')
    plot_effect_varpars(filens, filelabels, title=title,
                        outname=outname)
    

        
