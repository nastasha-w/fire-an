'''
plot phase diagrams with different weights
'''

import h5py
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import fire_an.explore.find_cpie_cat as fcpie
import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
from fire_an.ionrad.ion_utils import Linetable_PS20
import fire_an.simlists as sl 


def readpd(filen, rrange_rvir=(0.1, 1.)):
    with h5py.Filen(filen, 'r') as f:
        linpd = f['histogram/histogram'][:]
        if bool(f['histogram'].attrs['log']):
            linpd = 10**linpd
        rbins_rvir = f['axis_0/bins'][:]
        ri0 = np.where(np.isclose(rbins_rvir, rrange_rvir[0]))[0][0]
        ri1 = np.where(np.isclose(rbins_rvir, rrange_rvir[1]))[0][0]
        linpd = np.sum(linpd[ri0 : ri1, :, :], axis=0)

        tbins = rbins_rvir = f['axis_1/bins'][:]
        nHbins = rbins_rvir = f['axis_2/bins'][:]
    return {'linpd': linpd, 'Tbins': tbins, 'nHbins': nHbins}
        

def checkNe8pds_physmodel(physlabel, massset='m12', redshift=0.5,
                          rrange_rvir=(0.1, 1.)):
    ddir = '/projects/b1026/nastasha/hists/phasediagrams_all2/'
    fbase = ('hist_rcen_temperature_density_by_{wtstr}_'
             '{simname}_snap{snapnum}'
             '_bins1_v1.hdf5')
    filetemp = ddir + fbase
    
    mdir = '/projects/b1026/nastasha/imgs/pcie/test_minmaxTnHcuts/'
    outname = (f'phasediag_Ne8_Neon_mass_{physlabel}_{massset}'
               f'z{redshift:.1f}_{rrange_rvir[0]:.2f}_to_'
               f'{rrange_rvir[1]:.2f}_Rvir')
    outname = mdir + outname.replace('.', 'p') + '.pdf'

    zopts = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
    snapi = np.where(np.isclose(zopts, redshift))[0][0]
    snap_sr = sl.snaps_sr[snapi]
    snap_hr = sl.snaps_hr[snapi]
    snap_f2 = sl.snaps_f2md[snapi]
    sims_sr = sl.m13_hr_all2 + sl.m13_hr_all2
    sims_hr = sl.m12_hr_all2 + sl.m12_hr_all2
    sims_f2 = sl.m12_f2md

    allsims = sl.m12_sr_all2 + sl.m12_hr_all2 + sl.m12_f2md + \
              sl.m13_sr_all2 + sl.m13_sr_all2
    simnames = [sn for sn in allsims
                if (sn not in sl.buglist1) 
                    and sl.physlabel_from_simname(sn) == physlabel
                    and (sl.ic_from_simname(sn)).startswith(massset)]
    snaps = [snap_f2 if sn in sims_f2
             else snap_hr if sn in sims_hr
             else snap_sr if sn in sims_sr
             else None
             for sn in simnames]
    
    ion = 'Ne8'
    wts = ['gasmass', 'Neon', 'Ne8']
    wtlabels = {'gasmass': 'gas',
                'Neon': 'Neon',
                'Ne8': 'Ne VIII',
                'gasvol': 'Volme'}
    colors = tc.tol_cset('vibrant')
    wtcolors = {wt: col for wt, col in zip(wts, colors)}
    wtlinestyles = {wt: 'solid' for wt in wts}
    ibcolor = colors[len(wts)]
    iblinestyle = 'dashed'
    encllevels = [0.99, 0.9, 0.5]
    iblevels = [0.001, 0.01, 0.1]
    linewidths = [1.5, 1., 0.5]
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

    ionpds = {simname: readpd(filetemp.format(simname=simname, snapnum=snap,
                                              weight=ion), 
                              rrange_rvir=rrange_rvir)
              for simname, snap in zip(simnames, snaps)}
    vmax = max([np.max(ionpds[simname]['linpd']) for simname in simnames])
    vmax = np.log10(vmax)
    vmin = vmax - 10.

    npanels = len(simnames)
    ncols = 4
    nrows = (npanels - 1) // ncols + 1
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
    cax = fig.add_subplot(grid[:, -1])
    
    clabel = ('$\\log_{{10}} \\, \\mathrm{# ions}$')
    xlabel = ('$\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
              '\\; [\\mathrm{cm}^{-3}]$')
    ylabel = '$\\log_{10} \\, \\mathrm{T} \\; [\\mathrm{K}]$'
    cmap = pu.truncate_colormap('gist_yarg', minval=0., maxval=0.7)
    fontsize = 12

    lognHcut, logTcut, _ = fcpie.get_cie_pie_nHT(ion, redshift, 
                                                 useZ_log10sol=refZ)

    axes = []
    for axi, (simname, snap) in enumerate(zip(simnames, snaps)):
        xi = axi % ncols
        yi = axi // ncols
        doleft = xi == 0
        dobottom = axi >= npanels - ncols

        ax = fig.add_subplot(grid[yi, xi])
        axes.append(ax)
        ax.tick_params(which='both', direction='in', labelsize=fontsize - 1,
                       top=True, right=True, labelbottom=False,
                       labelleft=doleft)
        if dobottom:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        if doleft:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        
        iondat = ionpds[simname]
        toplot = np.log10(iondat['linpd'])
        extent = (iondat['nHbins'][0], iondat['nHbins'][1],
                  iondat['Tbins'][0], iondat['Tbins'][1])
        img = ax.imshow(toplot, extent=extent,
                        cmap=cmap, vmin=vmin, vmax=vmax, 
                        interpolation='nearest', origin='lower',
                        aspect='auto')
        pu.add_2dhist_contours(iondat['linpd'].T, 
                               [iondat['nHbins'], iondat['Tbins']],
                               histlegend=False, fraclevels=True, 
                               levels=encllevels, colors=wtcolors[ion],
                               linestyles=wtlinestyles[ion], 
                               linewidths=linewidths)
        for wt in wts[1:]:
            dat = readpd(filetemp.format(simname=simname, snapnum=snap,
                                         weight=wt), 
                         rrange_rvir=rrange_rvir)
            pu.add_2dhist_contours(dat['linpd'].T, 
                                   [dat['nHbins'], dat['Tbins']],
                                   histlegend=False, fraclevels=True, 
                                   levels=encllevels, colors=wtcolors[wt],
                                   linestyles=wtlinestyles[wt], 
                                   linewidths=linewidths)
        ax.contour(lognH, logTK, tab_T_nH, origin='lower',
                   levels=iblevels, linewidths=linewidths, colors=ibcolor,
                   linewdiths=linewidths, linestyles=iblinestyle)

        ax.text(0.05, 0.95, sl.ic_from_simname(simname), color='black',
                fontsize=fontsize, horizontalalignment='left',
                verticalalignment='top', transform=ax.transAxes)
        
    plt.colorbar(img, cax=cax, orientation='vertical', aspect=15.,
                 extend='min')
    cax.set_ylabel(clabel, fontsize=fontsize)

    for axi, ax in enumerate(axes):
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.plot(xlim, (logTcut,) * 2, linestyle='dashed', color='black')
        ax.plot((lognHcut,) * 2, (logTcut, ylim[1]), 
                linestyle='dashed', color='black')
        if axi == 0:
            ax.text(xlim[0], logTcut, 'PIE',
                    color='black', fontsize=fontsize,
                    horizontalalignment='left',
                    verticalalignment='top')
            ax.text(xlim[0], logTcut, 'C+PIE',
                    color='black', fontsize=fontsize,
                    horizontalalignment='left',
                    verticalalignment='bottom')
            ax.text(xlim[1], logTcut, 'CIE',
                    color='black', fontsize=fontsize,
                    horizontalalignment='right',
                    verticalalignment='bottom')
    
    handles = [mlines.Line2D((), (), linewidth=lw, 
                             label=f'{100. * cl:.0f}%',
                             color='gray', 
                             linestyle='solid')
               for cl, lw in zip(encllevels, linewidths)]
    legtitle = (f'encl. frac.')
    l0 = axes[0].legend(handles=handles, loc='lower right',
                        fontsize=fontsize - 2, 
                        title=legtitle,
                        title_fontsize=fontsize - 2,
                        handlelength=1.5)
    handles = [mlines.Line2D((), (), linewidth=lw, 
                             label=f'{cl:.0e}',
                             color=ibcolor, 
                             linestyle=iblinestyle)
               for cl, lw in zip(iblevels, linewidths)]
    legtitle = (f'{ion} frac.')
    handles = [mlines.Line2D((), (), linewidth=1., 
                             label=wtlabels[wt],
                             color=wtcolors[wt], 
                             linestyle='solid')
               for wt in wts]
    axes[0].legend(handles=handles, loc='upper right',
                   ontsize=fontsize - 2, 
                   title=legtitle,
                   title_fontsize=fontsize - 2,
                   handlelength=1.5)
    axes[0].add_artist(l0)
    legtitle = (f'encl. qty.')
    axes[1].legend(handles=handles, loc='lower right',
                   fontsize=fontsize - 2, 
                   title=legtitle,
                   title_fontsize=fontsize - 2,
                   handlelength=1.5)

    plt.savefig(outname, bbox_inches='tight')