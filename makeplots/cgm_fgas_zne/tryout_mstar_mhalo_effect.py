from tkinter import Y
import matplotlib.cm as cm
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np

import fire_an.makeplots.cgm_fgas_zne.readprop as rpr
import fire_an.makeplots.plot_utils as pu
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c

mdir = '/projects/b1026/nastasha/imgs/cgmprop/'

def plot_ne8mean_vs_mstar_mhalo(massset='m12'):
    df = rpr.readin_all_data(massset=massset, inclm12plus=True)
    df['meanNe8col'] = df['Ne8_numpart'] / (np.pi * df['Rvir_cm']**2)
    df['Mvir_Msun'] = df['Mvir_g'] / c.solar_mass
    df['Mstarcen'] = df['Mstarcen_g'] / c.solar_mass
    df = df.reset_index()

    simnames_sr = sl.m12_sr_all2 + sl.m13_hr_all2
    simnames_hr = sl.m12_hr_all2 + sl.m13_hr_all2 \
                  + sl.m12plus_f3nobh + sl.m12plus_f3nobh_lores
    simnames_f2 = sl.m12_f2md
    simnames_all = simnames_sr + simnames_hr + simnames_f2
    simnames_all = [sn for sn in simnames_all if sn not in sl.buglist2]
    snaps_hr = sl.snaps_hr
    snaps_sr = sl.snaps_sr
    snaps_f2md = sl.snaps_f2md
    
    vmin = np.log10(np.min(df['meanNe8col']))
    vmax = np.log10(np.max(df['meanNe8col']))
    __cmap = cm.get_cmap('viridis')
    _cmap = pu.truncate_colormap(__cmap, minval=0., maxval=1.)
    cmap = pu.paste_cmaps([__cmap], edges=[vmin, vmax])

    physmarkers = {'FIRE-2': 's',
                   'noBH': 'P',
                   'AGN-noCR': 'o',
                   'AGN-CR': 'd',
                   'noBH-m12+': 'X',
                   }
    panelsize = 2.5
    caxwidth = 0.2
    npanels = 7
    ncols = 4
    nrows = 2
    width_ratios = [panelsize] * ncols + [caxwidth]
    height_ratios = [panelsize] * nrows
    width = sum(width_ratios)
    height = sum(height_ratios)
    fontsize = 12

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(nrows=nrows, ncols=ncols + 1, 
                        width_ratios=width_ratios, 
                        height_ratios=height_ratios,
                        hspace=0., wspace=0.)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols])
            for i in range(7)]
    cax = fig.add_subplot(grid[:, ncols])

    for simname in simnames_all:
        f1 = df['simname'] == simname
        zv = np.log10(df.loc[f1, 'redshift'])
        zo = np.argsort(zv)
        xv = np.array(np.log10(df.loc[f1, 'Mvir_Msun']))[zo]
        yv = np.array(np.log10(df.loc[f1, 'Mstarcen']))[zo]
        _cv = np.array(np.log10(df.loc[f1, 'meanNe8col']))[zo]
        cv = cmap((_cv - vmin) / (vmax - vmin))
        
        physlabel = sl.physlabel_from_simname(simname)
        if physlabel == 'FIRE-2':
            axi = 0
        elif physlabel == 'noBH':
            axi = 1
        elif physlabel == 'noBH-m12+':
            if simname.endswith('r28000') or simname.endswith('r57000'):
                axi = 3
            else:
                axi = 2
        elif physlabel == 'AGN-noCR':
            axi = 4
        elif physlabel == 'AGN-CR':
            axi = 5
        _args = (xv, yv)
        _kwargs_ln = {'color': 'gray', 'linewidth': 0.5}
        _kwargs_sc = {'edgecolor': 'black',
                      's': 30., 'marker': physmarkers[physlabel]}
        
        axes[axi].plot(*_args, **_kwargs_ln)
        axes[6].plot(*_args, **_kwargs_ln)
        for _x, _y, _c in zip(xv, yv, cv):
            axes[axi].scatter((_x), (_y), facecolor=_c,
                              **_kwargs_sc)
            axes[6].scatter((_x), (_y), facecolor=_c,
                            **_kwargs_sc)      
    
    xlims = [ax.get_xlim() for ax in axes[:6]]
    xmin = min(xl[0] for xl in xlims)
    xmax = max(xl[1] for xl in xlims)
    ylims = [ax.get_ylim() for ax in axes[:6]]
    ymin = min(yl[0] for yl in ylims)
    ymax = max(yl[1] for yl in ylims)
    for axi, ax in enumerate(axes):
        doleft = axi % ncols == 0
        dobottom = axi >= npanels - ncols
        ax.tick_params(which='both', direction='in', top=True,
                       right=True, labelsize=fontsize - 1.,
                       labelbottom=dobottom, labelleft=doleft)
        ax.set_xlim((xmin, xmax))
        ax.set_ylim((ymin, ymax))
        if doleft:
            ax.set_ylabel('$\\log_{10} \\; \\mathrm{M}_{\\star, '
                          '\\mathrm{cen}} \\;'
                          '[\\mathrm{M}_{\\odot}]$',
                          fontsize=fontsize)
        if dobottom:
            ax.set_xlabel('$\\log_{10} \\; \\mathrm{M}_{'
                          '\\mathrm{vir}} \\;'
                          '[\\mathrm{M}_{\\odot}]$',
                          fontsize=fontsize)
        if axi == 0:
            axlabel = 'FIRE-2 noBH'
        elif axi == 1:
            axlabel = 'FIRE-3 noBH'
        elif axi == 2:
            axlabel = 'FIRE-3 noBH\nm12+ hi-res'
        elif axi == 3:
            axlabel = 'FIRE-3 noBH\nm12+ lo-res'
        elif axi == 4:
            axlabel = 'FIRE-3 AGN-noCR'
        elif axi == 5:
            axlabel = 'FIRE-3 AGN-CR'
        ax.text(0.05, 0.95, axlabel, fontsize=fontsize - 1,
                color='black', transform=ax.transAxes,
                horizontalalignment='left', verticalalignment='top')
    clabel = ('$\\log_{10} \\, \\langle\\mathrm{N}(\\mathrm{Ne\\,VIII})'
              '\\rangle \\; [\\mathrm{cm}^{-2}]$')
    pu.add_colorbar(ax=cax, vmin=vmin, vmax=vmax, cmap=cmap, 
                    fontsize=fontsize, clabel=clabel)
    plt.savefig(mdir + 'mstar_mvir_ne8av_m12.pdf', bbox_inches='tight')
        





