import h5py
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np

import fire_an.ionrad.get_cieranges as gcie
import fire_an.makeplots.tol_colors as tc
import fire_an.simlists as sl
import fire_an.utils.opts_locs as ol

profdatafn = ('/projects/b1026/nastasha/plotdata/'
              'radprof3d_nH_T_ZNe_by_vol_Ne8_opt1.hdf5')
mdir = '/projects/b1026/nastasha/imgs/3dprof/'

def plot_fewprof(simnames, snapnums, yranges={}, axsettitles=None,
                 outname=None):
    weights = ['gasvol', 'Ne8']
    wlabels = {'gasvol': 'Volume-wtd',
               'Ne8': 'Ne VIII-wtd'}
    cset = tc.tol_cset('vibrant')
    wcolors = {'gasvol': cset.black,
               'Ne8': cset.cyan}
    yqtys = ['temperature', 'hdens', 'NeonAbundance']
    ylabels = {'temperature': ('$\\log_{10} \\, \\mathrm{T}'
                               '\\; [\\mathrm{K}]$'),
               'hdens': ('$\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
                         '\\; [\\mathrm{cm}^{-3}]$'),
               'NeonAbundance': ('$\\log_{10} \\, \\mathrm{Z}_{\\mathrm{Ne}}'
                                 '\\; [\\mathrm{Z}_{\\mathrm{Ne}, \\odot}]$')
               }   
    ynorms = {'NeonAbundance': np.log10(ol.solar_abunds_ea['neon']),
              'temperature': 0.,
              'hdens': 0.,
              }
    percvals = [0.1, 0.5, 0.9]

    npsets = len(simnames)
    ncols = min(npsets, 4)
    nrows = len(yqtys) * ((npsets - 1) // ncols + 1)
    panelsize = 2.5
    fontsize = 12

    fig = plt.figure(figsize=(panelsize * ncols, panelsize * nrows))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows, hspace=0.0,
                        wspace=0.0)
    axsets = [[fig.add_subplot(grid[i // ncols + j, i % ncols]) 
               for j in range(len(yqtys))] for i in range(npsets)] 
    yminmax = {yqty: (np.inf, -np.inf) for yqty in yqtys}
    for seti, (axset, simn, snap) in enumerate(zip(axsets, simnames, snapnums)):
        for yi, (ax, yqty) in enumerate(zip(axset, yqtys)):
            doleft = seti % ncols == 0
            dobottom = (yi == len(yqtys) - 1) and (seti >= npsets - ncols)
            ax.tick_params(which='both', direction='in', labelsize=fontsize - 1,
                           top=True, right=True, labelbottom=dobottom,
                           labelleft=doleft)
            if dobottom:
                ax.set_xlabel('$\\mathrm{r}_{\\mathrm{cen}} \\,/\\,'
                              '\\mathrm{R}_{\\mathrm{vir}}$', 
                              fontsize=fontsize)
            if doleft:
                ax.set_ylabel(ylabels[yqty], 
                              fontsize=fontsize)
            ynorm = ynorms[yqty]
            for wi, weight in enumerate(weights):
                with h5py.File(profdatafn, 'r') as f:
                    grp = f[f'{simn}/snap_{snap}/{weight}_weighted/{yqty}']
                    rbins = grp['rbins'][:]
                    ylo = grp[f'perc_{percvals[0]:.2f}'][:]
                    ymid = grp[f'perc_{percvals[1]:.2f}'][:]
                    yhi = grp[f'perc_{percvals[2]:.2f}'][:]
                rcens = 0.5 * (rbins[:-1] + rbins[1:])
                
                _label = wlabels[weight] + ', median'
                ax.plot(rcens, ymid - ynorm, color=wcolors[weight], 
                        linestyle='solid', linewidth=2., label=_label)
                if wi == 0:
                    _label = wlabels[weight] + (f'; ${100 * percvals[0]:.0f}'
                                                ' \\endash'
                                                f'{100 * percvals[2]:.0f}$%')
                    ax.fill_between(rcens, ylo - ynorm, yhi - ynorm, 
                                    color=wcolors[weight], alpha=0.3,
                                    label=_label)
                else:
                    _label = wlabels[weight] + (f'; ${100 * percvals[0]:.0f}'
                                                ', '
                                                f'{100 * percvals[2]:.0f}$%')
                    ax.plot(rcens, ylo - ynorm, color=wcolors[weight], 
                            linestyle='dashed', linewidth=1.5, label=_label)
                    ax.plot(rcens, yhi - ynorm, color=wcolors[weight], 
                            linestyle='dashed', linewidth=1.5)
            if yqty == 'temperature':
                tcies = gcie.cieranges1['Ne8']
                for tcie in tcies:
                    ax.axhline(tcie, color='black', linestyle='dotted', linewidth=1.)

            ylim = ax.get_ylim()
            if yqty in yminmax:
                yminmax[yqty] = (min(yminmax[yqty][0], ylim[0]),
                                 max(yminmax[yqty][1], ylim[1]))
        if axsettitles is not None:
            axset[0].set_title(axsettitles[seti], fontsize=fontsize)
    
    for axset in axsets:
        for ax, yqty in zip(axset, yqtys):
            if yqty in yranges:
                yr = yranges[yqty]
            else:
                yr = yminmax[yqty]
            ax.set_ylim(yr)
    
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')
    
def plotsets_fewprof():
    sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    sims_f2md = sl.m12_f2md
    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    snaps_f2md = sl.snaps_f2md

    sims_all = sims_sr + sims_hr + sims_f2md
    ics_run = ['m12f', 'm13h113', 'm13h206'] # clean sample
    
    sortorder = {'FIRE-2': 0,
                 'noBH': 1,
                 'AGN-noCR': 2,
                 'AGN-CR': 3}
    def sortkey(simname):
        pl = sl.physlabel_from_simname(simname)
        return sortorder[pl]
    
    for ic in ics_run:
        sims_use = [simn for simn in sims_all if simn.startswith(ic)]
        sims_use.sort(key=sortkey)
        physlabels = [sl.plotlabel_from_physlabel[
                          sl.physlabel_from_simname(simn)]
                      for simn in sims_use]
        snaps = [(snaps_sr if simn in sims_sr 
                  else snaps_hr if simn in sims_hr
                  else snaps_f2md if simn in sims_f2md
                  else None)
                 for simn in sims_use]
        snaps = [[snaps[j][i] for j in range(len(sims_use))]
                  for i in range(6)]
        for zi, z in enumerate(np.arange(1.0, 0.45, -0.1)):
            _snaps = snaps[zi]
            outname = (f'rprof3d_T_nH_ZNe_by_Ne8_Vol_{ic}_z{z:.1f}')
            outname = mdir + outname.replace('.', 'p') + '.pdf'
            plot_fewprof(sims_use, _snaps, yranges={}, 
                         axsettitles=physlabels, outname=outname)
            
    



    