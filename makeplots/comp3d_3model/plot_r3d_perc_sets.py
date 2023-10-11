import h5py
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np

import fire_an.ionrad.get_cieranges as gcie
import fire_an.makeplots.tol_colors as tc
import fire_an.simlists as sl
import fire_an.utils.opts_locs as ol

profdatafn = ('/projects/b1026/nastasha/plotdata/'
              'radprof3d_nH_T_ZNe_by_vol_Ne8_opt2.hdf5')
mdir = '/projects/b1026/nastasha/imgs/3dprof/'

def plot_profset(simnames, snapnums, yranges={}, axsettitles=None,
                 outname=None):
    '''
    simnames, snapnums: matching lists of lists; inner lists are combined
    '''
    weights = ['gasvol', 'Ne8']
    wlabels = {'gasvol': 'V-wtd',
               'Ne8': 'Ne8-wtd'}
    cset = tc.tol_cset('vibrant')
    wcolors = {'gasvol': cset.black,
               'Ne8': cset.magenta}
    yqtys = ['temperature', 'hdens', 'NeonAbundance']
    ylabels = {'temperature': ('$\\log_{10} \\, \\mathrm{T}'
                               '\\; [\\mathrm{K}]$'),
               'hdens': ('$\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
                         '\\; [\\mathrm{cm}^{-3}]$'),
               'NeonAbundance': ('$\\log_{10} \\, \\mathrm{Z}_{\\mathrm{Ne}}'
                                 '\\; [\\mathrm{Z}_{\\mathrm{Ne}, \\odot}]$')
               }   
    ynorms = {'NeonAbundance': -2.9008431  ,
              'temperature': 0.,
              'hdens': 0.,
              }
    percvals = [0.1, 0.5, 0.9]

    npsets = len(simnames)
    if npsets <= 6:
        ncols = npsets
    else:
        # even things out a bit if more rows are needed
        nrows = (npsets - 1) // 6 + 1
        ncols = (npsets - 1) // nrows + 1
    nrows = len(yqtys) * ((npsets - 1) // ncols + 1)
    panelsize = 2.5
    fontsize = 12

    fig = plt.figure(figsize=(panelsize * ncols, panelsize * nrows))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows, hspace=0.0,
                        wspace=0.0)
    axsets = [[fig.add_subplot(grid[(i // ncols) * len(yqtys) + j, i % ncols])
               for j in range(len(yqtys))] for i in range(npsets)] 
    yminmax = {yqty: (np.inf, -np.inf) for yqty in yqtys}
    for seti, (axset, simsub, snapsub) in enumerate(zip(axsets, simnames, 
                                                        snapnums)):
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
                rbins = None
                ylo = []
                ymid = []
                yhi = []
                for simn, snap in zip(simsub, snapsub):
                    with h5py.File(profdatafn, 'r') as f:
                        grp = f[f'{simn}/snap_{snap}/{weight}_weighted/{yqty}']
                        _rbins = grp['rbins'][:]
                        if rbins is None:
                            rbins = _rbins
                        else:
                            if not np.allclose(rbins, _rbins):
                                raise RuntimeError('rbins mismatch within'
                                                   'panel set: '
                                                   f'{simsub}, {snapsub}')
                        ylo.append(grp[f'perc_{percvals[0]:.2f}'][:])
                        ymid.append(grp[f'perc_{percvals[1]:.2f}'][:])
                        yhi.append(grp[f'perc_{percvals[2]:.2f}'][:])
                    rcens = 0.5 * (rbins[:-1] + rbins[1:])
                ylom = np.median(np.array(ylo), axis=0)
                yhim = np.median(np.array(yhi), axis=0)
                ymidm = np.median(np.array(ymid), axis=0)
                yerr_lo = ymidm - np.min(np.array(ymid), axis=0)
                yerr_hi = np.max(np.array(ymid), axis=0) - ymidm
                _label = wlabels[weight] + ', med.'
                # errorevery needs a newer matplotlib version, I guess
                errsel = slice(wi, None, len(weights))
                ax.errorbar(rcens, ymidm - ynorm, yerr=np.zeros(len(rcens)),
                            color=wcolors[weight],
                            linestyle='solid', linewidth=2., label=_label)
                ax.errorbar(rcens[errsel], (ymidm - ynorm)[errsel], 
                            yerr=(yerr_lo[errsel], yerr_hi[errsel]),
                            color=wcolors[weight], linestyle='none', 
                            linewidth=2.)
                if wi == 0:
                    _label = wlabels[weight] + (f'; ${100 * percvals[0]:.0f}'
                                                ' \\endash'
                                                f'{100 * percvals[2]:.0f}$%')
                    ax.fill_between(rcens, ylom - ynorm, yhim - ynorm, 
                                    color=wcolors[weight], alpha=0.3,
                                    label=_label)
                else:
                    _label = wlabels[weight] + (f'; ${100 * percvals[0]:.0f}'
                                                ', '
                                                f'{100 * percvals[2]:.0f}$%')
                    ax.plot(rcens, ylom - ynorm, color=wcolors[weight], 
                            linestyle='dashed', linewidth=1.5, label=_label)
                    ax.plot(rcens, yhim - ynorm, color=wcolors[weight], 
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
            ax.grid(True)
    axsets[0][1].legend(fontsize=fontsize - 2., loc='upper right',
                        handlelength=1.5)
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotsets_physmodel():
    sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    sims_f2md = sl.m12_f2md
    sims_m12plus = sl.m12plus_f3nobh + sl.m12plus_f3nobh_lores
    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    snaps_f2md = sl.snaps_f2md

    sims_all = sims_sr + sims_hr + sims_f2md + sims_m12plus
    sims_all = [sn for sn in sims_all if sn not in sl.buglist2]
    physlabels = [sl.physlabel_from_simname(sn) for sn in sims_all]
    ics = [sl.ic_from_simname(sn) for sn in sims_all]
    simsets = {}
    for sn, ic, pm in zip(sims_all, ics, physlabels):
        catlabel = ic[:3] + '_' + pm
        if catlabel in simsets:
            simsets[catlabel].append(sn)
        else:
            simsets[catlabel] = [sn]
    #print(simsets)

    yranges = {'NeonAbundance': (-3.0, 0.5),
               'hdens': (-5.5, -0.5),
               'temperature': (5., 7.2)}
    for catlabel in simsets:
        simnames = simsets[catlabel]
        snaps = [(snaps_sr if simn in sims_sr 
                  else snaps_hr if simn in sims_hr
                  else snaps_f2md if simn in sims_f2md
                  else snaps_hr if simn in sims_m12plus
                  else None)
                 for simn in simnames]
        sims_use = [[sn] * 6 for sn in simnames]
        ics = [sl.ic_from_simname(sn) for sn in simnames]
        outname = (f'rprof3dset_T_nH_ZNe_by_Ne8_Vol_z0p5_to_1p0_{catlabel}')
        outname = mdir + outname.replace('.', 'p') + '.pdf'
        plot_profset(sims_use, snaps, yranges=yranges, axsettitles=ics,
                     outname=outname)
        