import h5py
import matplotlib.colors as mcolors
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib.transforms import Transform
import numpy as np

import fire_an.utils.math_utils as mu
import fire_an.makeplots.tol_colors as tc
import fire_an.simlists as sl

def getdata_hist(filen, yunit=1., xunit=1.):
    with h5py.File(filen, 'r') as f:
        hist = f['histogram/histogram'][:]
        islog = bool(f['histogram'].attrs['log'])
        if islog:
            hist = 10**hist
        rbins = f['axis_0/bins'][:] / xunit
        ybins = f['axis_1/bins'][:] / yunit   
        csmpath = 'Header/cosmopars'
        cosmopars = {key: val for key, val in f[csmpath].attrs.items()}
    basehist = hist / np.diff(rbins)[:, np.newaxis] \
                    / np.diff(ybins)[np.newaxis, :] \
                    / np.sum(hist)
    # for contour levels
    norm1hist = hist / np.sum(hist, axis=1)[:, np.newaxis]
    addorder = np.argsort(basehist, axis=1)
    backorder = np.argsort(addorder, axis=1)
    csum = np.cumsum(np.take_along_axis(norm1hist, addorder, 1), axis=1)
    fracvals_y = np.take_along_axis(csum, backorder, 1)
    return {'hist': hist, 'fracvals_y': fracvals_y, 
            'rbins': rbins, 'ybins': ybins, 
            'cosmopars': cosmopars}

def plot_vradtot_zphyscomp(filen_template, zfills, phicfills,
                           title=None, outname=None, ylabel=None,
                           physlabels=None):
    '''
    different redshifts in each panel, different panels for different
    phys. vars.
    same IC, weight between panels
    '''
    perclevels = [0.02, 0.1, 0.5, 0.9, 0.98]
    pstyles = [{'linestyle': 'dotted'},
               {'linestyle': 'dashed'},
               {'linestyle': 'solid'},
               {'linestyle': 'dashed'},
               {'linestyle': 'dotted'},
               ]
    fontsize = 12

    nz = len(zfills[0])
    _colors = tc.tol_cmap('rainbow_discrete', nz)
    zticks = np.linspace(0.5 / nz, 1. - 0.5 / nz, nz)
    colors = _colors(zticks)[::-1]
    # zip/list comprehension issues or something when using 
    # tol_cmap outputs directly
    colors = [mcolors.to_rgb(col) for col in colors]
    
    yunit = 1e5 # cm/s -> km/s
    data = [[getdata_hist(filen_template.format(**(zfill | phicfill)),
                          yunit=yunit, xunit=1.)
             for zfill in _zflist] 
            for _zflist, phicfill in zip(zfills, phicfills)]
    ncols = len(phicfills)
    panelsize = 2.5
    laxspace = 1.5
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] + [laxspace]
    hspace = 0.0
    wspace = 0.0
    height = sum(height_ratios) * (1. + hspace * sum(height_ratios) \
                * (len(height_ratios) - 1.)/ len(height_ratios))
    width = sum(width_ratios) * (1. + hspace * sum(width_ratios) \
                *  (len(width_ratios) - 1.)  / (len(width_ratios)))
    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols, nrows=2, wspace=wspace,
                        hspace=hspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[0, i]) for i in range(ncols)]
    lax = fig.add_subplot(grid[1, :])
    lax.axis('off')
    
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    xlabel = '$\\mathrm{r}_{\\mathrm{3D}} \\; [\\mathrm{R}_\\mathrm{vir}]$'
    for pi in range(ncols):
        zvals = []
        ax = axes[pi]
        doleft = pi == 0
        ax.tick_params(labelsize=fontsize - 1., which='both',
                       direction='in', top=True, right=True,
                       labelbottom=True, labelleft=doleft)
        if physlabels is not None:
            plab = physlabels[pi]
            ax.text(0.95, 0.95, plab, fontsize=fontsize,
                    horizontalalignment='right',
                    verticalalignment='top',
                    transform=ax.transAxes,
                    color='black')
        if doleft and ylabel is not None:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.axhline(0., linestyle='solid', color='gray', linewidth=1.)
        for zi in range(nz):
            _data = data[pi][zi]
            color = colors[zi]
            rbins = _data['rbins']
            rc = 0.5 * (rbins[:-1] + rbins[1:])
            ybins = _data['ybins']
            hist = _data['hist']
            zvals.append(_data['cosmopars']['z'])
            
            yv = mu.percentiles_from_histogram(hist, ybins, axis=1, 
                percentiles=np.array(perclevels))
            for pvi in range(len(perclevels)):
                ax.plot(rc, yv[pvi], color=color, **(pstyles[pvi]))
    ylims = [ax.get_ylim() for ax in axes]
    ymin = min([yl[0] for yl in ylims])
    ymin = max(ymin, -300.)
    ymax = max([yl[1] for yl in ylims])
    ymax = min(ymax, 500.)
    [ax.set_ylim((ymin, ymax)) for ax in axes]

    chandles = [mlines.Line2D((), (), color=color, label=f'$z={z:.1f}$')
                for color, z in zip(colors, zvals)]
    lshandles = [mlines.Line2D((), (), color='black', 
                               label=f'{pv * 100:.0f}%', **ls)
                 for ls, pv in zip(pstyles, perclevels)]
    lax.legend(handles=lshandles + chandles, fontsize=fontsize - 1., 
               ncol=int(np.floor(1.5 * ncols)),
               loc='upper center', bbox_to_anchor=(0.5, 0.65))
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotset_vradtot_zphyscomp(vtype='rad'):
    fdir = '/Users/nastasha/ciera/profiles/fire/vradtot_all2/'
    outdir = ('/Users/nastasha/ciera/projects_lead/fire3_ionabs/'
              'model3/vel3dcomp/')
    simnames = sl.m13_sr_all2 + sl.m13_hr_all2 \
               + sl.m12_sr_all2 + sl.m12_hr_all2
    sn_hr = sl.m13_hr_all2 + sl.m12_hr_all2
    sn_sr = sl.m13_sr_all2 + sl.m12_sr_all2
    icsets = []
    ylabel = (f'$v_{{\\mathrm{{{vtype}}}}}\\;'
              '[\\mathrm{km} \\, \\mathrm{s}^{-1}]$')
    for simname in simnames:
        ic = simname.split('_')[0]
        setind = np.where([np.all([_sn.startswith(ic) for _sn in icset])
                           for icset in icsets])[0]
        if len(setind) == 0:
            icsets.append([simname])
        else:
            icsets[setind[0]].append(simname)
    weights = ['gasmass', 'gasvol', 'Neon', 'Ne8', 'Mg10', 'O6']
    for icset in icsets:
        for weight in weights:
            ftemp = fdir + \
                    (f'hist_v{vtype}_by_{weight}_{{simname}}'
                    '_snap{snapnum}_bins1_v1_hvcen.hdf5')
            phicfills = [{'simname': simname} for simname in icset]
            zfills = [[{'snapnum': snap} for snap in sl.snaps_hr] 
                      if simname in sn_hr else 
                      [{'snapnum': snap} for snap in sl.snaps_sr]
                      if simname in sn_sr else 
                      None 
                      for simname in icset]
            physlabels = ['noAGN' if '_sdp1e10_' in simname
                          else 'AGN-CR' if 'MHDCRspec1' in simname
                          else 'AGN-noCR'
                          for simname in icset]
            ic = icset[0].split('_')[0]
            title = f'{ic} velocities, weighted by {weight}'
            outname = outdir + f'{ic}_zphyscomp_v{vtype}_by_{weight}_r3D.pdf'
            
            plot_vradtot_zphyscomp(ftemp, zfills, phicfills,
                                   title=title, outname=outname, 
                                   ylabel=ylabel, physlabels=physlabels)

def plotvradtot_zweightcomp(filen_template, zfills, wfills,
                            phicfills, wlabels=None,
                            physlabels=None, title=None,
                            outname=None, ylabel=None):
    '''
    different redshifts in each panel, different plots for
    different ics
    columns: weights
    rows: phys
    recommended: use 2 or 3 sample redshifts
    '''
    perclevels = [0.1, 0.9]
    pstyle = {'linestyle': 'dashed'}
    avstyle = {'linestyle': 'solid'}
    fontsize = 12

    nz = len(zfills[0])
    _colors = tc.tol_cmap('rainbow_discrete', nz)
    zticks = np.linspace(0.5 / nz, 1. - 0.5 / nz, nz)
    colors = _colors(zticks)[::-1]
    # zip/list comprehension issues or something when using 
    # tol_cmap outputs directly
    colors = [mcolors.to_rgb(col) for col in colors]
    
    yunit = 1e5 # cm/s -> km/s
    data = [[[getdata_hist(filen_template.format(**(zfill | phicfill | wfill)),
                          yunit=yunit, xunit=1.)
              for zfill in _zflist] 
             for wfill in wfills]
            for _zflist, phicfill in zip(zfills, phicfills)]
    ncols = len(wfills)
    nrows = len(phicfills)
    panelsize = 2.5
    laxspace = 1.5
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] * nrows + [laxspace]
    hspace = 0.0
    wspace = 0.0
    height = sum(height_ratios) * (1. + hspace * sum(height_ratios) \
                * (len(height_ratios) - 1.)/ len(height_ratios))
    width = sum(width_ratios) * (1. + hspace * sum(width_ratios) \
                *  (len(width_ratios) - 1.)  / (len(width_ratios)))
    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows + 1, wspace=wspace,
                        hspace=hspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [[fig.add_subplot(grid[j, i]) for i in range(ncols)]
            for j in range(nrows)]
    lax = fig.add_subplot(grid[nrows, :])
    lax.axis('off')
    
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    xlabel = '$\\mathrm{r}_{\\mathrm{3D}} \\; [\\mathrm{R}_\\mathrm{vir}]$'
    for wi in range(ncols):
        for pi in range(nrows):
            ax = axes[pi][wi]
            doleft = wi == 0
            ax.tick_params(labelsize=fontsize - 1., which='both',
                           direction='in', top=True, right=True,
                           labelbottom=True, labelleft=doleft)
            if physlabels is not None:
                plab = physlabels[pi]
                ax.text(0.95, 0.95, plab, fontsize=fontsize,
                        horizontalalignment='right',
                        verticalalignment='top',
                        transform=ax.transAxes,
                        color='black')
            if pi == 0 and wlabels is not None:
                ax.set_title(wlabels[wi], fontsize=fontsize)
            if doleft and ylabel is not None:
                ax.set_ylabel(ylabel, fontsize=fontsize)
            if pi == nrows - 1:
                ax.set_xlabel(xlabel, fontsize=fontsize)
            anyle0 = False
            zvals = []
            for zi in range(nz):
                _data = data[pi][wi][zi]
                color = colors[zi]
                rbins = _data['rbins']
                rc = 0.5 * (rbins[:-1] + rbins[1:])
                ybins = _data['ybins']
                hist = _data['hist']
                yc = 0.5 * (ybins[:-1] + ybins[1:])
                av = np.einsum('j,ij->i', yc, hist) / np.sum(hist, axis=1)
                #av = np.average(np.tile(yc[np.newaxis, :], (len(rc), 1)), 
                #                weights=hist, axis=1)
                zvals.append(_data['cosmopars']['z'])
                anyle0 |= ybins[0] < 0.
                
                yv = mu.percentiles_from_histogram(hist, ybins, axis=1, 
                    percentiles=np.array(perclevels))
                for pvi in range(len(perclevels)):
                    ax.plot(rc, yv[pvi], color=color, **pstyle)
                ax.plot(rc, av, color=color, **avstyle)
            if anyle0:
                ax.axhline(0., linestyle='solid', color='gray', linewidth=1.)
            
    ylims = [ax.get_ylim() for l1 in axes for ax in l1]
    ymin = min([yl[0] for yl in ylims])
    ymin = max(ymin, -300.)
    ymax = max([yl[1] for yl in ylims])
    ymax = min(ymax, 500.)
    [ax.set_ylim((ymin, ymax)) for l1 in axes for ax in l1]

    chandles = [mlines.Line2D((), (), color=color, label=f'$z={z:.1f}$')
                for color, z in zip(colors, zvals)]
    plab = ', '.join([f'{pv * 100:.0f}' for pv in perclevels]) + '%'
    lshandles = [mlines.Line2D((), (), color='black', 
                               label=plab, **pstyle),
                 mlines.Line2D((), (), color='black', 
                               label='wtd. mean', **avstyle),
                 ]
    lax.legend(handles=lshandles + chandles, fontsize=fontsize - 1., 
               ncol=int(np.floor(1.5 * ncols)),
               loc='upper center', bbox_to_anchor=(0.5, 0.65))
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotset_vradtot_zweightcomp(vtype='rad'):
    fdir = '/Users/nastasha/ciera/profiles/fire/vradtot_all2/'
    outdir = ('/Users/nastasha/ciera/projects_lead/fire3_ionabs/'
              'model3/vel3dcomp/')
    simnames = sl.m13_sr_all2 + sl.m13_hr_all2 \
               + sl.m12_sr_all2 + sl.m12_hr_all2
    sn_hr = sl.m13_hr_all2 + sl.m12_hr_all2
    sn_sr = sl.m13_sr_all2 + sl.m12_sr_all2
    icsets = []
    ylabel = (f'$v_{{\\mathrm{{{vtype}}}}}\\;'
              '[\\mathrm{km} \\, \\mathrm{s}^{-1}]$')
    for simname in simnames:
        ic = simname.split('_')[0]
        setind = np.where([np.all([_sn.startswith(ic) for _sn in icset])
                           for icset in icsets])[0]
        if len(setind) == 0:
            icsets.append([simname])
        else:
            icsets[setind[0]].append(simname)
    weights = ['gasmass', 'gasvol', 'Neon', 'Ne8', 'Mg10', 'O6']
    wlabs = ['gas mass', 'gas vol.', 'Neon', 'Ne8', 'Mg10', 'O6']
    snapsel = [0, 3, 5]
    snaps_hr = [sl.snaps_hr[i] for i in snapsel]
    snaps_sr = [sl.snaps_sr[i] for i in snapsel]
    for icset in icsets:
        wfills = [{'weight': weight} for weight in weights]
        ftemp = fdir + \
                (f'hist_v{vtype}_by_{{weight}}_{{simname}}'
                    '_snap{snapnum}_bins1_v1_hvcen.hdf5')
        phicfills = [{'simname': simname} for simname in icset]
        zfills = [[{'snapnum': snap} for snap in snaps_hr] 
                  if simname in sn_hr else 
                  [{'snapnum': snap} for snap in snaps_sr]
                  if simname in sn_sr else 
                  None 
                  for simname in icset]
        physlabels = ['noAGN' if '_sdp1e10_' in simname
                        else 'AGN-CR' if 'MHDCRspec1' in simname
                        else 'AGN-noCR'
                        for simname in icset]
        ic = icset[0].split('_')[0]
        title = f'{ic} {vtype}. velocities'
        outname = outdir + f'{ic}_physweightcomp_v{vtype}_by_r3D.pdf'
            
        plotvradtot_zweightcomp(ftemp, zfills, wfills,
                                phicfills, wlabels=wlabs,
                                physlabels=physlabels, title=title,
                                outname=outname, ylabel=ylabel)

def plotvradtot_icphyscomp(filen_template, icphfills):
    '''
    same weight and redshift in the plot, compare different ic/phys
    between. (new row for new phys)
    '''
