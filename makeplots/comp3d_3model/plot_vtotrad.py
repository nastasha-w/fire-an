import h5py
import matplotlib.colors as mcolors
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import fire_an.utils.math_utils as mu
import fire_an.makeplots.tol_colors as tc

def getdata_hist(filen, yunit=1., xunit=1.):
    with h5py.File(filen, 'r') as f:
        hist = f['histogram/histogram'][:]
        islog = bool(f['histogram'].attrs['log'])
        if islog:
            hist = 10**hist
        rbins = f['axis_0/bins'][:] / xunit
        ybins = f['axis_1/bins'][:] / yunit   
        csmpath = 'Header/inputpars/cosmopars'
        cosmopars = {key: val for key, val in f[csmpath].items()}
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
                           title=None, outname=None, ylabel=None):
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
    colors = _colors(zticks)
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
    laxspace = 1.0
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] + [laxspace]
    hspace = 0.0
    wspace = 0.0
    height = sum(height_ratios) \
             * (1. + hspace * sum(height_ratios) / (len(height_ratios) - 1.))
    width = sum(width_ratios) \
            * (1. + hspace * sum(width_ratios) / (len(width_ratios) - 1.))
    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols, nrows=2, wspace=wspace,
                        hspace=hspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = fig.add_subplot(grid[0, i] for i in range(ncols))
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
        if doleft and ylabel is not None:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.set_xlabel(xlabel, fontsize=fontsize)
        for zi in range(nz):
            _data = data[pi][zi]
            color = colors[zi]
            rbins = _data['rbins']
            rc = 0.5 * (rbins[:-1] + rbins[:1])
            ybins = _data['ybins']
            hist = _data['hist']
            zvals.append(_data['cosmopars']['z'])
            
            yv = mu.percentiles_from_histogram(hist, ybins, axis=-1, 
                percentiles=np.array(perclevels))
            for pvi in range(len(perclevels)):
                ax.plot(rc, yv[pvi], color=color, **(pstyles[pvi]))
    chandles = [mlines.line2D((), (), color=color, label=f'$z={z:.1f}$')
                for color, z in zip(colors, zvals)]
    lshandles = [mlines.line2D((), (), color='black', 
                               label=f'${pv * 100:.0f}%$', **ls)
                 for ls, pv in zip(pstyles, perclevels)]
    lax.legend(handles=lshandles + chandles, fontsize=fontsize - 1., 
               ncol=int(np.floor(1.5 * ncols)),
               loc='upper center', bbox_to_anchor=(0.5, 0.65))
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotvradtot_zweightcomp(filen_template, zfills, wfills):
    '''
    different redshifts in each panel, different plots for
    different ic/phys
    '''

def plotvradtot_icphyscomp(filen_template, icphfills):
    '''
    same weight and redshift in the plot, compare different ic/phys
    between. (new row for new phys)
    '''
