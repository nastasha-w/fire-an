import h5py
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsp
import numpy as np

import fire_an.makeplots.plot_utils as pu
import fire_an.utils.constants_and_units as c
import fire_an.utils.math_utils as mu

def get2dmap_r_vr(filen):
    with h5py.File(filen, 'r') as f:
        hist = f['histogram/histogram'][:]
        islog = bool(f['histogram/histogram'].attrs['log'])
        if islog:
            hist = 10**hist
        hist = np.sum(hist, axis=2)
        rbins_rvir = f['axis_0/bins'][:]
        vradbins_cmps = f['axis_1/bins'][:]
        cosmopars = {key: val for key, val 
                     in f['Header/cosmopars'].attrs.items()}
        mvir_g = f['Header/inputpars/halodata'].attrs['Mvir_g']
        rvir_cm = f['Header/inputpars/halodata'].attrs['Rvir_cm']
    vradbins_kmps = vradbins_cmps * 1e-5
    normedhist = hist / np.sum(hist) / np.diff(rbins_rvir)[:, np.newaxis] \
                 / np.diff(vradbins_kmps)[np.newaxis, :]
    lognormedhist = np.log10(normedhist)
    out = {'rbins_rvir': rbins_rvir,
           'vradbins_kmps': vradbins_kmps,
           'cosmopars': cosmopars,
           'mvir_g': mvir_g,
           'rvir_cm': rvir_cm,
           'lognormedhist': lognormedhist,
           'hist': hist}
    return out

def plot_r_vr_weights(filen_temp, weightfills, weightlabels=None,
                      title=None, outname=None):
    
    cmap = pu.truncate_colormap('gist_yarg', minval=0., maxval=0.7)
    data = [get2dmap_r_vr(filen_temp.format(**fill)) for fill in weightfills]
    vmax = max([np.max(datum['lognormedhist']) for datum in data])
    vmin = min([np.min(datum['lognormedhist']) for datum in data])
    vmin = max(vmin, vmax - 7.)
    extend = 'min'
    clabel = ('$\\log_{10} \\; \\partial^2 \\mathrm{weight} \\,/\\,'
              '\\Sigma \\, \\mathrm{weight}$'
              '\\,/\\, \\partial (r \\; [\\mathrm{R}_{\\mathrm{vir}}])'
              '\\,/\\, \\partial (v_{\\mathrm{r}}'
              '\\; [\\mathrm{km}\\, \\mathrm{s}^{-1}])')
    xlabel = '$r \\; [\\mathrm{R}_{\\mathrm{vir}}$'
    ylabel = '$v_{\\mathrm{r}}\\; [\\mathrm{km}\\, \\mathrm{s}^{-1}]$'
    fontsize = 12

    nw = len(weightfills)
    ncmax = 4
    panelsize = 2.5
    ncols = min(ncmax, nw)
    nrows = (nw - 1) // ncols + 1
    caxspace = 1.
    width_ratios = [panelsize] * ncols + [caxspace]
    height_ratios = [panelsize] * nrows
    cbar_orientation = 'vertical'
    _ncols = ncols + 1
    width = sum(width_ratios)
    height = sum(height_ratios)
    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=_ncols, nrows=nrows,
                        height_ratios=height_ratios,
                        width_ratios=width_ratios)
    axes = [fig.add_subplot(grid[i % ncols, i // ncols]) for i in range(nw)]
    cax = fig.add_subplot(grid[:, -1])

    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    for i, (ax, datum) in enumerate(zip(axes, data)):
        below = ncols <= nw - i
        left = i % ncols == 0
        ax.tick_params(which='both', direction='in', top=True,
                       right=True, labelbottom=below,
                       labelright=left)
        if left:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        if below:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        img = ax.pcolormesh(datum['rbins_rvir'], datum['vradbins_kmps'],
                            datum['lognormedhist'].T, cmap=cmap,
                            vmin=vmin, vmax=vmax, rasterized=True)
        ax.axvline(1., linestyle='dotted', color='black')
        vescvir = np.sqrt(c.gravity * datum['Mvir_g'] / datum['Rvir_cm']) \
                  * 1e-5
        ax.axhline(vescvir, linestyle='solid', color='black')
        ax.axhline(-vescvir, linestyle='solid', color='black')

        ylo, ymed, yhi = mu.percentiles_from_histogram(
            datum['hist'], datum['vradbins_kmps'], axis=-1, 
            percentiles=np.array([0.1, 0.5, 0.9]))
    
    plt.colorbar(img, cax=cax, extend=extend, orientation=cbar_orientation)
    cax.set_ylabel(clabel, fontsize=fontsize)





    



