import h5py
import matplotlib as mpl
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np

import fire_an.makeplots.plot_utils as pu
import fire_an.utils.math_utils as mu

def readin_vlos_vrad(filen, rrange_rvir=(0.1, 1.)):
    with h5py.File(filen, 'r') as f:
        hist = 10**f['histogram/histogram'][:]
        rbins_rvir = f['axis_0/bins'][:]
        vlos_bins = f['axis_1/bins'][:] * 1e-5
        vrad_bins = f['axis_2/bins'][:] * 1e-5
        irmin = np.where(np.isclose(rrange_rvir[0], rbins_rvir))
        irmax = np.where(np.isclose(rrange_rvir[-1], rbins_rvir))
        rsel = slice(irmin, irmax, None)
        hist = np.sum(hist[rsel, :, :], axis=0)
    out = {'hist': hist, 'vlos_bins': vlos_bins, 'vrad_bins': vrad_bins}
    return out

def combhist_vlos_vrad(filens, rrange_rvir=(0.1, 1.)):
    out = readin_vlos_vrad(filens[0])
    if len(filens) > 1:
        for filen in filens[1:]:
            o2 = readin_vlos_vrad(filen, rrange_rvir=rrange_rvir)
            e1 = [out['vlos_bins'], out['vrad_bins']]
            e2 = [o2['vlos_bins'], o2['vrad_bins']]
            ht, et = mu.combine_hists(out['hist'], o2['hist'], e1, e2,
                                      rtol=1e-5, atol=1e-8, add=True)
            out['hist'] = ht
            out['vlos_bins'] = et[0]
            out['vrad_bins'] = et[1]
    return out


def plot_vlos_vrad(filen_temp, fills_panel, fills_comb,
                   labels_panel=None, title=None, outname=None,
                   rrange_rvir=(0.1, 1.)):
    fontsize = 12
    npanels = len(fills_panel)
    ncmax = 3
    panelsize = 3.
    nc = min(ncmax, npanels)
    nr = (npanels - 1) // nc + 1
    wspace = 0.
    hspace = 0.
    width = nc * panelsize
    height = nr * panelsize
    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(nrows=nr, ncols=nc, hspace=hspace,
                        wspace=wspace)
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    axes = []
    xlabel = '$v_{\\mathrm{los}} \\; [\\mathrm{km}\\,\\mathrm{s}^{-1}]$'
    ylabel = '$v_{\\mathrm{rad}} \\; [\\mathrm{km}\\,\\mathrm{s}^{-1}]$'
    _cmap = mpl.cm.get_cmap('gist_gray')
    cmap = pu.truncate_colormap(_cmap, minval=0., maxval=0.7)
    clevels = (0.99, 0.9)
    for pi, pfill in enumerate(fills_panel):
        ax = fig.add_subplot(grid[pi // nc, pi % nc])
        axes.append(ax)
        left = pi % nc == 0
        bottom = pi >= npanels - nc
        ax.tick_params(labelsize=fontsize - 1, which='both', direction='in',
                       top=True, right=True, labelbottom=bottom, 
                       labelleft=left)
        if bottom:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        if left:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        if labels_panel is not None:
            _lab = labels_panel[pi]
            ax.text(0.5, 0.95, _lab, fontsize=fontsize, color='black',
                    horizontalalignment='center', verticalalignment='top',
                    transform=ax.transAxes)
        
        fills = [pfill.copy() for _ in range(len(fills_comb))]
        [fill.update(fc) for fill, fc in zip(fills, fills_comb)]
        filens = [filen_temp.format(**fill) for fill in fills]
        he = combhist_vlos_vrad(filens, rrange_rvir=rrange_rvir)
        e0 = he['vlos_bins']
        e1 = he['vrad_bins']
        _hist = he['hist']
        hist = np.log10(he['hist'] / np.sum(he['hist']))
        vmax = np.max(hist)
        vmin = vmax - 5.

        ax.pcolormesh(e0, e1, _hist.T, vmin=vmin, vmax=vmax, cmap=cmap,
                      rasterized=True)
        pu.add_2dhist_contours(ax, _hist, [e0, e1], histlegend=pi==0,
                               fraclevels=True, levels=clevels,
                               color='fuchsia', 
                               linestyles=['dashed', 'solid'])
    xlims = [ax.get_xlim() for ax in axes]
    xmin = np.min([xlim[0] for xlim in xlims])
    xmax = np.min([xlim[1] for xlim in xlims])
    [ax.set_xlim((xmin, xmax)) for ax in axes]
    ylims = [ax.get_ylim() for ax in axes]
    ymin = np.min([ylim[0] for ylim in ylims])
    ymax = np.min([ylim[1] for ylim in ylims])
    [ax.set_ylim((ymin, ymax)) for ax in axes]
    for axi, ax in enumerate(axes):
        plo = max(xmin, ymin)
        phi = min(xmax, ymax)
        ax.plot([plo, phi], [plo, phi], color='black', linestyle='dashed')
        if axi == 0:
            pr = phi - plo
            ax.text(phi - 0.05 * pr, phi - 0.05 * pr, 
                    '$v_{\\mathrm{r}} = v_{\\mathrm{los}}$', 
                    fontsize=fontsize, horizontalalignment='right',
                    verticalalignment='bottom')
        plo = -min(-xmin, ymax)
        phi = max(xmax, -ymin)
        ax.plot([plo, phi], [plo, phi], color='black', linestyle='dashed')
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

