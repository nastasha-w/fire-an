from re import L
import h5py
import matplotlib as mpl
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import fire_an.makeplots.plot_utils as pu
import fire_an.simlists as sl
import fire_an.utils.math_utils as mu

def get_model(simname):
    if '_sdp1e10_' in simname:
        return 'noBH'
    elif '_MHDCRspec1_' in simname:
        return 'AGN-CR'
    else:
        return 'AGN-noCR'

def readin_vlos_vrad(filen, rrange_rvir=(0.1, 1.)):
    with h5py.File(filen, 'r') as f:
        hist = 10**f['histogram/histogram'][:]
        rbins_rvir = f['axis_0/bins'][:]
        vlos_bins = f['axis_1/bins'][:] * 1e-5
        vrad_bins = f['axis_2/bins'][:] * 1e-5
        print(rbins_rvir)
        print(rrange_rvir)
        irmin = np.where(np.isclose(rrange_rvir[0], rbins_rvir))[0][0]
        irmax = np.where(np.isclose(rrange_rvir[-1], rbins_rvir))[0][0]
        rsel = slice(irmin, irmax, None)
        hist = np.sum(hist[rsel, :, :], axis=0)
    out = {'hist': hist, 'vlos_bins': vlos_bins, 'vrad_bins': vrad_bins}
    return out

def combhist_vlos_vrad(filens, rrange_rvir=(0.1, 1.)):
    out = readin_vlos_vrad(filens[0], rrange_rvir=rrange_rvir)
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
                   rrange_rvir=(0.1, 1.), maxv = 400.):
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
    _cmap = mpl.cm.get_cmap('gist_yarg')
    cmap = pu.truncate_colormap(_cmap, minval=0., maxval=0.7)
    clevels = (0.96, 0.8)
    ccolors = ['fuchsia'] * 2
    clinestyles = ['dashed', 'solid']
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
        vmin = vmax - 4.

        ax.pcolormesh(e0, e1, hist.T, vmin=vmin, vmax=vmax, cmap=cmap,
                      rasterized=True)
        pu.add_2dhist_contours(ax, _hist, [e0, e1], (0, 1), histlegend=False,
                               fraclevels=True, levels=clevels,
                               colors=ccolors, 
                               linestyles=clinestyles)
    xlims = [ax.get_xlim() for ax in axes]
    xmin = np.min([xlim[0] for xlim in xlims])
    xmin = max(xmin, -maxv)
    xmax = np.max([xlim[1] for xlim in xlims])
    xmax = min(xmax, maxv)
    [ax.set_xlim((xmin, xmax)) for ax in axes]
    ylims = [ax.get_ylim() for ax in axes]
    ymin = np.min([ylim[0] for ylim in ylims])
    ymin = max(ymin, -maxv)
    ymax = np.max([ylim[1] for ylim in ylims])
    ymax = min(ymax, maxv)
    [ax.set_ylim((ymin, ymax)) for ax in axes]
    for axi, ax in enumerate(axes):
        plo = max(xmin, ymin)
        phi = min(xmax, ymax)
        ax.plot([plo, phi], [plo, phi], color='black', linestyle='dotted',
                linewidth=1.)
        if axi == 0:
        #    pr = phi - plo
        #    ax.text(plo + 0.15 * pr, plo + 0.15 * pr, 
        #            '$v_{\\mathrm{r}} = v_{\\mathrm{los}}$', 
        #            fontsize=fontsize, horizontalalignment='left',
        #            verticalalignment='top')
            handles = [mlines.Line2D((), (), label=f'{level*100.:.0f}%',
                                     color=color, linestyle=linestyle)
                       for level, color, linestyle 
                       in zip(clevels, ccolors, clinestyles)]
            ax.legend(handles=handles, fontsize=fontsize - 2,
                      loc='lower center', bbox_to_anchor=(0.5, 0.0),
                      handlelength=1.5, ncol=2, columnspacing=1.,
                      handletextpad=0.5)
        plo = -min(-xmin, ymax)
        phi = max(xmax, -ymin)
        ax.plot([plo, phi], [-plo, -phi], color='black', linestyle='dotted',
                linewidth=1.)
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')


def plotset_vlos_vrad():
    fills_comb = [{'pax': 'x'}, {'pax': 'y'}, {'pax': 'z'}]
    ddir = '/projects/b1026/nastasha/hists/vdop_vrad_all2/'
    filen_temp = ('hist_vlos_{pax}ax_vrad_by_Ne8_{simname}'
                  '_snap{snapnum}_bins1_v1_hvcen.hdf5')
    odir = '/projects/b1026/nastasha/imgs/vlos_vr_hists/'

    simnames = sl.m12_sr_clean2 + sl.m13_sr_clean2 +\
               sl.m12_hr_clean2 + sl.m13_hr_clean2
    sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    snaps_hr = (min(sl.snaps_hr), max(sl.snaps_hr))
    snaps_sr = (min(sl.snaps_sr), max(sl.snaps_sr))
    snapindic = [1.0, 0.5]
    for sn in sl.buglist1:
        if sn in simnames:
            simnames.remove(sn)
    ics = [simname.split('_')[0] for simname in simnames]
    ics_clean = {ic for ic in ics if sum([ic == _ic for _ic in ics]) == 3}
    simsets = {ic: sorted([sn for sn in simnames if sn.split('_')[0] == ic])
               for ic in ics_clean}
    for ic in simsets:
        _sms = simsets[ic]
        _sns = [snaps_hr if _sm in sims_hr else
                snaps_sr if _sm in sims_sr else
                None 
                for _sm in _sms]
        _sms = [_sm for _sm in _sms for _ in range(len(snaps_hr))]
        _sns = [_sn for sub in _sns for _sn in sub]
        fills_panel = [{'simname': _sm, 'snapnum': _sn}
                       for _sm, _sn in zip(_sms, _sns)]
        maxv = 400. if ic.startswith('m12') else 600.
        for rrange_rvir in [(0.15, 1.), (0.15, 0.25), (0.45, 0.55), (0.9, 1.0)]:
            title = (f'{ic}, ${rrange_rvir[0]:.2f} \\endash'
                     f'{rrange_rvir[1]:.2f} \\,'
                     '\\mathrm{R}_{\\mathrm{vir}}$')
            labels_panel = [f'{get_model(_sm)}, z={snapindic[i%2]:.1f}' 
                            for i, _sm in enumerate(_sms)]
            rstr = f'{rrange_rvir[0]:.2f}_to_{rrange_rvir[1]:.2f}_Rvir'
            rstr = rstr.replace('.', 'p')
            outname = odir + f'hist_vdop_vrad_{ic}_z0p5_1p0_{rstr}.pdf'
            plot_vlos_vrad(ddir + filen_temp, fills_panel, fills_comb,
                           labels_panel=labels_panel, title=title, 
                           outname=outname, rrange_rvir=rrange_rvir,
                           maxv=maxv)
            break
        break
