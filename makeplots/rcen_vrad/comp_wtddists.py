
import h5py
from matplotlib import patheffects
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.math_utils as mu

def getphys(simname):
    phys = ('noBH' if '_sdp1e10_' in simname 
            else 'AGN-CR' if '_MHDCRspec1_' in simname 
            else 'AGN-noCR')
    return phys

def get_kindist(filen, rrange_rvir=(0.1, 1.0),
                vrrange=None,
                vrrange_units='kmps',
                dist_target='vr', vrbin_fac=4):
    '''
    when plotting vdists, those 5 km/s bins make
    for some seriously noisy lines.
    '''
    with h5py.File(filen, 'r') as f:
        hist = f['histogram/histogram'][:]
        islog = bool(f['histogram'].attrs['log'])
        if islog:
            hist = 10**hist
        mvir_g = f['Header/inputpars/halodata'].attrs['Mvir_g']
        rvir_cm = f['Header/inputpars/halodata'].attrs['Rvir_cm']
        vescvir_kmps = np.sqrt(c.gravity * mvir_g / rvir_cm) \
                       * 1e-5
        hsel = [slice(None, None, None)] * 3
        if dist_target == 'vr':
            sumaxes = (0, 2)
            binaxis = 1
            binconv = 1e-5
        elif dist_target == 'r':
            sumaxes = (1, 2)
            binaxis = 0
            binconv = 1.
        elif dist_target == 'T':
            if '_temperature_' not in filen:
                raise ValueError(f'{filen} does not contain T info')
            sumaxes = (0, 1)
            binaxis = 2
            binconv = 1.
        elif dist_target == 'nH':
            if '_hdens_' not in filen:
                raise ValueError(f'{filen} does not contain nH info')
            sumaxes = (0, 1)
            binaxis = 2
            binconv = 1.
        else:
            raise ValueError(f'Invalid dist_target {dist_target}')
        if rrange_rvir is not None:
            rbins_rvir = f['axis_0/bins'][:]
            rimin = np.where(np.isclose(rrange_rvir[0], rbins_rvir))[0][0]
            rimax = np.where(np.isclose(rrange_rvir[-1], rbins_rvir))[0][0]
            if rrange_rvir[0] == -np.inf:
                rimin = 0
            if rrange_rvir[-1] == np.inf:
                rimax = len(rbins_rvir)
            hsel[0] = slice(rimin, rimax, None)
        if vrrange is not None:
            if vrrange_units == 'kmps':
                vrrange = (vrrange[0] * 1e5, vrrange[-1] * 1e5)
            elif vrrange_units == 'vesc':
                vrrange = [vrrange[0] * vescvir_kmps * 1e5,
                           vrrange[-1] * vescvir_kmps * 1e5]
            vbins_cmps = f['axis_1/bins'][:]
            vimin = np.argmin(np.abs(vbins_cmps - vrrange[0]))
            vimax = np.argmin(np.abs(vbins_cmps - vrrange[-1]))
            if vrrange[0] == -np.inf:
                vimin = 0
            if vrrange[-1] == np.inf:
                vimax = len(vbins_cmps)
            hsel[1] = slice(vimin, vimax, None) 
        cosmopars = {key: val for key, val 
                     in f['Header/cosmopars'].attrs.items()}
        _hist = np.sum(hist[tuple(hsel)], axis=sumaxes)
        esel = slice(hsel[binaxis].start, 
                     hsel[binaxis].stop if hsel[binaxis].stop is None
                     else hsel[binaxis].stop + 1,
                     hsel[binaxis].step)
        _edges = f[f'axis_{binaxis}/bins'][esel]
        if dist_target == 'vr' and vrbin_fac != 1:
            rest = len(_hist) % vrbin_fac
            __edges = _edges[slice(None, None, vrbin_fac)]
            if rest != 0:
                _hist = np.append(_hist, np.zeros(vrbin_fac - rest))
                __edges = np.append(__edges, _edges[-1])
            _edges = __edges
            newshape = (len(_hist) // vrbin_fac, vrbin_fac)
            _hist = np.sum(_hist.reshape(newshape), axis=1)
        if np.any(np.logical_not(np.isfinite(_hist))):
            print(f'filen: inf/nan values in histogram')
            print(_hist)
        if np.any(np.logical_not(np.isfinite(_edges))):
            print(f'filen: inf/nan values in edges')
            print(_edges)
        _edges = _edges * binconv
    out = {'cosmopars': cosmopars,
           'vescvir_kmps': vescvir_kmps,
           'mvir_g': mvir_g,
           'rvir_cm': rvir_cm,
           'edges': _edges,
           'hist': _hist}
    return out

def fmtdcts(str, *args):
    dct = {}
    for arg in args:
        dct.update(arg)
    _str = str.format(**dct)
    return _str

def get_kindists(filen_temp, weights, simnames,
                 ssel=[0, 2, 5],
                 **kwargs):
    snaps_sr = [sl.snaps_sr[i] for i in ssel]
    snaps_hr = [sl.snaps_hr[i] for i in ssel]
    sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    sfills = [{'simname': simname} for simname in simnames]
    wfills = [{'weight': weight} for weight in weights]
    zfills = [[{'snapnum': snap} for snap in snaps_sr]
              if simname in sims_sr else 
              [{'snapnum': snap} for snap in snaps_hr]
              if simname in sims_hr else
              None
              for simname in simnames]
    out = [[[get_kindist(fmtdcts(filen_temp, sfill, wfill, zfill), **kwargs)
             for zfill in zlist]
            for sfill, zlist in zip(sfills, zfills)]
           for wfill in wfills]
    return out

def plot_1Ddists_zphysweightcomp(filen_temp, weights,
                                 simnames, 
                                 weightlabels=None, 
                                 rrange_rvir=(0.1, 1.0),
                                 vrrange=None,
                                 vrrange_units='kmps',
                                 outname=None,
                                 dist_target='vr',
                                 title=None):
    '''
    Each panel: different z -> different color,
        different phys. model -> different linestyle
    different panels: different weights
    just do the clean sample, 4 haloes
    '''
    if weightlabels is None:
        weightlabels = weights
    zvals = np.linspace(0.5, 1., 6)
    get_zcolors = pu.getzcolorfunc(zvals, ztol=1e-2)
    fontsize = 12
    physstyles = {'noBH': 'dashed',
                  'AGN-noCR': 'solid',
                  'AGN-CR': 'dotted'}
    linewidth = 1.5
    patheff = pu.getoutline(linewidth)

    if dist_target == 'vr':
        xlabel = '$v_{\\mathrm{rad}} \\; [\\mathrm{km}\\,\\mathrm{s}^{-1}]$'
        yl_add = '\\partial \\, v_{\\mathrm{rad}}'
    elif dist_target == 'r':
        xlabel =  '$r \\; [\\mathrm{R}_{\\mathrm{vir}}]$'
        yl_add = '\\partial \\, r'
    elif dist_target == 'T':
        xlabel = '$\\log_{10} \\, \\mathrm{T} \\; [\\mathrm{K}]$'
        yl_add = '\\partial \\, \\log_{10} \\,\\mathrm{T}'
    elif dist_target == 'nH':
        xlabel = ('$\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
                  ' \\; [\\mathrm{cm}^{-3}]$')
        yl_add = '\\partial \\, \\log_{10} \\,\\mathrm{n}_{\\mathrm{H}}'
    ylabel = ('$\\partial \\, \\mathrm{weight} \\,/\\, \\mathrm{tot.}'
              f' \\,/\\, {yl_add}$')

    data = get_kindists(filen_temp, weights, simnames,
                        rrange_rvir=rrange_rvir,
                        vrrange=vrrange,
                        vrrange_units=vrrange_units,
                        dist_target=dist_target)
    
    panelsize = 2.5
    ncmax = 4
    nw = len(weights)
    ncols = min(ncmax, nw)
    nrows = (nw - 1) // ncols + 1
    hspace = 0.
    wspace = 0.
    lspace = 2.
    height_ratios = [panelsize] * nrows + [lspace]
    width_ratios = [panelsize] * ncols
    width = sum(width_ratios)
    height = sum(height_ratios)
    ncol_legend = int(np.floor(1.7 * ncols))

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows + 1,
                        height_ratios=height_ratios,
                        width_ratios=width_ratios,
                        hspace=hspace, wspace=wspace)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols]) for i in range(nw)]
    lax = fig.add_subplot(grid[-1, :])
    lax.axis('off')

    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    
    ymaxmax = -np.inf
    xminmin = np.inf
    xmaxmax = -np.inf
    xmar = 0.05
    yrange = 2.5
    ymar_top = 0.07
    
    for wi, (wdata, wlabel) in enumerate(zip(data, weightlabels)):
        ax = axes[wi]
        bottom = wi >= nw - ncols
        left = wi % ncols == 0
        ax.tick_params(which='both', labelsize=fontsize - 1.,
                       direction='in', right=True, top=True,
                       labelleft=left, labelbottom=bottom)
        ax.grid(visible=True, which='major')
        ax.set_yscale('log')
        if left:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        if bottom:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.text(0.98, 0.98, wlabel, fontsize=fontsize,
                horizontalalignment='right', verticalalignment='top',
                transform=ax.transAxes, color='black')
        for (sdata, simname) in zip(wdata, simnames):
            phys = getphys(simname)
            linestyle = physstyles[phys]
            for datum in sdata:
                _hist = datum['hist']
                _edges = datum['edges']
                if np.sum(_hist) == 0.:
                    print('zero weight for ', wlabel, 
                          ' Rvir ranges: ', rrange_rvir,
                          ' simname: ', simname)
                    continue
                phist = _hist / np.sum(_hist) / np.diff(_edges)
                #phist = np.append(phist, [0.])
                zv = datum['cosmopars']['z']
                color = get_zcolors(zv)
                #ax.step(_edges, phist, where='post', color=color,
                #        linestyle=linestyle, linewidth=linewidth,
                #        path_effects=patheff)
                bcens = 0.5 * (_edges[:-1] + _edges[1:])
                ax.plot(bcens, phist, color=color,
                        linestyle=linestyle, linewidth=linewidth,
                        path_effects=patheff)
                _ymax = np.max(phist)
                ymaxmax = max(ymaxmax, _ymax)
                goodinds = np.where(phist >= _ymax * 10**(-yrange))[0]
                _xmin = _edges[goodinds[0]]
                _xmax = _edges[goodinds[-1] + 1]
                if np.isfinite(_xmin):
                    xminmin = min(xminmin, _xmin)
                if np.isfinite(_xmax):
                    xmaxmax = max(xmaxmax, _xmax)
    
    ymax = ymaxmax
    ymin = ymax * 10**(-yrange)
    print(ymax, ymin)
    ymax = ymax * 10**((np.log10(ymax) - np.log10(ymin)) * ymar_top)
    [ax.set_ylim((ymin, ymax)) for ax in axes]
    xmar = xmar * (xmaxmax - xminmin)
    xmin = xminmin - xmar
    xmax = xmaxmax + xmar
    [ax.set_xlim((xmin, xmax)) for ax in axes]

    handles1 = [mlines.Line2D((), (), color='gray', linestyle=physstyles[phys],
                              linewidth=linewidth, path_effects=patheff,
                              label=phys)
                for phys in physstyles]
    handles2 = [mlines.Line2D((), (), color=get_zcolors(zv), linestyle='solid',
                              linewidth=linewidth, path_effects=patheff,
                              label=f'z={zv:.1f}')
                for zv in zvals]
    lax.legend(handles=handles1 + handles2, fontsize=fontsize,
               ncol=ncol_legend, loc='upper center',
               bbox_to_anchor=(0.5, 0.65))
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')
            
def plot_vr_dists(plottype='vr_fixedr'):
    ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
    filetemp = ('hist_rcen_vcen_temperature_by_{weight}_{simname}'
                '_snap{snapnum}_bins1_v1_hvcen.hdf5')
    outdir = '/projects/b1026/nastasha/imgs/r_vr_hists/'
    weights = ['gasmass', 'gasvol', 'Oxygen', 'Neon',
               'O6', 'Ne8', 'O7', 'Ne9', 'O8', 'Ne10']
    weightlabels = ['Mass', 'Volume', 'Oxygen', 'Neon',
                    'O VI', 'Ne VIII', 'O VII', 'Ne IX', 'O VIII', 'Ne X']

    simnames = sl.m13_hr_clean2 + sl.m13_sr_clean2 +\
                [('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                    '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                    '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp2e-4_gacc31_fa0.5'),
                ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp1e10_gacc31_fa0.5'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp2e-4_gacc31_fa0.5'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp1e10_gacc31_fa0.5')]
    ics = []
    simsets = []
    for simname in simnames:
        ic = simname.split('_')[0]
        if ic in ics:
            ici = np.where([ic == _ic for _ic in ics])[0][0]
            simsets[ici].append(simname)
        else:
            ics.append(ic)
            simsets.append([simname])
    for ic, _simnames in zip(ics, simsets):
        for rrange_rvir in [(0.1, 1.), (0.15, 0.25), 
                            (0.45, 0.55), (0.9, 1.0)]:
            filen_temp = ddir + filetemp
            if plottype == 'vr_fixedr':
                title = (f'{ic}, ${rrange_rvir[0]:.2f} \\endash '
                        f'{rrange_rvir[-1]:.2f} \\,'
                        ' \\mathrm{R}_{\\mathrm{vir}}$')
                _outname = (f'dists_vr_{ic}_phys_z_wt_comp'
                            f'_{rrange_rvir[0]:.2f}'
                            f'_to_{rrange_rvir[0]:.2f}_Rvir')
                _outname = _outname.replace('.', 'p')
                outname = outdir +  _outname + '.pdf'
                plot_1Ddists_zphysweightcomp(filen_temp, weights,
                                            _simnames, 
                                            weightlabels=weightlabels, 
                                            rrange_rvir=rrange_rvir,
                                            vrrange=None,
                                            vrrange_units='kmps',
                                            outname=outname,
                                            dist_target='vr',
                                            title=title)
            elif plottype == 'rdist':
                title = (f'{ic}, ${rrange_rvir[0]:.2f} \\endash '
                        f'{rrange_rvir[-1]:.2f} \\,'
                        ' \\mathrm{R}_{\\mathrm{vir}}$')
                _outname = (f'dists_r_{ic}_phys_z_wt_comp'
                            f'_{rrange_rvir[0]:.2f}'
                            f'_to_{rrange_rvir[0]:.2f}_Rvir')
                _outname = _outname.replace('.', 'p')
                outname = outdir +  _outname + '.pdf'
                plot_1Ddists_zphysweightcomp(filen_temp, weights,
                                            _simnames, 
                                            weightlabels=weightlabels, 
                                            rrange_rvir=rrange_rvir,
                                            vrrange=None,
                                            vrrange_units='kmps',
                                            outname=outname,
                                            dist_target='r',
                                            title=title)
            elif plottype == 'Tdist_rranges':
                title = (f'{ic}, ${rrange_rvir[0]:.2f} \\endash '
                        f'{rrange_rvir[-1]:.2f} \\,'
                        ' \\mathrm{R}_{\\mathrm{vir}}$')
                _outname = (f'dists_T_{ic}_phys_z_wt_comp'
                            f'_{rrange_rvir[0]:.2f}'
                            f'_to_{rrange_rvir[0]:.2f}_Rvir')
                _outname = _outname.replace('.', 'p')
                outname = outdir +  _outname + '.pdf'
                plot_1Ddists_zphysweightcomp(filen_temp, weights,
                                            _simnames, 
                                            weightlabels=weightlabels, 
                                            rrange_rvir=rrange_rvir,
                                            vrrange=None,
                                            vrrange_units='kmps',
                                            outname=outname,
                                            dist_target='T',
                                            title=title)
            elif plottype == 'nHdist_rranges':
                filetemp = ('hist_rcen_vcen_hdens_by_{weight}_{simname}'
                            '_snap{snapnum}_bins1_v1_hvcen.hdf5')
                filen_temp = ddir + filetemp
                title = (f'{ic}, ${rrange_rvir[0]:.2f} \\endash '
                        f'{rrange_rvir[-1]:.2f} \\,'
                        ' \\mathrm{R}_{\\mathrm{vir}}$')
                _outname = (f'dists_nH_{ic}_phys_z_wt_comp'
                            f'_{rrange_rvir[0]:.2f}'
                            f'_to_{rrange_rvir[0]:.2f}_Rvir')
                _outname = _outname.replace('.', 'p')
                outname = outdir +  _outname + '.pdf'
                plot_1Ddists_zphysweightcomp(filen_temp, weights,
                                            _simnames, 
                                            weightlabels=weightlabels, 
                                            rrange_rvir=rrange_rvir,
                                            vrrange=None,
                                            vrrange_units='kmps',
                                            outname=outname,
                                            dist_target='nH',
                                            title=title)
                
            elif plottype == 'Tdist_vr_r_ranges':
                for vrrange in [(-np.inf, -0.5), (-0.5, 0.5), (0.5, np.inf)]:
                    vrrange_units = 'vesc'
                    title = (f'{ic}, ${rrange_rvir[0]:.2f} \\endash '
                            f'{rrange_rvir[-1]:.2f} \\,'
                            ' \\mathrm{R}_{\\mathrm{vir}}, '
                            f' {vrrange[0]:.1f} \\endash {vrrange[-1]:.1f}'
                            '\\, v_{\\mathrm{esc}}$')
                    _outname = (f'dists_T_{ic}_phys_z_wt_comp'
                                f'_{rrange_rvir[0]:.2f}'
                                f'_to_{rrange_rvir[-1]:.2f}_Rvir'
                                f'_vr_{vrrange[0]:.1f}_to_{vrrange[-1]:.1f}'
                                f'_{vrrange_units}')
                    _outname = _outname.replace('.', 'p')
                    outname = outdir +  _outname + '.pdf'
                    plot_1Ddists_zphysweightcomp(filen_temp, weights,
                                                _simnames, 
                                                weightlabels=weightlabels, 
                                                rrange_rvir=rrange_rvir,
                                                vrrange=vrrange,
                                                vrrange_units=vrrange_units,
                                                outname=outname,
                                                dist_target='T',
                                                title=title)
            elif plottype == 'nHdist_vr_r_ranges':
                for vrrange in [(-np.inf, -0.5), (-0.5, 0.5), (0.5, np.inf)]:
                    filetemp = ('hist_rcen_vcen_hdens_by_{weight}_{simname}'
                                '_snap{snapnum}_bins1_v1_hvcen.hdf5')
                    filen_temp = ddir + filetemp
                    vrrange_units = 'vesc'
                    title = (f'{ic}, ${rrange_rvir[0]:.2f} \\endash '
                            f'{rrange_rvir[-1]:.2f} \\,'
                            ' \\mathrm{R}_{\\mathrm{vir}}, '
                            f' {vrrange[0]:.1f} \\endash {vrrange[-1]:.1f}'
                            '\\, v_{\\mathrm{esc}}$')
                    _outname = (f'dists_nH_{ic}_phys_z_wt_comp'
                                f'_{rrange_rvir[0]:.2f}'
                                f'_to_{rrange_rvir[-1]:.2f}_Rvir'
                                f'_vr_{vrrange[0]:.1f}_to_{vrrange[-1]:.1f}'
                                f'_{vrrange_units}')
                    _outname = _outname.replace('.', 'p')
                    outname = outdir +  _outname + '.pdf'
                    plot_1Ddists_zphysweightcomp(filen_temp, weights,
                                                _simnames, 
                                                weightlabels=weightlabels, 
                                                rrange_rvir=rrange_rvir,
                                                vrrange=vrrange,
                                                vrrange_units=vrrange_units,
                                                outname=outname,
                                                dist_target='nH',
                                                title=title)

def comp_vdists_physmodels(filen_temp, simnames, weight='Volume',
                           title=None, outname=None,
                           rrange_rvir=(0.1, 1.), label_weightpart=None):
    
    ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
    _filen_temp = ('hist_rcen_vcen_temperature_by_{weight}_{simname}'
                   '_snap{snapnum}_bins1_v1_hvcen.hdf5')
    filen_temp = ddir + _filen_temp
    if label_weightpart is None:
        label_weightpart = '\\mathrm{weight}'
    xlabel = '$v_{\\mathrm{rad}} \\; [\\mathrm{km} \\, \\mathrm{s}^{-1}]$'
    ylabel = (f'$\\partial {label_weightpart} \\,/\\, '
              f'{label_weightpart}'
              '\\,/\\, \\partial v_{\\mathrm{rad}}'
              '\\;[\\mathrm{km}^{-1}\\mathrm{s}]$')

    data = get_kindists(filen_temp, [weight], simnames,
                        rrange_rvir=rrange_rvir,
                        vrrange=None,
                        vrrange_units='kmps',
                        dist_target='vr', vrbin_fac=4,
                        ssel=range(6))
    # out = {'cosmopars': cosmopars,
    #       'vescvir_kmps': vescvir_kmps,
    #       'mvir_g': mvir_g,
    #       'rvir_cm': rvir_cm,
    #       'edges': _edges,
    #       'hist': _hist}
    # out = [[[get_kindist(fmtdcts(filen_temp, sfill, wfill, zfill), **kwargs)
    #         for zfill in zlist]
    #        for sfill, zlist in zip(sfills, zfills)]
    #       for wfill in wfills]
    data = data[0]
    zs = list({np.round(datum['cosmopars']['z'], 2) 
               for l1 in data for datum in l1})
    zs.sort()
    get_zcolors = pu.getzcolorfunc(zs)
    
    panelsize = 2.5
    ncmax = 3
    ns = len(simnames)
    ncols = min(ncmax, ns)
    nrows = (ns - 1) // ncols + 1
    hspace = 0.
    wspace = 0.
    lspace = 2.
    height_ratios = [panelsize] * nrows + [lspace]
    width_ratios = [panelsize] * ncols
    width = sum(width_ratios)
    height = sum(height_ratios)
    ncol_legend = int(np.floor(1.7 * ncols))

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows + 1,
                        height_ratios=height_ratios,
                        width_ratios=width_ratios,
                        hspace=hspace, wspace=wspace)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols]) for i in range(ns)]
    lax = fig.add_subplot(grid[-1, :])
    lax.axis('off')

    fontsize = 12
    linewidth = 1.5
    patheff = pu.getoutline(linewidth)

    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    
    ymaxmax = -np.inf
    xminmin = np.inf
    xmaxmax = -np.inf
    xmar = 0.05
    yrange = 2.5
    ymar_top = 0.07
    
    for si, simname in enumerate(simnames):
        ax = axes[si]
        sdata = data[si]
        bottom = si >= ns - ncols
        left = si % ncols == 0
        ax.tick_params(which='both', labelsize=fontsize - 1.,
                       direction='in', right=True, top=True,
                       labelleft=left, labelbottom=bottom)
        ax.grid(visible=True, which='both')
        ax.set_yscale('log')
        if left:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        if bottom:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        ic = simname.split('_')[0]
        phys = getphys(simname)
        axlabel = f'{ic}, {phys}'
        ax.text(0.98, 0.98, axlabel, fontsize=fontsize,
                horizontalalignment='right', verticalalignment='top',
                transform=ax.transAxes, color='black')
        linestyle = 'solid'
        for datum in sdata:
            _hist = datum['hist']
            _edges = datum['edges']
            if np.sum(_hist) == 0.:
                print('zero weight for ', axlabel, 
                      ' Rvir ranges: ', rrange_rvir,
                      ' simname: ', simname)
                continue
            phist = _hist / np.sum(_hist) / np.diff(_edges)
            #phist = np.append(phist, [0.])
            zv = datum['cosmopars']['z']
            color = get_zcolors(zv)
            #ax.step(_edges, phist, where='post', color=color,
            #        linestyle=linestyle, linewidth=linewidth,
            #        path_effects=patheff)
            bcens = 0.5 * (_edges[:-1] + _edges[1:])
            ax.plot(bcens, phist, color=color,
                    linestyle=linestyle, linewidth=linewidth,
                    path_effects=patheff)
            _ymax = np.max(phist)
            ymaxmax = max(ymaxmax, _ymax)
            goodinds = np.where(phist >= _ymax * 10**(-yrange))[0]
            _xmin = _edges[goodinds[0]]
            _xmax = _edges[goodinds[-1] + 1]
            if np.isfinite(_xmin):
                xminmin = min(xminmin, _xmin)
            if np.isfinite(_xmax):
                xmaxmax = max(xmaxmax, _xmax)
    
    ymax = ymaxmax
    ymin = ymax * 10**(-yrange)
    print(ymax, ymin)
    ymax = ymax * 10**((np.log10(ymax) - np.log10(ymin)) * ymar_top)
    [ax.set_ylim((ymin, ymax)) for ax in axes]
    xmar = xmar * (xmaxmax - xminmin)
    xmin = xminmin - xmar
    xmax = xmaxmax + xmar
    [ax.set_xlim((xmin, xmax)) for ax in axes]

    #handles1 = [mlines.Line2D((), (), color='gray', linestyle=physstyles[phys],
    #                          linewidth=linewidth, path_effects=patheff,
    #                          label=phys)
    #            for phys in physstyles]
    handles2 = [mlines.Line2D((), (), color=get_zcolors(zv), linestyle='solid',
                              linewidth=linewidth, path_effects=patheff,
                              label=f'z={zv:.1f}')
                for zv in zs]
    lax.legend(handles=handles2, fontsize=fontsize,
               ncol=ncol_legend, loc='upper center',
               bbox_to_anchor=(0.5, 0.65))
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def compsets_vdists_physmodels():
    ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
    filetemp = ('hist_rcen_vcen_temperature_by_{weight}_{simname}'
                '_snap{snapnum}_bins1_v1_hvcen.hdf5')
    outdir = '/projects/b1026/nastasha/imgs/r_vr_hists/'

    simnames_m13 = sl.m13_hr_clean2 + sl.m13_sr_clean2
    simnames_m12 = [('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                    '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp2e-4_gacc31_fa0.5'),
                ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp1e10_gacc31_fa0.5'),
                ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                    '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp2e-4_gacc31_fa0.5'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp1e10_gacc31_fa0.5')]
    simnames_m13.sort()
    simnames_m12.sort()
    simsets = [simnames_m12, simnames_m13]
    simsetnames = ['m12', 'm13']
    for simset, simsetn in zip(simsets, simsetnames):
        for rrange in [(0.1, 1.), (0.15, 0.25), (0.45, 0.55), (0.9, 1.0)]:
            filen_temp = ddir + filetemp
            title = (f'{simsetn}, '
                     f'${rrange[0]:.2f} \\endash {rrange[-1]:.2f}'
                      ' \\, \\mathrm{R}_{\\mathrm{vir}}$')
            _outname = (f'vdists_physcomp_{simsetn}'
                        f'_{rrange[0]:.2f}'
                        f'_{rrange[-1]:.2f}_Rvir_v1')
            _outname = _outname.replace('.', 'p')
            outname = outdir + _outname + '.pdf'

            comp_vdists_physmodels(filen_temp, simset, weight='gasvol',
                                   title=title, outname=outname,
                                   rrange_rvir=rrange, 
                                   label_weightpart='\\mathrm{V}')
            
def comp_vperc_physmodels(filen_temp, simnames, weight='gasvol',
                          title=None, outname=None,
                          rrange_rvir=(0.1, 1.), label_weightpart=None,
                          simname_styles=None):
    
    ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
    _filen_temp = ('hist_rcen_vcen_temperature_by_{weight}_{simname}'
                   '_snap{snapnum}_bins1_v1_hvcen.hdf5')
    filen_temp = ddir + _filen_temp
    if label_weightpart is None:
        label_weightpart = '\\mathrm{weight}'
    ylabel = (f'$v_{{\\mathrm{{rad}}}} ({label_weightpart})'
              '\\; [\\mathrm{km} \\, \\mathrm{s}^{-1}]$')
    xlabel = '$z$'
    stats = ('mean', 'p10', 'p90')
    statlabels = ['mean', 'perc. 10', 'perc. 90']

    data = get_kindists(filen_temp, [weight], simnames,
                        rrange_rvir=rrange_rvir,
                        vrrange=None,
                        vrrange_units='kmps',
                        dist_target='vr', vrbin_fac=1,
                        ssel=range(6))
    data = data[0]
    zs = list({np.round(datum['cosmopars']['z'], 2) 
               for l1 in data for datum in l1})
    zs.sort()
    panelsize = 2.5
    nrows = len(stats)
    ncols = 3
    hspace = 0.
    wspace = 0.
    lspace = 2.
    height_ratios = [panelsize] * nrows + [lspace]
    width_ratios = [panelsize] * ncols
    width = sum(width_ratios)
    height = sum(height_ratios)
    ncol_legend = int(np.floor(1.7 * ncols))
    physcols = {'noBH': 0,
                'AGN-noCR': 1,
                'AGN-CR': 2}

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows + 1,
                        height_ratios=height_ratios,
                        width_ratios=width_ratios,
                        hspace=hspace, wspace=wspace)
    axes = [[fig.add_subplot(grid[i, j]) 
             for j in range(3)]
            for i in range(len(stats))]
    lax = fig.add_subplot(grid[-1, :])
    lax.axis('off')

    fontsize = 12
    linewidth = 1.0
    iccolors = sl.m12_iccolors
    iccolors.update(sl.m13_iccolors)

    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    
    for pi, phys in enumerate(physcols):
        for si, statl in enumerate(statlabels):
            ax = axes[si][pi]
            bottom = si == len(stats) - 1
            left = pi == 0
            ax.tick_params(which='both', labelsize=fontsize - 1.,
                           direction='in', right=True, top=True,
                           labelleft=left, labelbottom=bottom)
            ax.grid(visible=True, which='both')
            if left:
                ax.set_ylabel(ylabel, fontsize=fontsize)
            if bottom:
                ax.set_xlabel(xlabel, fontsize=fontsize)
            
            axlabel = f'{phys}, {statl}'
            ax.text(0.98, 0.98, axlabel, fontsize=fontsize,
                    horizontalalignment='right', verticalalignment='top',
                    transform=ax.transAxes, color='black')
    ics = []
    for smi, simname in enumerate(simnames):
        ic = simname.split('_')[0]
        ics.append(ic)
        phys = getphys(simname)
        coli = physcols[phys]
        sdata = data[smi]
        color = iccolors[ic]
        xs_toplot = [[] for _ in range(len(stats))]
        ys_toplot = [[] for _ in range(len(stats))]
        for datum in sdata:
            _hist = datum['hist']
            _edges = datum['edges']
            zv = datum['cosmopars']['z']
            if np.sum(_hist) == 0.:
                print('zero weight for ', axlabel, 
                      ' Rvir ranges: ', rrange_rvir,
                      ' simname: ', simname)
                continue
            for si, stat in enumerate(stats):
                if stat == 'mean':
                    bcens = 0.5 * (_edges[:-1] + _edges[1:])
                    av = np.average(bcens, weights=_hist)
                    ys_toplot[si].append(av)
                elif stat.startswith('p'):
                    pval = np.array([float(stat[1:]) / 100.])
                    perc = mu.percentiles_from_histogram(_hist, _edges, 
                                                         axis=0,
                                                         percentiles=pval)
                    ys_toplot[si].append(perc[0])
                xs_toplot[si].append(zv)
        for si, (xp, yp) in enumerate(zip(xs_toplot, ys_toplot)):
            xp = np.array(xp)
            yp = np.array(yp)
            xsort = np.argsort(xp)
            xp = xp[xsort]
            yp = yp[xsort]
            if simname_styles is None:
                ls = 'solid'
            else:
                ls = simname_styles[simname]
            axes[si][coli].plot(xp, yp, color=color,
                                linestyle=ls, linewidth=linewidth,
                                marker='o', markersize=5.)
    
    ylims = [[ax.get_ylim() for ax in l1] for l1 in axes ]
    ymin = [min([ylim[0] for ylim in l1]) for l1 in ylims]
    ymax = [max([ylim[1] for ylim in l1]) for l1 in ylims]
    [[ax.set_ylim((ymin[i], ymax[i])) for ax in l1] 
     for i, l1 in enumerate(axes)]
    xlims = [[ax.get_xlim() for ax in l1] for l1 in axes ]
    xmin = [min([xlim[0] for xlim in l1]) for l1 in xlims]
    xmax = [max([xlim[1] for xlim in l1]) for l1 in xlims]
    [[ax.set_xlim((xmin[i], xmax[i])) for ax in l1] 
     for i, l1 in enumerate(axes)]

    ics = list(set(ics))
    ics.sort()
    handles2 = [mlines.Line2D((), (), color=iccolors[ic], 
                              linestyle='solid',
                              linewidth=linewidth, marker='o',
                              markersize=5.,
                              label=ic)
                for ic in ics]
    lax.legend(handles=handles2, fontsize=fontsize,
               ncol=ncol_legend, loc='upper center',
               bbox_to_anchor=(0.5, 0.65))
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')
    
def compset_vperc_physmodels():
    ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
    filetemp = ('hist_rcen_vcen_temperature_by_{weight}_{simname}'
                '_snap{snapnum}_bins1_v1_hvcen.hdf5')
    outdir = '/projects/b1026/nastasha/imgs/r_vr_hists/'

    simnames_m13 = sl.m13_hr_clean2 + sl.m13_sr_clean2
    simnames_m12 = [('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                    '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp2e-4_gacc31_fa0.5'),
                ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp1e10_gacc31_fa0.5'),
                ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                    '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp2e-4_gacc31_fa0.5'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp1e10_gacc31_fa0.5')]
    simnames_m13.sort()
    simnames_m12.sort()
    simsets = [simnames_m12, simnames_m13]
    simsetnames = ['m12', 'm13']
    for simset, simsetn in zip(simsets, simsetnames):
        for rrange in [(0.1, 1.), (0.15, 0.25), (0.45, 0.55), (0.9, 1.0)]:
            filen_temp = ddir + filetemp
            title = (f'{simsetn}, '
                     f'${rrange[0]:.2f} \\endash {rrange[-1]:.2f}'
                      ' \\, \\mathrm{R}_{\\mathrm{vir}}$')
            _outname = (f'vperc_physcomp_{simsetn}'
                        f'_{rrange[0]:.2f}'
                        f'_{rrange[-1]:.2f}_Rvir_v1')
            _outname = _outname.replace('.', 'p')
            outname = outdir + _outname + '.pdf'

            comp_vperc_physmodels(filen_temp, simset, weight='gasvol',
                                  title=title, outname=outname,
                                  rrange_rvir=rrange, 
                                  label_weightpart='\\mathrm{V}')

def comp_vfrac_physmodels(filen_temp, simnames, weight='gasvol',
                          title=None, outname=None,
                          rrange_rvir=(0.1, 1.), 
                          label_weightpart=None,
                          vrranges=[(-np.inf, np.inf)],
                          vrranges_units=['kmps'],
                          simname_styles=None):
    
    ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
    _filen_temp = ('hist_rcen_vcen_temperature_by_{weight}_{simname}'
                   '_snap{snapnum}_bins1_v1_hvcen.hdf5')
    filen_temp = ddir + _filen_temp
    if label_weightpart is None:
        label_weightpart = '\\mathrm{weight}'
    ylabel = (f'${label_weightpart}\\, \\mathrm{{frac.}}$')
    xlabel = '$z$'
    allvkey = (-np.inf, np.inf)
    _vrranges = [allvkey] + vrranges
    _vrranges_units = ['kmps'] + vrranges_units
    data = {vrrange: get_kindists(filen_temp, [weight], simnames,
                                  rrange_rvir=rrange_rvir,
                                  vrrange=vrrange,
                                  vrrange_units=vrrange_units,
                                  dist_target='vr', vrbin_fac=1,
                                  ssel=range(6))[0]
            for vrrange, vrrange_units in zip(_vrranges, _vrranges_units)}
    zs = list({np.round(datum['cosmopars']['z'], 2) 
               for key in data for l1 in data[key] for datum in l1})
    zs.sort()
    panelsize = 2.5
    nrows = len(vrranges)
    ncols = 3
    hspace = 0.
    wspace = 0.
    lspace = 2.
    height_ratios = [panelsize] * nrows + [lspace]
    width_ratios = [panelsize] * ncols
    width = sum(width_ratios)
    height = sum(height_ratios)
    ncol_legend = int(np.floor(1.4 * ncols))
    physcols = {'noBH': 0,
                'AGN-noCR': 1,
                'AGN-CR': 2}

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows + 1,
                        height_ratios=height_ratios,
                        width_ratios=width_ratios,
                        hspace=hspace, wspace=wspace)
    axes = [[fig.add_subplot(grid[i, j]) 
             for j in range(3)]
            for i in range(len(vrranges))]
    lax = fig.add_subplot(grid[-1, :])
    lax.axis('off')

    fontsize = 12
    linewidth = 1.0
    iccolors = sl.m12_iccolors
    iccolors.update(sl.m13_iccolors)

    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    
    for pi, phys in enumerate(physcols):
        for vi, (vrrange, vrrange_units) in enumerate(zip(vrranges,
                                                          vrranges_units)):
            ax = axes[vi][pi]
            bottom = vi == len(vrranges) - 1
            left = pi == 0
            ax.tick_params(which='both', labelsize=fontsize - 1.,
                           direction='in', right=True, top=True,
                           labelleft=left, labelbottom=bottom)
            ax.grid(visible=True, which='both')
            if left:
                ax.set_ylabel(ylabel, fontsize=fontsize)
            if bottom:
                ax.set_xlabel(xlabel, fontsize=fontsize)
            if vrrange_units == 'kmps':
                upart = '\\mathrm{km} \\, / \\, \\mathrm{s}'
                v0 = f'{vrrange[0]:.0f}'
                v1 = f'{vrrange[1]:.0f}'
            elif vrrange_units == 'vesc':
                upart = 'v_{\\mathrm{esc}}'
                v0 = f'{vrrange[0]:.2f}'
                v1 = f'{vrrange[1]:.2f}'
            if np.isfinite(vrrange[0]) and np.isfinite(vrrange[1]):
                rpart = f'{v0}\\endash{v1}'
            elif np.isfinite(vrrange[0]):
                rpart = f'>{v0}'
            elif np.isfinite(vrrange[1]):
                rpart = f'<{v1}'
            vlabel = f'$v_{{\\mathrm{{rad}}}} {rpart} \\; {upart}$'
            axlabel = f'{vlabel}\n{phys}'
            ax.text(0.98, 0.98, axlabel, fontsize=fontsize,
                    horizontalalignment='right', verticalalignment='top',
                    transform=ax.transAxes, color='black')
    ics = []
    for smi, simname in enumerate(simnames):
        ic = simname.split('_')[0]
        ics.append(ic)
        phys = getphys(simname)
        coli = physcols[phys]
        color = iccolors[ic]
        xs_toplot = [[] for _ in range(len(vrranges))]
        ys_toplot = [[] for _ in range(len(vrranges))]
        for vi, vrrange in enumerate(vrranges):
            _data_sub = data[vrrange][smi]
            _data_all = data[allvkey][smi]
            for dtot, dsub in zip(_data_all, _data_sub):
                _sub = np.sum(dsub['hist'])
                _tot = np.sum(dtot['hist'])
                zv = dsub['cosmopars']['z']
                ys_toplot[vi].append(_sub / _tot)
                xs_toplot[vi].append(zv)
        for vi, (xp, yp) in enumerate(zip(xs_toplot, ys_toplot)):
            xp = np.array(xp)
            yp = np.array(yp)
            xsort = np.argsort(xp)
            xp = xp[xsort]
            yp = yp[xsort]
            if simname_styles is None:
                ls = 'solid'
            else:
                ls = simname_styles[simname]
            axes[vi][coli].plot(xp, yp, color=color,
                                linestyle=ls, linewidth=linewidth,
                                marker='o', markersize=5.)
    
    ylims = [[ax.get_ylim() for ax in l1] for l1 in axes ]
    ymin = [min([ylim[0] for ylim in l1]) for l1 in ylims]
    ymax = [max([ylim[1] for ylim in l1]) for l1 in ylims]
    yrs = [_max - _min for _max, _min in zip(ymax, ymin)]
    [[ax.set_ylim((ymin[i], ymax[i] + 0.15 * yrs[i])) for ax in l1]
     for i, l1 in enumerate(axes)]
    xlims = [[ax.get_xlim() for ax in l1] for l1 in axes ]
    xmin = [min([xlim[0] for xlim in l1]) for l1 in xlims]
    xmax = [max([xlim[1] for xlim in l1]) for l1 in xlims]
    [[ax.set_xlim((xmin[i], xmax[i])) for ax in l1] 
     for i, l1 in enumerate(axes)]

    ics = list(set(ics))
    ics.sort()
    handles2 = [mlines.Line2D((), (), color=iccolors[ic], 
                              linestyle='solid',
                              linewidth=linewidth, marker='o',
                              markersize=5.,
                              label=ic)
                for ic in ics]
    lax.legend(handles=handles2, fontsize=fontsize,
               ncol=ncol_legend, loc='upper center',
               bbox_to_anchor=(0.5, 0.65))
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def compset_vfrac_physmodels():
    ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
    filetemp = ('hist_rcen_vcen_temperature_by_{weight}_{simname}'
                '_snap{snapnum}_bins1_v1_hvcen.hdf5')
    outdir = '/projects/b1026/nastasha/imgs/r_vr_hists/'

    simnames_m13 = sl.m13_hr_clean2 + sl.m13_sr_clean2
    simnames_m12 = [('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                    '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp2e-4_gacc31_fa0.5'),
                ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp1e10_gacc31_fa0.5'),
                ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                    '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp2e-4_gacc31_fa0.5'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                    '_sdp1e10_gacc31_fa0.5')]
    simnames_m13.sort()
    simnames_m12.sort()
    simsets = [simnames_m12, simnames_m13]
    simsetnames = ['m12', 'm13']
    for simset, simsetn in zip(simsets, simsetnames):
        for rrange in [(0.1, 1.), (0.15, 0.25), (0.45, 0.55), (0.9, 1.0)]:
            filen_temp = ddir + filetemp
            title = (f'{simsetn}, '
                     f'${rrange[0]:.2f} \\endash {rrange[-1]:.2f}'
                      ' \\, \\mathrm{R}_{\\mathrm{vir}}$')
            _outname = (f'vpfrac_physcomp_{simsetn}'
                        f'_{rrange[0]:.2f}'
                        f'_{rrange[-1]:.2f}_Rvir_v1')
            _outname = _outname.replace('.', 'p')
            outname = outdir + _outname + '.pdf'
            if simsetn == 'm12':
                vrranges = [(-np.inf, -0.5), (-np.inf, -100.),
                            (0.5, np.inf), (100., np.inf)]
            elif simsetn == 'm13':
                vrranges = [(-np.inf, -0.5), (-np.inf, -150.),
                            (0.5, np.inf), (150., np.inf)]
            vrranges_units = (['vesc'] + ['kmps']) * 2

            comp_vfrac_physmodels(filen_temp, simset, weight='gasvol',
                                  title=title, outname=outname,
                                  rrange_rvir=rrange, 
                                  label_weightpart='\\mathrm{V}',
                                  vrranges=vrranges,
                                  vrranges_units=vrranges_units)

def compset_vfrac_physmodels_summaryversion():
    ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
    filetemp = ('hist_rcen_vcen_temperature_by_{weight}_{simname}'
                '_snap{snapnum}_bins1_v1_hvcen.hdf5')
    outdir = '/projects/b1026/nastasha/imgs/r_vr_hists/'

    simnames_m13 = sl.m13_hr_all2 + sl.m13_sr_all2
    simnames_m12 = sl.m12_hr_all2 + sl.m12_sr_all2
    for simname in sl.buglist1:
        if simname in simnames_m13:
            simnames_m13.remove(simname)
        if simname in simnames_m12:
            simnames_m12.remove(simname)
    sn_all = simnames_m12 + simnames_m13
    ics = []
    icsets = []
    for _sn in sn_all:
        ic = _sn.split('_')[0]
        if ic in ics:
            ici = np.where([ic == _ic for _ic in ics])[0][0]
            icsets[ici].append(_sn)
        else:
            ics.append(ic)
            icsets.append([_sn])
    icis_allphys = np.where([len(icset) == 3 for icset in icsets])
    ics_allphys = np.array(ics)[icis_allphys]
    simname_styles = {simname: ('solid' if np.any([simname.startswith(ic) 
                                                  for ic in ics_allphys])
                                 else 'dotted')
                      for simname in sn_all}
    simnames_m13.sort()
    simnames_m12.sort()
    simsets = [simnames_m12, simnames_m13]
    simsetnames = ['m12', 'm13']
    for simset, simsetn in zip(simsets, simsetnames):
        for rrange in [(0.1, 1.), (0.15, 0.25), (0.45, 0.55), (0.9, 1.0)]:
            filen_temp = ddir + filetemp
            #title = (f'{simsetn}, '
            #         f'${rrange[0]:.2f} \\endash {rrange[-1]:.2f}'
            #          ' \\, \\mathrm{R}_{\\mathrm{vir}}$')
            title = None
            _outname = (f's1_vpfrac_physcomp_{simsetn}'
                        f'_{rrange[0]:.2f}'
                        f'_{rrange[-1]:.2f}_Rvir')
            _outname = _outname.replace('.', 'p')
            outname = outdir + _outname + '.pdf'
            if simsetn == 'm12':
                vrranges = [(0.5, np.inf), (-np.inf, -0.5)]
            elif simsetn == 'm13':
                vrranges = [(0.5, np.inf), (-np.inf, -0.5)]
            vrranges_units = ['vesc'] * 2

            comp_vfrac_physmodels(filen_temp, simset, weight='gasvol',
                                  title=title, outname=outname,
                                  rrange_rvir=rrange, 
                                  label_weightpart='\\mathrm{V}',
                                  vrranges=vrranges,
                                  vrranges_units=vrranges_units,
                                  simname_styles=simname_styles)