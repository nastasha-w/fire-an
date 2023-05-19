
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
from pyparsing import line

def getphys(simname):
    phys = ('noBH' if '_sdp1e10_' in simname 
            else 'AGN-CR' if '_MHDCRspec1_' in simname 
            else 'AGN-noCR')
    return phys

def get_kindist(filen, rrange_rvir=(0.1, 1.0),
                vrrange=None,
                vrrange_units='kmps',
                dist_target='vr'):
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
        _hist = np.sum(hist[hsel], axes=sumaxes)
        _edges = f[f'axis_{binaxis}/bins'][:] * binconv
    out = {'cosmopars': cosmopars,
           'vescvir_kmps': vescvir_kmps,
           'mvir_g': mvir_g,
           'rvir_cm': rvir_cm,
           'edges': _edges,
           'hist': _hist}
    return out

def get_kindists(filen_temp, weights, simnames,
                 **kwargs):
    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    sfills = {'simname': simname for simname in simnames}
    wfills = {'weight': weight for weight in weights}
    zfills = [{'snapnum': snap for snap in snaps_sr}
              if simname in sims_sr else 
              {'snapnum': snap for snap in snaps_hr}
              if simname in sims_hr else
              None
              for simname in simnames]
    out = [[[get_kindist(filen_temp.format(**(sfill | wfill | zfill)),
                          **kwargs)
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
    get_zcolors = pu.getzcolorfunc(zvals, ztol=1e-3)
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
    ylabel = f'$\\partial \\, \\mathrm{{weight frac.}} \\,/\\, {yl_add}$'

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
    axes = [fig.add_subplot(grid[i % ncols, i // ncols]) for i in range(nw)]
    lax = fig.add_subplot(grid[-1, :])
    lax.axis('off')

    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    
    yminmax = np.inf
    xminmin = np.inf
    xmaxmax = -np.inf
    xmar = 0.05
    yrange = 4.
    ymar_top = 0.07
    
    for wi, (wdata, wlabel) in enumerate(zip(data, weightlabels)):
        ax = axes[wi]
        bottom = wi >= nw - ncols
        left = wi % ncols == 0
        ax.tick_params(which='both', labelsize=fontsize - 1.,
                       direction='in', right=True, top=True,
                       labelleft=left, labelbottom=bottom)
        ax.grid(visible=True, which='both')
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
                phist = _hist / np.sum(_hist) / np.diff(_edges)
                phist = np.append(phist, [0.])
                zv = datum['cosmopars']['z']
                color = get_zcolors(zv)
                ax.step(_edges, phist, where='post', color=color,
                        linestyle=linestyle, linewidth=linewidth,
                        patheffects=patheff)
                _ymax = np.max(phist)
                yminmax = min(yminmax, _ymax)
                goodinds = np.where(phist >= _ymax * 10**(-yrange))[0]
                _xmin = _edges[goodinds[0]]
                _xmax = _edges[goodinds[-1] + 1]
                xminmin = min(xminmin, _xmin)
                xmaxmax = max(xmaxmax, _xmax)
    
    ymax = yminmax
    ymin = ymax * 10**(-yrange)
    ymax = ymax + (ymax - ymin) * ymar_top
    [ax.set_ylim((ymin, ymax)) for ax in axes]
    xmar = xmar * (xmaxmax - xminmin)
    xmin = xminmin - xmar
    xmax = xmaxmax + xmar
    [ax.set_xlim((xmin, xmax)) for ax in axes]

    handles1 = [mlines.Line2D((), (), color='gray', linestyle=physstyles[phys],
                              linewidth=linewidth, patheffects=patheff,
                              label=phys)
                for phys in physstyles]
    handles2 = [mlines.Line2D((), (), color=get_zcolors(zv), linestyle='solid',
                              linewidth=linewidth, patheffects=patheff,
                              label=f'z={zv:.1f}')
                for zv in zvals]
    lax.legend(handles=handles1 + handles2, fontsize=fontsize,
               ncol=ncol_legend, legend_loc='upper center',
               bbox_to_anchor=(0.5, 0.65))
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')
    
            
def plot_vr_dists():
    ddir = '/projects/b1026/nastasha/hists/r_vr_clean2_nobug/'
    filetemp = ('hist_rcen_vcen_temperature_by_{{weight}}_{simname}'
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
        for rrange_rvir in [(0.1, 1.), (0.15, 0.15), 
                            (0.45, 0.55), (0.9, 1.0)]:
            filen_temp = ddir + filetemp
            title = (f'{ic}, ${rrange_rvir[0]:.2f} \\endash '
                     f'{rrange_rvir[-1]:.2f} \\,'
                     ' \\mathrm{R}_{\\mathrm{vir}}$')
            _outname = (f'dists_vr_{ic}_phys_z_wt_comp_{rrange_rvir[0]:.2f}'
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

        


