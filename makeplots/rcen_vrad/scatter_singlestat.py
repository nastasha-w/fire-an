from turtle import circle
import h5py 
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import fire_an.ionrad.get_cieranges as gcr
import fire_an.makeplots.plot_utils as pu
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.math_utils as mu


def getphys(simname):
    phys = ('noBH' if '_sdp1e10_' in simname 
            else 'AGN-CR' if '_MHDCRspec1_' in simname 
            else 'AGN-noCR')
    return phys

def fmtdcts(str, *args):
    dct = {}
    for arg in args:
        dct.update(arg)
    _str = str.format(**dct)
    return _str

def getvperc_rbins(filen, rrange_rvir=(0.1, 1.),
                   vperc=(0.1, 0.9)):
    with h5py.File(filen, 'r') as f:
        hist = f['histogram/histogram'][:]
        islog = bool(f['histogram'].attrs['log'])
        if islog:
            hist = 10**hist
        mvir_g = f['Header/inputpars/halodata'].attrs['Mvir_g']
        rvir_cm = f['Header/inputpars/halodata'].attrs['Rvir_cm']
        vescvir_kmps = np.sqrt(2. * c.gravity * mvir_g / rvir_cm) \
                       * 1e-5
        hsel = [slice(None, None, None)] * 3
        sumaxes = (0, 2)
        cosmopars = {key: val for key, val 
                     in f['Header/cosmopars'].attrs.items()}
        
        if rrange_rvir is not None:
            rbins_rvir = f['axis_0/bins'][:]
            rimin = np.where(np.isclose(rrange_rvir[0], rbins_rvir))[0][0]
            rimax = np.where(np.isclose(rrange_rvir[-1], rbins_rvir))[0][0]
            if rrange_rvir[0] == -np.inf:
                rimin = 0
            if rrange_rvir[-1] == np.inf:
                rimax = len(rbins_rvir)
            hsel[0] = slice(rimin, rimax, None)
        _hist = np.sum(hist[tuple(hsel)], axis=sumaxes)
        vedges = f['axis_1/bins'][:] * 1e-5 # cm/s -> km/s
        vcens = 0.5 * (vedges[:-1] + vedges[1:])
        mean = np.average(vcens, weights=_hist)
        perc = mu.percentiles_from_histogram(_hist, vedges, axis=0,
                                             percentiles=np.array(vperc))
    out = {'cosmopars': cosmopars,
           'vescvir_kmps': vescvir_kmps,
           'mvir_g': mvir_g,
           'rvir_cm': rvir_cm,
           'mid': mean,
           'lo': perc[0],
           'hi': perc[-1]}
    return out

def getthemperc_rbins(filen, rrange_rvir=(0.1, 1.),
                      vrrange=None, vrrange_units='kmps',
                      perc=np.array([0.1, 0.5, 0.9])):
    with h5py.File(filen, 'r') as f:
        hist = f['histogram/histogram'][:]
        islog = bool(f['histogram'].attrs['log'])
        if islog:
            hist = 10**hist
        mvir_g = f['Header/inputpars/halodata'].attrs['Mvir_g']
        rvir_cm = f['Header/inputpars/halodata'].attrs['Rvir_cm']
        vescvir_kmps = np.sqrt(2. * c.gravity * mvir_g / rvir_cm) \
                       * 1e-5
        hsel = [slice(None, None, None)] * 3
        sumaxes = (0, 1)
        cosmopars = {key: val for key, val 
                     in f['Header/cosmopars'].attrs.items()}
        
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
        _hist = np.sum(hist[tuple(hsel)], axis=sumaxes)
        qtedges = f['axis_2/bins'][:]
        perc = mu.percentiles_from_histogram(_hist, qtedges, axis=0,
                                             percentiles=np.array(perc))
    out = {'cosmopars': cosmopars,
           'vescvir_kmps': vescvir_kmps,
           'mvir_g': mvir_g,
           'rvir_cm': rvir_cm,
           'lo': perc[0],
           'mid': perc[1],
           'hi': perc[2]}
    return out
            
def plotdata_censcatter(datax, datay, xweightmap,
                        yweightmap, xlabel=None, ylabel=None,
                        ncols=None, fontsize=12,
                        wspace=0.0, hspace=0.0,
                        syncaxlims=True):
    '''
    plot data points against each other, with x/y points matched
    by ics/phys. model/snapshot. 

    datax/y: should be organised as 
    dict(simname: list(dict(weight: inner dict for each snapshot))
    the inner dict should contain keys 'cosmopars', 'lo', 'mid', 'hi'
    if 'lo' and 'hi' are omitted, no error bars are drawn.
    simnames, snapshots should be matched between datax and datay

    xweightmap, yweightmap: which weight key in datax and datay,
    respectively, to use in each subpanel (lists). 
    '''
    physmarkers = {'noBH': '*',
                   'AGN-noCR': 'o',
                   'AGN-CR': 'P'}
    markeredgewidth = 2
    size = 5
    fontsize = 12
    
    iccolors = sl.m12_iccolors 
    iccolors.update(sl.m13_iccolors)
    
    npanels = len(yweightmap)
    if ncols is None:
        ncmax = 4
        ncols = min(ncmax, npanels)
    nrows = (npanels - 1) // ncols + 1
    panelsize = 2.5
    laxspace = 1.5
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] * nrows + [laxspace]
    width = sum(width_ratios) \
            * (1. + wspace * sum(width_ratios) / (len(width_ratios) - 1)) 
    height = sum(height_ratios) \
            * (1. + hspace * sum(height_ratios) / (len(height_ratios) - 1))

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows + 1, hspace=hspace,
                        wspace=wspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols])
            for i in range(npanels)]
    lax = fig.add_subplot(grid[nrows, :])
    
    simnames = list(datay.keys())
    simnames.sort()
    ics = list({simname.split('_')[0] for simname in simnames})
    ics.sort()
    zs = {np.round(wdict[key]['cosmopars']['z'], 1) 
          for key1 in datay for wdict in datay[key1] for key in wdict}
    zs = list(zs)
    zs.sort()
    zcolors = pu.getzcolorfunc(zs)

    for axi, ax in enumerate(axes):
        bottom = axi >= npanels - ncols
        left = axi % ncols == 0
        ax.tick_params(which='both', direction='in', labelsize=fontsize - 1,
                       right=True, top=True, labelbottom=bottom,
                       labelleft=left)
        if bottom and xlabel is not None:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        if left and ylabel is not None:
            ax.set_ylabel(ylabel, fontsize=fontsize)
    lax.axis('off')
    if hspace == 0.:
        ltop = 0.60
    else:
        ltop = 1.
    if ncols == 4:
        ncol_leg = (1, 3, 2)
    elif ncols == 3:
        ncol_leg = (1, 2, 1)
    handles1 = [mlines.Line2D((), (), label=ic, markersize=size,
                              markerfacecolor='gray', marker='o',
                              markeredgecolor=iccolors[ic],
                              linestyle='none',
                              markeredgewidth=markeredgewidth)
                for ic in ics]
    leg1 = lax.legend(handles=handles1, fontsize=fontsize,
                      bbox_to_anchor=(0.0, ltop), loc='upper left',
                      handletextpad=0.4, columnspacing=1.0,
                      ncol=ncol_leg[0])
    handles2 = [mlines.Line2D((), (), label=f'z={z:.1f}', 
                              markersize=size,
                              markerfacecolor=zcolors(z), marker='o',
                              markeredgecolor='gray',
                              linestyle='none',
                              markeredgewidth=markeredgewidth)
                 for z in zs]
    leg2 = lax.legend(handles=handles2, fontsize=fontsize, 
                      ncol=ncol_leg[1],
                      bbox_to_anchor=(0.4, ltop), loc='upper center',
                      handletextpad=0.4, columnspacing=1.0)
    handles3 = [mlines.Line2D((), (), label=phys, 
                              markersize=size,
                              markerfacecolor='gray', 
                              marker=physmarkers[phys],
                              markeredgecolor='black',
                              linestyle='none',
                              markeredgewidth=markeredgewidth)
                for phys in physmarkers]
    leg3 = lax.legend(handles=handles3, fontsize=fontsize, 
                      ncol=ncol_leg[2],
                      bbox_to_anchor=(1.0, ltop), loc='upper right',
                      handletextpad=0.4, columnspacing=1.0)
    lax.add_artist(leg1)
    lax.add_artist(leg2)
    lax.add_artist(leg3)

    xmin = np.inf
    xmax = -np.inf
    ymin = np.inf
    ymax = -np.inf
    for simname in simnames:
        datax_zw = datax[simname]
        datay_zw = datay[simname]
        ic = simname.split('_')[0]
        phys = getphys(simname)
        marker = physmarkers[phys]
        edgecolor = iccolors[ic]
        for datax_w, datay_w in zip(datax_zw, datay_zw):
            for wi, (xweight, yweight) in enumerate(zip(
                    xweightmap, yweightmap)):
                ax = axes[wi]
                datumx = datax_w[xweight]
                datumy = datay_w[yweight]
                facecolor = zcolors(datumy['cosmopars']['z'])
                xmid = datumx['mid']
                ymid = datumy['mid']
                if 'lo' in datumx:
                    xlo = datumx['lo']
                else:
                    xlo = xmid
                if 'hi' in datumx:
                    xhi = datumx['hi']
                else:
                    xhi = xmid
                if 'lo' in datumy:
                    ylo = datumy['lo']
                else:
                    ylo = ymid
                if 'hi' in datumy:
                    yhi = datumy['hi']
                else:
                    yhi = ymid
                xerr = (xmid - xlo, xhi - xmid)
                yerr = (ymid - ylo, yhi - ymid)
                ax.errorbar([xmid], [ymid], xerr=[xerr], yerr=[yerr],
                            linestyle='none', ecolor=(0.5, 0.5, 0.5, 0.3),
                            marker=marker, markersize=size,
                            markerfacecolor=facecolor,
                            markeredgecolor=edgecolor)
                xmin = min(xmin, xmid)
                xmax = max(xmax, xmid)
                ymin = min(ymin, ymid)
                ymax = max(ymax, ymid)
    xran = xmax - xmin
    xmin = xmin - xran
    xmax = xmax + xran
    ymin = ymin - xran
    ymax = ymax + xran
    if syncaxlims:
        [ax.set_xlim((xmin, xmax)) for ax in axes]
        [ax.set_ylim((ymin, ymax)) for ax in axes]
    
    axdoc = {'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax,
             'fontsize': fontsize}
    return fig, axes, lax, axdoc

def plot_thermcomp_medscat(simnames, compqty='T',
                           haloweight='Volume',
                           rrange_rvir=(0.1, 1.),
                           vrrange=None, vrrange_units='kmps',
                           title=None, outname=None):
    ions = ['O6', 'Ne8', 'O7', 'Ne9', 'O8', 'Ne10']
    ionlabels = ['O VI', 'Ne VIII', 'O VII', 'Ne IX', 'O VIII', 'Ne X']
    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    sims_sr = sl.m13_sr_all2 + sl.m12_sr_all2
    sims_hr = sl.m13_hr_all2 + sl.m12_hr_all2
    ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
    filen_temp = ('hist_rcen_vcen_{compqty}_by_{weight}_{simname}'
                  '_snap{snapnum}_bins1_v1_hvcen.hdf5')
    if haloweight == 'Volume':
        fhaloweight = 'gasvol'
    elif haloweight == 'Mass':
        fhaloweight = 'gasmass'
    if compqty == 'T':
        fqty = 'temperature'
        flab = '\\mathrm{T}'
        funits = '\\mathrm{K}'
    elif compqty == 'nH':
        fqty = 'hdens'
        flab = '\\mathrm{n}_{\\mathrm{H}}'
        funits = '\\mathrm{cm}^{-3}'
    xlabel = (f'$\\left\\langle {flab} \\right\\rangle'
              f'_{{\\mathrm{{{haloweight}}}}}'
              f'\\; [{funits}]$')
    ylabel = (f'$\\left\\langle {flab} \\right\\rangle'
              '_{\\mathrm{ion}}'
              f'\\; [{funits}]$')
    showperc = (0.1, 0.5, 0.9)
    fontsize = 12
    
    snaplists = [snaps_sr if simname in sims_sr
                 else snaps_hr if simname in sims_hr
                 else None
                 for simname in simnames]
    datay = {simname: [{weight: 
                        getthemperc_rbins(ddir + filen_temp.format(
                        weight=weight, simname=simname, snapnum=snap,
                        compqty=fqty), 
                                          rrange_rvir=rrange_rvir,
                                          vrrange=vrrange, 
                                          vrrange_units=vrrange_units,
                                          perc=np.array(showperc))
              for weight in ions}
             for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)}
    datax = {simname: [{weight: 
                        getthemperc_rbins(ddir + filen_temp.format(
                        weight=weight, simname=simname, snapnum=snap,
                        compqty=fqty), 
                                          rrange_rvir=rrange_rvir,
                                          vrrange=vrrange, 
                                          vrrange_units=vrrange_units,
                                          perc=np.array(showperc))
              for weight in [fhaloweight]}
             for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)}
    xweightmap = [fhaloweight] * len(ions)
    yweightmap = ions
    fig, axes, lax, axdoc = plotdata_censcatter(datax, datay, xweightmap,
                                                yweightmap, xlabel=xlabel,
                                                ylabel=ylabel,
                                                ncols=None,
                                                fontsize=fontsize)
    xmin = axdoc['xmin']
    xmax = axdoc['xmax']
    ymin = axdoc['ymin']
    ymax = axdoc['ymax']
    eqp = (max(xmin, ymin), min(xmax, ymax))
    [ax.plot(eqp, eqp, color='black', linestyle='dotted', linewidth=1,
             zorder=-1) for ax in axes]
    fontsize = axdoc['fontsize']
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    for ioni, ionl in enumerate(ionlabels):
        axes[ioni].text(0.05, 0.95, ionl, fontsize=fontsize,
                horizontalalignment='left', verticalalignment='top',
                transform=axes[ioni].transAxes)
    
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotset_thermcomp_medscat(haloweight='Volume', qty='T'):
    rranges_rvir = [(0.1, 1.), (0.15, 0.25), (0.45, 0.55),
                    (0.9, 1.0)]
    simnames_m12 = [('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
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
    simnames_m13 = sl.m13_hr_clean2 + sl.m13_sr_clean2
    
    for rrange in rranges_rvir:
        for simset, simsetname in zip([simnames_m12, simnames_m13],
                                      ['m12', 'm13']):
            title = (f'{simsetname}, '
                     f'${rrange[0]:.2f} \\endash {rrange[-1]:.2f}'
                     ' \\, \\mathrm{R}_{\\mathrm{vir}}$')
            outdir = '/projects/b1026/nastasha/imgs/summary_plots/'
            _outname = (f'{qty}bias_ions_{simsetname}_{haloweight}'
                        f'_{rrange[0]:.2f}'
                        f'_{rrange[-1]:.2f}_Rvir_med_scatter'
                        f'_0.10_0.90_v1')
            _outname = _outname.replace('.', 'p')
            outname = outdir + _outname + '.pdf'
                           
            plot_thermcomp_medscat(simset, compqty=qty,
                                  haloweight=haloweight,
                                  rrange_rvir=rrange,
                                  vrrange=None, vrrange_units='kmps',
                                  title=title, outname=outname)
            
def plot_vmean_ion_vs_halo(simnames, haloweight='Volume',
                           errbperc=(0.1, 0.9),
                           rrange_rvir=(0.1, 1.),
                           title=None, outname=None):
    ions = ['O6', 'Ne8', 'O7', 'Ne9', 'O8', 'Ne10']
    ionlabels = ['O VI', 'Ne VIII', 'O VII', 'Ne IX', 'O VIII', 'Ne X']
    fontsize = 12
    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    sims_sr = sl.m13_sr_all2 + sl.m12_sr_all2
    sims_hr = sl.m13_hr_all2 + sl.m12_hr_all2
    ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
    filen_temp = ('hist_rcen_vcen_temperature_by_{weight}_{simname}'
                 '_snap{snapnum}_bins1_v1_hvcen.hdf5')
    if haloweight == 'Volume':
        fhaloweight = 'gasvol'
    elif haloweight == 'Mass':
        fhaloweight = 'gasmass'

    xlabel = ('$\\left\\langlev_{\\mathrm{rad}} \\right\\rangle'
              f'_{{\\mathrm{{{haloweight}}}}}'
              '\\; [\\mathrm{km} \\, \\mathrm{s}^{-1}]$')
    ylabel = ('$\\left\\langlev_{\\mathrm{rad}} \\right\\rangle'
              '_{\\mathrm{ion}}'
              '\\; [\\mathrm{km} \\, \\mathrm{s}^{-1}]$')

    snaplists = [snaps_sr if simname in sims_sr
                 else snaps_hr if simname in sims_hr
                 else None
                 for simname in simnames]
    datay = {simname: [{weight: getvperc_rbins(ddir + filen_temp.format(
                            weight=weight, simname=simname, snapnum=snap),
                            rrange_rvir=rrange_rvir, vperc=errbperc)
                        for weight in ions}
                       for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)}
    datax = {simname: [{weight: getvperc_rbins(ddir + filen_temp.format(
                            weight=weight, simname=simname, snapnum=snap),
                            rrange_rvir=rrange_rvir, vperc=errbperc)
                        for weight in [fhaloweight]}
                       for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)}
    yweightmap = ions
    xweightmap = [fhaloweight] * len(ions)
    fig, axes, lax, axdoc = plotdata_censcatter(datax, datay, xweightmap,
                                                yweightmap, xlabel=xlabel,
                                                ylabel=ylabel,
                                                ncols=None, fontsize=fontsize)
    xmin = axdoc['xmin']
    xmax = axdoc['xmax']
    ymin = axdoc['ymin']
    ymax = axdoc['ymax']
    eqp = (max(xmin, ymin), min(xmax, ymax))
    [ax.plot(eqp, eqp, color='black', linestyle='dotted', linewidth=1,
             zorder=-1) for ax in axes]
    fontsize = axdoc['fontsize']
    for ioni, ionl in enumerate(ionlabels):
        axes[ioni].text(0.05, 0.95, ionl, fontsize=fontsize,
                        horizontalalignment='left', verticalalignment='top',
                        transform=axes[ioni].transAxes)
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    
    eqp = (max(xmin, ymin), min(xmax, ymax))
    [ax.plot(eqp, eqp, color='black', linestyle='dotted', linewidth=1,
             zorder=-1) for ax in axes]

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotset_vmean_ion_vs_halo(haloweight='Volume'):
    errbperc = (0.1, 0.9)
    rranges_rvir = [(0.1, 1.), (0.15, 0.25), (0.45, 0.55),
                    (0.9, 1.0)]
    simnames_m12 = [('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
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
    simnames_m13 = sl.m13_hr_clean2 + sl.m13_sr_clean2
    
    for rrange in rranges_rvir:
        for simset, simsetname in zip([simnames_m12, simnames_m13],
                                      ['m12', 'm13']):
            title = (f'{simsetname}, '
                     f'${rrange[0]:.2f} \\endash {rrange[-1]:.2f}'
                     ' \\, \\mathrm{R}_{\\mathrm{vir}}$')
            outdir = '/projects/b1026/nastasha/imgs/summary_plots/'
            _outname = (f'vbias_ions_{haloweight}_{simsetname}'
                        f'_{rrange[0]:.2f}'
                        f'_{rrange[-1]:.2f}_Rvir_v_mean_scatter'
                        f'_{errbperc[0]:.2f}_{errbperc[-1]:.2f}_v1')
            _outname = _outname.replace('.', 'p')
            outname = outdir + _outname + '.pdf'
                           
            plot_vmean_ion_vs_halo(simset, haloweight=haloweight,
                                   errbperc=errbperc,
                                   rrange_rvir=rrange,
                                   title=title, 
                                   outname=outname)
            
def plot_thermcomp_vrange_medscat(simnames, vrranges, vrranges_units,
                                  compqty='T',
                                  rrange_rvir=(0.1, 1.),
                                  title=None, outname=None):
    wqty = 'gasvol'
    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    sims_sr = sl.m13_sr_all2 + sl.m12_sr_all2
    sims_hr = sl.m13_hr_all2 + sl.m12_hr_all2
    ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
    filen_temp = ('hist_rcen_vcen_{compqty}_by_{weight}_{simname}'
                  '_snap{snapnum}_bins1_v1_hvcen.hdf5')
    if compqty == 'T':
        fqty = 'temperature'
        flab = '\\log_{10} \\, \\mathrm{T}'
        funits = '\\mathrm{K}'
    elif compqty == 'nH':
        fqty = 'hdens'
        flab = '\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
        funits = '\\mathrm{cm}^{-3}'
    elif compqty == 'O':
        fqty = 'OxygenAbundance'
        flab = '\\log_{10} \\, \\mathrm{Z}\\mathrm{O}}'
        funits = ''
    elif compqty == 'Ne':
        fqty = 'NeonAbundance'
        flab = '\\log_{10} \\, \\mathrm{Z}\\mathrm{Ne}}'
        funits = ''

    xlabel = (f'$ {flab} \\; [{funits}]$' if funits != '' 
              else f'$ {flab}$')
    ylabel = ((f'$ {flab} \\; [{funits}],' if funits != '' 
               else f'$ {flab}$')
              + 'v_{\\mathrm{rad}}\\,\\mathrm{sel.}$')
    showperc = (0.1, 0.5, 0.9)
    fontsize = 12
    
    snaplists = [snaps_sr if simname in sims_sr
                 else snaps_hr if simname in sims_hr
                 else None
                 for simname in simnames]
    datay = {simname: [{vrrange: 
                        getthemperc_rbins(ddir + filen_temp.format(
                        weight=wqty, simname=simname, snapnum=snap,
                        compqty=fqty), 
                                          rrange_rvir=rrange_rvir,
                                          vrrange=vrrange, 
                                          vrrange_units=vrrange_units,
                                          perc=np.array(showperc))
              for vrrange, vrrange_units in zip(vrranges, vrranges_units)}
             for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)}
    allvkey = (-np.inf, np.inf)
    datax = {simname: [{allvkey: 
                        getthemperc_rbins(ddir + filen_temp.format(
                        weight=wqty, simname=simname, snapnum=snap,
                        compqty=fqty), 
                                          rrange_rvir=rrange_rvir,
                                          vrrange=None, 
                                          vrrange_units=None,
                                          perc=np.array(showperc))
              for _ in range(1)}
             for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)}
    yweightmap = vrranges
    xweightmap = [allvkey] * len(vrranges)
    fig, axes, lax, axdoc = plotdata_censcatter(datax, datay, xweightmap,
                                                yweightmap, xlabel=xlabel,
                                                ylabel=ylabel,
                                                ncols=None,
                                                fontsize=fontsize)
    xmin = axdoc['xmin']
    xmax = axdoc['xmax']
    ymin = axdoc['ymin']
    ymax = axdoc['ymax']
    eqp = (max(xmin, ymin), min(xmax, ymax))
    [ax.plot(eqp, eqp, color='black', linestyle='dotted', linewidth=1,
             zorder=-1) for ax in axes]
    fontsize = axdoc['fontsize']
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    for vi, (vrsel, vru) in enumerate(zip(vrranges, vrranges_units)):
        if vru == 'kmps':
            upart = '\\mathrm{km} \\, / \\, \\mathrm{s}'
            v0 = f'{vrsel[0]:.0f}'
            v1 = f'{vrsel[1]:.0f}'
        elif vru == 'vesc':
            upart = 'v_{\\mathrm{esc}}'
            v0 = f'{vrsel[0]:.2f}'
            v1 = f'{vrsel[1]:.2f}'
        if np.isfinite(vrsel[0]) and np.isfinite(vrsel[1]):
            rpart = f'{v0}\\endash{v1}'
        elif np.isfinite(vrsel[0]):
            rpart = f'>{v0}'
        elif np.isfinite(vrsel[1]):
            rpart = f'<{v1}'
        vlabel = f'$v_{{\\mathrm{{rad}}}} {rpart} \\; {upart}$'
        axes[vi].text(0.05, 0.95, vlabel, fontsize=fontsize,
                      horizontalalignment='left', verticalalignment='top',
                      transform=axes[vi].transAxes)
    
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotset_thermcomp_vrange_medscat(qty='T'):
    rranges_rvir = [(0.1, 1.), (0.15, 0.25), (0.45, 0.55),
                    (0.9, 1.0)]
    vrranges = [(50., np.inf), (100., np.inf), (150., np.inf),
                (0.5, np.inf),
                (-np.inf, -50.), (-np.inf, -100.), (-np.inf, -150.),
                (-np.inf, -0.5)]
    vrranges_units = (['kmps'] * 3 + ['vesc']) * 2
    simnames_m12 = [('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
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
    simnames_m13 = sl.m13_hr_clean2 + sl.m13_sr_clean2
    
    for rrange in rranges_rvir:
        for simset, simsetname in zip([simnames_m12, simnames_m13],
                                    ['m12', 'm13']):
            title = (f'{simsetname}, '
                    f'${rrange[0]:.2f} \\endash {rrange[-1]:.2f}'
                    ' \\, \\mathrm{R}_{\\mathrm{vir}}$')
            outdir = '/projects/b1026/nastasha/imgs/summary_plots/'
            _outname = (f'vbias_{simsetname}_{qty}_med_scatter'
                            '_0.10_0.90'
                        f'_{rrange[0]:.2f}'
                        f'_{rrange[-1]:.2f}_Rvir_'
                        f'vrrangevar_v1')
            _outname = _outname.replace('.', 'p')
            outname = outdir + _outname + '.pdf'
                        
            plot_thermcomp_vrange_medscat(simset, vrranges, 
                                          vrranges_units,
                                          compqty=qty,
                                          rrange_rvir=rrange,
                                          title=title, 
                                          outname=outname)

def summaryplot_ionvolcomp(simset='m12_clean2', rrange_rvir=(0.45, 0.55),
                           outname=None):
    ions = ['O6', 'Ne8', 'O7', 'O8']
    ionlabels = ['O VI', 'Ne VIII', 'O VII', 'O VIII']
    haloweight = 'gasvol'
    fqtys = ['temperature', 'hdens', 'OxygenAbundance', 'temperature']
    targetvals = ['T', 'nH', 'O', 'vr']
    showperc = (0.1, 0.5, 0.9)
    errbperc = (0.1, 0.9)

    buglist = sl.buglist1
    bugics = {simname.split('_')[0] for simname in buglist}
    bugics = list(bugics)
    if simset == 'm12_clean2':
        simnames = sl.m12_hr_clean2 + sl.m12_sr_clean2
        toscan = simnames.copy()
        for simname in toscan:
            if simname.split('_')[0] in bugics:
                simnames.remove(simname)
    elif simset == 'm12_all2':
        simnames = sl.m12_hr_all2 + sl.m12_sr_all2
        toscan = simnames.copy()
        for simname in toscan:
            if simname.split('_')[0] in buglist:
                simnames.remove(simname)
    elif simset == 'm13_clean2':
        simnames = sl.m13_hr_clean2 + sl.m13_sr_clean2
        toscan = simnames.copy()
        for simname in toscan:
            if simname.split('_')[0] in bugics:
                simnames.remove(simname)
    elif simset == 'm13_all2':
        simnames = sl.m13_hr_all2 + sl.m13_sr_all2
        toscan = simnames.copy()
        for simname in toscan:
            if simname.split('_')[0] in buglist:
                simnames.remove(simname)

    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    sims_sr = sl.m13_sr_all2 + sl.m12_sr_all2
    sims_hr = sl.m13_hr_all2 + sl.m12_hr_all2
    ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
    filen_temp = ('hist_rcen_vcen_{compqty}_by_{weight}_{simname}'
                  '_snap{snapnum}_bins1_v1_hvcen.hdf5')

    ylabels = [('$\\langle\\log_{10} \\, \\mathrm{T}'
                '\\rangle_{\\mathrm{ion}, \\mathrm{med.}}'
                '\\; [\\mathrm{K}]$'),
               ('$\\langle\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
                '\\rangle_{\\mathrm{ion}, \\mathrm{med.}}'
                '\\; [\\mathrm{cm}^{-3}]$'),
               ('$\\langle\\log_{10} \\, \\mathrm{Z}_{\\mathrm{O}}'
                '\\rangle_{\\mathrm{ion}, \\mathrm{med.}}$'),
               ('$\\langle v_{\\mathrm{rad}}'
                '\\rangle_{\\mathrm{ion}, \\mathrm{mean}}'
                '\\; [\\mathrm{km}\\, \\mathrm{s}^{-1}]$')]
    xlabels = [ylabel.replace('ion', 'V') for ylabel in ylabels]
    fontsize = 12
    
    snaplists = [snaps_sr if simname in sims_sr
                 else snaps_hr if simname in sims_hr
                 else None
                 for simname in simnames]
    datay = {simname: [{(weight, targetval): 
                        getthemperc_rbins(ddir + filen_temp.format(
                        weight=weight, simname=simname, snapnum=snap,
                        compqty=compqty), 
                                          rrange_rvir=rrange_rvir,
                                          vrrange=None, 
                                          vrrange_units=None,
                                          perc=np.array(showperc))
                        if targetval != 'vr' else
                        getvperc_rbins(ddir + filen_temp.format(
                        weight=weight, simname=simname, snapnum=snap,
                        compqty=compqty), 
                                       rrange_rvir=rrange_rvir,
                                       vperc=errbperc)
              for weight in ions 
              for compqty, targetval in zip(fqtys, targetvals)}
             for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)}
    datax = {simname: [{(haloweight, targetval): 
                        getthemperc_rbins(ddir + filen_temp.format(
                        weight=haloweight, simname=simname, snapnum=snap,
                        compqty=compqty), 
                                          rrange_rvir=rrange_rvir,
                                          vrrange=None, 
                                          vrrange_units=None,
                                          perc=np.array(showperc)) 
                        if targetval != 'vr' else
                        getvperc_rbins(ddir + filen_temp.format(
                        weight=haloweight, simname=simname, snapnum=snap,
                        compqty=compqty), 
                                       rrange_rvir=rrange_rvir,
                                       vperc=errbperc)
              for compqty, targetval in zip(fqtys, targetvals)}
             for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)}
    yweightmap = [(weight, targetval) 
                  for targetval in targetvals for weight in ions]
    xweightmap = [(haloweight, targetval) 
                  for targetval in targetvals for _ in ions]
    fig, axes, lax, axdoc = plotdata_censcatter(datax, datay, xweightmap,
                                                yweightmap, xlabel='',
                                                ylabel='',
                                                ncols=len(fqtys),
                                                fontsize=fontsize,
                                                syncaxlims=False,
                                                hspace=0.25)
    fontsize = axdoc['fontsize']
    nrows = len(fqtys)
    ncols = len(ions)
    for ri in range(nrows):
        axsel = slice(ri * ncols, (ri + 1) * ncols)
        xlims = [ax.get_xlim() for ax in axes[axsel]]
        xmin = min([xlim[0] for xlim in xlims])
        xmax = max([xlim[1] for xlim in xlims])
        [ax.set_xlim(*(xmin, xmax)) for ax in axes[axsel]]
        ylims = [ax.get_ylim() for ax in axes[axsel]]
        ymin = min([ylim[0] for ylim in ylims])
        ymax = max([ylim[1] for ylim in ylims])
        [ax.set_ylim(*(ymin, ymax)) for ax in axes[axsel]]
        eqp = (max(xmin, ymin), min(xmax, ymax))
        for ci in range(len(ions)):
            ax = axes[ri * ncols + ci]
            ax.plot(eqp, eqp, color='black', linestyle='dotted', 
                    linewidth=1, zorder=-1)
            ax.set_xlabel(xlabels[ri], fontsize=fontsize)
            ax.tick_params(labelbottom=True)
            if ci == 0:
                ax.set_ylabel(ylabels[ri], fontsize=fontsize)
            if ri == 0:
                ax.set_title(ionlabels[ci], fontsize=fontsize)
            if ri == 0:
                rng = gcr.cieranges1[ions[ci]]
                ax.axhline(rng[0], color='black', linestyle='dotted',
                           linewidth=1.)
                ax.axhline(rng[1], color='black', linestyle='dotted',
                           linewidth=1.)
            
                
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def runsummaryplots_ionvolcomp():
    simsets = ['m12_clean2', 'm12_all2', 'm13_clean2', 'm13_all2']
    rranges_rvir = [(0.15, 0.25), (0.45, 0.55), (0.9, 1.0)]
    outdir = '/projects/b1026/nastasha/imgs/summary_plots/'
    for simset in simsets:
        for rrange_rvir in rranges_rvir:
            outname = (f's1_vol_vs_ions_TnHZvr_{simset}'
                       f'_{rrange_rvir[0]:.2f}_to_{rrange_rvir[1]:.2f}'
                       '_rvir')
            outname = outdir + outname.replace('.', 'p') + '.pdf'
            summaryplot_ionvolcomp(simset=simset, rrange_rvir=rrange_rvir,
                                   outname=outname)

def summaryplot_inoutflowcomp(simset='m12_clean2', 
                              rrange_rvir=(0.45, 0.55),
                              outname=None):
    vranges_y = [(0.3, np.inf), (-np.inf, -0.3)]
    vranges_y_units = ['vesc', 'vesc']
    vrlabels = ['$v_{\\mathrm{rad}} > 0.3 v_{\\mathrm{esc}}$',
                '$v_{\\mathrm{rad}} < -0.3 v_{\\mathrm{esc}}$'
                ]
    haloweight = 'gasvol'
    fqtys = ['temperature', 'hdens', 'OxygenAbundance']
    targetvals = ['T', 'nH', 'O']
    showperc = (0.1, 0.5, 0.9)

    buglist = sl.buglist1
    bugics = {simname.split('_')[0] for simname in buglist}
    bugics = list(bugics)
    if simset == 'm12_clean2':
        simnames = sl.m12_hr_clean2 + sl.m12_sr_clean2
        toscan = simnames.copy()
        for simname in toscan:
            if simname.split('_')[0] in bugics:
                simnames.remove(simname)
    elif simset == 'm12_all2':
        simnames = sl.m12_hr_all2 + sl.m12_sr_all2
        toscan = simnames.copy()
        for simname in toscan:
            if simname.split('_')[0] in buglist:
                simnames.remove(simname)
    elif simset == 'm13_clean2':
        simnames = sl.m13_hr_clean2 + sl.m13_sr_clean2
        toscan = simnames.copy()
        for simname in toscan:
            if simname.split('_')[0] in bugics:
                simnames.remove(simname)
    elif simset == 'm13_all2':
        simnames = sl.m13_hr_all2 + sl.m13_sr_all2
        toscan = simnames.copy()
        for simname in toscan:
            if simname.split('_')[0] in buglist:
                simnames.remove(simname)

    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    sims_sr = sl.m13_sr_all2 + sl.m12_sr_all2
    sims_hr = sl.m13_hr_all2 + sl.m12_hr_all2
    ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
    filen_temp = ('hist_rcen_vcen_{compqty}_by_{weight}_{simname}'
                  '_snap{snapnum}_bins1_v1_hvcen.hdf5')

    xlabels = [('$\\langle\\log_{10} \\, \\mathrm{T}'
                '\\rangle_{\\mathrm{V}, \\mathrm{med.}}'
                '\\; [\\mathrm{K}]$'),
               ('$\\langle\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
                '\\rangle_{\\mathrm{V}, \\mathrm{med.}}'
                '\\; [\\mathrm{cm}^{-3}]$'),
               ('$\\langle\\log_{10} \\, \\mathrm{Z}_{\\mathrm{O}}'
                '\\rangle_{\\mathrm{V}, \\mathrm{med.}}$')]
    ylabels = [[xlabel.replace('V', 'V, out') for xlabel in xlabels],
               [xlabel.replace('V', 'V, in') for xlabel in xlabels]]
    fontsize = 12
    
    snaplists = [snaps_sr if simname in sims_sr
                 else snaps_hr if simname in sims_hr
                 else None
                 for simname in simnames]
    datay = {simname: [{(vrange, targetval): 
                        getthemperc_rbins(ddir + filen_temp.format(
                        weight='gasvol', simname=simname, snapnum=snap,
                        compqty=compqty), 
                                          rrange_rvir=rrange_rvir,
                                          vrrange=vrange, 
                                          vrrange_units=vunits,
                                          perc=np.array(showperc))
              for (vrange, vunits) in zip(vranges_y, vranges_y_units)
              for compqty, targetval in zip(fqtys, targetvals)}
             for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)}
    datax = {simname: [{('vall', targetval): 
                        getthemperc_rbins(ddir + filen_temp.format(
                        weight=haloweight, simname=simname, snapnum=snap,
                        compqty=compqty), 
                                          rrange_rvir=rrange_rvir,
                                          vrrange=None, 
                                          vrrange_units=None,
                                          perc=np.array(showperc))
              for compqty, targetval in zip(fqtys, targetvals)}
             for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)}
    yweightmap = [(vrange, targetval) 
                  for vrange in vranges_y for targetval in targetvals ]
    xweightmap = [('vall', targetval) 
                  for _ in vranges_y for targetval in targetvals]
    fig, axes, lax, axdoc = plotdata_censcatter(datax, datay, xweightmap,
                                                yweightmap, xlabel='',
                                                ylabel='',
                                                ncols=len(fqtys),
                                                fontsize=fontsize,
                                                syncaxlims=False,
                                                wspace=0.25)
    fontsize = axdoc['fontsize']
    nrows = len(vranges_y)
    ncols = len(targetvals)
    for ci in range(ncols):
        axsel = slice(ci, None, ncols)
        xlims = [ax.get_xlim() for ax in axes[axsel]]
        xmin = min([xlim[0] for xlim in xlims])
        xmax = max([xlim[1] for xlim in xlims])
        [ax.set_xlim(*(xmin, xmax)) for ax in axes[axsel]]
        ylims = [ax.get_ylim() for ax in axes[axsel]]
        ymin = min([ylim[0] for ylim in ylims])
        ymax = max([ylim[1] for ylim in ylims])
        [ax.set_ylim(*(ymin, ymax)) for ax in axes[axsel]]
        eqp = (max(xmin, ymin), min(xmax, ymax))
        for ri in range(nrows):
            ax = axes[ri * ncols + ci]
            ax.tick_params(labelleft=True)
            ax.plot(eqp, eqp, color='black', linestyle='dotted', 
                    linewidth=1, zorder=-1)
            if ri == nrows - 1:
                ax.set_xlabel(xlabels[ci], fontsize=fontsize)
                ax.tick_params(labelbottom=True)
            ax.set_ylabel(ylabels[ri][ci], fontsize=fontsize)
            if ci == 0:
                ax.text(0.05, 0.95, vrlabels[ri],
                        fontsize=fontsize, transform=ax.transAxes,
                        horizontalalignment='left',
                        verticalalignment='top')
                
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def runsummaryplots_inoutflowcomp():
    simsets = ['m12_clean2', 'm12_all2', 'm13_clean2', 'm13_all2']
    rranges_rvir = [(0.15, 0.25), (0.45, 0.55), (0.9, 1.0)]
    outdir = '/projects/b1026/nastasha/imgs/summary_plots/'
    for simset in simsets:
        for rrange_rvir in rranges_rvir:
            outname = (f's1_inflow_vs_outflow_TnHZ_{simset}'
                       f'_{rrange_rvir[0]:.2f}_to_{rrange_rvir[1]:.2f}'
                       '_rvir')
            outname = outdir + outname.replace('.', 'p') + '.pdf'
            summaryplot_inoutflowcomp(simset=simset, rrange_rvir=rrange_rvir,
                                      outname=outname)