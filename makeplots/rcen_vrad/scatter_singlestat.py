import h5py 
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

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
        vescvir_kmps = np.sqrt(c.gravity * mvir_g / rvir_cm) \
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
           'mean': mean,
           'perc': perc}
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
        vescvir_kmps = np.sqrt(c.gravity * mvir_g / rvir_cm) \
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
           'perc': perc}
    return out


def plot_vmean_ion_vs_halo(simnames, haloweight='Volume',
                           errbperc=(0.1, 0.9),
                           rrange_rvir=(0.1, 1.),
                           title=None, outname=None):
    ions = ['O6', 'Ne8', 'O7', 'Ne9', 'O8', 'Ne10']
    ionlabels = ['O VI', 'Ne VIII', 'O VII', 'Ne IX', 'O VIII', 'Ne X']
    physmarkers = {'noBH': '*',
                   'AGN-noCR': 'o',
                   'AGN-CR': 'P'}
    markeredgewidth = 2
    size = 5
    fontsize = 12
    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    sims_sr = sl.m13_sr_all2 + sl.m12_sr_all2
    sims_hr = sl.m13_hr_all2 + sl.m12_hr_all2
    ddir = '/projects/b1026/nastasha/hists/r_vr_clean2_nobug/'
    filen_temp = ('hist_rcen_vcen_temperature_by_{weight}_{simname}'
                 '_snap{snapnum}_bins1_v1_hvcen.hdf5')
    if haloweight == 'Volume':
        fhaloweight = 'gasvol'
    elif haloweight == 'Mass':
        fhaloweight = 'gasmass'
    iccolors = sl.m12_iccolors 
    iccolors.update(sl.m13_iccolors)

    xlabel = ('$\\left\\langlev_{\\mathrm{rad}} \\right\\rangle'
              f'_{{\\mathrm{{{haloweight}}}}}'
              '\\; [\\mathrm{km} \\, \\mathrm{s}^{-1}]$')
    ylabel = ('$\\left\\langlev_{\\mathrm{rad}} \\right\\rangle'
              '_{\\mathrm{ion}}'
              '\\; [\\mathrm{km} \\, \\mathrm{s}^{-1}]$')

    npanels = len(ions)
    ncmax = 4
    ncols = min(ncmax, npanels)
    nrows = (npanels - 1) // ncmax + 1
    panelsize = 2.5
    laxspace = 1.5
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] * nrows + [laxspace]
    width = sum(width_ratios)
    height = sum(height_ratios)

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows + 1, hspace=0.0,
                        wspace=0.0, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols])
            for i in range(npanels)]
    lax = fig.add_subplot(grid[nrows, :])
    
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    
    allweights = ions + [fhaloweight]
    snaplists = [snaps_sr if simname in sims_sr
                 else snaps_hr if simname in sims_hr
                 else None
                 for simname in simnames]
    data = [[{weight: getvperc_rbins(ddir + filen_temp.format(
                    weight=weight, simname=simname, snapnum=snap),
                             rrange_rvir=rrange_rvir, vperc=errbperc)
              for weight in allweights}
             for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)]
    ics = list({simname.split('_')[0] for simname in simnames})
    ics.sort()
    zs = {np.round(wdict[key]['cosmopars']['z'], 1) 
          for l1 in data for wdict in l1 for key in wdict}
    zs = list(zs)
    zs.sort()
    zcolors = pu.getzcolorfunc(zs)

    for axi, ax in enumerate(axes):
        bottom = axi >= npanels - ncols
        left = axi % ncols == 0
        ax.tick_params(which='both', direction='in', labelsize=fontsize - 1,
                       right=True, top=True, labelbottom=bottom,
                       labelleft=left)
        if bottom:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        if left:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.text(0.05, 0.95, ionlabels[axi], fontsize=fontsize,
                horizontalalignment='left', verticalalignment='top',
                transform=ax.transAxes)
    lax.axis('off')
    ltop = 0.60
    handles1 = [mlines.Line2D((), (), label=ic, markersize=size,
                              markerfacecolor='gray', marker='o',
                              markeredgecolor=iccolors[ic],
                              linestyle='none',
                              markeredgewidth=markeredgewidth)
                for ic in ics]
    leg1 = lax.legend(handles=handles1, fontsize=fontsize,
                      bbox_to_anchor=(0.0, ltop), loc='upper left',
                      handletextpad=0.4, columnspacing=1.0)
    handles2 = [mlines.Line2D((), (), label=f'z={z:.1f}', 
                              markersize=size,
                              markerfacecolor=zcolors(z), marker='o',
                              markeredgecolor='gray',
                              linestyle='none',
                              markeredgewidth=markeredgewidth)
                 for z in zs]
    leg2 = lax.legend(handles=handles2, fontsize=fontsize, ncol=3,
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
    leg3 = lax.legend(handles=handles3, fontsize=fontsize, ncol=2,
                      bbox_to_anchor=(1.0, ltop), loc='upper right',
                      handletextpad=0.4, columnspacing=1.0)
    lax.add_artist(leg1)
    lax.add_artist(leg2)
    lax.add_artist(leg3)

    xmin = np.inf
    xmax = -np.inf
    ymin = np.inf
    ymax = -np.inf
    for datazw, simname in zip(data, simnames):
        ic = simname.split('_')[0]
        phys = getphys(simname)
        marker = physmarkers[phys]
        edgecolor = iccolors[ic]
        for dataw in datazw:
            for ioni, ion in enumerate(ions):
                ax = axes[ioni]
                datumy = dataw[ion]
                datumx = dataw[fhaloweight]
                facecolor = zcolors(datumy['cosmopars']['z'])
                xmid = datumx['mean']
                xlo, xhi = datumx['perc']
                ymid = datumy['mean']
                ylo, yhi = datumy['perc']
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
    #vmar = 100.
    xran = xmax - xmin
    xmin = xmin - xran
    xmax = xmax + xran
    #yran = ymax - ymin
    ymin = ymin - xran
    ymax = ymax + xran
    #xlims = [ax.get_xlim() for ax in axes]
    #xmin = min([xlim[0] for xlim in xlims])
    #xmax = min([xlim[-1] for xlim in xlims])
    [ax.set_xlim((xmin, xmax)) for ax in axes]
    #ylims = [ax.get_ylim() for ax in axes]
    #ymin = min([ylim[0] for ylim in ylims])
    #ymax = min([ylim[-1] for ylim in ylims])
    [ax.set_ylim((ymin, ymax)) for ax in axes]
    eqp = (max(xmin, ymin), min(xmax, ymax))
    [ax.plot(eqp, eqp, color='black', linestyle='dotted', linewidth=1,
             zorder=-1) for ax in axes]

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotset_vmean_ion_vs_halo(haloweight='Volume'):
    errbperc = (0.1, 0.9)
    rranges_rvir = [(0.1, 1.), (0.15, 0.25), (0.45, 0.5),
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
            _outname = (f'vbias_ions_{haloweight}'
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
            
def plot_thermcomp_medscat(simnames,
                           compqty='T',
                           haloweight='Volume',
                           rrange_rvir=(0.1, 1.),
                           vrrange=None, vrrange_units='kmps',
                           title=None, outname=None):
    ions = ['O6', 'Ne8', 'O7', 'Ne9', 'O8', 'Ne10']
    ionlabels = ['O VI', 'Ne VIII', 'O VII', 'Ne IX', 'O VIII', 'Ne X']
    physmarkers = {'noBH': '*',
                   'AGN-noCR': 'o',
                   'AGN-CR': 'P'}
    markeredgewidth = 2
    size = 5
    fontsize = 12
    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    sims_sr = sl.m13_sr_all2 + sl.m12_sr_all2
    sims_hr = sl.m13_hr_all2 + sl.m12_hr_all2
    ddir = '/projects/b1026/nastasha/hists/r_vr_clean2_nobug/'
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
    
    iccolors = sl.m12_iccolors 
    iccolors.update(sl.m13_iccolors)

    xlabel = (f'$\\left\\langle {flab} \\right\\rangle'
              f'_{{\\mathrm{{{haloweight}}}}}'
              f'\\; [{funits}]$')
    ylabel = (f'$\\left\\langle {flab} \\right\\rangle'
              '_{\\mathrm{ion}}'
              f'\\; [{funits}]$')
    showperc = (0.1, 0.5, 0.9)
    npanels = len(ions)
    ncmax = 4
    ncols = min(ncmax, npanels)
    nrows = (npanels - 1) // ncmax + 1
    panelsize = 2.5
    laxspace = 1.5
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] * nrows + [laxspace]
    width = sum(width_ratios)
    height = sum(height_ratios)

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows + 1, hspace=0.0,
                        wspace=0.0, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols])
            for i in range(npanels)]
    lax = fig.add_subplot(grid[nrows, :])
    
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    
    allweights = ions + [fhaloweight]
    snaplists = [snaps_sr if simname in sims_sr
                 else snaps_hr if simname in sims_hr
                 else None
                 for simname in simnames]
    data = [[{weight: getthemperc_rbins(ddir + filen_temp.format(
                    weight=weight, simname=simname, snapnum=snap,
                    compqty=fqty), 
                                        rrange_rvir=rrange_rvir,
                                        vrrange=vrrange, 
                                        vrrange_units=vrrange_units,
                                        perc=np.array(showperc))
              for weight in allweights}
             for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)]
    ics = list({simname.split('_')[0] for simname in simnames})
    ics.sort()
    zs = {np.round(wdict[key]['cosmopars']['z'], 1) 
          for l1 in data for wdict in l1 for key in wdict}
    zs = list(zs)
    zs.sort()
    zcolors = pu.getzcolorfunc(zs)

    for axi, ax in enumerate(axes):
        bottom = axi >= npanels - ncols
        left = axi % ncols == 0
        ax.tick_params(which='both', direction='in', labelsize=fontsize - 1,
                       right=True, top=True, labelbottom=bottom,
                       labelleft=left)
        if bottom:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        if left:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.text(0.05, 0.95, ionlabels[axi], fontsize=fontsize,
                horizontalalignment='left', verticalalignment='top',
                transform=ax.transAxes)
    lax.axis('off')
    ltop = 0.60
    handles1 = [mlines.Line2D((), (), label=ic, markersize=size,
                              markerfacecolor='gray', marker='o',
                              markeredgecolor=iccolors[ic],
                              linestyle='none',
                              markeredgewidth=markeredgewidth)
                for ic in ics]
    leg1 = lax.legend(handles=handles1, fontsize=fontsize,
                      bbox_to_anchor=(0.0, ltop), loc='upper left',
                      handletextpad=0.4, columnspacing=1.0)
    handles2 = [mlines.Line2D((), (), label=f'z={z:.1f}', 
                              markersize=size,
                              markerfacecolor=zcolors(z), marker='o',
                              markeredgecolor='gray',
                              linestyle='none',
                              markeredgewidth=markeredgewidth)
                 for z in zs]
    leg2 = lax.legend(handles=handles2, fontsize=fontsize, ncol=3,
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
    leg3 = lax.legend(handles=handles3, fontsize=fontsize, ncol=2,
                      bbox_to_anchor=(1.0, ltop), loc='upper right',
                      handletextpad=0.4, columnspacing=1.0)
    lax.add_artist(leg1)
    lax.add_artist(leg2)
    lax.add_artist(leg3)

    xmin = np.inf
    xmax = -np.inf
    ymin = np.inf
    ymax = -np.inf
    for datazw, simname in zip(data, simnames):
        ic = simname.split('_')[0]
        phys = getphys(simname)
        marker = physmarkers[phys]
        edgecolor = iccolors[ic]
        for dataw in datazw:
            for ioni, ion in enumerate(ions):
                ax = axes[ioni]
                datumy = dataw[ion]
                datumx = dataw[fhaloweight]
                facecolor = zcolors(datumy['cosmopars']['z'])
                xlo, xmid, xhi = datumx['perc']
                ylo, ymid, yhi = datumy['perc']
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
    #vmar = 100.
    xran = xmax - xmin
    xmin = xmin - xran
    xmax = xmax + xran
    #yran = ymax - ymin
    ymin = ymin - xran
    ymax = ymax + xran
    #xlims = [ax.get_xlim() for ax in axes]
    #xmin = min([xlim[0] for xlim in xlims])
    #xmax = min([xlim[-1] for xlim in xlims])
    [ax.set_xlim((xmin, xmax)) for ax in axes]
    #ylims = [ax.get_ylim() for ax in axes]
    #ymin = min([ylim[0] for ylim in ylims])
    #ymax = min([ylim[-1] for ylim in ylims])
    [ax.set_ylim((ymin, ymax)) for ax in axes]
    eqp = (max(xmin, ymin), min(xmax, ymax))
    [ax.plot(eqp, eqp, color='black', linestyle='dotted', linewidth=1,
             zorder=-1) for ax in axes]

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')
        
    

    
    
    