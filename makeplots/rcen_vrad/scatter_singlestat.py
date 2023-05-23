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
        vedges = f['axis_1/bins'][:]
        vcens = 0.5 * (vedges[:-1] + vedges[1:])
        mean = np.average(vcens, weights=_hist)
        perc = mu.percentiles_from_histogram(_hist, vedges, 
                                                 percentiles=np.array(vperc))
    out = {'cosmopars': cosmopars,
           'vescvir_kmps': vescvir_kmps,
           'mvir_g': mvir_g,
           'rvir_cm': rvir_cm,
           'mean': mean,
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
    size = 20
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
    nrows = npanels // ncmax
    panelsize = 2.5
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] * nrows
    width = sum(width_ratios)
    height = sum(height_ratios)

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows, hspace=0.0,
                        wspace=0.0, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols])
            for i in range(npanels)]
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    
    allweights = ions + [fhaloweight]
    snaplists = [snaps_sr if simname in sims_sr
                 else snaps_hr if simname in sims_hr
                 else None
                 for simname in simnames]
    data = [[[getvperc_rbins(ddir + filen_temp.format(weight=weight, 
                                                      simname=simname,
                                                      snapnum=snap),
                             rrange_rvir=rrange_rvir, vperc=errbperc)
              for weight in allweights]
             for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)]
    ics = list({simname.split('_')[0] for simname in simnames})
    ics.sort()
    zs = {np.round(datum['cosmopars']['z'], 1) 
          for l1 in data for l2 in l1 for datum in l2}
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
        if axi == 0:
            handles = [mlines.Line2D((), (), label=ic, markersize=size,
                                     color='gray', markerstyle='o',
                                     markeredgecolor=iccolors[ic])
                       for ic in ics]
            ax.legend(handles, fontsize=fontsize)
        elif axi == 1:
            handles = [mlines.Line2D((), (), label=f'z={z:.1f}', 
                                     markersize=size,
                                     color=zcolors(z), markerstyle='o',
                                     markeredgecolor='gray')
                       for z in zs]
            ax.legend(handles, fontsize=fontsize)
        elif axi == 2:
            handles = [mlines.Line2D((), (), label=phys, 
                                     markersize=size,
                                     color='gray', 
                                     markerstyle=physmarkers[phys],
                                     markeredgecolor='black')
                       for phys in physmarkers]
    
    if outname is not None:
        plt.savefig(outname=outname)

        
    

    
    
    