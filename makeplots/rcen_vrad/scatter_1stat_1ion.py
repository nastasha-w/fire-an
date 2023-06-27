import matplotlib.pyplot as plt
import numpy as np

import fire_an.ionrad.get_cieranges as gcr
import fire_an.makeplots.rcen_vrad.scatter_singlestat as s3
import fire_an.simlists as sl


def plot_props_ion_vs_vol(ion, rranges_rvir_col, simnames,
                          outname=None, ionlabel=None,
                          coltitles=None):
    haloweight = 'gasvol'
    fqtys = ['temperature', 'hdens', 'NeonAbundance']
    targetvals = ['T', 'nH', 'N']
    showperc = (0.1, 0.5, 0.9)
    errbperc = (0.1, 0.9)
    if ionlabel is None:
        ionlabel = ion

    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    snaps_f2md = sl.snaps_f2md
    sims_sr = sl.m13_sr_all2 + sl.m12_sr_all2
    sims_hr = sl.m13_hr_all2 + sl.m12_hr_all2
    sims_f2md  = sl.snaps_f2md
    ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
    filen_temp = ('hist_rcen_vcen_{compqty}_by_{weight}_{simname}'
                  '_snap{snapnum}_bins1_v1_hvcen.hdf5')

    ylabels = [('$\\langle\\log_{10} \\, \\mathrm{T}'
                f'\\rangle_{{\\mathrm{{{ionlabel}}}, \\mathrm{{med.}}}}'
                '\\; [\\mathrm{K}]$'),
               ('$\\langle\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
                f'\\rangle_{{\\mathrm{{{ionlabel}}}, \\mathrm{{med.}}}}'
                '\\; [\\mathrm{cm}^{-3}]$'),
               ('$\\langle\\log_{10} \\, \\mathrm{Z}_{\\mathrm{N}}'
                f'\\rangle_{{\\mathrm{{{ionlabel}}}, \\mathrm{{med.}}}}$'),
                ]
    xlabels = [ylabel.replace(ionlabel, 'V') for ylabel in ylabels]
    fontsize = 12
    
    snaplists = [snaps_sr if simname in sims_sr
                 else snaps_hr if simname in sims_hr
                 else snaps_f2md if simname in sims_f2md
                 else None
                 for simname in simnames]
    datay = {simname: [{(rrange_rvir, targetval): 
                        s3.getthemperc_rbins(ddir + filen_temp.format(
                        weight=ion, simname=simname, snapnum=snap,
                        compqty=compqty), 
                                          rrange_rvir=rrange_rvir,
                                          vrrange=None, 
                                          vrrange_units=None,
                                          perc=np.array(showperc))
                        if targetval != 'vr' else
                        s3.getvperc_rbins(ddir + filen_temp.format(
                        weight=ion, simname=simname, snapnum=snap,
                        compqty=compqty), 
                                       rrange_rvir=rrange_rvir,
                                       vperc=errbperc)
              for compqty, targetval in zip(fqtys, targetvals)
              for rrange_rvir in rranges_rvir_col}
             for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)}
    datax = {simname: [{(rrange_rvir, targetval): 
                        s3.getthemperc_rbins(ddir + filen_temp.format(
                        weight=haloweight, simname=simname, snapnum=snap,
                        compqty=compqty), 
                                          rrange_rvir=rrange_rvir,
                                          vrrange=None, 
                                          vrrange_units=None,
                                          perc=np.array(showperc)) 
                        if targetval != 'vr' else
                        s3.getvperc_rbins(ddir + filen_temp.format(
                        weight=haloweight, simname=simname, snapnum=snap,
                        compqty=compqty), 
                                       rrange_rvir=rrange_rvir,
                                       vperc=errbperc)
              for compqty, targetval in zip(fqtys, targetvals)
              for rrange_rvir in rranges_rvir_col}
             for snap in snaplist]
            for snaplist, simname in zip(snaplists, simnames)}
    yweightmap = [(rrange_rvir, targetval) 
                  for targetval in targetvals
                  for rrange_rvir in rranges_rvir_col]
    xweightmap = [(rrange_rvir, targetval) 
                  for targetval in targetvals
                  for rrange_rvir in rranges_rvir_col]
    fig, axes, lax, axdoc = s3.plotdata_censcatter(datax, datay, xweightmap,
                                                   yweightmap, xlabel='',
                                                   ylabel='',
                                                   ncols=len(fqtys),
                                                   fontsize=fontsize,
                                                   syncaxlims=False,
                                                   hspace=0.28)
    fontsize = axdoc['fontsize']
    nrows = len(fqtys)
    ncols = len(rranges_rvir_col)
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
        for ci in range(len(rranges_rvir_col)):
            ax = axes[ri * ncols + ci]
            ax.plot(eqp, eqp, color='black', linestyle='dotted', 
                    linewidth=1, zorder=-1)
            ax.set_xlabel(xlabels[ri], fontsize=fontsize)
            ax.tick_params(labelbottom=True)
            if ci == 0:
                ax.set_ylabel(ylabels[ri], fontsize=fontsize)
            if ri == 0:
                if coltitles is not None:
                    ax.set_title(coltitles[ci], fontsize=fontsize)
            if ri == 0:
                rng = gcr.cieranges1[ion]
                ax.axhline(rng[0], color='black', linestyle='dotted',
                           linewidth=1.)
                ax.axhline(rng[1], color='black', linestyle='dotted',
                           linewidth=1.)
            if ri == 2 and ion == 'Ne8':
                logsolarmassfrac = -2.9008431 #copied from PS20 tables
                ax.axhline(logsolarmassfrac, 
                           color='black', linestyle='dotted',
                           linewidth=1.)
                ax.axvline(logsolarmassfrac, 
                           color='black', linestyle='dotted',
                           linewidth=1.)
            
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def runplotprop_ion_vs_vol(ion='Ne8', ionlabel='Ne VIII'):
    simsets = ['m12_clean2', 'm12_all2', 'm13_clean2', 'm13_all2']
    rranges_rvir = [(0.15, 0.25), (0.45, 0.55), (0.9, 1.0)]
    outdir = '/projects/b1026/nastasha/imgs/summary_plots/'
    coltitles = [f'${rrange[0]:.2}\\endash{rrange[1]:.2f}'
                 '\\,\\mathrm{R}_{\\mathrm{vir}}$'
                 for rrange in rranges_rvir]
    
    buglist = sl.buglist1
    bugics = {simname.split('_')[0] for simname in buglist}
    bugics = list(bugics)
    for simset in simsets:
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
        outname = (f's2_vol_vs_{ion}_TnHZ_{simset}'
                    'rrange_var')
        outname = outdir + outname.replace('.', 'p') + '.pdf'
        plot_props_ion_vs_vol(ion, rranges_rvir, simnames,
                                outname=outname, ionlabel=ionlabel,
                                coltitles=coltitles)
