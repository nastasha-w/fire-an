
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fire_an.analytic_halo.model_ionprof_pl as mip
import fire_an.makeplots.litcomp.obsdataread as odr
import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
import fire_an.mstar_mhalo.analytical as msmhan
import fire_an.mstar_mhalo.loader_smdpl_sfr as ldsmdpl
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu
import fire_an.utils.math_utils as mu


outdir = '/projects/b1026/nastasha/imgs/analytical/'

def plot_plmodel_datacomp():
    ion = 'Ne8'
    redshift_model = 0.75
    nsigma = 1.

    impactpars_kpc = np.linspace(5., 450., 50)
    logmvirs_msun = np.arange(10.7, 13.3, 0.3)
    fcgm = 0.9
    z_sol = 0.3
    redshift = 0.75
    plis = [0.20, 0.0, -0.20]

    data_bur = odr.readdata_b19(nsigma=nsigma)

    vmin_an = logmvirs_msun[0]
    vmax_an = logmvirs_msun[-1]
    vmin_obs = np.min(data_bur['logmvir_msun_lo'])
    vmax_obs = np.max(data_bur['logmvir_msun_hi'])
    vmin = min(vmin_an, vmin_obs)
    vmax = max(vmax_an, vmax_obs) 
    cmap = mcm.get_cmap('viridis')
    dc = np.average(np.diff(logmvirs_msun))
    bounds_disc = np.arange(vmin_an - 0.5 * dc, vmax_an + dc, dc)
    cmap_disc = pu.truncate_colormap(cmap, 
                                     (vmin_an - vmin) / (vmax - vmin),
                                     (vmax_an - vmin) / (vmax - vmin))
    norm_disc = mcolors.BoundaryNorm(bounds_disc, cmap.N)
    norm_cont = mcolors.Normalize(vmin=vmin, vmax=vmax)

    fig = plt.figure(figsize=(11.5, 3.))
    grid = gsp.GridSpec(ncols=6, nrows=1, wspace=0.0, 
                        width_ratios=(10., 10., 10., 1., 2., 1.))
    axes = [fig.add_subplot(grid[i]) for i in range(3)]
    cax = fig.add_subplot(grid[3])
    cax2 = fig.add_subplot(grid[5])
    fontsize = 12
    linewidth = 1.5
    path_effects = pu.getoutline(linewidth)
    
    for plii, pli in enumerate(plis):
        ax = axes[plii]
        doleft = plii == 0
        ax.set_xlabel('$\\mathrm{r}_{\\perp} \\; [\\mathrm{pkpc}]$',
                  fontsize=fontsize)
        if doleft:
            ax.set_ylabel('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\,VIII})'
                          '\\; [\\mathrm{cm}^{-2}]$',
                          fontsize=fontsize)
        ax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                       top=True, right=True, labelleft=doleft)
        ax.text(0.95, 0.95, f'$v_{{\\mathrm{{c}}}} \\propto r^{{{pli:.2f}}}$',
                transform=ax.transAxes, fontsize=fontsize,
                verticalalignment='top', horizontalalignment='right')

        for lmv in logmvirs_msun:
            hmod = mip.PLmodel(10**lmv, redshift, fcgm, z_sol, pli)
            cvs = hmod.coldensprof(ion, impactpars_kpc)
            ax.plot(impactpars_kpc, np.log10(cvs), 
                    color=cmap((lmv - vmin) / (vmax - vmin)),
                    linewidth=linewidth, path_effects=path_effects)
        
        for dbi in range(len(data_bur)):
            xv = data_bur['impact_parameter_kpc'][dbi]
            yv = data_bur['log_N_Ne8_pcm2'][dbi]
            isul = data_bur['log_N_Ne8_isUL'][dbi]
            yerr = data_bur['log_N_Ne8_pcm2_err'][dbi] if not isul else None
            cbest = data_bur['logmvir_msun_bestest'][dbi]
            clo = data_bur['logmvir_msun_lo'][dbi]
            chi = data_bur['logmvir_msun_hi'][dbi]

            marker = 'v' if isul else 'o'
            markersize = 10
            zobase = 5. - 1. * isul

            ax.errorbar([xv], [yv], yerr=yerr, 
                        linestyle='solid', elinewidth=1.5,
                        color='black', capsize=3,
                        zorder=zobase, fmt='none')
            ax.plot([xv], [yv], linestyle='None', marker=marker, 
                    markersize=markersize, 
                    markerfacecolor=cmap((clo - vmin) / (vmax - vmin)), 
                    markeredgecolor='black', 
                    zorder=zobase + 0.1, fillstyle='left')
            ax.plot([xv], [yv], linestyle='None', marker=marker, 
                    markersize=markersize, 
                    markerfacecolor=cmap((chi - vmin) / (vmax - vmin)), 
                    markeredgecolor='black', 
                    zorder=zobase + 0.1, fillstyle='right')
            ax.plot([xv], [yv], linestyle='None', marker=marker, 
                    markersize=markersize, 
                    markerfacecolor=cmap((cbest - vmin) / (vmax - vmin)), 
                    markeredgecolor='black', 
                    zorder=zobase + 0.2, fillstyle='bottom')

    ylims = [ax.get_ylim() for ax in axes]
    ymin = min([ylim[0] for ylim in ylims])
    ymin = max(ymin, 12.5)
    ymax = max([ylim[1] for ylim in ylims])
    [ax.set_ylim((ymin, ymax)) for ax in axes]

    scm1 = mcm.ScalarMappable(norm=norm_disc, cmap=cmap_disc)
    scm1.set_array(logmvirs_msun)
    cb1 = plt.colorbar(scm1, cax=cax,  orientation='vertical')
    scm2 = mcm.ScalarMappable(norm=norm_cont, cmap=cmap)
    scm2.set_array(logmvirs_msun)
    plt.colorbar(scm2, cax=cax2,  orientation='vertical')
    cax2.set_ylabel('$\\log_{10} \\, \\mathrm{M}_{\\mathrm{vir}}'
                   '\\; [\\mathrm{M}_{\\odot}]$',
                   fontsize=fontsize)
    cax2.tick_params(labelsize=fontsize - 1.)
    cax.tick_params(labelsize=fontsize - 1.)
    cb1.set_ticks(logmvirs_msun)

    outname = (f'prof_Ne8_analytical_pl_s19_z{redshift_model:.2f}'
               f'_vs_b19_loghalomass_{nsigma:.1f}sigma')
    outname = outname.replace('.', 'p')
    plt.savefig(outdir + outname + '.pdf', bbox_inches='tight')


def plot_plmodel_datacomp_Kvar():
    ion = 'Ne8'
    redshift_model = 0.75
    nsigmas = (1, 2)

    impactpars_kpc = np.linspace(5., 450., 50)
    logmvirs_msun = np.arange(11.3, 13.5, 0.3)
    fcgm = 1.0
    z_sol = 0.3
    redshift = 0.75
    plis_vc = [-0.1] #[0.0, -0.20]
    colors_vc = ['black'] #['black', 'blue']
    plis_k = [0.0, 2./3., 1.2]
    linestyles_k = ['dotted', 'dashed', 'solid']

    data_bur = odr.readdata_b19(nsigmas=nsigmas)
    
    panelsize = 2.5
    ncol_max = 4
    npanels = len(logmvirs_msun)
    ncols = min(npanels, ncol_max)
    nrows = (npanels - 1) // ncols + 1
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] * nrows
 
    fig = plt.figure(figsize=(sum(width_ratios), sum(height_ratios)))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows, wspace=0.0, 
                        hspace=0.0, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols]) 
            for i in range(npanels)]
    fontsize = 12
    linewidth = 1.5
    
    for mi, mvir in enumerate(logmvirs_msun):
        ax = axes[mi]
        doleft = mi % ncols == 0
        dobottom = npanels - mi <= ncols
        if dobottom:
            ax.set_xlabel('$\\mathrm{r}_{\\perp} \\; [\\mathrm{pkpc}]$',
                          fontsize=fontsize)
        if doleft:
            ax.set_ylabel('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\,VIII})'
                          '\\; [\\mathrm{cm}^{-2}]$',
                          fontsize=fontsize)
        ax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                       top=True, right=True, labelleft=doleft,
                       labelbottom=dobottom)
        axlabel = (f'$\\mathrm{{M}}_{{\\mathrm{{vir}}}} = 10^{{{mvir:.1f}}}'
                   '\\, \\mathrm{M}_{\\odot}$')
        ax.text(0.95, 0.95, axlabel,
                transform=ax.transAxes, fontsize=fontsize,
                verticalalignment='top', horizontalalignment='right')
        yvir_max = - np.inf
        for ls, pli_k in zip(linestyles_k, plis_k):
            for color, pli_vc in zip(colors_vc, plis_vc):
                print(10**mvir, redshift, fcgm, z_sol, pli_vc, pli_k)
                hmod = mip.PLmodel(10**mvir, redshift, fcgm, z_sol, pli_vc,
                                   pli_entropy=pli_k)
                cvs = hmod.coldensprof(ion, impactpars_kpc)
                ax.plot(impactpars_kpc, np.log10(cvs), 
                        color=color, linestyle=ls,
                        linewidth=linewidth)
                rvir = hmod.rvir_cgs / (c.cm_per_mpc * 1e-3)
                if rvir > impactpars_kpc[-1]:
                    continue
                yv = mu.linterpsolve(impactpars_kpc, np.log10(cvs), rvir)
                yvir_max = max(yvir_max, yv)
        if rvir <= impactpars_kpc[-1]:
            ax.plot([rvir, rvir], [11., yvir_max + 0.2], color='gray',
                    linestyle='solid', alpha=0.3, linewidth=3)
        ulsig0done = False
        ulsig1done = False
        detsig0done = False
        detsig1done = False
        for dbi in range(len(data_bur)):
            cloer = data_bur['logmvir_msun_loer'][dbi]
            chier = data_bur['logmvir_msun_hier'][dbi]
            if cloer > mvir or chier < mvir:
                continue

            xv = data_bur['impact_parameter_kpc'][dbi]
            yv = data_bur['log_N_Ne8_pcm2'][dbi]
            isul = data_bur['log_N_Ne8_isUL'][dbi]
            yerr = data_bur['log_N_Ne8_pcm2_err'][dbi] if not isul else None
            #cbest = data_bur['logmvir_msun_bestest'][dbi]
            clo = data_bur['logmvir_msun_lo'][dbi]
            chi = data_bur['logmvir_msun_hi'][dbi]
            
            issig0 = (clo <= mvir and chi >= mvir)
            _label = None
            if issig0:
                _color = 'black'
                if isul and not ulsig0done:
                    _label = ('UL, $\\Delta\\mathrm{M}'
                              f' < {nsigmas[0]}\\sigma$')
                    ulsig0done = True
                elif not isul and not detsig0done:
                    _label = ('det., $\\Delta\\mathrm{M}'
                              f' < {nsigmas[0]}\\sigma$')
                    detsig0done = True
            else:
                _color = 'gray'
                if isul and not ulsig1done:
                    _label = ('UL, $\\Delta\\mathrm{M}'
                              f' < {nsigmas[1]}\\sigma$')
                    ulsig1done = True
                elif not isul and not detsig1done:
                    _label = ('det., $\\Delta\\mathrm{M}'
                              f' < {nsigmas[1]}\\sigma$')
                    detsig1done = True
            marker = 'v' if isul else 'o'
            markersize = 5
            zobase = 5. - 1. * isul

            ax.errorbar([xv], [yv], yerr=yerr, 
                        linestyle='none', elinewidth=1.5,
                        color=_color, capsize=3,
                        zorder=zobase,
                        marker=marker, markersize=markersize,
                        markeredgecolor='black', markeredgewidth=1.0,
                        label=_label)
        if detsig1done and detsig0done and ulsig1done and ulsig0done:
            getlegax = mi
    ylims = [ax.get_ylim() for ax in axes]
    ymin = min([ylim[0] for ylim in ylims])
    ymin = max(ymin, 12.5)
    ymax = max([ylim[1] for ylim in ylims])
    [ax.set_ylim((ymin, ymax)) for ax in axes]

    handles1, _ = axes[getlegax].get_legend_handles_labels()
    axes[-1].legend(handles=handles1, fontsize=fontsize - 2,
                    loc='upper right', bbox_to_anchor=(1.0, 0.87),
                    handlelength=1.0, labelspacing=0.3,
                    handletextpad=0.4)
    if len(plis_k) > 1:
        handles2 = [mlines.Line2D(
                        (), (), color='black', linewidth=linewidth,
                        linestyle=ls,
                        label=f'$\\mathrm{{K}} \\propto r^{{{pli_k:.2f}}}$')
                    for pli_k, ls in zip(plis_k, linestyles_k)]
        axes[-2].legend(handles=handles2, fontsize=fontsize,
                        loc='upper right', bbox_to_anchor=(1.0, 0.87),
                        handlelength=1.5, labelspacing=0.1,
                        handletextpad=0.4)
    if len(plis_vc) > 1:
        handles3 = [mlines.Line2D(
                        (), (), color=color, linewidth=linewidth,
                        linestyle='solid',
                        label=('$\\mathrm{v}_{\\mathrm{c}} \\propto '
                               f'r^{{{pli_vc:.2f}}}$'))
                    for pli_vc, color in zip(plis_vc, colors_vc)]
        axes[-3].legend(handles=handles3, fontsize=fontsize,
                        loc='upper right', bbox_to_anchor=(1.0, 0.87),
                        handlelength=1.0, labelspacing=0.0,
                        handletextpad=0.4)
    
    pli_vc_str = 'pli_vc_' + \
                 '_'.join([f'{pli_vc:.2f}' for pli_vc in plis_vc])
    pli_k_str = 'pli_k_' + \
                 '_'.join([f'{pli_k:.2f}' for pli_k in plis_k])

    outname = (f'prof_Ne8_analytical_pl_s19_z{redshift_model:.2f}'
               f'_vs_b19_loghalomass_{pli_k_str}_{pli_vc_str}')
    outname = outname.replace('.', 'p')
    outname = outname.replace('-', 'm')
    plt.savefig(outdir + outname + '.pdf', bbox_inches='tight')

def plot_plmodel_datacomp_Kvar_fcgmvar():
    ion = 'Ne8'
    redshift_model = 0.75
    nsigmas = (1, 2)

    impactpars_kpc = np.linspace(5., 450., 50)
    logmvir_msun = 12.2
    fcgms = [1.0, 0.3, 0.1] #, 0.3]
    z_sols = [0.3, 0.3, 0.3] #, 1.0]
    _colors = tc.tol_cset('bright')
    colors_fZ = [_colors.blue, _colors.green, _colors.yellow, _colors.cyan]
    redshift = 0.75
    pli_vc = -0.1 #[0.0, -0.20]
    plis_vc = [pli_vc]
    #colors_vc = ['black'] #['black', 'blue']
    plis_k = [0.0, 2./3., 1.2]
    #linestyles_k = ['dotted', 'dashed', 'solid']

    data_bur = odr.readdata_b19(nsigmas=nsigmas)
    
    panelsize = 2.5
    ncol_max = 4
    npanels = len(plis_k)
    ncols = min(npanels, ncol_max)
    nrows = (npanels - 1) // ncols + 1
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] * nrows
 
    fig = plt.figure(figsize=(sum(width_ratios), sum(height_ratios)))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows, wspace=0.0, 
                        hspace=0.0, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols]) 
            for i in range(npanels)]
    fontsize = 12
    linewidth = 1.5
    
    for plii, pli_k in enumerate(plis_k):
        ax = axes[plii]
        doleft = plii % ncols == 0
        dobottom = npanels - plii <= ncols
        if dobottom:
            ax.set_xlabel('$\\mathrm{r}_{\\perp} \\; [\\mathrm{pkpc}]$',
                          fontsize=fontsize)
        if doleft:
            ax.set_ylabel('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\,VIII})'
                          '\\; [\\mathrm{cm}^{-2}]$',
                          fontsize=fontsize)
        ax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                       top=True, right=True, labelleft=doleft,
                       labelbottom=dobottom)
        axlabel = f'$\\mathrm{{K}} \\propto r^{{{pli_k:.2f}}}$'
        ax.text(0.08, 0.95, axlabel,
                transform=ax.transAxes, fontsize=fontsize,
                verticalalignment='top', horizontalalignment='left')
        yvir_max = - np.inf
        for color, fcgm, z_sol in zip(colors_fZ, fcgms, z_sols):
            hmod = mip.PLmodel(10**logmvir_msun, redshift, fcgm, z_sol, 
                                pli_vc, pli_entropy=pli_k)
            cvs = hmod.coldensprof(ion, impactpars_kpc)

            ax.plot(impactpars_kpc, np.log10(cvs), 
                    color=color, linestyle='solid',
                    linewidth=linewidth)
            rvir = hmod.rvir_cgs / (c.cm_per_mpc * 1e-3)
            if rvir > impactpars_kpc[-1]:
                continue
            yv = mu.linterpsolve(impactpars_kpc, np.log10(cvs), rvir)
            yvir_max = max(yvir_max, yv)
        if rvir <= impactpars_kpc[-1]:
            ax.plot([rvir, rvir], [11., yvir_max + 0.2], color='gray',
                    linestyle='solid', alpha=0.3, linewidth=3)
        ulsig0done = False
        ulsig1done = False
        detsig0done = False
        detsig1done = False
        for dbi in range(len(data_bur)):
            cloer = data_bur['logmvir_msun_loer'][dbi]
            chier = data_bur['logmvir_msun_hier'][dbi]
            if cloer > logmvir_msun or chier < logmvir_msun:
                continue

            xv = data_bur['impact_parameter_kpc'][dbi]
            yv = data_bur['log_N_Ne8_pcm2'][dbi]
            isul = data_bur['log_N_Ne8_isUL'][dbi]
            yerr = data_bur['log_N_Ne8_pcm2_err'][dbi] if not isul else None
            #cbest = data_bur['logmvir_msun_bestest'][dbi]
            clo = data_bur['logmvir_msun_lo'][dbi]
            chi = data_bur['logmvir_msun_hi'][dbi]
            
            issig0 = (clo <= logmvir_msun and chi >= logmvir_msun)
            _label = None
            if issig0:
                _color = 'black'
                if isul and not ulsig0done:
                    _label = ('UL, $\\Delta\\mathrm{M}'
                              f' < {nsigmas[0]}\\sigma$')
                    ulsig0done = True
                elif not isul and not detsig0done:
                    _label = ('det., $\\Delta\\mathrm{M}'
                              f' < {nsigmas[0]}\\sigma$')
                    detsig0done = True
            else:
                _color = 'gray'
                if isul and not ulsig1done:
                    _label = ('UL, $\\Delta\\mathrm{M}'
                              f' < {nsigmas[1]}\\sigma$')
                    ulsig1done = True
                elif not isul and not detsig1done:
                    _label = ('det., $\\Delta\\mathrm{M}'
                              f' < {nsigmas[1]}\\sigma$')
                    detsig1done = True
            marker = 'v' if isul else 'o'
            markersize = 5
            zobase = 5. - 1. * isul

            ax.errorbar([xv], [yv], yerr=yerr, 
                        linestyle='none', elinewidth=1.5,
                        color=_color, capsize=3,
                        zorder=zobase,
                        marker=marker, markersize=markersize,
                        markeredgecolor='black', markeredgewidth=1.0,
                        label=_label)
        if detsig1done and detsig0done and ulsig1done and ulsig0done:
            getlegax = plii
    ylims = [ax.get_ylim() for ax in axes]
    ymin = min([ylim[0] for ylim in ylims])
    ymin = max(ymin, 12.5)
    ymax = max([ylim[1] for ylim in ylims])
    #ymax = ymax + 0.5
    [ax.set_ylim((ymin, ymax)) for ax in axes]

    #handles1, _ = axes[getlegax].get_legend_handles_labels()
    #axes[-1].legend(handles=handles1, fontsize=fontsize - 2,
    #                loc='upper right', bbox_to_anchor=(1.0, 1.0),
    #                handlelength=1.0, labelspacing=0.3,
    #                handletextpad=0.4)
    #label=(f'${fcgm:.1f} \\,'
    #                                 '\\mathrm{f}_{\\mathrm{b, c}}$,\n'
    #                                 f'${z_sol:.1f} \\, '
    #                                 '\\mathrm{Z}_{\\odot}$')
    handles3 = [mlines.Line2D((), (), color=color, linewidth=2.,
                              label=('$\\mathrm{f}_{\\mathrm{CGM}} = '
                                     f'{fcgm:.1f}$'),
                              linestyle='solid')
                for color, fcgm, z_sol in zip(colors_fZ, fcgms, z_sols)]
    axes[0].legend(handles=handles3[:3], fontsize=fontsize - 2,
                   loc='upper right', bbox_to_anchor=(1.0, 1.0),
                   handlelength=1.0, labelspacing=0.0,
                   handletextpad=0.4)
    #axes[1].legend(handles=handles3[1:2], fontsize=fontsize - 2,
    #               loc='upper right', bbox_to_anchor=(1.0, 1.0),
    #               handlelength=1.0, labelspacing=0.0,
    #               handletextpad=0.4)
    #axes[2].legend(handles=handles3[2:], fontsize=fontsize - 2,
    #               loc='upper right', bbox_to_anchor=(1.0, 1.0),
    #               handlelength=1.0, labelspacing=0.0,
    #               handletextpad=0.4)
    
    pli_vc_str = 'pli_vc_' + \
                 '_'.join([f'{pli_vc:.2f}' for pli_vc in plis_vc])
    pli_k_str = 'pli_k_' + \
                 '_'.join([f'{pli_k:.2f}' for pli_k in plis_k])
    zsol_str = 'Zsol_' + \
                 '_'.join([f'{z_sol:.2f}' for z_sol in z_sols])

    outname = (f'prof_Ne8_analytical_pl_s19_z{redshift_model:.2f}'
               f'_vs_b19_mvir{logmvir_msun:.1f}'
               f'_{pli_k_str}_{pli_vc_str}_{zsol_str}_fcgm_var')
    outname = outname.replace('.', 'p')
    outname = outname.replace('-', 'm')
    plt.savefig(outdir + outname + '.pdf', bbox_inches='tight')

def plot_plmodel_datacomp_parvar():
    ion = 'Ne8'
    redshift_model = 0.75
    nsigmas = (1, 2)

    impactpars_kpc = np.linspace(5., 450., 50)
    logmvir_msun = 12.2
    fcgms = [1.0, 0.3, 0.1]
    fcgm_def = 0.3
    z_sols = [1.0, 0.3, 0.1]
    z_sol_def = 0.3
    _colors = tc.tol_cset('bright')
    colors_fZ = [_colors.blue, _colors.green, _colors.yellow]
    redshift = 0.75
    plis_vc = [0.0, -0.2, -0.5]
    pli_vc_def = -0.1 
    colors_pli_vc = [_colors.cyan, _colors.red, _colors.purple]
    plis_k = [0.0, 2./3., 1.2]

    data_bur = odr.readdata_b19(nsigmas=nsigmas)
    
    panelsize = 1.8
    ncols = 3
    nrows = 3
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] * nrows
    hspace = 0.4
    height = sum(height_ratios) * (1. + hspace * (nrows - 1) / nrows)
 
    fig = plt.figure(figsize=(sum(width_ratios), height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows, wspace=0.0, 
                        hspace=hspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [[fig.add_subplot(grid[i, j]) 
             for j in range(ncols) if i < 2 or j > 0]
            for i in range(nrows)]
    fontsize = 12
    linewidth = 1.5
    
    for ri in range(nrows):
        for ci in range(ncols):
            if ri == 2 and ci > 0:
                continue
            ax = axes[ri][ci]
            doleft = ci == 0
            dobottom = ri == 2 or (ri == 1 and ci != 1)
            if dobottom and ri == 2:
                ax.set_xlabel('$\\mathrm{r}_{\\perp} \\; [\\mathrm{pkpc}]$',
                              fontsize=fontsize)
            if doleft and ri == 1:
                ax.set_ylabel('$\\log_{10} \\, \\mathrm{N}('
                              '\\mathrm{Ne\\,VIII})'
                              '\\; [\\mathrm{cm}^{-2}]$',
                              fontsize=fontsize)
            ax.tick_params(which='both', direction='in', 
                           labelsize=fontsize - 1.,
                           top=True, right=True, labelleft=doleft,
                           labelbottom=dobottom)
            pli_k = plis_k[ci]
            if ri == 2:
                pli_k = plis_k[1]
            axlabel = f'$\\mathrm{{K}} \\propto r^{{{pli_k:.2f}}}$'
            ax.text(0.08, 0.95, axlabel,
                    transform=ax.transAxes, fontsize=fontsize,
                    verticalalignment='top', horizontalalignment='left')
            yvir_max = - np.inf
            
            labels_r = []
            colors_r = []
            cvs_r = []
            if ri == 0:
                pli_vc = pli_vc_def
                fcgm = fcgms[ci]
                z_sol = z_sol_def
                subtitle = (f'$\\mathrm{{Z}} = {z_sol:.1f} \\,'
                            '\\mathrm{Z}_{\\odot},\\;'
                            'v_{\\mathrm{c}} \\propto '
                            f'r^{{{pli_vc:.1f}}},\\;'
                            '\\mathrm{f}_{\\mathrm{CGM}}:$')
                print(pli_vc, pli_k)
   
                for fi, fcgm in enumerate(fcgms):
                    labels_r.append(f'${fcgm:.1f}$')
                    colors_r.append(colors_fZ[fi])
                    hmod = mip.PLmodel(10**logmvir_msun, redshift, fcgm, 
                                       z_sol, pli_vc, pli_entropy=pli_k)
                    cvs = hmod.coldensprof(ion, impactpars_kpc)
                    cvs_r.append(cvs)
                    rvir = hmod.rvir_cgs / (c.cm_per_mpc * 1e-3)
            elif ri == 1:
                pli_vc = pli_vc_def
                fcgm = fcgm_def
                z_sol = z_sols[ci]
                subtitle = ('$\\mathrm{f}_{\\mathrm{CGM}} '
                            f'= {fcgm:.1f},\\;'
                            'v_{\\mathrm{c}} \\propto '
                            f'r^{{{pli_vc:.1f}}},\\;'
                            '\\mathrm{Z} \\, / \\, \\mathrm{Z}_{\\odot}:$')
                print(pli_vc, pli_k)
   
                for zi, z_sol in enumerate(z_sols):
                    labels_r.append(f'${z_sol:.1f}$')
                    colors_r.append(colors_fZ[zi])
                    hmod = mip.PLmodel(10**logmvir_msun, redshift, fcgm, 
                                       z_sol, pli_vc, pli_entropy=pli_k)
                    cvs = hmod.coldensprof(ion, impactpars_kpc)
                    cvs_r.append(cvs)
                    rvir = hmod.rvir_cgs / (c.cm_per_mpc * 1e-3)
            elif ri == 2:
                if ci > 0: 
                    continue
                fcgm = fcgm_def
                z_sol = z_sol_def
                subtitle = ('$\\mathrm{f}_{\\mathrm{CGM}} '
                            f'= {fcgm:.1f}, \\;'
                            f'\\mathrm{{Z}} = {z_sol:.1f} \\,'
                            '\\mathrm{Z}_{\\odot} $'
                            )
                print(pli_vc, pli_k)
   
                for plii, pli_vc in enumerate(plis_vc):
                    labels_r.append('$v_{\\mathrm{c}} \\propto '
                                     f'r^{{{pli_vc:.1f}}} $')
                    colors_r.append(colors_pli_vc[plii])
                    hmod = mip.PLmodel(10**logmvir_msun, redshift, fcgm, 
                                       z_sol, pli_vc, pli_entropy=pli_k)
                    cvs = hmod.coldensprof(ion, impactpars_kpc)
                    cvs_r.append(cvs)
                    rvir = hmod.rvir_cgs / (c.cm_per_mpc * 1e-3)
            yvir_max = -np.inf
            for color, label, cvs in zip(colors_r, labels_r, cvs_r):
                ax.plot(impactpars_kpc, np.log10(cvs), 
                        color=color, linestyle='solid',
                        linewidth=linewidth, label=label)
                if rvir <= impactpars_kpc[-1]:
                    yv = mu.linterpsolve(impactpars_kpc, np.log10(cvs), rvir)
                    yvir_max = max(yvir_max, yv)
            if rvir <= impactpars_kpc[-1]:
                ax.plot([rvir, rvir], [11., yvir_max + 0.2], color='gray',
                        linestyle='solid', alpha=0.3, linewidth=3)
            if ri == 2:
                ax.set_title(subtitle, fontsize=fontsize - 1)
            elif ci == 0:
                ax.text(0.0, 1.02, subtitle, fontsize=fontsize - 1, 
                        horizontalalignment='left', 
                        verticalalignment='bottom',
                        transform=ax.transAxes)

            ulsig0done = False
            ulsig1done = False
            detsig0done = False
            detsig1done = False
            for dbi in range(len(data_bur)):
                cloer = data_bur['logmvir_msun_loer'][dbi]
                chier = data_bur['logmvir_msun_hier'][dbi]
                if cloer > logmvir_msun or chier < logmvir_msun:
                    continue

                xv = data_bur['impact_parameter_kpc'][dbi]
                yv = data_bur['log_N_Ne8_pcm2'][dbi]
                isul = data_bur['log_N_Ne8_isUL'][dbi]
                yerr = data_bur['log_N_Ne8_pcm2_err'][dbi] if not isul else None
                #cbest = data_bur['logmvir_msun_bestest'][dbi]
                clo = data_bur['logmvir_msun_lo'][dbi]
                chi = data_bur['logmvir_msun_hi'][dbi]
                
                issig0 = (clo <= logmvir_msun and chi >= logmvir_msun)
                _label = None
                if issig0:
                    _color = 'black'
                    if isul and not ulsig0done:
                        _label = ('UL, $\\Delta\\mathrm{M}'
                                f' < {nsigmas[0]}\\sigma$')
                        ulsig0done = True
                    elif not isul and not detsig0done:
                        _label = ('det., $\\Delta\\mathrm{M}'
                                f' < {nsigmas[0]}\\sigma$')
                        detsig0done = True
                else:
                    _color = 'gray'
                    if isul and not ulsig1done:
                        _label = ('UL, $\\Delta\\mathrm{M}'
                                f' < {nsigmas[1]}\\sigma$')
                        ulsig1done = True
                    elif not isul and not detsig1done:
                        _label = ('det., $\\Delta\\mathrm{M}'
                                f' < {nsigmas[1]}\\sigma$')
                        detsig1done = True
                marker = 'v' if isul else 'o'
                markersize = 5
                zobase = 5. - 1. * isul

                ax.errorbar([xv], [yv], yerr=yerr, 
                            linestyle='none', elinewidth=1.5,
                            color=_color, capsize=3,
                            zorder=zobase,
                            marker=marker, markersize=markersize,
                            markeredgecolor='black', markeredgewidth=1.0,
                            label=_label)
            if detsig1done and detsig0done and ulsig1done and ulsig0done:
                getlegax = (ri, ci)
    ylims = [ax.get_ylim() for sub in axes for ax in sub]
    ymin = min([ylim[0] for ylim in ylims])
    ymin = max(ymin, 12.5)
    ymax = max([ylim[1] for ylim in ylims])
    #ymax = ymax + 0.5
    [ax.set_ylim((ymin, ymax)) for sub in axes for ax in sub]

    handles0, labels0 = axes[0][0].get_legend_handles_labels()
    axes[0][1].legend(handles=handles0[:3], labels=labels0,
                      fontsize=fontsize - 1,
                      loc='lower left', bbox_to_anchor=(0.49, 0.97),
                      handlelength=1.0, labelspacing=0.0,
                      handletextpad=0.4, ncol=3, columnspacing=0.7,
                      frameon=False, borderpad=0.0)
    handles1, labels1 = axes[1][0].get_legend_handles_labels()
    axes[1][1].legend(handles=handles1[:3], labels=labels1,
                      fontsize=fontsize - 1,
                      loc='lower left', bbox_to_anchor=(0.52, 0.97),
                      handlelength=1.0, labelspacing=0.0,
                      handletextpad=0.4, ncol=3, columnspacing=0.7,
                      frameon=False, borderpad=0.0)
    handles2, labels2 = axes[2][0].get_legend_handles_labels()
    axes[2][1].legend(handles=handles2[:3], labels=labels2,
                      fontsize=fontsize - 1,
                      loc='lower left', bbox_to_anchor=(0.0, 0.0),
                      handlelength=1.0, labelspacing=0.0,
                      handletextpad=0.4)
    axes[2][1].axis('off')

    pli_vc_str = 'pli_vc_' + \
                 '_'.join([f'{pli_vc:.2f}' for pli_vc in plis_vc])
    pli_k_str = 'pli_k_' + \
                 '_'.join([f'{pli_k:.2f}' for pli_k in plis_k])
    zsol_str = 'Zsol_' + \
                 '_'.join([f'{z_sol:.2f}' for z_sol in z_sols])
    fcgm_str = 'fCGM_' + \
                 '_'.join([f'{fcgm:.2f}' for fcgm in fcgms])

    outname = (f'prof_Ne8_analytical_pl_s19_z{redshift_model:.2f}'
               f'_vs_b19_mvir{logmvir_msun:.1f}'
               f'_{pli_k_str}_{pli_vc_str}_{zsol_str}_{fcgm_str}')
    outname = outname.replace('.', 'p')
    outname = outname.replace('-', 'm')
    plt.savefig(outdir + outname + '.pdf', bbox_inches='tight')

