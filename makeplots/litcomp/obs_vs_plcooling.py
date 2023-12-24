import matplotlib.colors as mcolors
import matplotlib.cm as mcm
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fire_an.analytic_halo.model_ionprof_coolingpl as mic
import fire_an.makeplots.litcomp.obsdataread as odr
import fire_an.makeplots.tol_colors as tc
import fire_an.mstar_mhalo.loader_smdpl_sfr as ldsmdpl
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu
import fire_an.utils.math_utils as mu

outdir = '/projects/b1026/nastasha/imgs/analytical/'
oddir = '/projects/b1026/nastasha/extdata/'
q23filen = oddir + 'plotdata_q23_nsigmas_1_2.dat'
b19filen = oddir + 'plotdata_b19_nsigmas_1_2.dat'

def get_sfrs_mh(logmhs_msun, z=0.75, percentiles=(0.16, 0.5, 0.84)):
    histobj = ldsmdpl.SFRHMhists(np.array([z]))
    sfrs_tab, mhs_tab = histobj.getperc_sfrmh(z, mode='mhtosfr', 
                                              percvals=np.array(percentiles))
    out = [[mu.linterpsolve(mhs_tab, sfrs_tab_perc, mh_this) 
            for sfrs_tab_perc in sfrs_tab]
           for mh_this in logmhs_msun]
    return np.array(out)

def get_obsdata(obsdata=('B+19', 'Q+23')):
    nsigmas = (1, 2)
    odata = {}
    if 'B+19' in obsdata:
        data_bur = pd.read_csv(b19filen, sep='\t')
        cloer = data_bur['logmvir_msun_loer']
        chier = data_bur['logmvir_msun_hier']
        rp = data_bur['impact_parameter_kpc']
        col = data_bur['log_N_Ne8_pcm2']
        isul = data_bur['log_N_Ne8_isUL']
        yerr = np.array((data_bur['log_N_Ne8_pcm2_err'],
                         data_bur['log_N_Ne8_pcm2_err']))
        clo = data_bur['logmvir_msun_lo']
        chi = data_bur['logmvir_msun_hi']
        flag = data_bur['flagged_by_qu23']
        odata['B+19'] = {'cloer': cloer, 'chier': chier,
                         'clo': clo, 'chi': chi, 'isul': isul,
                         'rp': rp, 'col': col, 'yerr': yerr,
                         'flag': flag}
    if 'Q+23' in obsdata:
        data_qu = odr.getplotdata_cubs()
        cloer = data_qu['logmvir_msun_loer']
        chier = data_qu['logmvir_msun_hier']
        rp = data_qu['impactpar_kpc']
        col = data_qu['ne8col_logcm2']
        isul = data_qu['isul_ne8']
        yerr = np.array((data_qu['ne8col_2s_loerr_dex'],
                         data_qu['ne8col_2s_loerr_dex']))
        clo = data_qu['logmvir_msun_lo']
        chi = data_qu['logmvir_msun_hi']
        flag = np.zeros(rp.shape, dtype=bool)
        odata['Q+23'] = {'cloer': cloer, 'chier': chier,
                         'clo': clo, 'chi': chi, 'isul': isul,
                         'rp': rp, 'col': col, 'yerr': yerr,
                         'flag': flag}
    return odata

def addobsdata_panel(ax, odata, mvir, colors):
    nsigmas = (1, 2)
    for odkey in odata:
        _colors = colors[odkey]
        _odata = odata[odkey]
        cloer = _odata['cloer']
        chier = _odata['chier']
        clo = _odata['clo']
        chi = _odata['chi']
        isul = _odata['isul']
        rp = _odata['rp']
        col = _odata['col']
        yerr = _odata['yerr']
        flag = _odata['flag']

        f1 = np.logical_and(cloer <= mvir, chier >= mvir)
        issig0 = np.logical_and(clo <= mvir, chi >= mvir)
        ul_sig0 = np.logical_and(isul, f1)
        ul_sig0 &= issig0
        ul_sig1 = np.logical_and(isul, f1)
        ul_sig1 &= np.logical_not(issig0)
        msm = np.logical_and(np.logical_not(isul), np.logical_not(flag))
        msf = np.logical_and(np.logical_not(isul), flag)
        msm_sig0 = np.logical_and(msm, f1)
        msm_sig0 &= issig0
        msm_sig1 = np.logical_and(msm, f1)
        msm_sig1 &= np.logical_not(issig0)
        msf_sig0 = np.logical_and(msf, f1)
        msf_sig0 &= issig0
        msf_sig1 = np.logical_and(msf, f1)
        msf_sig1 &= np.logical_not(issig0)
        
        markersize = 5
        ax.errorbar(rp[ul_sig0], col[ul_sig0], None, 
                        linestyle='none', elinewidth=1.5,
                        color=_colors['1s'], capsize=3,
                        zorder=4.,
                        marker='v', markersize=markersize,
                        markeredgecolor='black', markeredgewidth=1.0,
                        label=('UL, $\\Delta\\mathrm{M}'
                            f' < {nsigmas[0]}\\sigma$'))
        ax.errorbar(rp[ul_sig1], col[ul_sig1], None, 
                        linestyle='none', elinewidth=1.5,
                        color=_colors['2s'], capsize=3,
                        zorder=4.,
                        marker='v', markersize=markersize,
                        markeredgecolor='black', markeredgewidth=1.0,
                        label=('UL, $\\Delta\\mathrm{M}'
                            f' < {nsigmas[1]}\\sigma$'))
        ax.errorbar(rp[msm_sig0], col[msm_sig0], yerr=yerr[:, msm_sig0], 
                        linestyle='none', elinewidth=1.5,
                        color=_colors['1s'], capsize=3,
                        zorder=5.,
                        marker='o', markersize=markersize,
                        markeredgecolor='black', markeredgewidth=1.0,
                        label=('det., $\\Delta\\mathrm{M}'
                            f' < {nsigmas[0]}\\sigma$'))
        ax.errorbar(rp[msm_sig1], col[msm_sig1], yerr=yerr[:, msm_sig1],
                        linestyle='none', elinewidth=1.5,
                        color=_colors['2s'], capsize=3,
                        zorder=5.,
                        marker='o', markersize=markersize,
                        markeredgecolor='black', markeredgewidth=1.0,
                        label=('det., $\\Delta\\mathrm{M}'
                            f' < {nsigmas[1]}\\sigma$'))
        ax.errorbar(rp[msf_sig0], col[msf_sig0], yerr=yerr[:, msf_sig0], 
                        linestyle='none', elinewidth=1.5,
                        color=_colors['1s'], capsize=3,
                        zorder=5.,
                        marker='o', markersize=markersize,
                        markeredgecolor=_colors['1s'], markeredgewidth=1.5,
                        markerfacecolor='none',
                        label=None)
        ax.errorbar(rp[msf_sig1], col[msf_sig1], yerr=yerr[:, msf_sig1],
                        linestyle='none', elinewidth=1.5,
                        color=_colors['2s'], capsize=3,
                        zorder=5.,
                        marker='o', markersize=markersize,
                        markeredgecolor=_colors['2s'], markeredgewidth=1.5,
                        markerfacecolor='none',
                        label=None)



def plot_plcoolingmodel_datacomp(obsdata=('Q+23', 'B+19')):
    ion = 'Ne8'
    redshift_model = 0.75
    nsigmas = (1, 2)

    impactpars_kpc = np.linspace(5., 450., 50)
    logmvirs_msun = np.arange(11.3, 13.5, 0.3)
    z_sol = 0.1
    redshift = 0.75
    plis_vc = [0.0, -0.1, -0.20]
    colors_vc = tc.tol_cset('vibrant')
    pmdots = [0.16, 0.5, 0.84]
    linestyles_pmdot = ['dotted', 'solid', 'dashed']
    
    if len(obsdata) == 1:
       _colors = {'1s': 'black', '2s': 'gray'}
       colors = {obsdata[0]: _colors}
    else:
        _colors = tc.tol_cset('high-contrast')
        c1 = mcolors.to_rgba(_colors.blue)
        c2 = mcolors.to_rgba(_colors.red)
        colors = {'B+19': {'1s': c1, '2s': c1[:3] + (0.5,)},
                  'Q+23': {'1s': c2, '2s': c2[:3] + (0.5,)}}
    odata = get_obsdata(obsdata=obsdata)

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
        sfrs = get_sfrs_mh([mvir], z=redshift, percentiles=pmdots)
        for sfr, ls in zip(sfrs[0, :], linestyles_pmdot):
            for color, pli_vc in zip(colors_vc, plis_vc):
                hmod = mic.PLmodel(10**mvir, redshift, 10**sfr, z_sol, pli_vc)
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
        addobsdata_panel(ax, odata, mvir, colors)

    ylims = [ax.get_ylim() for ax in axes]
    ymin = min([ylim[0] for ylim in ylims])
    ymin = max(ymin, 12.5)
    ymax = max([ylim[1] for ylim in ylims])
    [ax.set_ylim((ymin, ymax)) for ax in axes]
    [ax.set_xlim(0., 450.) for ax in axes]

    handles1, _ = axes[0].get_legend_handles_labels()
    axes[-1].legend(handles=handles1[:4], fontsize=fontsize - 2,
                    loc='upper right', bbox_to_anchor=(1.0, 0.86),
                    handlelength=1.0, labelspacing=0.3,
                    handletextpad=0.4)
    if len(pmdots) > 1:
        handles2 = [mlines.Line2D(
                        (), (), color='black', linewidth=linewidth,
                        linestyle=ls,
                        label=f'${pmdot*100:.0f}^{{\\mathrm{{th}}}}$% SFR')
                    for pmdot, ls in zip(pmdots, linestyles_pmdot)]
        axes[-2].legend(handles=handles2, fontsize=fontsize - 2,
                        loc='upper right', bbox_to_anchor=(1.0, 0.86),
                        handlelength=1.5, labelspacing=0.1,
                        handletextpad=0.4)
    if len(obsdata) > 1:
        handles2 = [mlines.Line2D((), (), linestyle='none', marker='s',
                                markersize=5, label=oset, 
                                markeredgecolor='black',
                                c=colors[oset]['1s'])
                    for oset in obsdata]
        if 'B+19' in obsdata:
            handles2 += [mlines.Line2D((), (), linestyle='none', marker='s',
                         markersize=5, label='B+19 (!)', 
                         markeredgecolor=colors['B+19']['1s'],
                         markeredgewidth=2.0, markerfacecolor='none',
                         c=colors['B+19']['1s'])]
        axes[-3].legend(handles=handles2,
                        fontsize=fontsize - 2,
                        loc='upper right', bbox_to_anchor=(1.0, 0.86),
                        handlelength=1.0, labelspacing=0.3,
                        handletextpad=0.4)
    if len(plis_vc) > 1:
        handles2 = [mlines.Line2D(
                        (), (), color=color, linewidth=linewidth,
                        linestyle='solid',
                        label='$v_{\\mathrm{c}} \\propto'
                              f' r^{{{pli:.1f}}}$')
                    for pli, color in zip(plis_vc, colors_vc)]
        axes[-4].legend(handles=handles2,
                        fontsize=fontsize - 2,
                        loc='upper right', bbox_to_anchor=(1.0, 0.86),
                        handlelength=1.0, labelspacing=0.0,
                        handletextpad=0.4)
        
    
    pli_vc_str = 'pli_vc_' + \
                 '_'.join([f'{pli_vc:.2f}' for pli_vc in plis_vc])
    pmdot_str = 'pmdot_' + \
                 '_'.join([f'{pmdot:.2f}' for pmdot in pmdots])

    ods = '_'.join(obsdata)
    outname = (f'prof_Ne8_analytical_coolingpl_s19_z{redshift_model:.2f}'
               f'_vs_{ods}_loghalomass_{pli_vc_str}_{pmdot_str}'
               f'_Zsol{z_sol:.2f}')
    outname = outname.replace('.', 'p')
    outname = outname.replace('-', 'm')
    plt.savefig(outdir + outname + '.pdf', bbox_inches='tight')

