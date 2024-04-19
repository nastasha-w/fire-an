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
oddir = '/projects/b1026/nastasha/extdata/'
q23filen = oddir + 'plotdata_q23_nsigmas_1_2.dat'
b19filen = oddir + 'plotdata_b19_nsigmas_1_2.dat'

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

def get_explabel(exponent, base='r', ndec_max=2):
    if exponent == 0.:
        out = '\\mathrm{\\;const.}' 
    elif np.round(exponent, ndec_max) == 1.:
        out = f'\\propto {{{base}}}'
    else:
        # get desired precision string
        fmtstr = f'{{exp:.{ndec_max}f}}'
        exppart = fmtstr.format(exp=exponent)
        # strip trailing zeros (and any decimal point left trailing)
        if '.' in exppart:
            while exppart.endswith('0'):
                exppart = exppart[:-1]
            if exppart.endswith('.'):
                exppart = exppart[:-1]
        out = f'\\propto {{{base}}}^{{{exppart}}}'
    return out


def plot_plmodel_datacomp_Kvar(obsdata=('B+19',)):
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
    plis_k = [0.0, 2./3., 1., 1.2]
    linestyles_k = ['dotted', 'dashed', 'solid', 'dashdot']
    
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
    if len(plis_k) > 1:
        handles2 = [mlines.Line2D(
                        (), (), color='black', linewidth=linewidth,
                        linestyle=ls,
                        label=f'$\\mathrm{{K}}'
                              f' {get_explabel(pli_k, base="r", ndec_max=2)}$')
                    for pli_k, ls in zip(plis_k, linestyles_k)]
        axes[-2].legend(handles=handles2, fontsize=fontsize - 1.,
                        loc='upper right', bbox_to_anchor=(0.93, 0.87),
                        handlelength=2., labelspacing=0.1,
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
                        loc='upper right', bbox_to_anchor=(1.0, 0.88),
                        handlelength=1.0, labelspacing=0.3,
                        handletextpad=0.4)
    
    pli_vc_str = 'pli_vc_' + \
                 '_'.join([f'{pli_vc:.2f}' for pli_vc in plis_vc])
    pli_k_str = 'pli_k_' + \
                 '_'.join([f'{pli_k:.2f}' for pli_k in plis_k])

    ods = '_'.join(obsdata)
    outname = (f'prof_Ne8_analytical_pl_s19_z{redshift_model:.2f}'
               f'_vs_{ods}_loghalomass_{pli_k_str}_{pli_vc_str}'
               '_nhnorm_0p1_to_1_Rvir')
    outname = outname.replace('.', 'p')
    outname = outname.replace('-', 'm')
    plt.savefig(outdir + outname + '.pdf', bbox_inches='tight')


def plot_plmodel_datacomp_parvar_old(obsdata=('Q+23', 'B+19')):
    ion = 'Ne8'
    redshift_model = 0.75
    nsigmas = (1, 2)

    impactpars_kpc = np.linspace(5., 450., 50)
    logmvir_msun = 12.2
    fcgms = [1.0, 0.3, 0.1]
    fcgm_def = 0.3
    #z_sols = [1.0, 0.3, 0.1]
    z_sol_def = 0.3
    _colors = tc.tol_cset('bright')
    colors_fZ = [_colors.blue, _colors.green, _colors.yellow]
    redshift = 0.75
    plis_vc = [0.0, -0.2, -0.5]
    pli_vc_def = -0.1 
    colors_pli_vc = [_colors.cyan, _colors.red, _colors.purple]
    plis_k = [0.0, 2./3., 1.2]
    
    if len(obsdata) == 1:
       _colors = {'1s': 'black', '2s': 'gray'}
       odcolors = {obsdata[0]: _colors}
    else:
        _colors = tc.tol_cset('high-contrast')
        c1 = mcolors.to_rgba(_colors.blue)
        c2 = mcolors.to_rgba(_colors.red)
        odcolors = {'B+19': {'1s': c1, '2s': c1[:3] + (0.5,)},
                    'Q+23': {'1s': c2, '2s': c2[:3] + (0.5,)}}
    odata = get_obsdata(obsdata=(obsdata))
    
    panelsize = 1.8
    ncols = 3
    nrows = 2
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] * nrows
    hspace = 0.4
    height = sum(height_ratios) * (1. + hspace * (nrows - 1) / nrows)
 
    fig = plt.figure(figsize=(sum(width_ratios), height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows, wspace=0.0, 
                        hspace=hspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [[fig.add_subplot(grid[i, j]) 
             for j in range(ncols) if i < 1 or j > 0]
            for i in range(nrows)]
    fontsize = 12
    linewidth = 1.5
    
    for ri in range(nrows):
        for ci in range(ncols):
            if ri == 1 and ci > 0:
                continue
            ax = axes[ri][ci]
            doleft = ci == 0
            dobottom = ri == 1 or (ri == 0 and ci != 1)
            if dobottom and ri == 1:
                ax.set_xlabel('$\\mathrm{r}_{\\perp} \\; [\\mathrm{pkpc}]$',
                              fontsize=fontsize)
            if doleft and ri == 0:
                ax.set_ylabel('$\\log_{10} \\, \\mathrm{N}('
                              '\\mathrm{Ne\\,VIII})'
                              '\\; [\\mathrm{cm}^{-2}]$',
                              fontsize=fontsize)
            ax.tick_params(which='both', direction='in', 
                           labelsize=fontsize - 1.,
                           top=True, right=True, labelleft=doleft,
                           labelbottom=dobottom)
            pli_k = plis_k[ci]
            if ri == 1:
                pli_k = plis_k[1]
            rslope = get_explabel(pli_k, base='r', ndec_max=2)
            axlabel = f'$\\mathrm{{K}} {rslope}$'
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
                    rslope = get_explabel(pli_vc, base='r', ndec_max=2)
                    labels_r.append(f'$v_{{\\mathrm{{c}}}} {rslope}$')
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
            if ri == 1:
                ax.set_title(subtitle, fontsize=fontsize - 1)
            elif ci == 0:
                ax.text(0.0, 1.02, subtitle, fontsize=fontsize - 1, 
                        horizontalalignment='left', 
                        verticalalignment='bottom',
                        transform=ax.transAxes)

            addobsdata_panel(ax, odata, logmvir_msun, odcolors)

    ylims = [ax.get_ylim() for sub in axes for ax in sub]
    ymin = min([ylim[0] for ylim in ylims])
    ymin = max(ymin, 12.5)
    ymax = max([ylim[1] for ylim in ylims])
    #ymax = ymax + 0.5
    [ax.set_ylim((ymin, ymax)) for sub in axes for ax in sub]
    [ax.set_xlim(0., 450.) for sub in axes for ax in sub]

    handles0, labels0 = axes[0][0].get_legend_handles_labels()
    axes[0][1].legend(handles=handles0[:3], labels=labels0,
                      fontsize=fontsize - 1,
                      loc='lower left', bbox_to_anchor=(0.49, 0.97),
                      handlelength=1.0, labelspacing=0.0,
                      handletextpad=0.4, ncol=3, columnspacing=0.7,
                      frameon=False, borderpad=0.0)
    handles2, labels2 = axes[1][0].get_legend_handles_labels()
    axes[1][1].legend(handles=handles2[:3], labels=labels2,
                      fontsize=fontsize - 1,
                      loc='lower left', bbox_to_anchor=(0.0, 0.0),
                      handlelength=1.0, labelspacing=0.0,
                      handletextpad=0.4)
    axes[1][1].axis('off')

    pli_vc_str = 'pli_vc_' + \
                 '_'.join([f'{pli_vc:.2f}' for pli_vc in plis_vc])
    pli_k_str = 'pli_k_' + \
                 '_'.join([f'{pli_k:.2f}' for pli_k in plis_k])
    #zsol_str = 'Zsol_' + \
    #             '_'.join([f'{z_sol:.2f}' for z_sol in z_sols])
    fcgm_str = 'fCGM_' + \
                 '_'.join([f'{fcgm:.2f}' for fcgm in fcgms])
    ods = '_'.join(obsdata)
    outname = (f'prof_Ne8_analytical_pl_s19_z{redshift_model:.2f}'
               f'_vs_{ods}_mvir{logmvir_msun:.1f}'
               f'_{pli_k_str}_{pli_vc_str}_Zsol_0.30_{fcgm_str}'
               '_nhnorm_0p1_to_1_Rvir')
    outname = outname.replace('.', 'p')
    outname = outname.replace('-', 'm')
    plt.savefig(outdir + outname + '.pdf', bbox_inches='tight')

def plot_plmodel_datacomp_parvar(obsdata=('Q+23', 'B+19')):
    ion = 'Ne8'
    redshift_model = 0.75
    nsigmas = (1, 2)

    impactpars_kpc = np.linspace(5., 465., 50)
    logmvir_msun = 12.2
    fcgms = [1.0, 0.3, 0.1]
    fcgm_def = 0.3
    #z_sols = [1.0, 0.3, 0.1]
    z_sol_def = 0.3
    _colors = tc.tol_cset('bright')
    colors_fZ = [_colors.blue, _colors.green, _colors.yellow]
    redshift = 0.75
    plis_vc = [0.0, -0.2, -0.5]
    pli_vc_def = -0.1 
    colors_pli_vc = [_colors.cyan, _colors.red, _colors.purple]
    plis_k = [0.0, 2./3., 1., 1.2]
    
    if len(obsdata) == 1:
       _colors = {'1s': 'black', '2s': 'gray'}
       odcolors = {obsdata[0]: _colors}
    else:
        _colors = tc.tol_cset('high-contrast')
        c1 = mcolors.to_rgba(_colors.blue)
        c2 = mcolors.to_rgba(_colors.red)
        odcolors = {'B+19': {'1s': c1, '2s': c1[:3] + (0.5,)},
                    'Q+23': {'1s': c2, '2s': c2[:3] + (0.5,)}}
    odata = get_obsdata(obsdata=(obsdata))
    
    panelsize = 1.8
    ncols = len(plis_k)
    nrows = 2
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] * nrows
    hspace = 0.3
    height = sum(height_ratios) * (1. + hspace * (nrows - 1) / nrows)
 
    fig = plt.figure(figsize=(sum(width_ratios), height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows, wspace=0.0, 
                        hspace=hspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [[fig.add_subplot(grid[i, j]) 
             for j in range(ncols) if i < 1 or j > 0]
             for i in range(nrows)]
    fontsize = 12
    linewidth = 1.5
    
    for ri in range(nrows):
        for ci in range(ncols):
            if ri == 1 and ci >= len(plis_k) - 2:
                continue
            ax = axes[ri][ci]
            doleft = ci == 0
            dobottom = ri == 1 or (ri == 0 and 
                                   (ci == 0 or ci == len(plis_k) - 1))
            if dobottom and ri == 1:
                ax.set_xlabel('$\\mathrm{r}_{\\perp} \\; [\\mathrm{pkpc}]$',
                              fontsize=fontsize)
            if doleft and ri == 0:
                ax.set_ylabel('$\\log_{10} \\, \\mathrm{N}('
                              '\\mathrm{Ne\\,VIII})'
                              '\\; [\\mathrm{cm}^{-2}]$',
                              fontsize=fontsize)
            ax.tick_params(which='both', direction='in', 
                           labelsize=fontsize - 1.,
                           top=True, right=True, labelleft=doleft,
                           labelbottom=dobottom)
            pli_k = plis_k[ci]
            if ri == 1:
                pli_k = plis_k[ci + 1]
            rslope = get_explabel(pli_k, base='r', ndec_max=2)
            axlabel = f'$\\mathrm{{K}} {rslope}$'
            ax.text(0.92, 0.92, axlabel,
                    transform=ax.transAxes, fontsize=fontsize,
                    verticalalignment='top', horizontalalignment='right')
            yvir_max = -np.inf
            
            labels_r = []
            colors_r = []
            cvs_r = []
            if ri == 0:
                pli_vc = pli_vc_def
                #fcgm = fcgms[ci]
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
                if ci >= len(plis_k) - 2: 
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
                    rslope = get_explabel(pli_vc, base='r', ndec_max=2)
                    labels_r.append(f'$v_{{\\mathrm{{c}}}} {rslope}$')
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
            if ri == 1 and ci == 0:
                #ax.set_title(subtitle, fontsize=fontsize - 1)
                ax.text(1.0, 1.02, subtitle, fontsize=fontsize - 1.,
                       transform=ax.transAxes,
                       horizontalalignment='center',
                       verticalalignment='bottom')
            elif ri == 0 and ci == 0:
                ax.text(0.65, 1.02, subtitle, fontsize=fontsize - 1, 
                        horizontalalignment='left', 
                        verticalalignment='bottom',
                        transform=ax.transAxes)

            addobsdata_panel(ax, odata, logmvir_msun, odcolors)

    ylims = [ax.get_ylim() for sub in axes for ax in sub]
    ymin = min([ylim[0] for ylim in ylims])
    ymin = max(ymin, 12.5)
    ymax = max([ylim[1] for ylim in ylims])
    #ymax = ymax + 0.5
    [ax.set_ylim((ymin, ymax)) for sub in axes for ax in sub]
    # 470 instead of 450 so x labels don't overlap
    [ax.set_xlim(0., 470.) for sub in axes for ax in sub]

    handles0, labels0 = axes[0][0].get_legend_handles_labels()
    axes[0][1].legend(handles=handles0[:3], labels=labels0,
                      fontsize=fontsize - 1,
                      loc='lower left', bbox_to_anchor=(1.11, 0.97),
                      handlelength=1.0, labelspacing=0.0,
                      handletextpad=0.4, ncol=3, columnspacing=0.7,
                      frameon=False, borderpad=0.0)
    handles2, labels2 = axes[1][0].get_legend_handles_labels()
    ci_leg = len(plis_k) - 2
    axes[1][ci_leg].legend(handles=handles2[:3], labels=labels2,
                           fontsize=fontsize - 1,
                           loc='lower left', bbox_to_anchor=(0.0, 0.0),
                           handlelength=1.0, labelspacing=0.0,
                           handletextpad=0.4)
    axes[1][ci_leg].axis('off')

    pli_vc_str = 'pli_vc_' + \
                 '_'.join([f'{pli_vc:.2f}' for pli_vc in plis_vc])
    pli_k_str = 'pli_k_' + \
                 '_'.join([f'{pli_k:.2f}' for pli_k in plis_k])
    #zsol_str = 'Zsol_' + \
    #             '_'.join([f'{z_sol:.2f}' for z_sol in z_sols])
    fcgm_str = 'fCGM_' + \
                 '_'.join([f'{fcgm:.2f}' for fcgm in fcgms])
    ods = '_'.join(obsdata)
    outname = (f'prof_Ne8_analytical_pl_s19_z{redshift_model:.2f}'
               f'_vs_{ods}_mvir{logmvir_msun:.1f}'
               f'_{pli_k_str}_{pli_vc_str}_Zsol_0.30_{fcgm_str}'
               '_nhnorm_0p1_to_1_Rvir')
    outname = outname.replace('.', 'p')
    outname = outname.replace('-', 'm')
    plt.savefig(outdir + outname + '.pdf', bbox_inches='tight')


def plot_plmodel_datacomp_parvar_talk(obsdata=('Q+23', 'B+19')):
    ion = 'Ne8'
    redshift_model = 0.75

    impactpars_kpc = np.linspace(5., 450., 50)
    logmvir_msun = 12.2
    fcgms = [1.0, 0.3, 0.1]
    fcgm_def = 0.3
    #z_sols = [1.0, 0.3, 0.1]
    z_sol_def = 0.3
    _colors = tc.tol_cset('bright')
    colors_fZ = [_colors.blue, _colors.green, _colors.yellow]
    #plis_vc = [0.0, -0.2, -0.5]
    pli_vc_def = -0.1 
    #colors_pli_vc = [_colors.cyan, _colors.red, _colors.purple]
    plis_k = [0.0, 2./3., 1.2]
    
    if len(obsdata) == 1:
       _colors = {'1s': 'black', '2s': 'gray'}
       odcolors = {obsdata[0]: _colors}
    else:
        _colors = tc.tol_cset('high-contrast')
        c1 = mcolors.to_rgba(_colors.blue)
        c2 = mcolors.to_rgba(_colors.red)
        odcolors = {'B+19': {'1s': c1, '2s': c1[:3] + (0.5,)},
                    'Q+23': {'1s': c2, '2s': c2[:3] + (0.5,)}}
    odata = get_obsdata(obsdata=(obsdata))
    
    panelsize = 1.8
    ncols = 3
    nrows = 1
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] * nrows
    hspace = 0.4
    height = sum(height_ratios) * (1. + hspace * (nrows - 1) / nrows)
 
    fig = plt.figure(figsize=(sum(width_ratios), height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows, wspace=0.0, 
                        hspace=hspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [[fig.add_subplot(grid[i, j]) 
             for j in range(ncols) if i < 1 or j > 0]
            for i in range(nrows)]
    fontsize = 12
    linewidth = 1.5
    
    for ri in range(nrows):
        for ci in range(ncols):
            ax = axes[ri][ci]
            doleft = ci == 0
            dobottom = True #ri == 1 or (ri == 0 and ci != 1)
            if dobottom:
                ax.set_xlabel('$\\mathrm{r}_{\\perp} \\; [\\mathrm{pkpc}]$',
                              fontsize=fontsize)
            if doleft and ri == 0:
                ax.set_ylabel('$\\log_{10} \\, \\mathrm{N}('
                              '\\mathrm{Ne\\,VIII})'
                              '\\; [\\mathrm{cm}^{-2}]$',
                              fontsize=fontsize)
            ax.tick_params(which='both', direction='in', 
                           labelsize=fontsize - 1.,
                           top=True, right=True, labelleft=doleft,
                           labelbottom=dobottom)
            pli_k = plis_k[ci]
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
                    hmod = mip.PLmodel(10**logmvir_msun, redshift_model, fcgm,
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
            if ci == 0 and ri == 0:
                ax.text(0.0, 1.02, subtitle, fontsize=fontsize - 1, 
                        horizontalalignment='left', 
                        verticalalignment='bottom',
                        transform=ax.transAxes)

            addobsdata_panel(ax, odata, logmvir_msun, odcolors)

    ylims = [ax.get_ylim() for sub in axes for ax in sub]
    ymin = min([ylim[0] for ylim in ylims])
    ymin = max(ymin, 12.5)
    ymax = max([ylim[1] for ylim in ylims])
    #ymax = ymax + 0.5
    [ax.set_ylim((ymin, ymax)) for sub in axes for ax in sub]
    [ax.set_xlim(0., 450.) for sub in axes for ax in sub]

    handles0, labels0 = axes[0][0].get_legend_handles_labels()
    axes[0][1].legend(handles=handles0[:3], labels=labels0,
                      fontsize=fontsize - 1,
                      loc='lower left', bbox_to_anchor=(0.49, 0.97),
                      handlelength=1.0, labelspacing=0.0,
                      handletextpad=0.4, ncol=3, columnspacing=0.7,
                      frameon=False, borderpad=0.0)
    #handles2, labels2 = axes[1][0].get_legend_handles_labels()
    #axes[1][1].legend(handles=handles2[:3], labels=labels2,
    #                  fontsize=fontsize - 1,
    #                  loc='lower left', bbox_to_anchor=(0.0, 0.0),
    #                  handlelength=1.0, labelspacing=0.0,
    #                  handletextpad=0.4)
    #axes[1][1].axis('off')

    #pli_vc_str = 'pli_vc_' + \
    #             '_'.join([f'{pli_vc:.2f}' for pli_vc in plis_vc])
    pli_k_str = 'pli_k_' + \
                 '_'.join([f'{pli_k:.2f}' for pli_k in plis_k])
    #zsol_str = 'Zsol_' + \
    #             '_'.join([f'{z_sol:.2f}' for z_sol in z_sols])
    fcgm_str = 'fCGM_' + \
                 '_'.join([f'{fcgm:.2f}' for fcgm in fcgms])
    ods = '_'.join(obsdata)
    outname = (f'prof_Ne8_analytical_pl_s19_z{redshift_model:.2f}'
               f'_vs_{ods}_mvir{logmvir_msun:.1f}'
               f'_{pli_k_str}_Zsol_0.30_{fcgm_str}')
    outname = outname.replace('.', 'p')
    outname = outname.replace('-', 'm')
    plt.savefig(outdir + outname + '.pdf', bbox_inches='tight')