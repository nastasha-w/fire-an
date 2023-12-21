import h5py
import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import fire_an.analytic_halo.model_ionprof_js as mip
import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
import fire_an.utils.constants_and_units as c
import fire_an.utils.opts_locs as ol

def readin_hm_data(filen, tomatch, ion, mdotp_target=0.5):
    with h5py.File(mip.outdir_profiles + filen) as f:
        grpns_all = list(f.keys())
        grpns_use = [grpn for grpn in grpns_all if grpn.startswith(tomatch)]
        toplot = {}
        for grpn in grpns_use:
            grp = f[grpn]
            if bool(grp.attrs['failed']):
                continue
            lmvir = grp.attrs['logMvir_Msun_BN98']
            mdotp = grp.attrs['mdot_percentile_at_Mvir']
            if not np.isclose(mdotp, mdotp_target):
                continue
            _sub = {}
            sgrp = grp[f'coldens_{ion}']
            _sub['b_kpc'] = sgrp['impactpar_cm'][:] / (c.cm_per_mpc * 1e-3)
            _sub['logN_cm2'] = np.log10(sgrp['coldens_cm2'][:])
            _sub['fcgm'] = grp.attrs['fCGM']
            _sub['mdotp'] = mdotp
            _sub['mdot_msunpyr'] = grp.attrs['mdot_MsunperYr']
            toplot[lmvir] = _sub
    return toplot


def plot_profiles_hm(filen, redshift, zsol, plind, outname=None):
    tomatch = (f'z{redshift:.2f}_Zsolar{zsol:.2e}_vcplind{plind:.2f}')
    data = readin_hm_data(filen, tomatch, 'Ne8', mdotp_target=0.5)
    lmvs = np.sort(list(data.keys()))
    print(lmvs)
    vmin = lmvs[0]
    vmax = lmvs[-1]
    carray = [lmvs[0]] + [0.5 * (lmvs[i] + lmvs[i + 1]) 
                            for i in range(len(lmvs) - 1)] + [lmvs[1]]
    carray = np.array(carray)
    cmap = mcm.get_cmap('viridis')

    fig = plt.figure(figsize=(5.5, 5.))
    grid = gsp.GridSpec(ncols=2, nrows=1, wspace=0.1, 
                        width_ratios=(10., 1.))
    ax = fig.add_subplot(grid[0])
    cax = fig.add_subplot(grid[1])
    fontsize = 12
    linewidth = 1.5
    path_effects = pu.getoutline(linewidth)

    for lmv in lmvs:
        xv = data[lmv]['b_kpc']
        yv = data[lmv]['logN_cm2']
        ax.plot(xv, yv, color=cmap((lmv - vmin) / (vmax - vmin)),
                linewidth=linewidth, path_effects=path_effects)
    ax.set_xlabel('$\\mathrm{r}_{\\perp} \\; [\\mathrm{pkpc}]$',
                  fontsize=fontsize)
    ax.set_ylabel('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\,VIII})'
                  '\\; [\\mathrm{cm}^{-2}]$',
                  fontsize=fontsize)
    ax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                   top=True, right=True)
    
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    scm = mcm.ScalarMappable(norm=norm, cmap=cmap)
    scm.set_array(carray)
    plt.colorbar(scm, cax=cax, orientation='vertical')
    cax.set_ylabel('$\\log_{10} \\, \\mathrm{M}_{\\mathrm{vir}}'
                   '\\; [\\mathrm{M}_{\\odot}]$',
                   fontsize=fontsize)
    cax.tick_params(labelsize=fontsize - 1.)
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def check_fcgm_sfr(filen, outname=None):
    redshift = 0.75
    zsol = 0.3
    plind = 0.0
    mdotp_targets = (0.16, 0.5, 0.84)
    ls_mdotp = ('dotted', 'solid', 'dashed')
    tomatch = (f'z{redshift:.2f}_Zsolar{zsol:.2e}_vcplind{plind:.2f}')
    data = {}
    for mdp_tar in mdotp_targets:
        data[mdp_tar] = readin_hm_data(filen, tomatch, 'Ne8', 
                                       mdotp_target=mdp_tar)

    fig = plt.figure(figsize=(7.5, 3.5))
    grid = gsp.GridSpec(ncols=2, nrows=1, wspace=0.35, 
                        width_ratios=(1., 1.))
    fcgmax = fig.add_subplot(grid[0])
    sfrax = fig.add_subplot(grid[1])
    fontsize = 12
    linewidth = 1.5

    for mdp, ls in zip(mdotp_targets, ls_mdotp):
        label = f'{mdp*100:.0f}$^{{\\mathrm{{th}}}}$% SFR'
        xvs = []
        fcgms = []
        sfrs = []
        print(mdp)
        print(list(data[mdp].keys()))
        for lmv in data[mdp].keys():
            xvs.append(lmv)
            fcgm = data[mdp][lmv]['fcgm'] \
                * mip.cosmo_base['omegam'] / mip.cosmo_base['omegab']
            lsfr = np.log10(data[mdp][lmv]['mdot_msunpyr']) 
            fcgms.append(fcgm)
            sfrs.append(lsfr)
        xvs = np.array(xvs)
        po = np.argsort(xvs)
        fcgmax.plot(xvs[po], np.array(fcgms)[po], color='black',
                    linewidth=linewidth, linestyle=ls,
                    label=label)
        sfrax.plot(xvs[po], np.array(sfrs)[po], color='black',
                   linewidth=linewidth, linestyle=ls,
                   label=label)
    
    xlabel = ('$\\log_{10} \\, \\mathrm{M}_{\\mathrm{vir, BN98}}'
              ' \\; [\\mathrm{M}_{\\odot}]$')
    fcgmax.set_xlabel(xlabel, fontsize=fontsize)
    sfrax.set_xlabel(xlabel, fontsize=fontsize)
    fcgmax.set_ylabel('$\\mathrm{f}_{\\mathrm{CGM}}$',
                      fontsize=fontsize)
    sfrax.set_ylabel('$\\log_{10} \\, \\mathrm{SFR} \\; '
                     '[\\mathrm{M}_{\\odot} \\mathrm{yr}^{-1}]$',
                      fontsize=fontsize)
    fcgmax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                       top=True, right=True)
    sfrax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                      top=True, right=True)
    sfrax.legend(fontsize=fontsize - 1.)
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')


def check_fcgm_sfr_vars(filen, outname=None):
    redshift = 0.75
    zsols = [0.1, 0.3, 1.0]
    plinds = [0.0, -0.1, -0.2]
    mdotp_targets = (0.16, 0.5, 0.84)
    ls_zsol = ('dotted', 'solid', 'dashed')
    cl_plind = tc.tol_cset('vibrant')
    lw_mdotp = [1.0, 2.0, 1.4]

    fig = plt.figure(figsize=(7.5, 3.5))
    grid = gsp.GridSpec(ncols=2, nrows=1, wspace=0.35, 
                        width_ratios=(1., 1.))
    fcgmax = fig.add_subplot(grid[0])
    sfrax = fig.add_subplot(grid[1])
    fontsize = 12
    
    for zsol, ls in zip(zsols, ls_zsol):
        for plind, cl in zip(plinds, cl_plind):
            for mdp_tar, lw in zip(mdotp_targets, lw_mdotp):
                tomatch = (f'z{redshift:.2f}_Zsolar{zsol:.2e}'
                           f'_vcplind{plind:.2f}_')
                data = readin_hm_data(filen, tomatch, 'Ne8', 
                                      mdotp_target=mdp_tar)

                xvs = []
                fcgms = []
                sfrs = []
                for lmv in data.keys():
                    xvs.append(lmv)
                    fcgm = data[lmv]['fcgm'] \
                        * mip.cosmo_base['omegam'] / mip.cosmo_base['omegab']
                    lsfr = np.log10(data[lmv]['mdot_msunpyr']) 
                    fcgms.append(fcgm)
                    sfrs.append(lsfr)
                xvs = np.array(xvs)
                po = np.argsort(xvs)
                fcgmax.plot(xvs[po], np.array(fcgms)[po], color=cl,
                            linewidth=lw, linestyle=ls)
                sfrax.plot(xvs[po], np.array(sfrs)[po], color=cl,
                           linewidth=lw, linestyle=ls)
    
    xlabel = ('$\\log_{10} \\, \\mathrm{M}_{\\mathrm{vir, BN98}}'
              ' \\; [\\mathrm{M}_{\\odot}]$')
    fcgmax.set_xlabel(xlabel, fontsize=fontsize)
    sfrax.set_xlabel(xlabel, fontsize=fontsize)
    fcgmax.set_ylabel('$\\mathrm{f}_{\\mathrm{CGM}}$',
                      fontsize=fontsize)
    sfrax.set_ylabel('$\\log_{10} \\, \\mathrm{SFR} \\; '
                     '[\\mathrm{M}_{\\odot} \\mathrm{yr}^{-1}]$',
                      fontsize=fontsize)
    fcgmax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                       top=True, right=True)
    sfrax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                      top=True, right=True)
    sfrax.legend(fontsize=fontsize - 1.)

    mdp_han = [mlines.Line2D((), (), color='black', linestyle='solid',
                             linewidth=lw,
                             label=(f'{mdp*100:.0f}$^{{\\mathrm{{th}}}}$'
                                     '% SFR'))
               for lw, mdp in zip(lw_mdotp, mdotp_targets)]
    zsol_han = [mlines.Line2D((), (), color='black', linestyle=ls,
                              linewidth=1.5,
                              label=(f'{zsol:.1f} $\\mathrm{{Z}}'
                                     '_{\\odot}$'))
                for ls, zsol in zip(ls_zsol, zsols)]
    pli_han = [mlines.Line2D((), (), color=cl, linestyle='solid',
                             linewidth=1.5,
                             label=('$v_{\\mathrm{c}} \\propto '
                                    f'r^{{{plind:.1f}}}$'))
                for cl, plind in zip(cl_plind, plinds)]
    fcgmax.legend(handles=mdp_han + zsol_han, fontsize=fontsize - 1)
    sfrax.legend(handles=pli_han, fontsize=fontsize - 1)
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')


def plot_coldensprof_vars(filen, outname=None):
    lmvs_tar = np.arange(11.0, 13.65, 0.2)

    redshift = 0.75
    zsols = [0.1, 0.3, 1.0]
    plinds = [0.0, -0.1, -0.2]
    mdotp_targets = (0.16, 0.5, 0.84)
    ls_zsol = ('dotted', 'solid', 'dashed')
    cl_plind = tc.tol_cset('vibrant')
    lw_mdotp = [1.0, 2.0, 1.4]
    
    npanels = len(lmvs_tar)
    panelsize = 2.5
    ncols_max = 4
    ncols = min(ncols_max, npanels)
    nrows = (npanels - 1) // ncols + 1
    hspace = 0.2
    wspace = 0.2
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] * nrows
    width = sum(width_ratios) * (1. + wspace * (ncols - 1.) / ncols)
    height = sum(height_ratios) * (1. + hspace * (nrows - 1.) / nrows)

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows, wspace=wspace, 
                        hspace=hspace,
                        width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols]) 
            for i in range(npanels)]
    fontsize = 12
    
    for zsol, ls in zip(zsols, ls_zsol):
        for plind, cl in zip(plinds, cl_plind):
            for mdp_tar, lw in zip(mdotp_targets, lw_mdotp):
                tomatch = (f'z{redshift:.2f}_Zsolar{zsol:.2e}'
                           f'_vcplind{plind:.2f}_')
                data = readin_hm_data(filen, tomatch, 'Ne8', 
                                      mdotp_target=mdp_tar)
                lmvs = [key for key in data 
                        if np.any(np.isclose(key, lmvs_tar))]
                for lmv in lmvs:
                    axi = np.where(np.isclose(lmv, lmvs_tar))[0][0]
                    xv = data[lmv]['b_kpc']
                    yv = data[lmv]['logN_cm2']
                    axes[axi].plot(xv, yv, color=cl, linewidth=lw, 
                                   linestyle=ls)
    for axi in range(npanels):
        ax = axes[axi]
        doleft = axi % ncols == 0
        dobottom = axi >= npanels - ncols
        ax.tick_params(which='both', direction='in', right=True,
                       top=True)
        if doleft:
            ax.set_ylabel('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\,VIII})'
                          '\\; [\\mathrm{cm}^{-2}]$', fontsize=fontsize)
        if dobottom:
            ax.set_xlabel('$\\mathrm{r}_{\\perp} \\; [\\mathrm{pkpc}]$',
                          fontsize=fontsize)
        axtitle = ('$\\mathrm{M}_{\\mathrm{vir}} ='
                   f'10^{{{lmvs_tar[axi]:.1f}}}'
                   '\\, \\mathrm{M}_{\\odot}$')
        ax.text(0.95, 0.95, axtitle, fontsize=fontsize, 
                transform=ax.transAxes, horizontalalignment='right',
                verticalalignment='top')
    
    mdp_han = [mlines.Line2D((), (), color='black', linestyle='solid',
                             linewidth=lw,
                             label=(f'{mdp*100:.0f}$^{{\\mathrm{{th}}}}$'
                                     '% SFR'))
               for lw, mdp in zip(lw_mdotp, mdotp_targets)]
    zsol_han = [mlines.Line2D((), (), color='black', linestyle=ls,
                              linewidth=1.5,
                              label=(f'{zsol:.1f} $\\mathrm{{Z}}'
                                     '_{\\odot}$'))
                for ls, zsol in zip(ls_zsol, zsols)]
    pli_han = [mlines.Line2D((), (), color=cl, linestyle='solid',
                             linewidth=1.5,
                             label=('$v_{\\mathrm{c}} \\propto '
                                    f'r^{{{plind:.1f}}}$'))
                for cl, plind in zip(cl_plind, plinds)]
    if npanels >= 3:
        axes[0].legend(handles=mdp_han, fontsize=fontsize - 1)
        axes[1].legend(handles=zsol_han, fontsize=fontsize - 1)
        axes[2].legend(handles=pli_han, fontsize=fontsize - 1)
    elif npanels == 2:
        axes[0].legend(handles=mdp_han + zsol_han, 
                       fontsize=fontsize - 1)
        axes[1].legend(handles=pli_han, fontsize=fontsize - 1)
    elif npanels == 1:
        axes[0].legend(handles=mdp_han + zsol_han + pli_han, 
                       fontsize=fontsize - 1)
    
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

        

                




    