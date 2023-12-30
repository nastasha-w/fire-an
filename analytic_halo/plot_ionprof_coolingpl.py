
import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import fire_an.analytic_halo.model_ionprof_coolingpl as mic
import fire_an.ionrad.ion_utils as iu
import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
import fire_an.mstar_mhalo.loader_smdpl_sfr as ldsmdpl
import fire_an.utils.constants_and_units as c
import fire_an.utils.math_utils as mu

def get_sfrs_mh(logmhs_msun, z=0.75, percentiles=(0.16, 0.5, 0.84)):
    histobj = ldsmdpl.SFRHMhists(np.array([z]))
    sfrs_tab, mhs_tab = histobj.getperc_sfrmh(z, mode='mhtosfr', 
                                              percvals=np.array(percentiles))
    out = [[mu.linterpsolve(mhs_tab, sfrs_tab_perc, mh_this) 
            for sfrs_tab_perc in sfrs_tab]
           for mh_this in logmhs_msun]
    return np.array(out)


def check_fcgm_sfr_vars(outname=None):
    redshift = 0.75
    zsols = [0.1, 0.3, 1.0]
    plinds = [0.0, -0.1, -0.2]
    mdotp_targets = (0.16, 0.5, 0.84)
    ls_zsol = ('dotted', 'solid', 'dashed')
    cl_plind = tc.tol_cset('vibrant')
    lw_mdotp = [1.0, 2.5, 1.7]
    mvirs_logmsun = np.arange(11.0, 13.65, 0.1)

    fig = plt.figure(figsize=(11.5, 3.5))
    grid = gsp.GridSpec(ncols=3, nrows=1, wspace=0.40, 
                        width_ratios=(1., 1., 1.))
    fcgmax = fig.add_subplot(grid[0])
    sfrax = fig.add_subplot(grid[1])
    coolax = fig.add_subplot(grid[2])
    fontsize = 12
    
    for zsol, ls in zip(zsols, ls_zsol):
        for plind, cl in zip(plinds, cl_plind):
            for mdp_tar, lw in zip(mdotp_targets, lw_mdotp):
                sfrs_logmsun = get_sfrs_mh(mvirs_logmsun, z=redshift, 
                                            percentiles=(mdp_tar,))
                models_cur = {mv: mic.PLmodel(mv, redshift, sf, 
                                              zsol, plind)
                              for mv, sf in zip(10**mvirs_logmsun, 
                                                10**sfrs_logmsun[:, 0])}

                fcgms = []
                coolrates = []
                sfrs = sfrs_logmsun
                for mv in models_cur.keys():
                    model = models_cur[mv]
                    redges = np.linspace(0.1 * model.rvir_cgs, model.rvir_cgs,
                                         1000)
                    rcens = 0.5 * (redges[:-1] + redges[1:])
                    dr = np.diff(redges)
                    nHs = model.nH_cm3(rcens)
                    fcgm = np.sum(4. * np.pi * nHs * rcens**2 * dr)
                    fcgm = (fcgm / model.hmassfrac) * c.atomw_H * c.u \
                           / c.solar_mass
                    fcgm = fcgm / mv * mic.cosmopars_base_fire['omegam'] \
                           / mic.cosmopars_base_fire['omegab']
                    fcgms.append(fcgm)
                    coolrates.append(np.log10(model._Lambda_cgs))
                xvs = mvirs_logmsun
                po = np.argsort(xvs)
                fcgmax.plot(xvs[po], np.array(fcgms)[po], color=cl,
                            linewidth=lw, linestyle=ls)
                sfrax.plot(xvs[po], np.array(sfrs)[po], color=cl,
                           linewidth=lw, linestyle=ls)
                coolax.plot(xvs[po], np.array(coolrates)[po], color=cl,
                           linewidth=lw, linestyle=ls)
    
    xlabel = ('$\\log_{10} \\, \\mathrm{M}_{\\mathrm{vir, BN98}}'
              ' \\; [\\mathrm{M}_{\\odot}]$')
    fcgmax.set_xlabel(xlabel, fontsize=fontsize)
    sfrax.set_xlabel(xlabel, fontsize=fontsize)
    coolax.set_xlabel(xlabel, fontsize=fontsize)
    fcgmax.set_ylabel('$\\mathrm{f}_{\\mathrm{CGM}}$',
                      fontsize=fontsize)
    sfrax.set_ylabel('$\\log_{10} \\, \\mathrm{SFR} \\; '
                     '[\\mathrm{M}_{\\odot} \\mathrm{yr}^{-1}]$',
                      fontsize=fontsize)
    coolax.set_ylabel('$\\log_{10} \\, \\Lambda \\; '
                      '[\\mathrm{erg}\\,\\mathrm{cm}^{3}\\mathrm{s}^{-1}]$',
                      fontsize=fontsize)
    fcgmax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                       top=True, right=True)
    sfrax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                      top=True, right=True)
    coolax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                       top=True, right=True)

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
    fcgmax.legend(handles=mdp_han, fontsize=fontsize - 1,
                  handlelength=1.2, labelspacing=0.3, ncol=1,
                  handletextpad=0.4)
    coolax.legend(handles=zsol_han, fontsize=fontsize - 1,
                  handlelength=1.4, labelspacing=0.3, ncol=1,
                  handletextpad=0.4)
    sfrax.legend(handles=pli_han, fontsize=fontsize - 1,
                 handlelength=1., labelspacing=0.2,
                 handletextpad=0.4)
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')


def plot_TnH_vars(outname=None, xrange=(-6.5, -2.), yrange=(5.2, 7.5)):
    redshift = 0.75
    zsols = [0.1, 0.3, 1.0]
    plinds = [0.0, -0.1, -0.2]
    mdotp_targets = (0.16, 0.5, 0.84)
    ls_zsol = ('dotted', 'solid', 'dashed')
    cl_mdotp = tc.tol_cset('vibrant')
    lw_plind = [1.0, 2.5, 1.7]
    mvirs_logmsun = np.arange(11.0, 13.65, 0.3)
    
    iontab = iu.Linetable_PS20('Ne8', redshift, 
                               emission=False, lintable=True)
    iontab.findiontable()
    iZ = np.argmin(np.abs(iontab.logZsol - zsols[1]))
    table_T_nH = iontab.iontable_T_Z_nH[:, iZ, :]
    table_lTk = iontab.logTK
    table_lncc = iontab.lognHcm3
    table_extent = (table_lncc[0], table_lncc[-1],
                    table_lTk[0], table_lTk[-1])
    _cmap = mcm.get_cmap('gist_yarg')
    cmap = pu.truncate_colormap(_cmap, minval=0., maxval=0.7)
    vmax = np.log10(np.max(table_T_nH ))
    vmin = vmax - 5.

    panelsize = 2.5
    ncols_max = 4
    caxwidth = 0.5
    npanels = len(mvirs_logmsun)
    ncols = min(npanels, ncols_max)
    nrows = (npanels - 1) // ncols + 1
    width_ratios = [panelsize] * ncols + [caxwidth]
    height_ratios = [panelsize] * nrows
    wspace = 0.
    hspace = 0.
    width = sum(width_ratios) * (1. + (ncols - 1.) / ncols * wspace)
    height = sum(height_ratios) * (1. + (nrows - 1.) / nrows * hspace)
    
    fig = plt.figure(figsize=(width, height)) 
    grid = gsp.GridSpec(ncols=ncols + 1, nrows=nrows, wspace=wspace,
                        hspace=hspace, height_ratios=height_ratios, 
                        width_ratios=width_ratios)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols]) 
            for i in range(npanels)]
    cax = fig.add_subplot(grid[slice(None, min(1, nrows), None), ncols])
    fontsize = 12
    
    for axi, lmv in enumerate(mvirs_logmsun):
        ax = axes[axi]
        dobottom = axi >= npanels - ncols
        doleft = axi % ncols == 0
        title = ('$\\mathrm{M}_{\\mathrm{vir}} \\,'
                 f'= 10^{{{lmv:.1f}}} \\mathrm{{M}}_{{\\odot}}$')
        ax.text(0.95, 0.95, title, fontsize=fontsize,
                horizontalalignment='right', verticalalignment='top',
                transform=ax.transAxes)
        ax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                       top=True, right=True, labelbottom=dobottom,
                       labelleft=doleft)
        if dobottom:
            ax.set_xlabel('$\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
                          '\\; [\\mathrm{cm}^{-3}]$', fontsize=fontsize)
        if doleft:
            ax.set_ylabel('$\\log_{10} \\, \\mathrm{T}'
                          '\\; [\\mathrm{K}]$', fontsize=fontsize)
        img = ax.imshow(np.log10(table_T_nH), extent=table_extent, 
                        origin='lower',
                        interpolation='nearest', cmap=cmap, vmin=vmin,
                        vmax=vmax, aspect='auto')

        for zsol, ls in zip(zsols, ls_zsol):
            for plind, lw in zip(plinds, lw_plind):
                for mdp_tar, cl in zip(mdotp_targets, cl_mdotp):
                    sfr_logmsun = get_sfrs_mh([lmv], z=redshift, 
                                              percentiles=(mdp_tar,))[0][0]
                    model = mic.PLmodel(10**lmv, redshift, 10**sfr_logmsun, 
                                        zsol, plind)
                    rvs = np.linspace(0.1 * model.rvir_cgs, model.rvir_cgs,
                                      20)
                    xvs = np.log10(model.nH_cm3(rvs))
                    yvs = np.log10(model.t_K(rvs))

                    ax.plot(xvs, yvs, color=cl, linewidth=lw, linestyle=ls)
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
    plt.colorbar(img, cax=cax, extend='min', aspect=1. / 15.,
                 orientation='vertical')
    cax.set_ylabel('$\\log_{10}\\, \\mathrm{Ne\\,VIII} \\,/\\, \\mathrm{Ne}$',
                  fontsize=fontsize)
    cax.tick_params(which='both', labelsize=fontsize - 1.)
    
    
    mdp_han = [mlines.Line2D((), (), color=cl, linestyle='solid',
                             linewidth=2.,
                             label=(f'{mdp*100:.0f}$^{{\\mathrm{{th}}}}$'
                                     '% SFR'))
               for cl, mdp in zip(cl_mdotp, mdotp_targets)]
    zsol_han = [mlines.Line2D((), (), color='black', linestyle=ls,
                              linewidth=1.5,
                              label=(f'{zsol:.1f} $\\mathrm{{Z}}'
                                     '_{\\odot}$'))
                for ls, zsol in zip(ls_zsol, zsols)]
    pli_han = [mlines.Line2D((), (), color='black', linestyle='solid',
                             linewidth=lw,
                             label=('$v_{\\mathrm{c}} \\propto '
                                    f'r^{{{plind:.1f}}}$'))
                for lw, plind in zip(lw_plind, plinds)]
    axes[0].legend(handles=mdp_han, fontsize=fontsize - 1,
                   handlelength=1.2, labelspacing=0.3, ncol=1,
                   handletextpad=0.4, loc='upper left',
                   bbox_to_anchor=(0.05, 0.85))
    axes[1].legend(handles=zsol_han, fontsize=fontsize - 1,
                   handlelength=1.4, labelspacing=0.3, ncol=1,
                   handletextpad=0.4, loc='upper left',
                   bbox_to_anchor=(0.05, 0.85))
    axes[2].legend(handles=pli_han, fontsize=fontsize - 1,
                   handlelength=1., labelspacing=0.2,
                   handletextpad=0.4, loc='upper left',
                   bbox_to_anchor=(0.05, 0.85))
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plot_TnHionrac_prof(outname=None, logmvir_msun=11.3):
    redshift = 0.75
    zsols = [0.1, 0.3, 1.0]
    plinds = [0.0, -0.1, -0.2]
    mdotp_targets = (0.16, 0.5, 0.84)
    ls_zsol = ('dotted', 'solid', 'dashed')
    cl_mdotp = tc.tol_cset('vibrant')
    
    iontab = iu.Linetable_PS20('Ne8', redshift, 
                               emission=False, lintable=True)
    iontab.findiontable()

    panelsize = 2.5
    ncols_max = 3
    caxwidth = 0.5
    npanels = 9
    ncols = min(npanels, ncols_max)
    nrows = (npanels - 1) // ncols + 1
    width_ratios = [panelsize] * ncols + [caxwidth]
    height_ratios = [panelsize] * nrows
    wspace = 0.
    hspace = 0.
    width = sum(width_ratios) * (1. + (ncols - 1.) / ncols * wspace)
    height = sum(height_ratios) * (1. + (nrows - 1.) / nrows * hspace)
    
    fig = plt.figure(figsize=(width, height)) 
    grid = gsp.GridSpec(ncols=ncols + 1, nrows=nrows, wspace=wspace,
                        hspace=hspace, height_ratios=height_ratios, 
                        width_ratios=width_ratios)
    axes = [[fig.add_subplot(grid[i, j])
             for j in range(ncols)]
            for i in range(nrows)]
            
    fontsize = 12

    title = ('$\\mathrm{M}_{\\mathrm{vir}} \\,'
             f'= 10^{{{logmvir_msun:.1f}}} \\mathrm{{M}}_{{\\odot}}$')
    fig.suptitle(title, fontsize=fontsize)
    
    for coli, pli in enumerate(plinds):
        for rowi, pltv in enumerate(['T', 'nH', 'ionfrac']):
            ax = axes[rowi][coli]
            dobottom = rowi == ncols - 1
            doleft = coli == 0
            ax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                        top=True, right=True, labelbottom=dobottom,
                        labelleft=doleft)
            if rowi == 0:
                axtitle = f'$v_{{\\mathrm{{c}}}} \\propto r^{{{pli:.1f}}}$'
                ax.set_title(axtitle, fontsize=fontsize)
            if dobottom:
                ax.set_xlabel('$\ \\mathrm{r}_{\\mathrm{3D}}'
                            '\\; [\\mathrm{R}_{\\mathrm{vir}}]$',
                            fontsize=fontsize)
            if doleft:
                if pltv == 'T':
                    label = ('$\\log_{10} \\, \\mathrm{T}'
                             '\\; [\\mathrm{K}]$')
                elif pltv == 'nH':
                    label = ('$\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
                             '\\; [\\mathrm{cm}^{-3}]$')
                elif pltv == 'ionfrac':
                    label = ('$\\log_{10} \\, \\mathrm{Ne\\,VIII}'
                             '\\,/\\, \\mathrm{Ne}$')
                ax.set_ylabel(label, fontsize=fontsize)
    
    for coli, pli in enumerate(plinds):
        for zsol, ls in zip(zsols, ls_zsol):
            for mdp_tar, cl in zip(mdotp_targets, cl_mdotp):
                sfr_logmsun = get_sfrs_mh([logmvir_msun], z=redshift, 
                                              percentiles=(mdp_tar,))[0][0]
                model = mic.PLmodel(10**logmvir_msun, redshift, 
                                    10**sfr_logmsun, 
                                    zsol, pli)
                rvs = np.linspace(0.1 * model.rvir_cgs, model.rvir_cgs,
                                20)
                dens = np.log10(model.nH_cm3(rvs))
                temp = np.log10(model.t_K(rvs))
                zvs = zsol * iontab.solarZ * np.ones(dens.shape)
                dct_T_Z_nH = {'logT': temp,
                              'logZ': np.log10(zvs),
                              'lognH': dens}
                ionf = iontab.find_ionbal(dct_T_Z_nH, log=True)
                for rowi, pltv in enumerate(['T', 'nH', 'ionfrac']):
                    if pltv == 'T':
                        axes[rowi][coli].plot(rvs / model.rvir_cgs, 
                                              temp, color=cl, 
                                              linewidth=1.5, linestyle=ls)
                    elif pltv == 'nH':
                        axes[rowi][coli].plot(rvs / model.rvir_cgs, 
                                              dens, color=cl, 
                                              linewidth=1.5, linestyle=ls)
                    elif pltv == 'ionfrac':
                        axes[rowi][coli].plot(rvs / model.rvir_cgs, 
                                              ionf, color=cl, 
                                              linewidth=1.5, linestyle=ls)

    for rowi in range(nrows):
        axrow = axes[rowi]
        ylims = [axrow[coli].get_ylim() for coli in range(ncols)]
        ymin = min([y[0] for y in ylims])
        ymax = max([y[-1] for y in ylims])
        [axrow[coli].set_ylim(ymin, ymax) for coli in range(ncols)]
    
    mdp_han = [mlines.Line2D((), (), color=cl, linestyle='solid',
                             linewidth=2.,
                             label=(f'{mdp*100:.0f}$^{{\\mathrm{{th}}}}$'
                                     '% SFR'))
               for cl, mdp in zip(cl_mdotp, mdotp_targets)]
    zsol_han = [mlines.Line2D((), (), color='black', linestyle=ls,
                              linewidth=1.5,
                              label=(f'{zsol:.1f} $\\mathrm{{Z}}'
                                     '_{\\odot}$'))
                for ls, zsol in zip(ls_zsol, zsols)]
    #pli_han = [mlines.Line2D((), (), color='black', linestyle='solid',
    #                         linewidth=lw,
    #                         label=('$v_{\\mathrm{c}} \\propto '
    #                                f'r^{{{plind:.1f}}}$'))
    #            for lw, plind in zip(lw_plind, plinds)]
    axes[0][0].legend(handles=mdp_han, fontsize=fontsize - 1,
                   handlelength=1.2, labelspacing=0.3, ncol=1,
                   handletextpad=0.4, loc='lower left',
                   bbox_to_anchor=(0.05, 0.05))
    axes[0][1].legend(handles=zsol_han, fontsize=fontsize - 1,
                   handlelength=1.4, labelspacing=0.3, ncol=1,
                   handletextpad=0.4, loc='lower left',
                   bbox_to_anchor=(0.05, 0.05))
    #axes[2].legend(handles=pli_han, fontsize=fontsize - 1,
    #               handlelength=1., labelspacing=0.2,
    #               handletextpad=0.4, loc='upper left',
    #               bbox_to_anchor=(0.05, 0.85))
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

