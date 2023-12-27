
import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import fire_an.analytic_halo.model_ionprof_coolingpl as mic
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


def plot_TnH_vars(outname=None):
    redshift = 0.75
    zsols = [0.1, 0.3, 1.0]
    plinds = [0.0, -0.1, -0.2]
    mdotp_targets = (0.16, 0.5, 0.84)
    ls_zsol = ('dotted', 'solid', 'dashed')
    cl_plind = tc.tol_cset('vibrant')
    lw_mdotp = [1.0, 2.5, 1.7]
    mvirs_logmsun = np.arange(11.0, 13.65, 0.3)
    vmin = mvirs_logmsun[0]
    vmax = mvirs_logmsun[-1]

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