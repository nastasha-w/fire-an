
from re import M
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
import matplotlib.gridspec as gsp
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fire_an.analytic_halo.model_ionprof_pl as mip
import fire_an.makeplots.plot_utils as pu
import fire_an.mstar_mhalo.analytical as msmhan
import fire_an.mstar_mhalo.loader_smdpl_sfr as ldsmdpl
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu
import fire_an.utils.math_utils as mu


outdir = '/projects/b1026/nastasha/imgs/analytical/'
oddir = '/projects/b1026/nastasha/extdata/'
ofilen = oddir + 'data_burchett_etal_2019_table1.txt'


def calcmhalodist(logmstar_msun, sigmalogmstar, redshift):
    histobj = ldsmdpl.SMHMhists(np.array([redshift]), binsize=0.1)
    msbins_hist = histobj.getbins(redshift, mode='ms')
    _Pcbin_ms = msmhan.cumulgauss((msbins_hist - logmstar_msun) 
                                  / sigmalogmstar)
    _Pbin_ms = np.diff(_Pcbin_ms)
    mhp, _mhbins = histobj.matrixconv(_Pbin_ms, redshift, 
                                      mode='mstomh')    
    return _mhbins, mhp

def readdata_b19(nsigma=2):
    data_bur = pd.read_csv(ofilen, comment='#', sep='\t')
    #cosmopars_bur = {'h': 0.677, 'omegam': 0.31, 'omegalambda': 0.69}
    sig2t = msmhan.cumulgauss(nsigma) - msmhan.cumulgauss(-nsigma)
    cumulP_lo = 0.5 - 0.5 * sig2t
    cumulP_hi = 0.5 + 0.5 * sig2t

    logmstars_msun = data_bur['log_Mstar_Msun']
    sigmalogmstars = data_bur['log_Mstar_Msun_err']
    #isul = data_bur['log_N_Ne8_isUL']
    reshifts = data_bur['zgal']

    logmvir_msun_bestest = []
    logmvir_msun_lo = []
    logmvir_msun_hi = []
    for lgmstar_msun, slgmstar, _z in zip(logmstars_msun, sigmalogmstars,
                                          reshifts):
        mhbins, mhP = calcmhalodist(lgmstar_msun, slgmstar, _z)
        mhpd = mhP / np.diff(mhbins)
        mhcens = 0.5 * (mhbins[:-1] + mhbins[1:])
        logmvir_msun_bestest.append(mhcens[np.argmax(mhpd)])
        mlo = mu.linterpsolve(np.cumsum(mhP), mhbins[1:], cumulP_lo)
        logmvir_msun_lo.append(mlo)
        mhi = mu.linterpsolve(np.cumsum(mhP), mhbins[1:], cumulP_hi)
        logmvir_msun_hi.append(mhi)
    logmvir_msun_bestest = np.array(logmvir_msun_bestest)
    logmvir_msun_lo = np.array(logmvir_msun_lo)
    logmvir_msun_hi = np.array(logmvir_msun_hi)

    data_bur = data_bur.assign(logmvir_msun_bestest=logmvir_msun_bestest,
                               logmvir_msun_lo=logmvir_msun_lo,
                               logmvir_msun_hi=logmvir_msun_hi)
    return data_bur
    

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

    data_bur = readdata_b19(nsigma=nsigma)

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