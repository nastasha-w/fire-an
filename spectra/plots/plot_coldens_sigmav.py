
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
import fire_an.spectra.plots.select_obsdata as so

# TODO: actually do the data comparison for CASBaH
# Zhijie: CUBS sigmas are rel. to tau-weighted velocity
def datacomp(simfilen: str,
             outfilen: str | None = None):
    fmtsim = {'color': 'gray',
              'marker': 'o',
              's': 3.,
              'alpha': 0.5}
    cs = tc.tol_cset('bright')
    fmtcubs = {'color': cs[0],
               'marker': 's',
               'markersize': 5.,
               'alpha': 1.}
    fmtcasbah = {'color': cs[1],
                 'marker': 'p',
                 'markersize': 5.,
                 'alpha': 1.}
    simdata = pd.read_csv(simfilen, sep='\t')
    fig = plt.figure(figsize=(5., 5.))
    ax = fig.add_subplot()
    fontsize = 12
    
    ax.set_xlabel('$\\log_{10} \\, \\mathrm{N} \\; [\\mathrm{cm}^{-2}]$',
                  fontsize=fontsize)
    ax.set_ylabel('$\\sigma\\left(v - \\langle v \\rangle_{\\tau}\\right) \\; '
                  '[\\mathrm{km}/\\mathrm{s}]$',
                  fontsize=fontsize)
    ax.tick_params(direction='in', which='both', labelsize=fontsize - 1.,
                   left=True, top=True)
    xp = simdata['Ntot_logcm2']
    yp = simdata['sigmav_kmps']
    ax.scatter(xp, yp,
               label='FIRE-2', rasterized=True, **fmtsim)
    binsize = 0.1
    xlo = np.floor(np.min(xp) / binsize) * binsize
    xhi = np.ceil(np.max(xp) / binsize) * binsize
    xbins = np.arange(xlo, xhi + 0.5 * binsize, binsize)
    xcen = 0.5 * (xbins[:-1] + xbins[1:])
    yperc, _, xsel = pu.get_perc_and_points(xp, yp, xbins,
                           percentiles=(10., 50., 90.),
                           mincount_x=20,
                           getoutliers_y=False, getmincounts_x=False,
                           x_extremes_only=False)
    ax.plot(xcen[xsel], (yperc[1:-1, 0])[xsel],
            color='black', linestyle='dashed',
            linewidth=1.2, label='FIRE-2: 10, 90%')
    ax.plot(xcen[xsel], (yperc[1:-1, 2])[xsel], 
            color='black', linestyle='dashed',
            linewidth=1.2, label=None)
    ax.plot(xcen[xsel], (yperc[1:-1, 1])[xsel], 
            color='black', linestyle='solid',
            linewidth=1.5, label='FIRE-2: median')
    
    cubssys = so.selectsystems(simdata['Mvir_logMsun'],
                               np.array([0.5, 1.0]),
                               'cubs',
                               Mhmargin_dex=0.2,
                               zmargin=0.05,
                               method='bestest_inrange')
    isul_sigma = np.isnan(cubssys['sigmav [km/s]'])
    noul_sigma = np.logical_not(isul_sigma)
    cubssys_uls = cubssys[isul_sigma]
    cubssys_dets = cubssys[noul_sigma]
    ax.errorbar(cubssys_dets['Nline [logcm2]'],
                cubssys_dets['sigmav [km/s]'],
                xerr=(cubssys_dets['err_Nline_minus [dex]'],
                      cubssys_dets['err_Nline_plus [dex]']),
                yerr=(cubssys_dets['err_sigmav_minus [km/s]'],
                      cubssys_dets['err_sigmav_plus [km/s]']),
                linestyle='none',
                label='CUBS',
                **fmtcubs)
    ax.errorbar(cubssys_uls['Nline [logcm2]'],
                cubssys_uls['err_sigmav_plus [km/s]'],
                xerr=(cubssys_uls['err_Nline_minus [dex]'],
                      cubssys_uls['err_Nline_plus [dex]']),
                linestyle='none',
                **fmtcubs)
    for i in cubssys_uls.index:
        # arrow head width/length seem to be in data units
        ax.arrow(cubssys_uls.at[i, 'Nline [logcm2]'],
                cubssys_uls.at[i, 'err_sigmav_plus [km/s]'],
                0, -15., head_width=0.05, head_length=3.,
                color=fmtcubs['color'], alpha=fmtcubs['alpha'])
    ax.set_xlim(12.4, 14.8)
    ax.set_ylim(0., 250.)
    ax.legend(fontsize=fontsize - 1., handlelength=1.5)

    if outfilen is not None:
        plt.savefig(outfilen, bbox_inches='tight')


def plotpanel_colby(ax, cax, df, cqty, fontsize=12):
    if cqty == 'Mvir':
        ccol = 'Mvir_logMsun'
        clabel = ('$\\log_{10} \\, \\mathrm{M}_{\\mathrm{vir}}'
                  '\\; [\\mathrm{M}_{\odot}]$')
    elif cqty == 'Mstar':
        ccol = 'Mstar_logMsun'
        clabel = ('$\\log_{10} \\, \\mathrm{M}_{\\star}'
                  '\\; [\\mathrm{M}_{\odot}]$')
    elif cqty == 'ipar':
        ccol = 'impactpar_kpc'
        clabel = ('$\mathrm{r}_{\\perp}'
                  '\\; [\\mathrm{kpc}]$')
    elif cqty == 'sigmav':
        ccol = 'sigmav_kmps'
        clabel = ('$\\sigma(v - \\langle v \\rangle_{\\tau})'
                  '\\; [\\mathrm{km} / \\mathrm{s}]$')
    elif cqty == 'sigmav_gal':
        ccol = 'sigmav_galv_kmps'
        clabel = ('$\\sigma(v - v_{\\mathrm{gal}})'
                  '\\; [\\mathrm{km} / \\mathrm{s}]$')
    cmapn = 'rainbow'
    cdat = df[ccol]
    vmin = np.min(cdat)
    vmax = np.max(cdat)
    cmap = pu.paste_cmaps([cmapn], edges=[vmin, vmax])
    cp = cmap((cdat - vmin) / (vmax - vmin))
    pu.add_colorbar(cax, vmin=vmin, vmax=vmax, cmap=cmap,
                    clabel=clabel, fontsize=fontsize,
                    extend='neither', orientation='vertical')

    xp = df['Ntot_logcm2']
    yp = df['sigmav_kmps']
    ax.scatter(xp, yp, c=cp, marker='o', s=3., alpha=0.5,
               label='FIRE-2', rasterized=True, )
    binsize = 0.1
    xlo = np.floor(np.min(xp) / binsize) * binsize
    xhi = np.ceil(np.max(xp) / binsize) * binsize
    xbins = np.arange(xlo, xhi + 0.5 * binsize, binsize)
    xcen = 0.5 * (xbins[:-1] + xbins[1:])
    yperc, _, xsel = pu.get_perc_and_points(xp, yp, xbins,
                           percentiles=(10., 50., 90.),
                           mincount_x=20,
                           getoutliers_y=False, getmincounts_x=False,
                           x_extremes_only=False)
    print(xcen.shape, yperc.shape)
    ax.plot(xcen[xsel], (yperc[1:-1, 0])[xsel],
            color='black', linestyle='dashed',
            linewidth=1.2, label='FIRE-2: 10, 90%')
    ax.plot(xcen[xsel], (yperc[1:-1, 2])[xsel], 
            color='black', linestyle='dashed',
            linewidth=1.2, label=None)
    ax.plot(xcen[xsel], (yperc[1:-1, 1])[xsel], 
            color='black', linestyle='solid',
            linewidth=1.5, label='FIRE-2: median')
    
def plot3rdfactor_N_sigma(simdatafilen: str,
                          outfilen: str | None = None):
    fontsize = 12
    
    fig = plt.figure(figsize=(7.5, 6.))
    grid1 = gsp.GridSpec(nrows=2, ncols=2,
                         wspace=0.4, hspace=0.05,
                        )
    grid = [gsp.GridSpecFromSubplotSpec(
                ncols=2, nrows=1, subplot_spec=grid1[i // 2, i %2],
                hspace=0., wspace=0.05, width_ratios=[5., 0.5])
            for i in range(4)]
    axes = [fig.add_subplot(grid[i][j]) 
            for i in range(4) for j in range(2)]

    simdata = pd.read_csv(simdatafilen, sep='\t')
    plotpanel_colby(axes[0], axes[1], simdata, 'Mvir', fontsize=fontsize)
    plotpanel_colby(axes[2], axes[3], simdata, 'Mstar', fontsize=fontsize)
    plotpanel_colby(axes[4], axes[5], simdata, 'ipar', fontsize=fontsize)
    plotpanel_colby(axes[6], axes[7], simdata, 'sigmavgal', 
                    fontsize=fontsize)
    axes[0].legend(fontsize=fontsize - 1.)

    for axi, ax in enumerate(axes):
        if axi % 2 == 1: # caxes
            axes[axi].tick_params(which='both', labelsize=fontsize - 1.)
            continue
        _axi = axi // 2
        dobottom = _axi >= 2
        doleft = _axi % 2 == 0
        if dobottom:
             ax.set_xlabel('$\\log_{10} \\, \\mathrm{N} \\;'
                           ' [\\mathrm{cm}^{-2}]$',
                           fontsize=fontsize)
        if doleft:
            ax.set_ylabel('$\\sigma(v - \\langle v \\rangle_{\\tau}) \\;'
                          '[\\mathrm{km}/\\mathrm{s}]$',
                          fontsize=fontsize)
        ax.tick_params(direction='in', which='both',
                       labelsize=fontsize - 1.,
                       left=True, top=True,
                       labelleft=doleft, labelbottom=dobottom)
    if outfilen is not None:
        plt.savefig(outfilen, bbox_inches='tight')