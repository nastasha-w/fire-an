import h5py
import matplotlib.gridspec as gsp
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fire_an.makeplots.litcomp.obsdataread as odr
import fire_an.makeplots.litcomp.obs_datasel as ods
import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
import fire_an.simlists as sl
import fire_an.utils.math_utils as mu


proffilen = ('/projects/b1026/nastasha/plotdata/'
             'coldens_radprof_Ne8_opt2.hdf5')
proffilen_f3x = ('/projects/b1026/nastasha/plotdata/'
                 'coldens_radprof_Ne8_opt3.hdf5')
mdir = '/projects/b1026/nastasha/imgs/datacomp/'
oddir = '/projects/b1026/nastasha/extdata/'
q23filen = oddir + 'plotdata_q23_nsigmas_1_2.dat'
b19filen = oddir + 'plotdata_b19_nsigmas_1_2.dat'

def runningpercentiles(xvals, yvals, yperc=0.5, npoints=5):
    xp = []
    yp = []
    _xvals = np.array(xvals)
    _yvals = np.array(yvals)
    xorder = np.argsort(_xvals)
    #print('------ start ------')
    #print('yperc: ', yperc)
    #print(_xvals[xorder])
    #print(_yvals[xorder])
    for i in range(len(xvals) - npoints + 1):
        sel = slice(i, i + npoints, None)
        #print('x in: ', (_xvals[xorder])[sel])
        #print('y in: ', (_yvals[xorder])[sel])
        # median of x points
        xp.append(np.quantile((_xvals[xorder])[sel], 0.5))
        yp.append(np.quantile((_yvals[xorder])[sel], yperc))
        #print('x, y calc: ', xp[-1], ', ', yp[-1])
    #print('------ end ------')
    return np.array(xp), np.array(yp)

def plot_obscomp(massset='m12', obssample='B+19', zr='z0.5-1.0',
                 ricut_pkpc=450., sample='main',
                 percentilesul=False, npoints_percul=(5, 10)):

    percs_shading = ['0.1', '0.9']
    perc_mid = '0.5'
    samplestr = ''
    sample_ds = sample
    firefilen = proffilen
    if massset == 'm12':
        if sample == 'main':
            samplestr = '_main'
            physmodels = ['FIRE-2', 'noBH', 'AGN-noCR', 'AGN-CR']
        elif sample == 'inclm12plus':
            samplestr =  '_inclm12plus'
            physmodels = ['FIRE-2', 'noBH', 'noBH-m12+', 'AGN-noCR', 'AGN-CR']
            sample_ds = 'main'
        elif sample == 'm12_f3nobh_comp':
            samplestr =  '_m12_f3nobh_comp'
            physmodels = ['noBH', 'noBH-m12+']
            sample_ds = 'm12_f3nobh_comp'
        elif sample == 'fire2':
            samplestr = '_fire2'
            physmodels = ['FIRE-2']
        elif sample == 'f3xtest':
            samplestr = '_f3xtest'
            physmodels = ['FIRE-2', 'noBH', 'FIRE-3x-constpterm', 
                          'FIRE-3x-scmodules']
            sample = 'f3xtest'
            firefilen = proffilen_f3x
    elif massset == 'm13':
        physmodels = ['noBH', 'AGN-noCR', 'AGN-CR']

    # get obs. data
    mass_minmax, z_minmax = ods.get_M_z_boxes_fire(sample=sample_ds)
    if obssample == 'B+19':
        obsdata = pd.read_csv(b19filen, sep='\t')
        mh_obs = obsdata['logmvir_msun_bestest'].copy()
        isul_obs = obsdata['log_N_Ne8_isUL'].copy()
        ipar_obs = obsdata['impact_parameter_kpc'].copy() 
        cd_obs = obsdata['log_N_Ne8_pcm2'].copy()
        z_obs = obsdata['zgal'].copy()
        cderr_obs = np.array((obsdata['log_N_Ne8_pcm2_err'].copy(),
                              obsdata['log_N_Ne8_pcm2_err'].copy()))
        flagged = obsdata['flagged_by_qu23'].copy()
    elif  obssample == 'Q+23':
        obsdata = pd.read_csv(q23filen, sep='\t')
        cd_obs = obsdata['ne8col_logcm2'].copy()
        cdmeas = np.logical_not(np.isnan(cd_obs))
        cd_obs = cd_obs[cdmeas]
        mh_obs = obsdata['logmvir_msun_bestest'].copy()[cdmeas]
        isul_obs = obsdata['isul_ne8'].copy()[cdmeas]
        ipar_obs = obsdata['impactpar_kpc'].copy()[cdmeas] 
        z_obs = obsdata['z_gal'].copy()[cdmeas]
        cderr_obs = np.array((obsdata['ne8col_2s_loerr_dex'].copy()[cdmeas],
                              obsdata['ne8col_2s_hierr_dex'].copy()[cdmeas]))
        flagged = np.zeros(isul_obs.shape, dtype=bool)
    mainsel = np.logical_and(mh_obs >= mass_minmax[0][massset][0],
                             mh_obs <= mass_minmax[0][massset][1])
    mainsel &= np.logical_and(z_obs >= z_minmax[0][massset][0],
                              z_obs <= z_minmax[0][massset][1])
    restsel = np.logical_and(mh_obs >= mass_minmax[1][massset][0],
                             mh_obs <= mass_minmax[1][massset][1])
    restsel &= np.logical_and(z_obs >= z_minmax[1][massset][0],
                              z_obs <= z_minmax[1][massset][1])
    restsel &= np.logical_not(mainsel)
    notul = np.logical_not(isul_obs)
    color_main = 'black'
    color_rest = 'gray'
    dsel_main_noul = np.logical_and(mainsel, notul)
    dsel_rest_noul = np.logical_and(restsel, notul)
    dsel_main_ul = np.logical_and(mainsel, isul_obs)
    dsel_rest_ul = np.logical_and(restsel, isul_obs)
    
    fontsize = 12    
    panelsize = 2.5
    npanels = len(physmodels)
    ncols = npanels
    nrows = 1
    figheight = nrows * panelsize
    figwidth = ncols * panelsize

    fig = plt.figure(figsize=(figwidth, figheight))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows, hspace=0.0,
                        wspace=0.0)
    axes = [fig.add_subplot(grid[0, i]) for i in range(npanels)]
    figtitle = massset + ' halos'
    # avoid overlap with panel titles
    fig.suptitle(figtitle, fontsize=fontsize + 1, y=0.977, 
                 verticalalignment='bottom')

    for axi, (physmodel, ax) in enumerate(zip(physmodels, axes)):
        doleft = axi % ncols == 0
        ax.tick_params(which='both', direction='in', top=True,
                       right=True, labelbottom=True, labelleft=doleft) 
        if doleft:
            ax.set_ylabel('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\,VIII})'
                          '\\; [\\mathrm{cm}^{-2}]$', fontsize=fontsize)
        ax.set_xlabel('$\\mathrm{r}_{\\perp} \\; [\\mathrm{pkpc}]$',
                      fontsize=fontsize)
        ax.set_title(sl.plotlabel_from_physlabel[physmodel],
                     fontsize=fontsize)

        # plot FIRE data
        grpn = massset + '_' + physmodel + '_' + zr
        with h5py.File(firefilen, 'r') as f:
            grp = f[grpn]
            rbins = grp['rbins'][:]
            yv_mid = grp[f'perc-{perc_mid}'][:]
            yv_lo = grp[f'perc-{percs_shading[0]}'][:]
            yv_hi = grp[f'perc-{percs_shading[1]}'][:]
        rcen = 0.5 * (rbins[:-1] + rbins[1:])
        ax.plot(rcen, yv_mid, linewidth=1.5, linestyle='solid',
                color='black', label='FIRE med.')
        ax.fill_between(rcen, yv_lo, yv_hi,
                        color='black', alpha=0.2)
        # plot obs. data
        labelsdone = False
        for dsel_ul, dsel_noul, color in [(dsel_main_ul, dsel_main_noul, 
                                           color_main),
                                          (dsel_rest_ul, dsel_rest_noul, 
                                           color_rest)]:
            if not labelsdone:
                ullabel = obssample + ' UL'
                noul_mainlabel = obssample
                if sum(flagged) > 0:
                    noul_flaggedlabel = obssample + ' (!)'
                else:
                    noul_flaggedlabel = None
            else:
                ullabel = None
                noul_mainlabel = None
                noul_flaggedlabel = None
            dsel_noul_main = np.logical_and(dsel_noul, 
                                            np.logical_not(flagged))
            dsel_noul_flag = np.logical_and(dsel_noul, flagged)
            ax.errorbar(ipar_obs[dsel_noul_main], cd_obs[dsel_noul_main],
                        yerr=cderr_obs[:, dsel_noul_main], 
                        linestyle='None', elinewidth=2., marker='o', 
                        markersize=7, color=color, capsize=3,
                        label=noul_mainlabel, zorder=5,
                        markeredgecolor='black', ecolor='black')
            ax.errorbar(ipar_obs[dsel_noul_flag], cd_obs[dsel_noul_flag],
                        yerr=cderr_obs[:, dsel_noul_flag], 
                        linestyle='None', elinewidth=2., marker='o', 
                        markersize=7, markerfacecolor='none', 
                        markeredgecolor=color,
                        capsize=3, markeredgewidth=1.5,
                        label=noul_flaggedlabel, zorder=5,
                        ecolor='black')
            ax.scatter(ipar_obs[dsel_ul], cd_obs[dsel_ul],
                        linestyle='None', marker='v', 
                        s=30, facecolors='none', edgecolors=color, 
                        label=ullabel, zorder=5, linewidths=1.5)
            labelsdone = True
        if percentilesul:
            xv_goodmatch = (ipar_obs[np.logical_or(dsel_main_ul,
                                                   dsel_main_noul)]).copy()
            #xv_all = ipar_obs.copy() # includes stuff at wrong mass, z, ipar
            yv_goodmatch = (cd_obs[np.logical_or(dsel_main_ul,
                                                 dsel_main_noul)]).copy()
            #yv_all = cd_obs.copy()
            xp_goodmed, yp_goodmed = \
                runningpercentiles(xv_goodmatch, yv_goodmatch, 
                                   yperc=float(perc_mid), 
                                   npoints=npoints_percul[0])
            xp_goodhi, yp_goodhi = \
                runningpercentiles(xv_goodmatch, yv_goodmatch, 
                                   yperc=float(percs_shading[1]), 
                                   npoints=npoints_percul[1])
            #xp_allmed, yp_allmed = \
            #    runningpercentiles(xv_all, yv_all, 
            #                       yperc=float(perc_mid), 
            #                       npoints=npoints_percul[0])
            #xp_allhi, yp_allhi = \
            #    runningpercentiles(xv_all, yv_all, 
            #                       yperc=float(percs_shading[1]), 
            #                       npoints=npoints_percul[1])
            ax.plot(xp_goodmed, yp_goodmed, linestyle='dashed', 
                    color='red', zorder=10)
            ax.plot(xp_goodhi, yp_goodhi, linestyle='dotted', 
                    color='red', zorder=10)
            #ax.plot(xp_allmed, yp_allmed, linestyle='dashed', 
            #        color='pink', zorder=10)
            #ax.plot(xp_allhi, yp_allhi, linestyle='dotted', 
            #        color='pink', zorder=10)
    # sync ax limits
    ylims = [ax.get_ylim() for ax in axes]
    ymax = max([yl[1] for yl in ylims])
    ymin = min([yl[0] for yl in ylims])
    ymin = max(ymin, 12.5)
    [ax.set_ylim((ymin, ymax)) for ax in axes]
    xlims = [ax.get_xlim() for ax in axes]
    xmin = min([xl[0] for xl in xlims])
    xmax = max([xl[1] for xl in xlims])
    xmax = min(xmax, ricut_pkpc)
    xmin = max(xmin, 0.)
    [ax.set_xlim((xmin, xmax)) for ax in axes]
    
    handles, _ = ax.get_legend_handles_labels()
    leg0 = axes[0].legend(handles=handles[:1],
                          fontsize=fontsize - 3., loc='lower left',
                          handlelength=1., ncol=1, handletextpad=0.3,
                          columnspacing=1.0, labelspacing=0.3,
                          borderaxespad=0.2)
    axes[0].legend(handles=handles[1:],
                   fontsize=fontsize - 3., loc='upper right',
                   handlelength=1., ncol=1, handletextpad=0.3,
                   columnspacing=1.0, labelspacing=0.3,
                   borderaxespad=0.2)
    axes[0].add_artist(leg0)
    if percentilesul:
        handles = [mlines.Line2D((), (), color=color_main,
                                 linestyle='dashed', 
                                 label='run. perc. ' + perc_mid 
                                        + f', {npoints_percul[0]} pts.'),
                   mlines.Line2D((), (), color=color_main,
                                 linestyle='dotted', 
                                 label='run. perc. ' + percs_shading[1] 
                                        + f', {npoints_percul[1]} pts.'),
                    mlines.Line2D((), (), color='red',
                                  linestyle='dashdot', 
                                  label='good Mvir match'),
                    mlines.Line2D((), (), color='pink',
                                  linestyle='dashdot', 
                                  label='ok Mvir match')]
        axes[1].legend(handles=handles[:2], 
                       fontsize=fontsize - 3., loc='upper right',
                       handlelength=1., ncol=1, handletextpad=0.3,
                       columnspacing=1.0, labelspacing=0.3,
                       borderaxespad=0.2)
        axes[2].legend(handles=handles[2:], 
                       fontsize=fontsize - 3., loc='upper right',
                       handlelength=1., ncol=1, handletextpad=0.3,
                       columnspacing=1.0, labelspacing=0.3,
                       borderaxespad=0.2)
        perculstr = (f'_percUL_running_{perc_mid}_{npoints_percul[0]}pts'
                     f'_{percs_shading[1]}_{npoints_percul[1]}pts')
        perculstr = perculstr.replace('.', 'p')
    else:
        perculstr = ''
    outname = mdir + (f'coldenscomp_Ne8_{obssample}_vs_{massset}_at_{zr}'
                      f'_opt2{samplestr}{perculstr}.pdf')
    plt.savefig(outname, bbox_inches='tight')

def plot_obscomp_percul(massset='m12', physmodel='FIRE-2',
                        npoints_percul=((8, 15), (10, 20), (15, 30)),
                        ricut_pkpc=450., zr='z0.5-1.0'):
    '''
    nice version of this plot made with the default settings
    '''
    if not hasattr(npoints_percul[0], '__len__'):
        npoints_percul = (npoints_percul,) * 3

    percs_shading = ['0.1', '0.9']
    perc_mid = '0.5'
    sample_ds = 'main' # used for M-z obs. selection
    firefilen = proffilen

    # get obs. data
    mass_minmax, z_minmax = ods.get_M_z_boxes_fire(sample=sample_ds)
    obsdata_b19 = pd.read_csv(b19filen, sep='\t')
    mh_b19 = obsdata_b19['logmvir_msun_bestest'].copy()
    isul_b19 = obsdata_b19['log_N_Ne8_isUL'].copy()
    ipar_b19 = obsdata_b19['impact_parameter_kpc'].copy() 
    cd_b19 = obsdata_b19['log_N_Ne8_pcm2'].copy()
    z_b19 = obsdata_b19['zgal'].copy()
    cderr_b19 = np.array((obsdata_b19['log_N_Ne8_pcm2_err'].copy(),
                            obsdata_b19['log_N_Ne8_pcm2_err'].copy()))
    flagged_b19 = obsdata_b19['flagged_by_qu23'].copy()

    obsdata_q23 = pd.read_csv(q23filen, sep='\t')
    cd_q23 = obsdata_q23['ne8col_logcm2'].copy()
    cdmeas_q23 = np.logical_not(np.isnan(cd_q23))
    cd_q23 = cd_q23[cdmeas_q23]
    mh_q23 = obsdata_q23['logmvir_msun_bestest'].copy()[cdmeas_q23]
    isul_q23 = obsdata_q23['isul_ne8'].copy()[cdmeas_q23]
    ipar_q23 = obsdata_q23['impactpar_kpc'].copy()[cdmeas_q23] 
    z_q23 = obsdata_q23['z_gal'].copy()[cdmeas_q23]
    cderr_q23 = np.array((
        obsdata_q23['ne8col_2s_loerr_dex'].copy()[cdmeas_q23],
        obsdata_q23['ne8col_2s_hierr_dex'].copy()[cdmeas_q23]))
    flagged_q23 = np.zeros(isul_q23.shape, dtype=bool)

    isq23 = np.zeros((len(cd_q23) + len(cd_b19),), dtype=bool)
    isq23[:len(cd_q23)] = True
    isb19 = np.zeros((len(cd_q23) + len(cd_b19),), dtype=bool)
    isb19[len(cd_q23):] = True
    mh_obs = np.append(mh_q23, mh_b19)
    z_obs = np.append(z_q23, z_b19)
    ipar_obs = np.append(ipar_q23, ipar_b19)
    cd_obs = np.append(cd_q23, cd_b19)
    cderr_obs = np.append(cderr_q23, cderr_b19, axis=1)
    isul_obs = np.append(isul_q23, isul_b19)
    flagged = np.append(flagged_q23, flagged_b19)

    mainsel = np.logical_and(mh_obs >= mass_minmax[0][massset][0],
                             mh_obs <= mass_minmax[0][massset][1])
    mainsel &= np.logical_and(z_obs >= z_minmax[0][massset][0],
                              z_obs <= z_minmax[0][massset][1])
    restsel = np.logical_and(mh_obs >= mass_minmax[1][massset][0],
                             mh_obs <= mass_minmax[1][massset][1])
    restsel &= np.logical_and(z_obs >= z_minmax[1][massset][0],
                              z_obs <= z_minmax[1][massset][1])
    restsel &= np.logical_not(mainsel)
    notul = np.logical_not(isul_obs)
    
    # copied from the analytical model comparison plot
    _colors = tc.tol_cset('high-contrast')
    c1 = mcolors.to_rgba(_colors.blue)
    c2 = mcolors.to_rgba(_colors.red)
    colors = {'B+19': {'1s': c1, '2s': c1[:3] + (0.5,)},
              'Q+23': {'1s': c2, '2s': c2[:3] + (0.5,)}}
    dsel_main_noul = np.logical_and(mainsel, notul)
    dsel_rest_noul = np.logical_and(restsel, notul)
    dsel_main_ul = np.logical_and(mainsel, isul_obs)
    dsel_rest_ul = np.logical_and(restsel, isul_obs)
    
    fontsize = 12    
    panelsize = 2.5
    npanels = 3
    ncols = npanels
    nrows = 1
    figheight = nrows * panelsize
    figwidth = ncols * panelsize

    fig = plt.figure(figsize=(figwidth, figheight))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows, hspace=0.0,
                        wspace=0.0)
    axes = [fig.add_subplot(grid[0, i]) for i in range(npanels)]
    figtitle = f'{massset} halos, {sl.plotlabel_from_physlabel[physmodel]}'
    fig.suptitle(figtitle, fontsize=fontsize)

    for axi, ax in enumerate(axes):
        doleft = axi % ncols == 0
        ax.tick_params(which='both', direction='in', top=True,
                       right=True, labelbottom=True, labelleft=doleft) 
        if doleft:
            ax.set_ylabel('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\,VIII})'
                          '\\; [\\mathrm{cm}^{-2}]$', fontsize=fontsize)
        ax.set_xlabel('$\\mathrm{r}_{\\perp} \\; [\\mathrm{pkpc}]$',
                      fontsize=fontsize)

        # plot FIRE data
        grpn = massset + '_' + physmodel + '_' + zr
        with h5py.File(firefilen, 'r') as f:
            grp = f[grpn]
            rbins = grp['rbins'][:]
            yv_mid = grp[f'perc-{perc_mid}'][:]
            yv_lo = grp[f'perc-{percs_shading[0]}'][:]
            yv_hi = grp[f'perc-{percs_shading[1]}'][:]
        rcen = 0.5 * (rbins[:-1] + rbins[1:])
        ax.plot(rcen, yv_mid, linewidth=1.5, linestyle='solid',
                color='black')
        ax.fill_between(rcen, yv_lo, yv_hi,
                        color='black', alpha=0.2)
    # plot obs. data
    labelsdone = False
    for dsel_ul, dsel_noul, ckey in [(dsel_main_ul, dsel_main_noul, '1s'),
                                     (dsel_rest_ul, dsel_rest_noul, '2s')]:

        dsel_noul_main = np.logical_and(dsel_noul, 
                                        np.logical_not(flagged))
        dsel_noul_flag = np.logical_and(dsel_noul, flagged)
        for obssel, obsset in ((isb19, 'B+19'), (isq23, 'Q+23')):
            if obsset == 'B+19':
                axis = [0, 2]
            elif obsset == 'Q+23':
                axis = [1, 2]
            for axi in axis:
                _msel = np.logical_and(dsel_noul_main, obssel)
                _fsel = np.logical_and(dsel_noul_flag, obssel)
                _usel = np.logical_and(dsel_ul, obssel)
                if axi == 2 and not labelsdone:
                    ullabel = obsset + ' UL'
                    noul_mainlabel = obsset
                    if sum(_fsel) > 0:
                        noul_flaggedlabel = obsset + ' (!)'
                    else:
                        noul_flaggedlabel = None
                else:
                    ullabel = None
                    noul_mainlabel = None
                    noul_flaggedlabel = None
                axes[axi].errorbar(ipar_obs[_msel], cd_obs[_msel],
                                   yerr=cderr_obs[:, _msel], 
                                   linestyle='None', elinewidth=2.,
                                   marker='o', 
                                   markersize=6, color=colors[obsset][ckey],
                                   capsize=3, zorder=5,
                                   markeredgecolor='black',
                                   ecolor='black', label=noul_mainlabel)
                axes[axi].errorbar(ipar_obs[_fsel], cd_obs[_fsel],
                                   yerr=cderr_obs[:, _fsel], 
                                   linestyle='None', elinewidth=2., 
                                   marker='o', 
                                   markersize=6, markerfacecolor='none', 
                                   markeredgecolor=colors[obsset][ckey],
                                   capsize=3, markeredgewidth=1.5,
                                   zorder=5, ecolor='black',
                                   label=noul_flaggedlabel)
                axes[axi].scatter(ipar_obs[_usel], cd_obs[_usel],
                                  linestyle='None', marker='v', 
                                  s=12, facecolors='none',
                                  edgecolors=colors[obsset][ckey], 
                                  zorder=4, linewidths=1.2,
                                  label=ullabel)
        labelsdone = True
    allsel = np.logical_or(mainsel, restsel)
    allsel = np.logical_and(allsel, ipar_obs <= ricut_pkpc)
    for axi, (ax, obssets) in \
            enumerate(zip(axes, [('B+19',), ('Q+23',), ('B+19', 'Q+23')])):
        if obssets == ('B+19',):
            obssel = np.logical_and(allsel, isb19)
        elif obssets == ('Q+23',):
            obssel = np.logical_and(allsel, isq23)
        elif obssets == ('B+19', 'Q+23'):
            obssel = np.logical_or(isb19, isq23) 
            obssel = np.logical_and(allsel, obssel)
        print(np.sum(obssel))
        _np_percul = npoints_percul[axi]
        xv_all = ipar_obs[obssel].copy()
        yv_all = cd_obs[obssel].copy()
        xp_allmed, yp_allmed = \
            runningpercentiles(xv_all, yv_all, 
                                yperc=float(perc_mid), 
                                npoints=_np_percul[0])
        xp_allhi, yp_allhi = \
            runningpercentiles(xv_all, yv_all, 
                                yperc=float(percs_shading[1]), 
                                npoints=_np_percul[1])
        obsperccolor = (0.35, 0.35, 0.35)
        pe = pu.getoutline(2.)
        ax.plot(xp_allmed, yp_allmed, linestyle='solid', 
                color=obsperccolor, linewidth=2., zorder=6,
                path_effects=pe)
        pe = pu.getoutline(1.3)
        ax.plot(xp_allhi, yp_allhi, linestyle='solid', 
                color=obsperccolor, linewidth=1.3, zorder=6,
                path_effects=pe)
        # add arrows indicating these are upper limits
        offset_x = 15. # kpc (x units)
        length = 0.18 # dex (y units)
        dx = 0.
        dy = -1. * length
        hwidth = 12.
        hlength = 0.07
        akwa = {'facecolor': obsperccolor, 
                'edgecolor': 'black', 
                'length_includes_head': True,
                'head_length': hlength,
                'head_width': hwidth,
                'zorder': 5.9,
                'width': 0.004}
        for xp, yp in [(xp_allmed, yp_allmed), (xp_allhi, yp_allhi)]:
            x0 = xp[0] + offset_x
            y0 = mu.linterpsolve(xp, yp, x0)
            ax.arrow(x0, y0, dx, dy, **akwa)
            x1 = xp[-1] - offset_x
            y1 = mu.linterpsolve(xp, yp, x1)
            ax.arrow(x1, y1, dx, dy, **akwa)
            x2 = 0.5 * (xp[0] + xp[-1])
            y2 = mu.linterpsolve(xp, yp, x2)
            ax.arrow(x2, y2, dx, dy, **akwa)
    # sync ax limits
    ylims = [ax.get_ylim() for ax in axes]
    ymax = max([yl[1] for yl in ylims])
    ymin = min([yl[0] for yl in ylims])
    ymin = max(ymin, 12.5)
    [ax.set_ylim((ymin, ymax)) for ax in axes]
    xlims = [ax.get_xlim() for ax in axes]
    xmin = min([xl[0] for xl in xlims])
    xmax = max([xl[1] for xl in xlims])
    xmax = min(xmax, ricut_pkpc)
    xmin = max(xmin, 0.)
    [ax.set_xlim((xmin, xmax)) for ax in axes]
    
    handles_obs, _ = axes[2].get_legend_handles_labels()
    b19_hsel = [0, 2, 3]
    q23_hsel = [1, 4]
    leg0 = axes[0].legend(handles=[handles_obs[hi] for hi in b19_hsel],
                   fontsize=fontsize - 3., loc='upper right',
                   handlelength=1., ncol=1, handletextpad=0.3,
                   columnspacing=1.0, labelspacing=0.3,
                   borderaxespad=0.2)
    axes[1].legend(handles=[handles_obs[hi] for hi in q23_hsel],
                   fontsize=fontsize - 3., loc='upper right',
                   handlelength=1., ncol=1, handletextpad=0.3,
                   columnspacing=1.0, labelspacing=0.3,
                   borderaxespad=0.2)
    # + f', {npoints_percul[2][0]} pts.'
    # + f', {npoints_percul[2][1]} pts.'
    handles = [mlines.Line2D((), (), color=(0.35, 0.35, 0.35),
                             linestyle='solid', linewidth=1.3,
                             path_effects=pu.getoutline(1.3),
                             label=f'{float(percs_shading[1]) * 100: .0f}% UL'
                             ),
               mlines.Line2D((), (), color=(0.35, 0.35, 0.35),
                             linestyle='solid', linewidth=2.5,
                             path_effects=pu.getoutline(2.),
                             label=f'{float(perc_mid) * 100: .0f}% UL'
                             )]
    axes[2].legend(handles=handles, 
                   fontsize=fontsize - 3., loc='upper right',
                   handlelength=1.5, ncol=1, handletextpad=0.3,
                   columnspacing=1.0, labelspacing=0.3,
                   borderaxespad=0.2)
    ## copied legend defaults from docs
    #bbox_text = {'facecolor': 'white',
    #             'edgecolor': '0.8',
    #             'alpha': 0.8,
    #             'boxstyle': 'round'}
    #note0 = (f'{float(perc_mid) * 100: .0f}%: ' 
    #         f'{npoints_percul[0][0]} pts., '
    #         f'{float(percs_shading[1]) * 100: .0f}%: '
    #         f'{npoints_percul[0][1]} pts.')
    #axes[0].text(0.035, 0.025, note0, color='black',
    #             fontsize=fontsize - 3., transform=axes[0].transAxes,
    #             horizontalalignment='left', verticalalignment='bottom',
    #             bbox=bbox_text)
    #note1 = (f'{float(perc_mid) * 100: .0f}%: ' 
    #         f'{npoints_percul[1][0]} pts., '
    #         f'{float(percs_shading[1]) * 100: .0f}%: '
    #         f'{npoints_percul[1][1]} pts.')
    #axes[1].text(0.035, 0.025, note1, color='black',
    #             fontsize=fontsize - 3., transform=axes[1].transAxes,
    #             horizontalalignment='left', verticalalignment='bottom',
    #             bbox=bbox_text)
    handlesl = [mlines.Line2D((), (), linewidth=1.5, linestyle='solid',
                              color='black', 
                              label= sl.plotlabel_from_physlabel[physmodel]
                                     + ' med.')]
    axes[0].legend(handles=handlesl,
                          fontsize=fontsize - 3., loc='lower left',
                          handlelength=1., ncol=1, handletextpad=0.3,
                          columnspacing=1.0, labelspacing=0.3,
                          borderaxespad=0.2)
    axes[0].add_artist(leg0)

    nplist0 = [pts[0] for pts in npoints_percul]
    ptsstr0 = '_'.join([str(pts) for pts in nplist0])
    nplist1 = [pts[1] for pts in npoints_percul]
    ptsstr1 = '_'.join([str(pts) for pts in nplist1])
    perculstr = (f'_percUL_running_{perc_mid}_{ptsstr0}pts'
                 f'_{percs_shading[1]}_{ptsstr1}pts')
    outname = (f'coldenscomp_Ne8_obs_vs_{massset}_at_{zr}'
               f'_opt2_{physmodel}{perculstr}')
    outname = mdir + outname.replace('.', 'p') + '.pdf'
    plt.savefig(outname, bbox_inches='tight')

def runplots_obscomp():
    ricut_pkpc = 450.
    plot_obscomp(massset='m12', obssample='B+19', zr='z0.5-1.0',
                 ricut_pkpc=ricut_pkpc, sample='main')
    plot_obscomp(massset='m13', obssample='B+19', zr='z0.5-1.0',
                 ricut_pkpc=ricut_pkpc, sample='main')
    plot_obscomp(massset='m12', obssample='Q+23', zr='z0.5-1.0',
                 ricut_pkpc=ricut_pkpc, sample='main')
    plot_obscomp(massset='m13', obssample='Q+23', zr='z0.5-1.0',
                 ricut_pkpc=ricut_pkpc, sample='main')
    plot_obscomp(massset='m12', obssample='Q+23', zr='z0.5-0.7',
                 ricut_pkpc=ricut_pkpc, sample='main')
    plot_obscomp(massset='m13', obssample='Q+23', zr='z0.5-0.7',
                 ricut_pkpc=ricut_pkpc, sample='main')
    # make sure mass selection is consisent
    ods.plotMz_obs_fire_2panel(ricut_pkpc=ricut_pkpc, sample='main') 

def runplots_appendix():
    ricut_pkpc = 450.
    plot_obscomp(massset='m12', obssample='B+19', zr='z0.5-1.0',
                 ricut_pkpc=ricut_pkpc, sample='m12_f3nobh_comp')
    plot_obscomp(massset='m12', obssample='Q+23', zr='z0.5-1.0',
                 ricut_pkpc=ricut_pkpc, sample='m12_f3nobh_comp')
    ods.plotMz_obs_fire_2panel(ricut_pkpc=ricut_pkpc, 
                               sample='m12_f3nobh_comp') 
    
def runplots_runningperc_ul():
    plot_obscomp(massset='m12', obssample='B+19', zr='z0.5-1.0',
                 ricut_pkpc=450., sample='main',
                 percentilesul=True, npoints_percul=(5, 10))
    plot_obscomp(massset='m12', obssample='B+19', zr='z0.5-1.0',
                 ricut_pkpc=450., sample='main',
                 percentilesul=True, npoints_percul=(7, 12))
    plot_obscomp(massset='m12', obssample='B+19', zr='z0.5-1.0',
                 ricut_pkpc=450., sample='main',
                 percentilesul=True, npoints_percul=(10, 15))
    plot_obscomp(massset='m12', obssample='Q+23', zr='z0.5-1.0',
                 ricut_pkpc=450., sample='main',
                 percentilesul=True, npoints_percul=(5, 10))
    plot_obscomp(massset='m12', obssample='Q+23', zr='z0.5-1.0',
                 ricut_pkpc=450., sample='main',
                 percentilesul=True, npoints_percul=(7, 12))
    plot_obscomp(massset='m12', obssample='Q+23', zr='z0.5-1.0',
                 ricut_pkpc=450., sample='main',
                 percentilesul=True, npoints_percul=(10, 15))