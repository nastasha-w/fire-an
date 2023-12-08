import h5py
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fire_an.makeplots.litcomp.obsdataread as odr
import fire_an.makeplots.litcomp.obs_datasel as ods
import fire_an.simlists as sl


proffilen = ('/projects/b1026/nastasha/plotdata/'
             'coldens_radprof_Ne8_opt2.hdf5')
mdir = '/projects/b1026/nastasha/imgs/datacomp/'
oddir = '/projects/b1026/nastasha/extdata/'
q23filen = oddir + 'plotdata_q23_nsigmas_1_2.dat'
b19filen = oddir + 'plotdata_b19_nsigmas_1_2.dat'

def plot_obscomp(massset='m12', obssample='B+19', zr='z0.5-1.0',
                 ricut_pkpc=450., sample='main'):

    percs_shading = ['0.1', '0.9']
    perc_mid = '0.5'
    samplestr = ''
    sample_ds = sample
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
        with h5py.File(proffilen, 'r') as f:
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
    
    # sync ax limits
    ylims = [ax.get_ylim() for ax in axes]
    ymax = min([yl[1] for yl in ylims])
    ymin = max([yl[0] for yl in ylims])
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

    outname = mdir + (f'coldenscomp_Ne8_{obssample}_vs_{massset}_at_{zr}'
                      f'_opt2{samplestr}.pdf')
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