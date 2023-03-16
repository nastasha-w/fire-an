import h5py
import numpy as np

import matplotlib.collections as mcol
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.gridspec as gsp
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt

import fire_an.makeplots.get_2dprof as gpr
import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
import fire_an.utils.constants_and_units as c


def compare_profiles_physmodels(modelfilelists,
                                legendlablists,
                                physmodel_labels,
                                rbins_pkpc, title=None, 
                                ylabel=None, ax=None,
                                outname=None, domean=True,
                                doscatter=True, colors=None):
    '''
    modelfilelists, legendlablists: lists of lists, shapes should match
    inner index: physics model differences
    outer index: different versions with that physics model 
                 (e.g., different ics, redshifts)
    '''
    fontsize = 12
    if ax is None:
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5.5, 5.))
        if title is not None:
            fig.suptitle(title, fontsize=fontsize)
        if ylabel is not None:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.set_xlabel('$\\mathrm{r}_{\\perp}$ [pkpc]', fontsize=fontsize)
    
    if len(modelfilelists[0]) <= 4:
        lss = ['solid', 'dashed', 'dashdot', 'dashdot']
    else:
        lss = ['solid'] * len(modelfilelists[0])
    if domean:
        lw_med = 1.
        lw_av = 2.5
    else:
        lw_med = 2.
        lw_av = 0.5
    alpha_range = 0.3
    if colors is None:
        colors = tc.tol_cset('bright')

    rcens = 0.5 * (rbins_pkpc[1:] + rbins_pkpc[:-1])
    #print(modelfilelists)
    for vi, (modelfiles, leglabels) in \
             enumerate(zip(modelfilelists, legendlablists)):
        #print(modelfiles)
        for pi, (modelf, ls, leglab) in \
                enumerate(zip(modelfiles, lss, leglabels)):
            #print(modelf)
            with h5py.File(modelf, 'r') as f:  
                rvir_bh_pkpc = f['Header/inputpars/halodata'].attrs['Rvir_cm']
                rvir_bh_pkpc /= (c.cm_per_mpc * 1e-3)
        
            av, med, p90, p10 = gpr.get_profile_massmap(
                modelf, rbins_pkpc, rbin_units='pkpc', 
                profiles=['av-lin', 'perc-0.5', 'perc-0.9', 'perc-0.1'])
            
            color = colors[len(modelfiles) * vi + pi]
            if domean:
                ax.plot(rcens, av, color=color, linestyle=ls, 
                        linewidth=lw_av, label=leglab)
            ax.plot(rcens, med, color=color, linestyle=ls, 
                    linewidth=lw_med, 
                    label=None if domean else leglab)
            if doscatter:
                ax.fill_between(rcens, p10, p90, color=color, 
                                alpha=alpha_range, linestyle=ls,
                                linewidth=0.5)
            # indicate Rvir on the curves
            yp_bh_med = pu.linterpsolve(rcens, med, rvir_bh_pkpc)
            if domean:
                yp_bh_av = pu.linterpsolve(rcens, av, rvir_bh_pkpc)
                ax.scatter([rvir_bh_pkpc] * 2, [yp_bh_av, yp_bh_med],
                            marker='o', c=color, s=60)
            else:
                ax.scatter([rvir_bh_pkpc], [yp_bh_med],
                            marker='o', c=color, s=60)

    handles3, labels = ax.get_legend_handles_labels()
    handles1 = [mlines.Line2D((), (), linewidth=lw_med, linestyle='solid',
                            label='med.', color='black'),
                mlines.Line2D((), (), linestyle=None, marker='o',
                            label='$\\mathrm{R}_{\\mathrm{vir}}$', 
                            color='black', markersize=5)
                ]
    if domean:
        handles1 = [mlines.Line2D((), (), linewidth=lw_av, linestyle='solid',
                                  label='mean', color='black')] + handles1 
    if doscatter:
        handles1 = handles1 + \
                   [mpatch.Patch(label='perc. 10-90', linewidth=0.5, 
                                 color='black', alpha=alpha_range)]
    if len(modelfilelists) > 1: # more than one line per phys model
        handles2 = [mlines.Line2D((), (), linewidth=1.7, linestyle=ls,
                                  label=lab, color='black')\
                    for ls, lab in zip(lss, physmodel_labels)]
    else:
        handles2 = []
    ax.legend(handles=handles1 + handles2 + handles3, fontsize=fontsize)
    ax.set_xscale('log')

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def compare_profilesets_physmodel(mapset='clean_set1', var='zev'):

    if mapset == 'clean_set1':
        outdir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/profiles_clean_set1/'
        
        mdir = '/Users/nastasha/ciera/sim_maps/fire/clean_set1/'
        tpl_m12_noBH = 'coldens_{qt}_{ic}_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5_snap258_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5'
        tpl_m12_AGN_CR = 'coldens_{qt}_{ic}_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000_snap50_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5'
        tpl_m12_AGN_noCR = 'coldens_{qt}_{ic}_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5_snap258_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5'
        tpl_m13_noBH = 'coldens_{qt}_{ic}_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5_snap50_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5' 
        tpl_m13_AGN_CR = 'coldens_{qt}_{ic}_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000_snap50_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5' 
        tpl_m13_AGN_noCR_h113 = 'coldens_{qt}_{ic}_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5_snap258_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5' 
        tpl_m13_AGN_noCR_h206 = 'coldens_{qt}_{ic}_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31_fa0.5_snap258_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5' 
        model_labels = ['no BH', 'AGN, no CR', 'AGN + CR']
        ics_m12 = ['m12f'] # bug: 'm12m'
        ics_m13 = ['m13h206', 'm13h113']
        qts = ['gas-mass', 'O6', 'Ne8', 'Mg10', 'H1']
        h1_imstr = '_ionfrac-fromsim'
        
        fnss_m12 = [[[mdir + tpl_m12_noBH.format(ic=ic, qt=ion, 
                          im=h1_imstr if ion == 'H1' else ''),
                      mdir + tpl_m12_AGN_noCR.format(ic=ic, qt=ion,
                          im=h1_imstr if ion == 'H1' else ''),
                      mdir + tpl_m12_AGN_CR.format(ic=ic, qt=ion, 
                          im=h1_imstr if ion == 'H1' else ''),
                      ]] \
                    for ic in ics_m12 for ion in qts]
        fnss_m13 = [[[mdir + tpl_m13_noBH.format(ic=ic, qt=ion, 
                          im=h1_imstr if ion == 'H1' else ''),
                      mdir + tpl_m13_AGN_noCR_h206.format(ic=ic, 
                          qt=ion, im=h1_imstr if ion == 'H1' else '')\
                         if 'h206' in ic else
                         mdir + tpl_m13_AGN_noCR_h113.format(ic=ic, 
                             qt=ion, im=h1_imstr if ion == 'H1' else ''),
                      mdir + tpl_m13_AGN_CR.format(ic=ic, qt=ion, 
                          im=h1_imstr if ion == 'H1' else ''),
                     ]] \
                     for ic in ics_m13 for ion in qts]
        fnsss = fnss_m12 + fnss_m13
        model_labelss = [[model_labels]] * len(fnsss)
        physmodel_labelss = [model_labels] * len(fnsss)

        rbinss = [np.linspace(0., 300., 50) if 'm12' in ic else \
                  np.linspace(0., 700., 50) if 'm13' in ic else \
                  None \
                  for ic in ics_m12 + ics_m13 for qt in qts]
        mlabel = ('$\\log_{10} \\, \\Sigma_{\\mathrm{gas}}'+\
                  ' \\, [\\mathrm{g} \\, \\mathrm{cm}^2]$')
        ilabel = ('$\\log_{{10}} \\, \\mathrm{{N}}(\\mathrm{{{ion}}})'+\
                  ' \\, [\\mathrm{{cm}}^2]$')
        ylabels = [mlabel if 'mass' in qt else\
                   ilabel.format(ion=qt) \
                   for ic in  ics_m12 + ics_m13 for qt in qts]
        outtemp = 'rprofcomp_clean_set1_{ic}_{qt}_z0p5_noBH_AGNCR_AGNnoCR.pdf'
        outnames = [outdir + outtemp.format(ic=ic, qt=qt) \
                    for ic in ics_m12 + ics_m13 for qt in qts]

        ttpl_m12 = ('noBH-noCR: m7e3, AGN-noCR: m7e3, sdp2e-4,'
                    ' AGN-CR: m6e4, sdp1e-4, fcr1e-3_vw3000') 
        ttpl_m13h113 = ('noBH-noCR: m3e5, AGN-noCR: m3e4, sdp1e-4,'
                        ' AGN-CR: m3e5, sdp1e-4, fcr1e-3_vw3000') 
        ttpl_m13h206 = ('noBH-noCR: m3e5, AGN-noCR: m3e4, sdp3e-4,'
                        ' AGN-CR: m3e5, sdp1e-4, fcr1e-3_vw3000') 
        title_tpl= '{ion} profiles, {ic} z=0.5, {ttpl}'
        titles = [title_tpl.format(ic=ic, ion=ion, 
                                   ttpl=ttpl_m12 if 'm12' in ic else 
                                        ttpl_m13h206 if 'm13h206' in ic else\
                                        ttpl_m13h113)\
                  for ic in ics_m12 + ics_m13 for ion in qts]
        
        
        _outname = 'rprof_noBH_AGN_AGNCR_set1_{ic}_z0p5_gas_{ion}_nomean.pdf'
        outnames = [_outname.format(ic=ic, ion=ion) \
                    for ic in ics_m12 + ics_m13 for ion in qts]
        domean = False
    if mapset == 'clean_set1-2':
        outdir = ('/Users/nastasha/ciera/projects_lead/fire3_ionabs/'
                  'profiles_clean/')
        mdir1 = '/Users/nastasha/ciera/sim_maps/fire/clean_set1/'
        mdir2 = '/Users/nastasha/ciera/sim_maps/fire/clean_set2/'
        tpl_m12_noBH =          ('coldens_{qt}_{ic}_m7e3_MHD_fire3'
                                 '_fireBH_Sep182021_hr_crdiffc690_sdp1e10'
                                 '_gacc31_fa0.5_snap{snap}'
                                 '_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5')
        tpl_m12_AGN_CR =        ('coldens_{qt}_{ic}_m6e4_MHDCRspec1_fire3'
                                 '_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
                                 '_gacc31_fa0.5_fcr1e-3_vw3000_snap{snap}'
                                 '_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5')
        tpl_m12_AGN_noCR =      ('coldens_{qt}_{ic}_m7e3_MHD_fire3'
                                 '_fireBH_Sep182021_hr_crdiffc690_sdp2e-4'
                                 '_gacc31_fa0.5_snap{snap}'
                                 '_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5')
        tpl_m13_noBH =          ('coldens_{qt}_{ic}_m3e5_MHD_fire3'
                                 '_fireBH_Sep182021_crdiffc690_sdp1e10'
                                 '_gacc31_fa0.5_snap{snap}'
                                 '_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5') 
        tpl_m13_AGN_CR =        ('coldens_{qt}_{ic}_m3e5_MHDCRspec1_fire3'
                                 '_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
                                 '_gacc31_fa0.5_fcr1e-3_vw3000_snap{snap}'
                                 '_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5') 
        tpl_m13_AGN_noCR_h113 = ('coldens_{qt}_{ic}_m3e4_MHD_fire3'
                                 '_fireBH_Sep182021_hr_crdiffc690_sdp1e-4'
                                 '_gacc31_fa0.5_snap{snap}'
                                 '_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5') 
        tpl_m13_AGN_noCR_h206 = ('coldens_{qt}_{ic}_m3e4_MHD_fire3_fireBH'
                                 '_Sep182021_hr_crdiffc690_sdp3e-4'
                                 '_gacc31_fa0.5_snap{snap}'
                                 '_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5') 
        model_labels = ['no BH', 'AGN, no CR', 'AGN + CR']
        ics_m12 = ['m12f'] # bug: 'm12m'
        ics_m13 = ['m13h206', 'm13h113']
        qts = ['gas-mass', 'O6', 'Ne8', 'Mg10', 'H1']
        h1_imstr = '_ionfrac-fromsim'
        snaps_sr = [50, 49, 48, 47, 46, 45]
        snaps_hr =  [258, 240, 224, 210, 197, 186]
        redshifts = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

        mlabel = ('$\\log_{10} \\, \\Sigma_{\\mathrm{gas}}'+\
                  ' \\, [\\mathrm{g} \\, \\mathrm{cm}^2]$')
        ilabel = ('$\\log_{{10}} \\, \\mathrm{{N}}(\\mathrm{{{ion}}})'+\
                  ' \\, [\\mathrm{{cm}}^2]$')
        
        if var == 'model':
            ttpl_m12 = ('noBH: m7e3, AGN-noCR: m7e3, sdp2e-4,'
                       ' AGN-CR: m6e4, sdp1e-4, fcr1e-3_vw3000') 
            ttpl_m13h113 = ('noBH: m3e5, AGN-noCR: m3e4, sdp1e-4,'
                            ' AGN-CR: m3e5, sdp1e-4, fcr1e-3_vw3000') 
            ttpl_m13h206 = ('noBH: m3e5, AGN-noCR: m3e4, sdp3e-4,'
                            ' AGN-CR: m3e5, sdp1e-4, fcr1e-3_vw3000') 

            fnss_m12 = [[[tpl_m12_noBH.format(ic=ic, qt=qt, 
                              im=h1_imstr if qt == 'H1' else '', 
                              snap=snhr),
                          tpl_m12_AGN_noCR.format(ic=ic, qt=qt, 
                              im=h1_imstr if qt == 'H1' else '', 
                              snap=snhr),
                          tpl_m12_AGN_CR.format(ic=ic, qt=qt, 
                              im=h1_imstr if qt == 'H1' else '', 
                              snap=snsr),
                          ]] \
                          for ic in ics_m12 for qt in qts \
                          for snhr, snsr in zip(snaps_hr, snaps_sr)]
            fnss_m13 = [[[tpl_m13_noBH.format(ic=ic, qt=qt,
                              im=h1_imstr if qt == 'H1' else '',
                              snap=snsr),
                          tpl_m13_AGN_noCR_h206.format(ic=ic,
                              qt=qt, im=h1_imstr if qt == 'H1' else '',
                              snap=snhr)\
                            if 'h206' in ic else
                            tpl_m13_AGN_noCR_h113.format(ic=ic, 
                                qt=qt, im=h1_imstr if qt == 'H1' else '',
                                snap=snhr),
                          tpl_m13_AGN_CR.format(ic=ic, qt=qt, 
                                im=h1_imstr if qt == 'H1' else '', 
                                snap=snsr),
                         ]] \
                         for ic in ics_m13 for qt in qts\
                         for snhr, snsr in zip(snaps_hr, snaps_sr)]
            fnsss = fnss_m12 + fnss_m13
            fnsss = [[[mdir1 + fn if 'snap50' in fn or 'snap258' in fn else\
                       mdir2 + fn for fn in lst] for lst in _] for _ in fnsss]

            model_labelss = [[model_labels]] * len(fnsss)
            physmodel_labelss = [model_labels] * len(fnsss)

            rbinss = [np.linspace(0., 300., 50) if 'm12' in ic else \
                      np.linspace(0., 600., 50) if 'm13' in ic else \
                      None \
                      for ic in ics_m12 + ics_m13 for qt in qts\
                      for snhr, snsr in zip(snaps_hr, snaps_sr)]
            ylabels = [mlabel if 'mass' in qt else\
                       ilabel.format(ion=qt) \
                       for ic in  ics_m12 + ics_m13 for qt in qts\
                       for snhr, snsr in zip(snaps_hr, snaps_sr)]
            title_tpl= '{ion} profiles, {ic} z={z:.1f}, {ttpl}'
            titles = [title_tpl.format(ic=ic, ion=ion, z=z,
                                       ttpl=ttpl_m12 if 'm12' in ic else 
                                            ttpl_m13h206 if 'm13h206' in ic else\
                                            ttpl_m13h113)\
                      for ic in ics_m12 + ics_m13 for ion in qts\
                      for z in redshifts]
            
            outtemp = ('rprofcomp_clean_set1_set2_{ic}_{qt}'
                       '_z{z}_noBH_AGNCR_AGNnoCR.pdf')
            outnames = [outtemp.format(ic=ic, qt=qt,
                            z=('{:.1f}'.format(z)).replace('.', 'p')) \
                        for ic in ics_m12 + ics_m13 for qt in qts \
                        for z in redshifts]
            domean = False

    for fnss, model_labels, physmodel_labels, title, outname, \
            ylabel, rbins_pkpc in \
            zip(fnsss, model_labelss, physmodel_labelss, titles, outnames, 
                ylabels, rbinss):    
        outname = outdir + outname

        if title is not None:
            if len(title) > 70:
                splitopts = np.where([char == ',' for char in title])[0]
                optchoice = np.argmin(np.abs(splitopts - 0.5 * len(title)))
                splitind = splitopts[optchoice]
                title = (title[:splitind + 1]).strip() + '\n' +\
                        (title[splitind + 1:]).strip()
        compare_profiles_physmodels(fnss, model_labels,
                                    physmodel_labels,
                                    rbins_pkpc, title=title, 
                                    ylabel=ylabel, ax=None,
                                    outname=outname, domean=domean)

def compare_physmodel_z_panels(modelfiles_per_panel,
                               legendlablists_per_panel,
                               physmodel_labels_per_panel,
                               titles_per_panel,
                               ensemble_legend_labels,
                               rbins_pkpc_per_panel,
                               figtitle=None, ylabel=None,
                               outname=None):
    '''
    It's assumed modelfiles_per_panel and legendlablists_per_panel,
    and physmodel_labels_per_panel are lists of lists of strings.
    ensemble_legend_labels is a list of strings.
    '''
    nmodels = len(modelfiles_per_panel)
    ncols = min(3, nmodels + 1)
    nrows = (nmodels + 1) // 3 + 1
    panelsize = 2.5
    hspace = 0.25
    wspace = 0.25
    height_ratios = [panelsize] * nrows
    width_ratios = [panelsize] * ncols
    height = sum(height_ratios) * (1. + hspace * (ncols - 1.))
    width = sum(width_ratios) * (1. + wspace * (nrows - 1.))
    
    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(nrows=nrows, ncols=ncols, 
                        hspace=hspace, wspace=wspace, 
                        width_ratios=width_ratios, 
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols]) \
            for i in range(nmodels + 1)]
    li = ncols - (nmodels + 1) % ncols
    lax = fig.add_subplot(grid[nrows - 1, -1 * li:])

    fontsize = 12
    nlines = len(modelfiles_per_panel[0])
    _colors = tc.tol_cmap('rainbow_discrete', nlines)
    colors = _colors(np.linspace(1. / nlines, 1. - 1. / nlines, nlines))
    # zip/list comprehension issues or something when using 
    # tol_cmap outputs directly
    colors = [mcolors.to_rgb(col) for col in colors]
    modelcolors = tc.tol_cset('bright')

    if figtitle is not None:
        fig.suptitle(figtitle, fontsize=fontsize)

    for axi, (modelfiles, legendlablists, physmodel_labels, title, 
              rbins_pkpc, ax) in \
            enumerate(zip(modelfiles_per_panel, legendlablists_per_panel, 
                      physmodel_labels_per_panel, titles_per_panel, 
                      rbins_pkpc_per_panel, axes)):
        if axi % ncols == 0:
            _ylabel = ylabel
        else:
            _ylabel = None
        compare_profiles_physmodels([modelfiles],
                                    [legendlablists],
                                    physmodel_labels,
                                    rbins_pkpc, title=None, 
                                    ylabel=_ylabel, ax=ax,
                                    outname=None, domean=False,
                                    doscatter=False, colors=colors)
        ax.set_title(title, fontsize=fontsize)
        # legend in its own panel
        if axi == 0:
            handles2, labels2 = ax.get_legend_handles_labels()
        ax.get_legend().remove() # remove annoying empty rectangle

        if axi >= nmodels + 1 - ncols:
            ax.set_xlabel('$\\mathrm{r}_{\\perp} \\, [\\mathrm{pkpc}]$',
                          fontsize=fontsize)
        else:
            ax.set_xlabel(None)
        
        av, med, p90, p10 = gpr.get_profile_massmap(
                modelfiles, rbins_pkpc, rbin_units='pkpc', 
                profiles=['av-lin', 'perc-0.5', 'perc-0.9', 'perc-0.1'])
        rcens = 0.5 * (rbins_pkpc[:-1] + rbins_pkpc[1:])
        ax.plot(rcens, med, color='black', linewidth=2, linestyle='solid')
        ax.fill_between(rcens, p10, p90, color='black', alpha=0.3)

        c_last = modelcolors[axi]
        axes[nmodels].plot(rcens, med, color=c_last, linewidth=2, 
                           linestyle='solid', 
                           label=ensemble_legend_labels[axi])
        axes[nmodels].fill_between(rcens, p10, p90, color=c_last, alpha=0.3)
    
    if nmodels % ncols == 0:
        axes[nmodels].set_ylabel(ylabel, fontsize=fontsize)
    axes[nmodels].set_xlabel('$\\mathrm{r}_{\\perp} \\, [\\mathrm{pkpc}]$',
                             fontsize=fontsize)
    axes[nmodels].tick_params(direction='in', which='both', 
                              labelsize=fontsize - 1)
    axes[nmodels].legend(fontsize=fontsize, bbox_to_anchor=(1.05, 0.00),
                         loc='lower left')
    axes[nmodels].set_xscale('log')
    axes[nmodels].set_title('ensemble profiles')
    
    # sync y ranges
    ylims = [ax.get_ylim() for ax in axes]
    ymin = min([ylim[0] for ylim in ylims])
    if 'gas-mass' not in modelfiles_per_panel[0][0]:
        ymin = max(ymin, 10.)
    ymax = max([ylim[1] for ylim in ylims])    
    [ax.set_ylim((ymin, ymax)) for ax in axes]

    lax.axis('off')
    line = [[(0, 0)]]
    lc = mcol.LineCollection(line * len(colors), linestyle='solid', 
                             linewidth=2., colors=colors)
    handles1 = [lc,
                mlines.Line2D((), (), linestyle='solid', linewidth=2,
                              label='all med.', color='black'),
                mpatch.Patch(label='all perc. 10-90', linewidth=0.5, 
                             color='black', alpha=0.3),
                mlines.Line2D((), (), linestyle=None, marker='o',
                              label='$\\mathrm{R}_{\\mathrm{vir}}$', 
                              color='gray', markersize=5)
                ]
    labels1 = ['median', 'all med.', 'all perc. 10-90', 
               '$\\mathrm{R}_{\\mathrm{vir}}$']
    lax.legend(handles1 + handles2, labels1 + labels2, 
               fontsize=fontsize, ncols=2,
               handler_map={type(lc): pu.HandlerDashedLines()},
               bbox_to_anchor=(0.95, 0.95), loc='upper right')

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def compare_physmodel_z_mapsets(mapset='clean_set1-2'):

    if mapset == 'clean_set1-2':
        outdir = ('/Users/nastasha/ciera/projects_lead/fire3_ionabs/'
                  'profiles_clean/')
        mdir1 = '/Users/nastasha/ciera/sim_maps/fire/clean_set1/'
        mdir2 = '/Users/nastasha/ciera/sim_maps/fire/clean_set2/'
        tpl_m12_noBH =          ('coldens_{qt}_{ic}_m7e3_MHD_fire3'
                                 '_fireBH_Sep182021_hr_crdiffc690_sdp1e10'
                                 '_gacc31_fa0.5_snap{snap}'
                                 '_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5')
        tpl_m12_AGN_CR =        ('coldens_{qt}_{ic}_m6e4_MHDCRspec1_fire3'
                                 '_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
                                 '_gacc31_fa0.5_fcr1e-3_vw3000_snap{snap}'
                                 '_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5')
        tpl_m12_AGN_noCR =      ('coldens_{qt}_{ic}_m7e3_MHD_fire3'
                                 '_fireBH_Sep182021_hr_crdiffc690_sdp2e-4'
                                 '_gacc31_fa0.5_snap{snap}'
                                 '_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5')
        tpl_m13_noBH =          ('coldens_{qt}_{ic}_m3e5_MHD_fire3'
                                 '_fireBH_Sep182021_crdiffc690_sdp1e10'
                                 '_gacc31_fa0.5_snap{snap}'
                                 '_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5') 
        tpl_m13_AGN_CR =        ('coldens_{qt}_{ic}_m3e5_MHDCRspec1_fire3'
                                 '_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
                                 '_gacc31_fa0.5_fcr1e-3_vw3000_snap{snap}'
                                 '_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5') 
        tpl_m13_AGN_noCR_h113 = ('coldens_{qt}_{ic}_m3e4_MHD_fire3'
                                 '_fireBH_Sep182021_hr_crdiffc690_sdp1e-4'
                                 '_gacc31_fa0.5_snap{snap}'
                                 '_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5') 
        tpl_m13_AGN_noCR_h206 = ('coldens_{qt}_{ic}_m3e4_MHD_fire3_fireBH'
                                 '_Sep182021_hr_crdiffc690_sdp3e-4'
                                 '_gacc31_fa0.5_snap{snap}'
                                 '_shrink-sph-cen_BN98_2rvir{im}_v2.hdf5') 
        model_labels = ['no BH', 'AGN, no CR', 'AGN + CR']
        ics_m12 = ['m12f'] # bug: 'm12m'
        ics_m13 = ['m13h206', 'm13h113']
        qts = ['gas-mass', 'O6', 'Ne8', 'Mg10', 'H1']
        h1_imstr = '_ionfrac-fromsim'
        # low to high z -> colored blue to red
        snaps_sr = [50, 49, 48, 47, 46, 45]
        snaps_hr =  [258, 240, 224, 210, 197, 186]
        redshifts = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

        mlabel = ('$\\log_{10} \\, \\Sigma_{\\mathrm{gas}}'+\
                  ' \\, [\\mathrm{g} \\, \\mathrm{cm}^2]$')
        ilabel = ('$\\log_{{10}} \\, \\mathrm{{N}}(\\mathrm{{{ion}}})'+\
                  ' \\, [\\mathrm{{cm}}^2]$')
        
        ics = ics_m12 + ics_m13
        ensemble_legend_labels = ['all noBH', 'all AGN-noCR', 'all AGN-CR']
        titles_per_panel_m12  = ['noBH: m7e3', 'AGN-noCR: m7e3, sdp2e-4',
                                 'AGN-CR:\nm6e4, sdp1e-4, fcr1e-3_vw3000'] 
        titles_per_panel_m13h113 = ['noBH: m3e5', 'AGN-noCR: m3e4, sdp1e-4',
                                    'AGN-CR:\nm3e5, sdp1e-4, fcr1e-3_vw3000'] 
        titles_per_panel_m13h206 = ['noBH: m3e5', 'AGN-noCR: m3e4, sdp3e-4',
                                    'AGN-CR:\nm3e5, sdp1e-4, fcr1e-3_vw3000']

        fnss_m12 = [[[tpl_m12_noBH.format(ic=ic, qt=qt, 
                            im=h1_imstr if qt == 'H1' else '', 
                            snap=snhr)\
                      for snhr in snaps_hr],
                     [tpl_m12_AGN_noCR.format(ic=ic, qt=qt, 
                            im=h1_imstr if qt == 'H1' else '', 
                            snap=snhr)\
                      for snhr in snaps_hr],
                     [tpl_m12_AGN_CR.format(ic=ic, qt=qt, 
                            im=h1_imstr if qt == 'H1' else '', 
                            snap=snsr)\
                      for snsr in snaps_sr],
                     ] \
                    for ic in ics_m12 for qt in qts]
        fnss_m13 = [[[tpl_m13_noBH.format(ic=ic, qt=qt,
                            im=h1_imstr if qt == 'H1' else '',
                            snap=snsr)\
                      for snsr in snaps_sr],
                     [tpl_m13_AGN_noCR_h206.format(ic=ic,
                            qt=qt, im=h1_imstr if qt == 'H1' else '',
                            snap=snhr)\
                      if 'h206' in ic else \
                      tpl_m13_AGN_noCR_h113.format(ic=ic, 
                            qt=qt, im=h1_imstr if qt == 'H1' else '',
                            snap=snhr) \
                      for snhr in snaps_hr],
                     [tpl_m13_AGN_CR.format(ic=ic, qt=qt, 
                            im=h1_imstr if qt == 'H1' else '', 
                            snap=snsr)\
                      for snsr in snaps_sr],
                     ] \
                     for ic in ics_m13 for qt in qts]
        fnsss = fnss_m12 + fnss_m13
        fnsss = [[[mdir1 + fn if 'snap50' in fn or 'snap258' in fn else\
                   mdir2 + fn for fn in lst] for lst in _] for _ in fnsss]

        legendlabss = [[['z={:.1f}'.format(z) for z in redshifts]\
                        if i == 0 else [None] * len(redshifts) \
                        for i in range(len(fnsss[0]))] \
                       for ic in ics for qt in qts]
        physmodel_labelss = [[['z={:.1f}'.format(z) for z in redshifts]
                              for i in range(len(fnsss[0]))] \
                             for ic in ics_m12 + ics_m13 for qt in qts]
        
        rbinss = [[np.linspace(0., 300., 50) if 'm12' in ic else \
                   np.linspace(0., 600., 50) if 'm13' in ic else \
                   None] * len(fnsss[0]) \
                  for ic in ics for qt in qts]
        ylabels = [mlabel if 'mass' in qt else\
                   ilabel.format(ion=qt) \
                   for ic in  ics for qt in qts]
        title_tpl= '{qt} profiles, {ic}'
        figtitles = [title_tpl.format(ic=ic, qt=qt)\
                      for ic in ics for qt in qts]
        ensemble_legend_labelss = [ensemble_legend_labels \
                                   for ic in ics for qt in qts]
        ptitles_list = [titles_per_panel_m12 if 'm12' in ic else\
                        titles_per_panel_m13h113 if ic == 'm13h113' else\
                        titles_per_panel_m13h206 if ic == 'm13h206' else\
                        None\
                        for ic in ics for qt in qts]
        outtemp = ('rprofcomp_clean_set1_set2_{ic}_{qt}'
                   '_z0p5-1_noBH_AGNCR_AGNnoCR.pdf')
        outnames = [outtemp.format(ic=ic, qt=qt) \
                    for ic in ics for qt in qts]

    for fnss, physmodel_labels, title, outname, ptitles,\
            ylabel, rbins_pkpc, legendlabs, enslabs in \
            zip(fnsss, physmodel_labelss, figtitles, outnames, ptitles_list,
                ylabels, rbinss, legendlabss, ensemble_legend_labelss):    
        outname = outdir + outname

        if title is not None:
            if len(title) > 70:
                splitopts = np.where([char == ',' for char in title])[0]
                optchoice = np.argmin(np.abs(splitopts - 0.5 * len(title)))
                splitind = splitopts[optchoice]
                title = (title[:splitind + 1]).strip() + '\n' +\
                        (title[splitind + 1:]).strip()
        compare_physmodel_z_panels(fnss,
                                  legendlabs,
                                  physmodel_labels,
                                  ptitles,
                                  enslabs,
                                  rbins_pkpc,
                                  figtitle=title, ylabel=ylabel,
                                  outname=outname)