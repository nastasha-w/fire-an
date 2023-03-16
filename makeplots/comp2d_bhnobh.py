
import h5py
import numpy as np

import matplotlib.lines as mlines
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt

import fire_an.makeplots.get_2dprof as gpr
import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
import fire_an.utils.constants_and_units as c

def compare_profiles_BHnoBH(bhfiles, nobhfiles,
                            legendlab_bhfs, legendlab_nobhfs,
                            rbins_pkpc, title=None, 
                            ylabel=None, ax=None,
                            outname=None):
    fontsize = 12
    if ax is None:
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5.5, 5.))
        if title is not None:
            fig.suptitle(title, fontsize=fontsize)
        if ylabel is not None:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.set_xlabel('$\\mathrm{r}_{\\perp}$ [pkpc]', fontsize=fontsize)

    ls_bh = 'solid'
    ls_nobh = 'dashed'
    lw_med = 1.
    lw_av = 2.5
    alpha_range = 0.3
    colors = tc.tol_cset('bright')

    rcens = 0.5 * (rbins_pkpc[1:] + rbins_pkpc[:-1])
    
    for fi, (bhfn, nobhfn, bhlab, nobhlab) in \
            enumerate(zip(bhfiles, nobhfiles, 
                          legendlab_bhfs, legendlab_nobhfs)):
        with h5py.File(bhfn, 'r') as f:  
            rvir_bh_pkpc = f['Header/inputpars/halodata'].attrs['Rvir_cm']
            rvir_bh_pkpc /= (c.cm_per_mpc * 1e-3)
        with h5py.File(nobhfn, 'r') as f:  
            rvir_nobh_pkpc = f['Header/inputpars/halodata'].attrs['Rvir_cm']
            rvir_nobh_pkpc /= (c.cm_per_mpc * 1e-3)
        
        bh_av, bh_med, bh_90, bh_10 = gpr.get_profile_massmap(
            bhfn, rbins_pkpc, rbin_units='pkpc', 
            profiles=['av-lin', 'perc-0.5', 'perc-0.9', 'perc-0.1'])
        
        nb_av, nb_med, nb_90, nb_10 = gpr.get_profile_massmap(
            nobhfn, rbins_pkpc, rbin_units='pkpc', 
            profiles=['av-lin', 'perc-0.5', 'perc-0.9', 'perc-0.1'])
        
        color_bh = colors[2 * fi]
        color_nobh = colors[2 * fi + 1]
        ax.plot(rcens, bh_av, color=color_bh, linestyle=ls_bh, 
                linewidth=lw_av, label=bhlab)
        ax.plot(rcens, nb_av, color=color_nobh, linestyle=ls_nobh, 
                linewidth=lw_av, label=nobhlab)
        ax.plot(rcens, bh_med, color=color_bh, linestyle=ls_bh, 
                linewidth=lw_med)
        ax.plot(rcens, nb_med, color=color_nobh, linestyle=ls_nobh, 
                linewidth=lw_med)
        ax.fill_between(rcens, bh_10, bh_90, color=color_bh, 
                        alpha=alpha_range, linestyle=ls_bh,
                        linewidth=0.5)
        ax.fill_between(rcens, nb_10, nb_90, color=color_nobh, 
                        alpha=alpha_range, linestyle=ls_nobh,
                        linewidth=0.5)
        # indicate Rvir on the curves
        yp_bh_av = pu.linterpsolve(rcens, bh_av, rvir_bh_pkpc)
        yp_bh_med = pu.linterpsolve(rcens, bh_med, rvir_bh_pkpc)
        ax.scatter([rvir_bh_pkpc] * 2, [yp_bh_av, yp_bh_med], marker='o',
                   c=color_bh, s=15)
        yp_nb_av = pu.linterpsolve(rcens, nb_av, rvir_nobh_pkpc)
        yp_nb_med = pu.linterpsolve(rcens, nb_med, rvir_nobh_pkpc)
        ax.scatter([rvir_nobh_pkpc] * 2, [yp_nb_av, yp_nb_med], marker='o',
                   c=color_nobh, s=15)

    handles2, labels = ax.get_legend_handles_labels()
    handles1 = [mlines.Line2D((), (), linewidth=lw_med, linestyle=ls_bh,
                              label='med., w/ BH', color='black'),
                mlines.Line2D((), (), linewidth=lw_med, linestyle=ls_nobh,
                              label='med., no BH', color='black'),
                mlines.Line2D((), (), linewidth=lw_av, linestyle=ls_bh,
                              label='mean, w/ BH', color='black'),
                mlines.Line2D((), (), linewidth=lw_av, linestyle=ls_nobh,
                              label='mean, no BH', color='black'),
                mpatch.Patch(label='perc. 10-90', linewidth=0.5, 
                             color='black', alpha=alpha_range),
                mlines.Line2D((), (), linestyle=None, marker='o',
                              label='$\\mathrm{R}_{\\mathrm{vir}}$', 
                              color='black', markersize=5)
                ]
    ax.legend(handles=handles1 + handles2, fontsize=fontsize)
    ax.set_xscale('log')

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')
        
def compare_profilesets(mapset, compmode):

    if mapset == 4:
        outdir = ('/Users/nastasha/ciera/projects_lead/'
                  'fire3_ionabs/profiles_BH_noBH/')
        qts = ['gas-mass', 'Mg10', 'Ne8', 'O6', 'N5', 'Mg2']
        snaps = [500, 258, 186]
        bhmodes = {'nobh': 'sdp1e10',
                   'bh': 'sdp2e-4'}
        ics = ['m12f', 'm12m']
        
        fdir = '/Users/nastasha/ciera/sim_maps/fire/set4_BH_noBH/'
        filetemp = ('coldens_{qt}_{ic}_m7e3_MHD_fire3_fireBH_Sep182021'
                    '_hr_crdiffc690_{bh}_gacc31_fa0.5_snap{snap}_'
                    'shrink-sph-cen_BN98_2rvir_v2.hdf5')
        
        if compmode == 'indiv':
            filesets_bh = [[fdir + filetemp.format(bh=bhmodes['bh'], 
                                                   snap=snap,
                                                   ic=ic, qt=qt)] \
                            for ic in ics for snap in snaps for qt in qts]
            filesets_nobh = [[fdir + filetemp.format(bh=bhmodes['nobh'], 
                                                     snap=snap,
                                                     ic=ic, qt=qt)] \
                            for ic in ics for snap in snaps for qt in qts]
            ttemp = ('set 4, {qt}, snapshot {snap}, {ic}_m7e3,\n'
                     'MHD_fire3_fireBH_Sep18202_hr_crdiffc690, gacc31_fa0.5'
                     )
            titles = [ttemp.format(snap=snap, ic=ic, qt=qt)\
                      for ic in ics for snap in snaps for qt in qts]
            lls_bh = [['BH']] * len(filesets_bh)
            lls_nobh = [['no BH']] * len(filesets_bh)
            rbinss = [np.linspace(0., 300., 50) if 'm12' in ic else \
                      np.linspace(0., 600., 50) if 'm13' in ic else \
                      None \
                      for ic in ics for snap in snaps for qt in qts]
            mlabel = ('$\\log_{10} \\, \\Sigma_{\\mathrm{gas}}'
                      ' \\, [\\mathrm{g} \\, \\mathrm{cm}^2]$')
            ilabel = ('$\\log_{{10}} \\, \\mathrm{{N}}(\\mathrm{{{ion}}})'
                      ' \\, [\\mathrm{{cm}}^2]$')
            ylabels = [mlabel if 'mass' in qt else\
                       ilabel.format(ion=qt) \
                       for ic in ics for snap in snaps for qt in qts]
            outtemp = 'rprofcomp_set4_{ic}_{qt}_snap{snap}_BHnoBH.pdf'
            outnames = [outdir + outtemp.format(ic=ic, qt=qt, snap=snap) \
                       for ic in ics for snap in snaps for qt in qts]

    if mapset in [5, 6]:
        outdir = ('/Users/nastasha/ciera/projects_lead/'
                  'fire3_ionabs/profiles_BH_noBH/')
        qts = ['gas-mass', 'Mg10', 'Ne8', 'O6', 'N5', 'Mg2']
        if mapset == 5:
           snaps = [60, 50, 45]
           ics = ['m13h007', 'm13h206']
        else:
            snaps = [50, 45]
            ics = ['m13h002']        
        fdir = '/Users/nastasha/ciera/sim_maps/fire/set{}_BH_noBH/'
        fdir = fdir.format(mapset)
        filetemp_nobh = ('coldens_{qt}_{ic}_m3e5_MHD_fire3_fireBH_Sep182021'
                         '_crdiffc690_sdp1e10_gacc31_fa0.5_snap{snap}_shrink'
                         '-sph-cen_BN98_2rvir_v2.hdf5')
        filetemp_bh = ('coldens_{qt}_{ic}_m3e5_MHDCRspec1_fire3_fireBH'
                       '_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5'
                       '_fcr1e-3_vw3000_snap{snap}_shrink-sph-cen_BN98'
                       '_2rvir_v2.hdf5')
        
        if compmode == 'indiv':
            filesets_bh = [[fdir + filetemp_bh.format(snap=snap,
                                                      ic=ic, qt=qt)] \
                            for ic in ics for snap in snaps for qt in qts]
            filesets_nobh = [[fdir + filetemp_nobh.format(snap=snap,
                                                          ic=ic, qt=qt)] \
                            for ic in ics for snap in snaps for qt in qts]
            ttemp = ('set {set_}, {qt}, snapshot {snap}, {ic}_m3e5,\n'
                     'no BH: MHD_fire3_fireBH_Sep182021_crdiffc690,'
                     ' gacc31_fa0.5\n'
                     'BH: MHDCRspec1, fireCR1_Oct252021_crdiffc1_sdp1e-4,'
                     ' fcr1e-3_vw3000'
                     )
            titles = [ttemp.format(set_=mapset, snap=snap, ic=ic, qt=qt)\
                      for ic in ics for snap in snaps for qt in qts]
            lls_bh = [['BH']] * len(filesets_bh)
            lls_nobh = [['no BH']] * len(filesets_bh)
            rbinss = [np.linspace(0., 300., 50) if 'm12' in ic else \
                      np.linspace(0., 700., 50) if 'm13' in ic else \
                      None \
                      for ic in ics for snap in snaps for qt in qts]
            mlabel = ('$\\log_{10} \\, \\Sigma_{\\mathrm{gas}}'
                      ' \\, [\\mathrm{g} \\, \\mathrm{cm}^2]$')
            ilabel = ('$\\log_{{10}} \\, \\mathrm{{N}}(\\mathrm{{{ion}}})'
                      ' \\, [\\mathrm{{cm}}^2]$')
            ylabels = [mlabel if 'mass' in qt else\
                       ilabel.format(ion=qt) \
                       for ic in ics for snap in snaps for qt in qts]
            outtemp = 'rprofcomp_set4_{ic}_{qt}_snap{snap}_BHnoBH.pdf'
            outnames = [outdir + outtemp.format(ic=ic, qt=qt, snap=snap) \
                        for ic in ics for snap in snaps for qt in qts]

    for fs_bh, fs_nobh, ll_bh, ll_nobh, \
            rbins_pkpc, title, ylabel, outname in \
            zip(filesets_bh, filesets_nobh, lls_bh, lls_nobh,
                rbinss, titles, ylabels, outnames):
        compare_profiles_BHnoBH(fs_bh, fs_nobh, ll_bh, ll_nobh,
                                rbins_pkpc, title=title, 
                                ylabel=ylabel, ax=None,
                                outname=outname)