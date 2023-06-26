'''
First version.
'''

import h5py
import numpy as np
import os
import pandas as pd

import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt

import fire_an.makeplots.get_2dprof as gpr
import fire_an.makeplots.tol_colors as tc
import fire_an.makeplots.plot_utils as pu
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu
import fire_an.utils.opts_locs as ol


def plot_radprof_m12i_CR_comp(smallrange=True):
    '''
    Rough comparison to Ji, Chan, et al. (2020)
    '''
    if smallrange:
        rbins = np.linspace(0., 0.8, 16)
        yranges = {'si4': (10.8, 15.1), 'n5': (10.9, 14.7), 
                   'o6': (13., 15.3), 'ne8': (12.8, 15.1)}
        rcens = 0.5 * (rbins[:-1] + rbins[1:])
    else:
        rbins = np.append([-3.5], np.linspace(-1., 0.3, 14))
        rbins = 10**rbins
        rcens = 0.5 * (rbins[:-1] + rbins[1:])
        yranges = {'si4': (6., 15.1), 'n5': (6., 14.7), 
                   'o6': (11., 15.3), 'ne8': (10.4, 14.7)}
        
    axions = {'si4': 0, 'n5': 1, 'o6': 2, 'ne8': 3}
    snapcolors = {277: 'green', 600: 'black'}
    snaplabels = {277: 'z=1', 600: 'z=0'}
    ions = list(axions.keys())
    snapshots = list(snaplabels.keys())
    fdir = '/Users/nastasha/ciera/tests/fire_start/map_tests/'
    filens = ['coldens_n5_m12i_noAGNfb_CR-diff-coeff-690_FIRE-2_snap277_shrink-sph-cen_BN98_2rvir_v1.hdf5',
              'coldens_n5_m12i_noAGNfb_CR-diff-coeff-690_FIRE-2_snap600_shrink-sph-cen_BN98_2rvir_v1.hdf5',
              'coldens_ne8_m12i_noAGNfb_CR-diff-coeff-690_FIRE-2_snap277_shrink-sph-cen_BN98_2rvir_v1.hdf5',
              'coldens_ne8_m12i_noAGNfb_CR-diff-coeff-690_FIRE-2_snap600_shrink-sph-cen_BN98_2rvir_v1.hdf5',
              'coldens_o6_m12i_noAGNfb_CR-diff-coeff-690_FIRE-2_snap277_shrink-sph-cen_BN98_2rvir_v1.hdf5',
              'coldens_o6_m12i_noAGNfb_CR-diff-coeff-690_FIRE-2_snap600_shrink-sph-cen_BN98_2rvir_v1.hdf5',
              'coldens_si4_m12i_noAGNfb_CR-diff-coeff-690_FIRE-2_snap277_shrink-sph-cen_BN98_2rvir_v1.hdf5',
              'coldens_si4_m12i_noAGNfb_CR-diff-coeff-690_FIRE-2_snap600_shrink-sph-cen_BN98_2rvir_v1.hdf5',
              ]
    filens = [fdir + filen for filen in filens]

    fig = plt.figure(figsize=(5.5, 5.))
    grid = gsp.GridSpec(nrows=2, ncols=2, hspace=0.3, wspace=0.5)
    axes = [fig.add_subplot(grid[i // 2, i % 2]) for i in range(4)]
    fontsize = 12
    
    title = 'm12i with CRs, linear average and full\nrange column densities, z-projection'
    fig.suptitle(title, fontsize=fontsize)

    for filen in filens:   
        rd, cd = gpr.get_rval_massmap(filen, units='Rvir')
        rinds = np.searchsorted(rbins, rd) - 1
        cd_by_bin = [cd[rinds == i] for i in range(len(rbins))]
        cd_av = np.log10(np.array([np.average(10**cds) for cds in cd_by_bin]))
        cd_range = np.array([[np.min(cds), np.max(cds)] for cds in cd_by_bin])
        cd_range[cd_range == -np.inf] = -100.

        ion = ions[np.where(['coldens_{}'.format(ion) in filen \
                   for ion in ions])[0][0]] 
        snapnum = snapshots[np.where(['snap{}'.format(snap) in filen \
                             for snap in snapshots])[0][0]] 
        print(ion, snapnum)
        ax = axes[axions[ion]]
        snaplabel = snaplabels[snapnum]
        color = snapcolors[snapnum]

        ax.plot(rcens, cd_av[:-1], color=color, label=snaplabel)
        ax.fill_between(rcens, cd_range[:-1, 0], cd_range[:-1, 1], color=color,
                        alpha=0.3)
    for ion in ions:
        print(ion)
        print(axions[ion])
        ax = axes[axions[ion]]
        ax.set_xlabel('$r_{\perp} \\; [\\mathrm{R}_{\\mathrm{vir}}]$', 
                      fontsize=fontsize)
        ax.set_ylabel('$\\log_{10} \\, \\mathrm{N} \\; [\\mathrm{cm}^{-2}]$',
                      fontsize=fontsize)
        ax.tick_params(labelsize=fontsize - 1, direction='in', which='both')
        ax.set_ylim(yranges[ion])
        if not smallrange:
            ax.set_xscale('log')
        ax.text(0.05, 0.05, ion, fontsize=fontsize,
                horizontalalignment='left', verticalalignment='bottom',
                transform=ax.transAxes)

    axes[axions[ions[0]]].legend(fontsize=fontsize)
    outname = fdir + 'radprof_coldens_ji-chan-etal-2020_comp_m12i' + \
                     '_noAGNfb_CR-diff-coeff-690_FIRE-2'
    outname = outname + '_smallrad.pdf' if smallrange else \
              outname + '_largerad.pdf'
    plt.savefig(outname, bbox_inches='tight')

def plotMz_burchett_etal_2019_cleansample():
    datadir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/'
    dfilen = datadir + 'data_burchett_etal_2019_table1.txt'
    #TODO CHECK: R200c or R200m!
    # assuming impact parameters are physical/proper kpc
    #TODO ask: table footnote f says 2 systems for one Ne VIII absorber
    #          but only one line with that footnote and N(Ne VIII) value
    data_bur = pd.read_csv(dfilen, comment='#', sep='\t')
    ## calculate halo masses
    # from Burchett et al. (2019):
    cosmopars_bur = {'h': 0.677, 'omegam': 0.31, 'omegalambda': 0.69}
    def hmfunc(x):
        csm = cosmopars_bur.copy()
        csm.update({'z': x.zgal, 'a': 1. / (1. + x.zgal)})
        mv = cu.mvir_from_rvir(x.rvir_kpc * 1e-3 * c.cm_per_mpc, 
                               csm, meandef='200m')
        return mv / c.solar_mass
    data_bur = data_bur.assign(Mvir_Msun=lambda x: hmfunc(x))
    
    ## FIRE data
    snapfiles = [('m13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                  '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                 ('m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr'
                  '_crdiffc690_sdp1e-4_gacc31_fa0.5'),
                 ('m13h113_m3e5_MHD_fire3_fireBH_Sep182021'
                  '_crdiffc690_sdp1e10_gacc31_fa0.5'),
                 ('m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                  '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                 ('m13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr'
                  '_crdiffc690_sdp3e-4_gacc31_fa0.5'),
                 ('m13h206_m3e5_MHD_fire3_fireBH_Sep182021'
                  '_crdiffc690_sdp1e10_gacc31_fa0.5'),
                 ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                  '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                 ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr'
                  '_crdiffc690_sdp2e-4_gacc31_fa0.5'),
                 ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr'
                  '_crdiffc690_sdp1e10_gacc31_fa0.5'),
                 ]
    snaplabels = ['m13h113 AGN-CR', 'm13h113 AGN-noCR', 'm13h113 noBH',
                  'm13h206 AGN-CR', 'm13h206 AGN-noCR', 'm13h206 noBH',
                  'm12f AGN-CR', 'm12f AGN-noCR', 'm12f noBH',
                  ]
    firedataf = ol.filen_halocenrvir
    firemasses = []
    fireredshifts = []
    fireradii = []
    firecens = []
    with h5py.File(firedataf, 'r') as f:
        meandef = '200m'
        for sfn in snapfiles:
            _lsm = []
            _lsz = []
            _lsr = []
            _lfc = []
            smgrp = f[sfn]
            sngrpns = [key for key in smgrp.keys() if key.startswith('snap_')]
            for sngrpn in sngrpns:
                sngrp = smgrp[sngrpn]
                _lsz.append(sngrp['cosmopars'].attrs['z'])
                _lsm.append(sngrp[f'cen0/Rvir_{meandef}'].attrs['Mvir_g'])
                _lsr.append(sngrp[f'cen0/Rvir_{meandef}'].attrs['Rvir_cm'])
                _lfc.append(np.array([sngrp[f'cen0'].attrs[f'{_ax}c_cm']
                                      for _ax in ['X', 'Y', 'Z']]))
            zo = np.argsort(_lsz)
            _lsz = np.array(_lsz)[zo]
            _lsm = np.array(_lsm)[zo]
            _lsr = np.array(_lsr)[zo]
            _lfc = np.array(_lfc)[zo]
            _lsm /= c.solar_mass
            _lsr /= (c.cm_per_mpc * 1e-3)
            _lfc /= (c.cm_per_mpc * 1e-3)
            firemasses.append(_lsm)
            fireredshifts.append(_lsz)
            fireradii.append(_lsr)
            firecens.append(_lfc)
    cset = tc.tol_cset('bright')
    pcolors = {'AGN-CR': cset.green,
               'AGN-noCR': cset.red,
               'noBH': cset.blue,
               }
    hstyles = {'m13h113': 'solid',
               'm13h206': 'dashdot',
               'm12f': 'dashed',
               }
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5.5, 5.))
    fontsize = 12
    # Burchett et al. data
    z_bur = data_bur['zgal']
    m_bur = data_bur['Mvir_Msun']
    isul = data_bur['log_N_Ne8_isUL']
    ax.scatter(z_bur[isul], m_bur[isul], marker='o', linestyle='None',
               edgecolor='black', facecolor='none', label='Bur.+19 (UL)',
               s=40)
    ax.scatter(z_bur[np.logical_not(isul)], m_bur[np.logical_not(isul)], 
               marker='o', linestyle='None', edgecolor='black', 
               facecolor='black', label='Bur.+19', s=40)
    #for snapn, zs, ms, rs, cs in zip(snapfiles, fireredshifts, firemasses,
    #                                 fireradii, firecens):
    #    print()
    #    print(f'{snapn}:')
    #    print(f'z:\t ', '\t'.join([f'{x:.2f}' for x in zs]))
    #    print(f'M:\t ', '\t'.join([f'{x:.4e}' for x in ms]))
    #    print(f'R:\t ', '\t'.join([f'{x:.2f}' for x in rs]))
    #    print(f'C:\t ', '\t'.join([', '.join([f'{x:.2f}' for x in y]) 
    #                               for y in cs]))
    # FIRE data
    mass_minmax = {'m13': (np.inf, -np.inf),
                   'm12': (np.inf, -np.inf)}
    z_minmax = {'m13': (np.inf, -np.inf),
                'm12': (np.inf, -np.inf)}
    zmar = 0.05
    mmar = 0.2
    for snaplabel, firemass, fireredshift \
            in zip(snaplabels, firemasses, fireredshifts):
        hlabel, plabel = snaplabel.split(' ')
        ax.plot(fireredshift, firemass, color=pcolors[plabel],
                linestyle=hstyles[hlabel], linewidth=1.5,
                marker='o', markersize=4)
        for key in mass_minmax:
            if snaplabel.startswith(key):
                _prev = mass_minmax[key]
                mmin = min(_prev[0], np.min(firemass))
                mmax = max(_prev[1], np.max(firemass))
                mass_minmax[key] = (mmin, mmax)
                _prev = z_minmax[key]
                zmin = min(_prev[0], np.min(fireredshift))
                zmax = max(_prev[1], np.max(fireredshift))
                z_minmax[key] = (zmin, zmax)
    mass_minmax = {key: (10**(np.log10(mass_minmax[key][0]) - mmar),
                         10**(np.log10(mass_minmax[key][1]) + mmar))
                   for key in mass_minmax}
    z_minmax = {key: (z_minmax[key][0] - zmar, z_minmax[key][1] + zmar)
                for key in z_minmax}
    print(mass_minmax)
    print(z_minmax)
    for key in mass_minmax:
        line1 = [mass_minmax[key][0], mass_minmax[key][1], 
                 mass_minmax[key][1], mass_minmax[key][0], 
                 mass_minmax[key][0]]
        line0 = [z_minmax[key][0], z_minmax[key][0], 
                 z_minmax[key][1], z_minmax[key][1], 
                 z_minmax[key][0]]
        ax.plot(line0, line1, linestyle='solid', color=cset.grey,
                linewidth=2)
        ax.text(z_minmax[key][1] - 0.02, mass_minmax[key][1] * 0.95,
                key, fontsize=fontsize, color='gray',
                horizontalalignment='right', verticalalignment='top')

    ax.set_yscale('log')
    ax.set_ylabel(f'$\\mathrm{{M}}_{{\\mathrm{{{meandef}}}}} \\;'
                   ' [\\mathrm{M}_{\\odot}]$', fontsize=fontsize)
    ax.set_xlabel('redshift', fontsize=fontsize)
    ax.tick_params(labelsize=fontsize - 1., direction='in', which='both')
    
    _handles, _ = ax.get_legend_handles_labels()
    handles1 = [mlines.Line2D((), (), linewidth=1.5, linestyle=hstyles[key],
                              label=key, color='black') \
                for key in hstyles]
    handles2 = [mlines.Line2D((), (), linewidth=1.5, color=pcolors[key],
                              label=key, linestyle='solid', marker='o',
                              markersize=4)\
                for key in pcolors]
    handles = _handles + handles1 + handles2
    ax.legend(handles=handles, fontsize=fontsize - 1, handlelength=2.0)

    outdir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/datacomp/'
    outfilen = outdir + 'mass_z_selection_clean.pdf'
    plt.savefig(outfilen, bbox_inches='tight')

def plotMz_burchett_etal_2019_model3():
    datadir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/'
    dfilen = datadir + 'data_burchett_etal_2019_table1.txt'
    #TODO CHECK: R200c or R200m!
    # assuming impact parameters are physical/proper kpc
    #TODO ask: table footnote f says 2 systems for one Ne VIII absorber
    #          but only one line with that footnote and N(Ne VIII) value
    data_bur = pd.read_csv(dfilen, comment='#', sep='\t')
    ## calculate halo masses
    # from Burchett et al. (2019):
    cosmopars_bur = {'h': 0.677, 'omegam': 0.31, 'omegalambda': 0.69}
    def hmfunc(x):
        csm = cosmopars_bur.copy()
        csm.update({'z': x.zgal, 'a': 1. / (1. + x.zgal)})
        mv = cu.mvir_from_rvir(x.rvir_kpc * 1e-3 * c.cm_per_mpc, 
                               csm, meandef='200m')
        return mv / c.solar_mass
    data_bur = data_bur.assign(Mvir_Msun=lambda x: hmfunc(x))
    
    ## FIRE data
    m13_nobh = sl.m13_nobh_clean1 + sl.m13_nobh_rest1
    m13_agnnocr = sl.m13_agnnocr_clean1 + sl.m13_agnnocr_rest1
    m13_agncr = sl.m13_agncr_clean1 + sl.m13_agncr_rest1
    m12_nobh = sl.m12_nobh_clean1 + sl.m12_nobh_rest1
    m12_agnnocr = sl.m12_agnnocr_clean1 + sl.m12_agnnocr_rest1
    m12_agncr = sl.m12_agncr_clean1 + sl.m12_agncr_rest1
    snapfiles = m13_nobh + m13_agnnocr + m13_agncr + \
                m12_nobh + m12_agnnocr + m12_agncr
    nobh = m12_nobh + m13_nobh
    agnnocr = m12_agnnocr + m13_agnnocr
    agncr = m12_agncr + m13_agncr
    ics = [filen.split('_')[0] for filen in snapfiles]
    snaplabels = [ic + ' noBH' if filen in nobh else
                  ic + ' AGN-noCR' if filen in agnnocr else
                  ic + ' AGN-CR' if filen in agncr else
                  None
                  for filen, ic in zip(snapfiles, ics)]
    firedataf = ol.filen_halocenrvir
    firemasses = []
    fireredshifts = []
    fireradii = []
    firecens = []
    with h5py.File(firedataf, 'r') as f:
        meandef = '200m'
        for sfn in snapfiles:
            _lsm = []
            _lsz = []
            _lsr = []
            _lfc = []
            smgrp = f[sfn]
            sngrpns = [key for key in smgrp.keys() if key.startswith('snap_')]
            for sngrpn in sngrpns:
                sngrp = smgrp[sngrpn]
                _lsz.append(sngrp['cosmopars'].attrs['z'])
                _lsm.append(sngrp[f'cen0/Rvir_{meandef}'].attrs['Mvir_g'])
                _lsr.append(sngrp[f'cen0/Rvir_{meandef}'].attrs['Rvir_cm'])
                _lfc.append(np.array([sngrp[f'cen0'].attrs[f'{_ax}c_cm']
                                      for _ax in ['X', 'Y', 'Z']]))
            zo = np.argsort(_lsz)
            _lsz = np.array(_lsz)[zo]
            _lsm = np.array(_lsm)[zo]
            _lsr = np.array(_lsr)[zo]
            _lfc = np.array(_lfc)[zo]
            _lsm /= c.solar_mass
            _lsr /= (c.cm_per_mpc * 1e-3)
            _lfc /= (c.cm_per_mpc * 1e-3)
            firemasses.append(_lsm)
            fireredshifts.append(_lsz)
            fireradii.append(_lsr)
            firecens.append(_lfc)

    fig = plt.figure(figsize=(5.5, 7.))
    grid = gsp.GridSpec(nrows=2, ncols=1,
                        height_ratios=[5., 2.], hspace=0.25)
    ax = fig.add_subplot(grid[0, 0])
    lax = fig.add_subplot(grid[1, 0])
    fontsize = 12
    # Burchett et al. data
    z_bur = data_bur['zgal']
    m_bur = data_bur['Mvir_Msun']
    isul = data_bur['log_N_Ne8_isUL']
    ax.scatter(z_bur[isul], m_bur[isul], marker='o', linestyle='None',
               edgecolor='black', facecolor='none', label='Bur.+19 (UL)',
               s=40, zorder=5)
    ax.scatter(z_bur[np.logical_not(isul)], m_bur[np.logical_not(isul)], 
               marker='o', linestyle='None', edgecolor='black', 
               facecolor='black', label='Bur.+19', s=40, 
               zorder=5)
    
    mass_minmax = {'m13': (np.inf, -np.inf),
                   'm12': (np.inf, -np.inf)}
    z_minmax = {'m13': (np.inf, -np.inf),
                'm12': (np.inf, -np.inf)}
    zmar = 0.05
    mmar = 0.2
    m13ics_used = set()
    m12ics_used = set()
    physmodels_used = set()
    for snapfile, snaplabel, firemass, fireredshift \
            in zip(snapfiles[::-1], snaplabels[::-1], firemasses[::-1],
                   fireredshifts[::-1]):
        hlabel, plabel = snaplabel.split(' ')
        if hlabel.startswith('m13'):
            color = sl.m13_iccolors[hlabel]
            m13ics_used.add(hlabel)
        elif hlabel.startswith('m12'):
            color = sl.m12_iccolors[hlabel]
            m12ics_used.add(hlabel)
        linestyle = sl.physlinestyles[plabel]
        physmodels_used.add(plabel)
        if snapfile in sl.buglist1:
            marker = 'x'
            ms = 5
        else:
            marker = 'o'
            ms = 3
        ax.plot(fireredshift, firemass, color=color,
                linestyle=linestyle, linewidth=1.5,
                marker=marker, markersize=ms)
        for key in mass_minmax:
            if snaplabel.startswith(key):
                _prev = mass_minmax[key]
                mmin = min(_prev[0], np.min(firemass))
                mmax = max(_prev[1], np.max(firemass))
                mass_minmax[key] = (mmin, mmax)
                _prev = z_minmax[key]
                zmin = min(_prev[0], np.min(fireredshift))
                zmax = max(_prev[1], np.max(fireredshift))
                z_minmax[key] = (zmin, zmax)
    mass_minmax = {key: (10**(np.log10(mass_minmax[key][0]) - mmar),
                         10**(np.log10(mass_minmax[key][1]) + mmar))
                   for key in mass_minmax}
    z_minmax = {key: (z_minmax[key][0] - zmar, z_minmax[key][1] + zmar)
                for key in z_minmax}
    print(mass_minmax)
    print(z_minmax)
    for key in mass_minmax:
        line1 = [mass_minmax[key][0], mass_minmax[key][1], 
                 mass_minmax[key][1], mass_minmax[key][0], 
                 mass_minmax[key][0]]
        line0 = [z_minmax[key][0], z_minmax[key][0], 
                 z_minmax[key][1], z_minmax[key][1], 
                 z_minmax[key][0]]
        ax.plot(line0, line1, linestyle='solid', color='gray',
                linewidth=2)
        ax.text(z_minmax[key][1] - 0.02, mass_minmax[key][1] * 0.95,
                key, fontsize=fontsize, color='gray',
                horizontalalignment='right', verticalalignment='top')

    ax.set_yscale('log')
    ax.set_ylabel(f'$\\mathrm{{M}}_{{\\mathrm{{{meandef}}}}} \\;'
                   ' [\\mathrm{M}_{\\odot}]$', fontsize=fontsize)
    ax.set_xlabel('redshift', fontsize=fontsize)
    ax.tick_params(labelsize=fontsize - 1., direction='in', which='both')
    
    m13ics_used = sorted(list(m13ics_used))
    m12ics_used = sorted(list(m12ics_used))
    physmodels_used = sorted(list(physmodels_used))
    _handles, _ = ax.get_legend_handles_labels()
    handles1 = [mlines.Line2D((), (), linewidth=1.5, 
                              linestyle=sl.physlinestyles[key],
                              label=key, color='black') \
                for key in physmodels_used]
    handles2 = [mlines.Line2D((), (), linewidth=1.5, 
                              color=sl.m12_iccolors[key],
                              label=key, linestyle='solid', marker='o',
                              markersize=4)\
                for key in m12ics_used]
    handles3 = [mlines.Line2D((), (), linewidth=1.5, 
                              color=sl.m13_iccolors[key],
                              label=key, linestyle='solid', marker='o',
                              markersize=4)\
                for key in m13ics_used]
    handles = _handles + handles1
    ax.legend(handles=handles, fontsize=fontsize - 1, handlelength=2.0)
    lhandles = handles2 + handles3
    lax.axis('off')
    lax.legend(handles=lhandles, loc='upper center',
               bbox_to_anchor=(0.5, 1.0), fontsize=fontsize, ncols=3)
    outdir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/datacomp/'
    outfilen = outdir + 'mass_z_selection_model3.pdf'
    plt.savefig(outfilen, bbox_inches='tight')


def datacomp_ne8_burchett19_clean(dset='m12'):
    '''
    dset: ['m12', 'm13']
        which halo/data selection to use
    '''
    fontsize = 12
    xlabel = '$\\mathrm{r}_{\\perp}$ [pkpc]'
    ylabel = ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\, VIII}) \\;'
              '[\\mathrm{cm}^{-2}]$')
    axtitles = ['noBH', 'AGN-noCR', 'AGN-CR']
    phystab = {'noBH': lambda x: '_sdp1e10_' in x,
               'AGN-CR': lambda x: '_MHDCRspec1_' in x,
               'AGN-noCR': lambda x: ('_sdp1e10_' not in x 
                                      and '_MHDCRspec1_' not in x),
              }
    ffilentemp = ('{dir}/coldens_Ne8_{simn}_snap{snap}'
                 '_shrink-sph-cen_BN98_2rvir_v2.hdf5')
    fdir_opts = ['/Users/nastasha/ciera/sim_maps/fire/clean_set1/',
                 '/Users/nastasha/ciera/sim_maps/fire/clean_set2/']
    oddir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/'
    ofilen = oddir + 'data_burchett_etal_2019_table1.txt'
    datacomprange_m200m_msun = {'m13': (2437844520477.7627, 
                                        13425998015000.441),
                                'm12': (403416630932.0638, 
                                        1709940889606.5674)}
    datacomprange_z = {'m13': (0.4488065752755633, 1.0500000000106098), 
                       'm12': (0.44880657526818074, 1.0500000000006244)}
    if dset == 'm12':
        rbins_pkpc = np.linspace(0., 300., 50)
        simnames = [('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr'
                     '_crdiffc690_sdp2e-4_gacc31_fa0.5'),
                    ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr'
                     '_crdiffc690_sdp1e10_gacc31_fa0.5'),
                   ]
        simnames_sup = []
        snaps_sr = [45, 46, 47, 48, 49, 50]
        snaps_hr = [186, 197, 210, 224, 240, 258]
        sims_sr = ['m12f_m6e4']
        sims_hr = ['m12f_m7e3']
    elif dset == 'm13':
        rbins_pkpc = np.linspace(0., 600., 50)
        simnames = [('m13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                    '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ('m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr'
                     '_crdiffc690_sdp1e-4_gacc31_fa0.5'),
                    ('m13h113_m3e5_MHD_fire3_fireBH_Sep182021'
                     '_crdiffc690_sdp1e10_gacc31_fa0.5'),
                    ('m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ('m13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr'
                     '_crdiffc690_sdp3e-4_gacc31_fa0.5'),
                    ('m13h206_m3e5_MHD_fire3_fireBH_Sep182021'
                     '_crdiffc690_sdp1e10_gacc31_fa0.5'),
                    ]
        simnames_sup = []
        snaps_sr = [45, 46, 47, 48, 49, 50]
        snaps_hr = [186, 197, 210, 224, 240, 258]
        sims_sr = ['m13h113_m3e5', 'm13h206_m3e5']
        sims_hr = ['m13h113_m3e4', 'm13h206_m3e4']
    lw_main = 2.
    lw_sup = 1.
    alpha_range = 0.3
    colors = tc.tol_cset('muted')
    linestyles = ['solid', 'dashed', 'dashdot', 'dotted']
    icusedlist = []
    rcens = 0.5 * (rbins_pkpc[1:] + rbins_pkpc[:-1])

    panelsize = 3.
    numcols = len(axtitles)
    numrows = 1
    legheight = 0.5
    wspace = 0.
    width_ratios = [panelsize] * numcols
    hspace = 0.25
    height_ratios = [panelsize] * numrows + [legheight]
    width = sum(width_ratios) * (1. + 1. / (len(width_ratios) - 1.) * wspace)
    height = sum(height_ratios) \
             * (1. + 1. / (len(height_ratios) - 1.) * hspace)
    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(nrows=numrows + 1, ncols=numcols, hspace=hspace, 
                        wspace=wspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[0, i]) for i in range(numcols)]
    lax = fig.add_subplot(grid[1, :])

    title = 'Burchett et al. (2019) data vs. FIRE-3 z=0.5-1.0'
    fig.suptitle(title, fontsize=fontsize)

    #print(modelfilelists)
    for simi, simn in enumerate(simnames + simnames_sup):
        if np.any([simn.startswith(st) for st in sims_sr]):
            snapnums = snaps_sr
        elif np.any([simn.startswith(st) for st in sims_hr]):
            snapnums = snaps_hr
        else:
            msg = (f'No snap list for simname {simn}; options:'
                   f'{sims_hr},\n{sims_sr}')
            raise RuntimeError(msg)
        simlabel = simn.split('_')[0]
        filens = [[ffilentemp.format(dir=fdir, simn=simn, snap=snap)
                   for fdir in fdir_opts 
                   if os.path.isfile(ffilentemp.format(dir=fdir, 
                                                       simn=simn, 
                                                       snap=snap))
                   ][0]
                  for snap in snapnums]
        rcens = 0.5 * (rbins_pkpc[:-1] + rbins_pkpc[1:])
        plo, pmed, phi = gpr.get_profile_massmap(filens, rbins_pkpc,
                                                 rbin_units='pkpc',
                                                 profiles=['perc-0.1', 
                                                           'perc-0.5', 
                                                           'perc-0.9'])
        if simlabel in icusedlist:
            _label = None
            ici = np.where([simlabel == _ic for _ic in icusedlist])[0][0]
        else:
            _label = simlabel
            icusedlist.append(simlabel)
            ici = len(icusedlist) - 1 
        for _plab in axtitles:
            if phystab[_plab](simn):
                plab = _plab
                break
        ismain = simn in simnames
        lw = lw_main if ismain else lw_sup
        ax = axes[np.where([plab == axt for axt in axtitles])[0][0]] 
        color = colors[ici % len(colors)]
        ls = linestyles[ici // len(colors)]
        ax.plot(rcens, pmed, color=color, linestyle=ls, linewidth=lw, 
                label=_label)
        if ismain:
            ax.fill_between(rcens, plo, phi, color=color, 
                            alpha=alpha_range, linestyle=ls,
                            linewidth=0.5)
            
    data_bur = pd.read_csv(ofilen, comment='#', sep='\t')
    cosmopars_bur = {'h': 0.677, 'omegam': 0.31, 'omegalambda': 0.69}
    def hmfunc(x):
        csm = cosmopars_bur.copy()
        csm.update({'z': x.zgal, 'a': 1. / (1. + x.zgal)})
        mv = cu.mvir_from_rvir(x.rvir_kpc * 1e-3 * c.cm_per_mpc, 
                               csm, meandef='200m')
        return mv / c.solar_mass
    data_bur = data_bur.assign(Mvir_Msun=lambda x: hmfunc(x))
    minmax_m200m_msun = datacomprange_m200m_msun[dset]
    minmax_z = datacomprange_z[dset]
    data_bur = data_bur[data_bur['Mvir_Msun'] >= minmax_m200m_msun[0]]
    data_bur = data_bur[data_bur['Mvir_Msun'] <= minmax_m200m_msun[1]]
    data_bur = data_bur[data_bur['zgal'] >= minmax_z[0]]
    data_bur = data_bur[data_bur['zgal'] <= minmax_z[1]]

    for axi, (ax, axtitle) in enumerate(zip(axes, axtitles)):
        ax.set_title(axtitle, fontsize=fontsize)
        ax.set_xlabel(xlabel, fontsize=fontsize)
        if axi == 0:
            ax.set_ylabel(ylabel, fontsize=fontsize)
            _label = 'Burchett+19'
            _ullabel = 'Burchett+19 (UL)'
        else:
            _label = None
            _ullabel = None
        ax.tick_params(labelsize=fontsize - 1., direction='in', which='both',
                       top=True, right=True, labelleft=axi == 0)
        isul = data_bur['log_N_Ne8_isUL'].copy()
        notul = np.logical_not(isul)
        ax.errorbar(data_bur['impact_parameter_kpc'][notul], 
                    data_bur['log_N_Ne8_pcm2'][notul],
                    yerr=data_bur['log_N_Ne8_pcm2_err'][notul], 
                    linestyle='None', elinewidth=1.5, marker='o', 
                    markersize=3, color='black', capsize=3,
                    label=_label)
        ax.scatter(data_bur['impact_parameter_kpc'][isul], 
                   data_bur['log_N_Ne8_pcm2'][isul],
                   linestyle='None', marker='v', 
                   s=10, facecolors='none', edgecolors='black', 
                   label=_ullabel)
    ylims = [ax.get_ylim() for ax in axes]
    ymin = min([yl[0] for yl in ylims])
    ymax = min([yl[1] for yl in ylims])
    [ax.set_ylim((ymin, ymax)) for ax in axes]

    hlist = []
    for ax in axes:
        _h, _ = ax.get_legend_handles_labels()
        hlist = hlist + _h
    handles1 = [mlines.Line2D((), (), linewidth=lw_main, linestyle='solid',
                            label='FIRE median', color='black'),
                mpatch.Patch(label='FIRE perc. 10-90', linewidth=0.5, 
                             color='black', alpha=alpha_range)
                ]
    lax.axis('off')
    lax.legend(handles=handles1 + hlist, fontsize=fontsize, ncols=numcols,
               loc='upper center')
    outdir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/datacomp/'
    outname = outdir + f'N_Ne8comp_cleanset1-2_{dset}.pdf'
    plt.savefig(outname, bbox_inches='tight')


def datacomp_ne8_burchett19(mset='m12', dset='model3'):
    '''
    mset: ['m12', 'm13']
        which halo mass selection to use
    dset: ['model3', 'clean']
        which set of simulations (ICs, phys, redshift) to use
    '''
    fontsize = 12
    xlabel = '$\\mathrm{r}_{\\perp}$ [pkpc]'
    ylabel = ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\, VIII}) \\;'
              '[\\mathrm{cm}^{-2}]$')
    axtitles = ['noBH', 'AGN-noCR', 'AGN-CR']
    phystab = {'noBH': lambda x: '_sdp1e10_' in x,
               'AGN-CR': lambda x: '_MHDCRspec1_' in x,
               'AGN-noCR': lambda x: ('_sdp1e10_' not in x 
                                      and '_MHDCRspec1_' not in x),
              }
    if dset == 'clean':
        # laptop
        ffilentemp = ('{dir}/coldens_Ne8_{simn}_snap{snap}'
                    '_shrink-sph-cen_BN98_2rvir_v2.hdf5')
        fdir_opts = ['/Users/nastasha/ciera/sim_maps/fire/clean_set1/',
                    '/Users/nastasha/ciera/sim_maps/fire/clean_set2/']
        oddir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/'
        ofilen = oddir + 'data_burchett_etal_2019_table1.txt'
        datacomprange_m200m_msun = {'m13': (2437844520477.7627, 
                                            13425998015000.441),
                                    'm12': (403416630932.0638, 
                                            1709940889606.5674)}
        datacomprange_z = {'m13': (0.4488065752755633, 1.0500000000106098),
                           'm12': (0.44880657526818074, 1.0500000000006244)}
        otherfills = [{}]
        if dset == 'm12':
            rbins_pkpc = np.linspace(0., 300., 50)
            simnames = [('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                        '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                        ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr'
                        '_crdiffc690_sdp2e-4_gacc31_fa0.5'),
                        ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr'
                        '_crdiffc690_sdp1e10_gacc31_fa0.5'),
                    ]
            simnames_sup = []
            snaps_sr = [45, 46, 47, 48, 49, 50]
            snaps_hr = [186, 197, 210, 224, 240, 258]
            sims_sr = ['m12f_m6e4']
            sims_hr = ['m12f_m7e3']
        elif dset == 'm13':
            rbins_pkpc = np.linspace(0., 600., 50)
            simnames = [('m13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                        '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                        ('m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr'
                        '_crdiffc690_sdp1e-4_gacc31_fa0.5'),
                        ('m13h113_m3e5_MHD_fire3_fireBH_Sep182021'
                        '_crdiffc690_sdp1e10_gacc31_fa0.5'),
                        ('m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                        '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                        ('m13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr'
                        '_crdiffc690_sdp3e-4_gacc31_fa0.5'),
                        ('m13h206_m3e5_MHD_fire3_fireBH_Sep182021'
                        '_crdiffc690_sdp1e10_gacc31_fa0.5'),
                        ]
            simnames_sup = []
            snaps_sr = [45, 46, 47, 48, 49, 50]
            snaps_hr = [186, 197, 210, 224, 240, 258]
            sims_sr = ['m13h113_m3e5', 'm13h206_m3e5']
            sims_hr = ['m13h113_m3e4', 'm13h206_m3e4']
    elif dset == 'model3':
        ffilentemp = ('{dir}/coldens_Ne8_{simn}_snap{snap}'
                    '_shrink-sph-cen_BN98_2rvir_{pax}-proj_v3.hdf5')
        fdir_opts = ['/projects/b1026/nastasha/maps/set3_model3/',
                    ]
        oddir = '/projects/b1026/nastasha/extdata/'
        ofilen = oddir + 'data_burchett_etal_2019_table1.txt'
        otherfills = [{'pax': 'x'}, {'pax': 'y'}, {'pax': 'z'}]
        datacomprange_m200m_msun = {'m13': (2437844520477.7627, 
                                            42578790794863.43), 
                                    'm12': (175827359425.83875, 
                                            1940067864780.1567),
                                    }
        datacomprange_z = {'m13': (0.4488065752755633, 1.0500000000106098),
                           'm12': (0.4488065752678155, 1.0500000000105174),
                           }
        if mset == 'm12':
            rbins_pkpc = np.linspace(0., 450., 50)
            snaps_sr = sl.snaplists['m12_sr']
            snaps_hr = sl.snaplists['m12_hr']
            simnames = sl.m12_sr_clean1 + sl.m12_hr_clean1
            simnames_sup = sl.m12_sr_rest1 + sl.m12_hr_rest1
            sims_sr = sl.m12_sr_all1
            sims_hr = sl.m12_hr_all1
        elif mset == 'm13':
            rbins_pkpc = np.linspace(0., 700., 50)
            snaps_sr = sl.snaplists['m13_sr']
            snaps_hr = sl.snaplists['m13_hr']
            simnames = sl.m13_sr_clean1 + sl.m13_hr_clean1
            simnames_sup = sl.m13_sr_rest1 + sl.m13_hr_rest1
            sims_sr = sl.m13_sr_all1
            sims_hr = sl.m13_hr_all1
    lw_main = 2.
    lw_sup = 1.
    alpha_range = 0.3
    colors = sl.m13_iccolors.copy()
    colors.update(sl.m12_iccolors)
    icusedlist = []
    rcens = 0.5 * (rbins_pkpc[1:] + rbins_pkpc[:-1])

    panelsize = 3.
    numcols = len(axtitles)
    numrows = 1
    legheight = 0.5
    wspace = 0.
    width_ratios = [panelsize] * numcols
    hspace = 0.25
    height_ratios = [panelsize] * numrows + [legheight]
    width = sum(width_ratios) * (1. + 1. / (len(width_ratios) - 1.) * wspace)
    height = sum(height_ratios) \
             * (1. + 1. / (len(height_ratios) - 1.) * hspace)
    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(nrows=numrows + 1, ncols=numcols, hspace=hspace, 
                        wspace=wspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[0, i]) for i in range(numcols)]
    lax = fig.add_subplot(grid[1, :])

    title = 'Burchett et al. (2019) data vs. FIRE-3 z=0.5-1.0'
    fig.suptitle(title, fontsize=fontsize)

    #print(modelfilelists)
    for simi, simn in enumerate(simnames + simnames_sup):
        if np.any([simn.startswith(st) for st in sims_sr]):
            snapnums = snaps_sr
        elif np.any([simn.startswith(st) for st in sims_hr]):
            snapnums = snaps_hr
        else:
            msg = (f'No snap list for simname {simn}; options:'
                   f'{sims_hr},\n{sims_sr}')
            raise RuntimeError(msg)
        simlabel = simn.split('_')[0]
        filens = [[ffilentemp.format(dir=fdir, simn=simn, snap=snap,
                                     **ofill)
                   for fdir in fdir_opts 
                   if os.path.isfile(ffilentemp.format(dir=fdir, 
                                                       simn=simn, 
                                                       snap=snap,
                                                       **ofill))
                   ][0]
                  for snap in snapnums for ofill in otherfills]
        plo, pmed, phi = gpr.get_profile_massmap(filens, rbins_pkpc,
                                                 rbin_units='pkpc',
                                                 profiles=['perc-0.1', 
                                                           'perc-0.5', 
                                                           'perc-0.9'])
        if simlabel in icusedlist:
            _label = None
            ici = np.where([simlabel == _ic for _ic in icusedlist])[0][0]
        else:
            _label = simlabel
            icusedlist.append(simlabel)
            ici = len(icusedlist) - 1 
        for _plab in axtitles:
            if phystab[_plab](simn):
                plab = _plab
                break
        ismain = simn in simnames
        lw = lw_main if ismain else lw_sup
        ax = axes[np.where([plab == axt for axt in axtitles])[0][0]] 
        color = colors[simlabel]
        ls = 'solid'
        ax.plot(rcens, pmed, color=color, linestyle=ls, linewidth=lw, 
                label=_label, path_effects=pu.getoutline(lw))
        if ismain:
            ax.fill_between(rcens, plo, phi, color=color, 
                            alpha=alpha_range, linestyle=ls,
                            linewidth=0.5)
            
    data_bur = pd.read_csv(ofilen, comment='#', sep='\t')
    cosmopars_bur = {'h': 0.677, 'omegam': 0.31, 'omegalambda': 0.69}
    def hmfunc(x):
        csm = cosmopars_bur.copy()
        csm.update({'z': x.zgal, 'a': 1. / (1. + x.zgal)})
        mv = cu.mvir_from_rvir(x.rvir_kpc * 1e-3 * c.cm_per_mpc, 
                               csm, meandef='200m')
        return mv / c.solar_mass
    data_bur = data_bur.assign(Mvir_Msun=lambda x: hmfunc(x))
    minmax_m200m_msun = datacomprange_m200m_msun[mset]
    minmax_z = datacomprange_z[mset]
    data_bur = data_bur[data_bur['Mvir_Msun'] >= minmax_m200m_msun[0]]
    data_bur = data_bur[data_bur['Mvir_Msun'] <= minmax_m200m_msun[1]]
    data_bur = data_bur[data_bur['zgal'] >= minmax_z[0]]
    data_bur = data_bur[data_bur['zgal'] <= minmax_z[1]]

    for axi, (ax, axtitle) in enumerate(zip(axes, axtitles)):
        ax.set_title(axtitle, fontsize=fontsize)
        ax.set_xlabel(xlabel, fontsize=fontsize)
        if axi == 0:
            ax.set_ylabel(ylabel, fontsize=fontsize)
            _label = 'Burchett+19'
            _ullabel = 'Burchett+19 (UL)'
        else:
            _label = None
            _ullabel = None
        ax.tick_params(labelsize=fontsize - 1., direction='in', which='both',
                       top=True, right=True, labelleft=axi == 0)
        isul = data_bur['log_N_Ne8_isUL'].copy()
        notul = np.logical_not(isul)
        ax.errorbar(data_bur['impact_parameter_kpc'][notul], 
                    data_bur['log_N_Ne8_pcm2'][notul],
                    yerr=data_bur['log_N_Ne8_pcm2_err'][notul], 
                    linestyle='None', elinewidth=1.5, marker='o', 
                    markersize=3, color='black', capsize=3,
                    label=_label, zorder=5)
        ax.scatter(data_bur['impact_parameter_kpc'][isul], 
                   data_bur['log_N_Ne8_pcm2'][isul],
                   linestyle='None', marker='v', 
                   s=10, facecolors='none', edgecolors='black', 
                   label=_ullabel, zorder=5)
    ylims = [ax.get_ylim() for ax in axes]
    ymin = min([yl[0] for yl in ylims])
    ymin = max(ymin, 11.1)
    ymax = min([yl[1] for yl in ylims])
    [ax.set_ylim((ymin, ymax)) for ax in axes]

    hlist = []
    for ax in axes:
        _h, _ = ax.get_legend_handles_labels()
        hlist = hlist + _h
    handles1 = [mlines.Line2D((), (), linewidth=lw_main, linestyle='solid',
                            label='FIRE median', color='black'),
                mpatch.Patch(label='FIRE perc. 10-90', linewidth=0.5, 
                             color='black', alpha=alpha_range)
                ]
    lax.axis('off')
    lax.legend(handles=handles1 + hlist, fontsize=fontsize, ncol=numcols,
               loc='upper center')
    outdir = '/projects/b1026/nastasha/imgs/datacomp/'
    outname = outdir + f'N_Ne8comp_cleanset1-2_{mset}_{dset}.pdf'
    plt.savefig(outname, bbox_inches='tight')