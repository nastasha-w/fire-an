

import numpy as np

import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt

import makeplots.get_2dprof as gpr


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