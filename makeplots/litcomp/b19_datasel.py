'''
Burchett et al. (2018) data sample comparison selection
'''

import h5py
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fire_an.mainfunc.cengalprop as cgp
import fire_an.mainfunc.haloprop as hp
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu
import fire_an.utils.opts_locs as ol
import fire_an.makeplots.litcomp.b19_vs_analytical as bva

def readin_halodata(simnames, meandef='200m', zmin=0.45, zmax=1.05):
    firedataf = ol.filen_halocenrvir
    firemasses = []
    fireredshifts = []
    fireradii = []
    firecens = []
    with h5py.File(firedataf, 'r') as f:
        for sfn in simnames:
            _lsm = []
            _lsz = []
            _lsr = []
            _lfc = []
            smgrp = f[sfn]
            sngrpns = [key for key in smgrp.keys() if key.startswith('snap_')]
            for sngrpn in sngrpns:
                sngrp = smgrp[sngrpn]
                _zv = sngrp['cosmopars'].attrs['z']
                if _zv < zmin or _zv > zmax:
                    continue
                _lsz.append(_zv)
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
    return firemasses, fireredshifts, fireradii, firecens

def readin_cengaldata(simnames):
    datafn = ol.dir_halodata + 'pvcengal.hdf5'
    masses = []
    zs = []
    with h5py.File(datafn, 'r') as f:
        for sfn in simnames:
            _lsm = []
            _lsz = []
            smgrp = f[sfn]
            sngrpns = [key for key in smgrp.keys() if key.startswith('snap_')]
            for sngrpn in sngrpns:
                sngrp = smgrp[sngrpn]
                csmpath = 'pv0/doc/halodata_doc_dict/cosmopars_dict'
                _lsz.append(sngrp[csmpath].attrs['z'])
                _lsm.append(sngrp['pv0/doc'].attrs['mstar_gal_g'])
            zo = np.argsort(_lsz)
            _lsz = np.array(_lsz)[zo]
            _lsm = np.array(_lsm)[zo]
            _lsm /= c.solar_mass
            _lsm = np.log10(_lsm)
            masses.append(_lsm)
            zs.append(_lsz)
    return masses, zs

def plotMz_burchett_etal_2019(hset='clean', masscomp='halo_recalc'):
    '''
    hset: 'clean' or 'all'
    masscomp: 'halo', 'halo_recalc', or 'stellar'
    '''
    #datadir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/'
    datadir = '/projects/b1026/nastasha/extdata/'
    dfilen = datadir + 'data_burchett_etal_2019_table1.txt'
    if masscomp in ['halo', 'stellar']:
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
    elif masscomp in ['halo_recalc']:
        data_bur = bva.readdata_b19(nsigma=1)
    
    ## FIRE data
    m13_nobh = sl.m13_nobh_clean2 + sl.m13_nobh_rest2
    m13_agnnocr = sl.m13_agnnocr_clean2 + sl.m13_agnnocr_rest2
    m13_agncr = sl.m13_agncr_clean2 + sl.m13_agncr_rest2
    m12_nobh = sl.m12_nobh_clean2 + sl.m12_nobh_rest2
    m12_agnnocr = sl.m12_agnnocr_clean2 + sl.m12_agnnocr_rest2
    m12_agncr = sl.m12_agncr_clean2 + sl.m12_agncr_rest2
    m12_f2md = sl.m12_f2md
    _snapfiles = m13_nobh + m13_agnnocr + m13_agncr + \
                 m12_nobh + m12_agnnocr + m12_agncr + \
                 m12_f2md
    for sn in sl.buglist1:
        if sn in _snapfiles:
            _snapfiles.remove(sn)
    _ics = [sl.ic_from_simname(filen) for filen in _snapfiles]
    if hset == 'all':
        snapfiles = _snapfiles 
        ics = _ics
    elif hset == 'clean':
        _ics = np.array(_ics)
        _snapfiles = np.array(_snapfiles)
        icsel = np.array([sum([_ic == ic for _ic in _ics]) == 4 for 
                          ic in _ics])
        ics = _ics[icsel]
        snapfiles = _snapfiles[icsel]

    #nobh = m12_nobh + m13_nobh
    #agnnocr = m12_agnnocr + m13_agnnocr
    #agncr = m12_agncr + m13_agncr
    #f2md = m12_f2md
    snaplabels = [ic + ' ' + sl.physlabel_from_simname(filen)
                  for filen, ic in zip(snapfiles, ics)]

    if masscomp == 'halo':
        meandef = '200m'
        firemasses, fireredshifts, fireradii, firecens = \
            readin_halodata(snapfiles, meandef=meandef,
                            zmin=0.45, zmax=1.05)
    elif masscomp == 'stellar':
        firemasses, fireredshifts= \
            readin_cengaldata(snapfiles)
    elif masscomp == 'halo_recalc':
        meandef = 'BN98'
        firemasses, fireredshifts, fireradii, firecens = \
            readin_halodata(snapfiles, meandef=meandef,
                            zmin=0.45, zmax=1.05)
        firemasses = np.log10(firemasses)

    fig = plt.figure(figsize=(5.5, 7.))
    grid = gsp.GridSpec(nrows=2, ncols=1,
                        height_ratios=[5., 2.], hspace=0.25)
    ax = fig.add_subplot(grid[0, 0])
    lax = fig.add_subplot(grid[1, 0])
    fontsize = 12
    # Burchett et al. data
    if masscomp == 'halo':
        z_bur = data_bur['zgal']
        m_bur = data_bur['Mvir_Msun']
        isul = data_bur['log_N_Ne8_isUL']
        ax.scatter(z_bur[isul], m_bur[isul], marker='o', linestyle='None',
                edgecolor='black', facecolor='none', label='Bur.+19\n(UL)',
                s=40, zorder=5)
        ax.scatter(z_bur[np.logical_not(isul)], m_bur[np.logical_not(isul)], 
                marker='o', linestyle='None', edgecolor='black', 
                facecolor='black', label='Bur.+19', s=40, 
                zorder=5)
    elif masscomp == 'halo_recalc':
        z_bur = data_bur['zgal']
        m_bur = data_bur['logmvir_msun_bestest']
        m_bur_err = np.array([data_bur['logmvir_msun_bestest'] 
                                 - data_bur['logmvir_msun_lo'],
                              data_bur['logmvir_msun_hi']
                                 - data_bur['logmvir_msun_bestest']])
        isul = data_bur['log_N_Ne8_isUL']
        noul = np.logical_not(isul)
        ax.errorbar(z_bur[noul], m_bur[noul], yerr=m_bur_err[:, noul], 
                    linestyle='None', elinewidth=1.5, marker='o', 
                    markersize=7, color='black', capsize=3,
                    zorder=5, label='Bur.+19')
        ax.errorbar(z_bur[isul], m_bur[isul], yerr=m_bur_err[:, isul], 
                    color='black',
                    linestyle='None', elinewidth=1.5, marker='o', 
                    markersize=7, markeredgecolor='black', capsize=3,
                    markerfacecolor='none', zorder=5, label='Bur.+19\n(UL)')
    elif masscomp == 'stellar':
        isul = data_bur['log_N_Ne8_isUL']
        noul = np.logical_not(isul)
        z_bur = data_bur['zgal']
        m_bur = data_bur['log_Mstar_Msun']
        m_bur_err = data_bur['log_Mstar_Msun_err']
        ax.errorbar(z_bur[noul], m_bur[noul], yerr=m_bur_err[noul], 
                    linestyle='None', elinewidth=1.5, marker='o', 
                    markersize=4, color='black', capsize=3,
                    zorder=5, label='Bur.+19')
        ax.errorbar(z_bur[isul], m_bur[isul], yerr=m_bur_err[isul], 
                    color='black',
                    linestyle='None', elinewidth=1.5, marker='o', 
                    markersize=4, markeredgecolor='black', capsize=3,
                    markerfacecolor='none', zorder=5, label='Bur.+19\n(UL)')

    mass_minmax = {'m13': (np.inf, -np.inf),
                   'm12': (np.inf, -np.inf)}
    z_minmax = {'m13': (np.inf, -np.inf),
                'm12': (np.inf, -np.inf)}
    zmars = (0.05, 0.1)
    mmars = (0.2, 0.4)
    m13ics_used = set()
    m12ics_used = set()
    physmodels_used = set()
    for snapfile, snaplabel, firemass, fireredshift \
            in zip(snapfiles[::-1], snaplabels[::-1], firemasses[::-1],
                   fireredshifts[::-1]):
        if snapfile in sl.buglist1:
            continue
        hlabel, plabel = snaplabel.split(' ')
        if hlabel.startswith('m13'):
            color = sl.m13_iccolors[hlabel]
            m13ics_used.add(hlabel)
        elif hlabel.startswith('m12'):
            color = sl.m12_iccolors[hlabel]
            m12ics_used.add(hlabel)
        
        linestyle = sl.physlinestyles[plabel]
        physmodels_used.add(plabel)

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
    if masscomp in ['halo']:
        mass_minmax = [{key: (10**(np.log10(mass_minmax[key][0]) - mmar),
                              10**(np.log10(mass_minmax[key][1]) + mmar))
                        for key in mass_minmax}
                       for mmar in mmars]
    elif masscomp in ['halo_recalc']:
        mass_minmax = [{key: (mass_minmax[key][0] - mmar,
                              mass_minmax[key][1] + mmar)
                        for key in mass_minmax}
                       for mmar in mmars]
    elif masscomp == 'stellar':
        mass_minmax = [{key: (mass_minmax[key][0] - mmar,
                              mass_minmax[key][1] + mmar)
                        for key in mass_minmax}
                       for mmar in mmars]
    z_minmax = [{key: (z_minmax[key][0] - zmar, z_minmax[key][1] + zmar)
                 for key in z_minmax}
                for zmar in zmars]
    print(mass_minmax)
    print(z_minmax)
    for mi, ls in enumerate(['solid', 'dashed']):
        for key in mass_minmax[mi]:
            line1 = [mass_minmax[mi][key][0], mass_minmax[mi][key][1], 
                     mass_minmax[mi][key][1], mass_minmax[mi][key][0], 
                     mass_minmax[mi][key][0]]
            line0 = [z_minmax[mi][key][0], z_minmax[mi][key][0], 
                     z_minmax[mi][key][1], z_minmax[mi][key][1], 
                     z_minmax[mi][key][0]]
            ax.plot(line0, line1, linestyle=ls, color='gray',
                    linewidth=2, alpha=0.5)
    for key in ['m12', 'm13']:
        if key == 'm12':
            i0 = 0
            i1 = 0
        elif key == 'm13':
            i0 = 1
            i1 = 1
        if masscomp in ['halo']:
            top = mass_minmax[i0][key][i1] * 0.95
        elif masscomp in ['stellar', 'halo_recalc']:
            top = mass_minmax[i0][key][i1] - 0.05
        ax.text(z_minmax[0][key][1] - 0.02, top,
                key, fontsize=fontsize, color='gray',
                horizontalalignment='right', verticalalignment='top')

    if masscomp in ['halo']:
        ax.set_yscale('log')
        ax.set_ylabel(f'$\\mathrm{{M}}_{{\\mathrm{{{meandef}}}}} \\;'
                    ' [\\mathrm{M}_{\\odot}]$', fontsize=fontsize)
    elif masscomp in ['halo_recalc']:
        ax.set_ylabel(f'$\\mathrm{{M}}_{{\\mathrm{{vir}}}} \\;'
                    ' [\\mathrm{M}_{\\odot}]$', fontsize=fontsize)
    elif masscomp == 'stellar':
        ax.set_ylabel('$\\log_{10}\\, \\mathrm{M}_{*} \\;'
                      ' [\\mathrm{M}_{\\odot}]$',
                       fontsize=fontsize)
    ax.set_xlabel('redshift', fontsize=fontsize)
    ax.tick_params(labelsize=fontsize - 1., direction='in', which='both')
    xl = ax.get_xlim()
    ax.set_xlim((xl[0], 1.53))
    
    m13ics_used = sorted(list(m13ics_used))
    m12ics_used = sorted(list(m12ics_used))
    physmodels_used = sorted(list(physmodels_used))
    _handles, _ = ax.get_legend_handles_labels()
    _labels1 = {pm: (sl.plotlabel_from_physlabel[pm]).replace(' ', '\n')
                for pm in physmodels_used}
    handles1 = [mlines.Line2D((), (), linewidth=1.5, 
                              linestyle=sl.physlinestyles[key],
                              label=_labels1[key],
                              color='black') \
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
               bbox_to_anchor=(0.5, 1.0), fontsize=fontsize, ncol=3)
    #outdir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/datacomp/'
    outdir = '/projects/b1026/nastasha/imgs/datacomp/'
    outfilen = outdir + f'{masscomp}mass_z_selection_model3_{hset}.pdf'
    plt.savefig(outfilen, bbox_inches='tight')
    return mass_minmax, z_minmax

def plot_smhm_burchett19():
    #datadir = '/projects/b1026/nastasha/extdata/'
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
                               csm, meandef='200c')
        return mv / c.solar_mass
    data_bur = data_bur.assign(Mvir_Msun=lambda x: hmfunc(x))
    z_bur = data_bur['zgal']
    ms_bur = data_bur['log_Mstar_Msun']
    ms_bur_err = data_bur['log_Mstar_Msun_err']
    mh_bur = data_bur['Mvir_Msun']

    fig = plt.figure(figsize=(5.5, 5.))
    grid = gsp.GridSpec(ncols=2, nrows=1, wspace=0.1, width_ratios=[10., 1.])
    ax = fig.add_subplot(grid[0, 0])
    cax = fig.add_subplot(grid[0, 1])
    
    cmap = cm.get_cmap('jet')
    zmin = np.min(z_bur)
    zmax = np.max(z_bur)
    def cmapcall(z):
        return cmap((z - zmin)/ (zmax - zmin))
    caxnorm = mcolors.Normalize(vmin=zmin, vmax=zmax)
    for i in range(len(z_bur)):
        ax.errorbar(ms_bur[i], np.log10(mh_bur[i]), xerr=ms_bur_err[i], 
                    linestyle='None', elinewidth=1.5, marker='o', 
                    markersize=6, markerfacecolor=cmapcall(z_bur[i]), 
                    capsize=3, markeredgecolor='black',
                    ecolor='black')
    
    caxmap = cm.ScalarMappable(cmap=cmap, norm=caxnorm)
    plt.colorbar(caxmap, label='z', cax=cax)
    ax.set_xlabel('$\\log_{10}\\, \\mathrm{M}_{*}\\; [\\mathrm{M}_{\\odot}]$')
    ax.set_ylabel('$\\log_{10}\\, \\mathrm{M}_{\\mathrm{vir}} \\;'
                  ' [\\mathrm{M}_{\\odot}]$')
    #outdir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/datacomp/'
    outdir = '/projects/b1026/nastasha/imgs/datacomp/'
    outfilen = outdir + f'stellar_mass_halo_mass_burchett_etal_2019.pdf'
    plt.savefig(outfilen, bbox_inches='tight')
