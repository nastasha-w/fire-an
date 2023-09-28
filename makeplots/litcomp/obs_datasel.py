import h5py
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fire_an.makeplots.litcomp.b19_vs_analytical as bva
import fire_an.makeplots.litcomp.cubs7_qu_etal_dataread as cubsdr
import fire_an.makeplots.tol_colors as tc
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.opts_locs as ol

def readin_halodata(simnames, meandef='200m', zmin=0.45, zmax=1.05):
    firedataf = ol.filen_halocenrvir
    firemasses = []
    fireredshifts = []
    fireradii = []
    _simnames = []
    #firecens = []
    with h5py.File(firedataf, 'r') as f:
        for sfn in simnames:
            _lsm = []
            _lsz = []
            _lsr = []
            _lsn = []
            #_lfc = []
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
                #_lfc.append(np.array([sngrp[f'cen0'].attrs[f'{_ax}c_cm']
                #                      for _ax in ['X', 'Y', 'Z']]))
            zo = np.argsort(_lsz)
            _lsz = np.array(_lsz)[zo]
            _lsm = np.array(_lsm)[zo]
            _lsr = np.array(_lsr)[zo]
            #_lfc = np.array(_lfc)[zo]
            _lsm /= c.solar_mass
            _lsr /= (c.cm_per_mpc * 1e-3)
            #_lfc /= (c.cm_per_mpc * 1e-3)
            firemasses.append(_lsm)
            fireredshifts.append(_lsz)
            fireradii.append(_lsr)
            _simnames.append(sfn)
            #firecens.append(_lfc)
    return firemasses, fireredshifts, fireradii, _simnames

def plotMz_obs_fire(obsdata=('Q+23', 'B+19')):
    '''
    All halo masses calculated using the UM methods, always halo mass
    selection, boundaries from all non-bug ICs, but only plot one 
    phys. model for each IC for legibility.
    '''
    plotdata_obs = {}
    _colors = tc.tol_cset('high-contrast')
    if len(obsdata) == 2:
        obscolors = [_colors.blue, _colors.red]
    elif len(obsdata) == 1:
        obscolors = ['gray']
    if 'B+19' in obsdata:
        data_bur = bva.readdata_b19(nsigmas=1)
        z_bur = data_bur['zgal']
        m_bur = data_bur['logmvir_msun_bestest']
        m_bur_err = np.array([data_bur['logmvir_msun_bestest'] 
                                 - data_bur['logmvir_msun_lo'],
                              data_bur['logmvir_msun_hi']
                                 - data_bur['logmvir_msun_bestest']])
        isul_bur = data_bur['log_N_Ne8_isUL']
        noul_bur = np.logical_not(isul_bur)
        plotdata_obs['B+19'] = {'z': z_bur,
                               'mh': m_bur,
                               'mherr': m_bur_err,
                               'isul': isul_bur,
                               'noul': noul_bur}
    if 'Q+23' in obsdata:
        data_qu = cubsdr.getplotdata_cubs()
        z_qu = data_qu['z_gal']
        m_qu = data_qu['logmvir_msun_bestest']
        m_qu_err = np.array([data_qu['logmvir_msun_bestest'] 
                                 - data_qu['logmvir_msun_lo'],
                              data_qu['logmvir_msun_hi']
                                 - data_qu['logmvir_msun_bestest']])
        isul_qu = data_qu['isul_ne8']
        noul_qu = np.logical_not(isul_qu)
        # can't compare to missing data
        f1 = np.logical_not(np.isnan(m_qu))
        plotdata_obs['Q+23'] = {'z': z_qu[f1],
                               'mh': m_qu[f1],
                               'mherr': m_qu_err[:, f1],
                               'isul': isul_qu[f1],
                               'noul': noul_qu[f1]}

    ## FIRE data
    simnames = sl.m12_hr_all2 + sl.m12_sr_all2 + sl.m12_f2md \
               + sl.m13_hr_all2 + sl.m13_sr_all2 
    for sn in sl.buglist1:
        if sn in simnames:
            simnames.remove(sn)
    meandef = 'BN98'
    firemasses, fireredshifts, fireradii, firesimnames = \
        readin_halodata(simnames, meandef=meandef,
                        zmin=0.45, zmax=1.05)
    firemasses = np.log10(firemasses)

    fig = plt.figure(figsize=(5.5, 5.))
    ax = fig.add_subplot(1, 1, 1)
    fontsize = 12
   
    for obslabel, color in zip(obsdata, obscolors):     
        _data = plotdata_obs[obslabel]
        noul_label = obslabel + ' (det.)'
        ul_label = obslabel + ' (UL)'   
        noul = _data['noul']
        isul = _data['isul']
        ax.errorbar(_data['z'][noul], _data['mh'][noul],
                    yerr=_data['mherr'][:, noul], 
                    linestyle='None', elinewidth=1.5, marker='o', 
                    markersize=7, color=color, capsize=3,
                    zorder=5, label=noul_label, alpha=1.)
        ax.errorbar(_data['z'][isul], _data['mh'][isul], 
                    yerr=_data['mherr'][:, isul], 
                    color=color,
                    linestyle='None', elinewidth=1.5, marker='o', 
                    markersize=7, markeredgecolor=color, capsize=3,
                    markerfacecolor='none', zorder=5, label=ul_label,
                    alpha=0.42)

    mass_minmax = {'m13': (np.inf, -np.inf),
                   'm12': (np.inf, -np.inf)}
    z_minmax = {'m13': (np.inf, -np.inf),
                'm12': (np.inf, -np.inf)}
    zmars = (0.05, 0.1)
    mmars = (0.2, 0.4)
    ics_used = set()
    labeldone = False
    for simname, firemass, fireredshift \
            in zip(firesimnames, firemasses, fireredshifts):
        linestyle = 'solid'
        marker = 'o'
        ms = 3
        ic = sl.ic_from_simname(simname)
        if ic not in ics_used: # one curve per IC
            if not labeldone: 
                label = 'FIRE'
            else:
                label = None
            ax.plot(fireredshift, firemass, color='black',
                    linestyle=linestyle, linewidth=1.5,
                    marker=marker, markersize=ms,
                    label=label)
            labeldone = True
        for key in mass_minmax:
            if ic.startswith(key):
                print(key, ic, simname, firemass)
                _prev = mass_minmax[key]
                mmin = min(_prev[0], np.min(firemass))
                mmax = max(_prev[1], np.max(firemass))
                mass_minmax[key] = (mmin, mmax)
                _prev = z_minmax[key]
                zmin = min(_prev[0], np.min(fireredshift))
                zmax = max(_prev[1], np.max(fireredshift))
                z_minmax[key] = (zmin, zmax)
        ics_used.add(ic)

    mass_minmax = [{key: (mass_minmax[key][0] - mmar,
                            mass_minmax[key][1] + mmar)
                    for key in mass_minmax}
                    for mmar in mmars]
    z_minmax = [{key: (z_minmax[key][0] - zmar, z_minmax[key][1] + zmar)
                 for key in z_minmax}
                for zmar in zmars]
    print(mass_minmax)
    print(z_minmax)
    for mi, ls in enumerate(['solid']): # only plot one box
        for key in mass_minmax[mi]:
            line1 = [mass_minmax[mi][key][0], mass_minmax[mi][key][1], 
                     mass_minmax[mi][key][1], mass_minmax[mi][key][0], 
                     mass_minmax[mi][key][0]]
            line0 = [z_minmax[mi][key][0], z_minmax[mi][key][0], 
                     z_minmax[mi][key][1], z_minmax[mi][key][1], 
                     z_minmax[mi][key][0]]
            ax.plot(line0, line1, linestyle=ls, color='gray',
                    linewidth=2, alpha=0.5)

    ax.set_ylabel('$\\mathrm{M}_{\\mathrm{vir}} \\;'
                  ' [\\mathrm{M}_{\\odot}]$', fontsize=fontsize)
    ax.set_xlabel('redshift', fontsize=fontsize)
    ax.tick_params(labelsize=fontsize - 1., direction='in', which='both')
    xl = ax.get_xlim()
    if 'B+19' in obsdata:
        pass
        #ax.set_xlim((xl[0], 1.53))
    else:
        ax.set_xlim((xl[0], 0.87))
    
    for key in ['m12', 'm13']:
        if key == 'm12':
            i0 = 0
            i1 = 0
        elif key == 'm13':
            i0 = 1
            i1 = 1
        right = z_minmax[0][key][1] - 0.02
        xlim = ax.get_xlim()
        if right > xlim[1]:
            right = min(right, xlim[1] - 0.03)
            if key == 'm12':
                i1 = 1
        top = mass_minmax[i0][key][i1] - 0.05
        ax.text(right, top,
                key, fontsize=fontsize, color='gray',
                horizontalalignment='right', verticalalignment='top')
        
    #_handles, _ = ax.get_legend_handles_labels()
    #handles1 = [mlines.Line2D((), (), linewidth=1.5, 
    #                          linestyle='solid',
    #                          label='FIRE',
    #                          color='black')]
    #handles = _handles #+ handles1
    ax.legend(fontsize=fontsize - 1, handlelength=1.5)
    
    outdir = '/projects/b1026/nastasha/imgs/datacomp/'
    obss = '_'.join(obsdata)
    outfilen = outdir + ('recalc_halomass_z_selection_all2_simplified'
                         f'_{obss}.pdf')
    plt.savefig(outfilen, bbox_inches='tight')
    return mass_minmax, z_minmax