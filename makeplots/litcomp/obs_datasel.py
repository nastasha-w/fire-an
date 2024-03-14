import h5py
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fire_an.makeplots.litcomp.obsdataread as odr
import fire_an.makeplots.tol_colors as tc
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.opts_locs as ol

oddir = '/projects/b1026/nastasha/extdata/'
q23filen = oddir + 'plotdata_q23_nsigmas_1_2.dat'
b19filen = oddir + 'plotdata_b19_nsigmas_1_2.dat'

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
        data_bur = odr.readdata_b19(nsigmas=1)
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
        data_qu = odr.getplotdata_cubs()
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
    for sn in sl.buglist2:
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
            if not np.isfinite(mass_minmax[mi][key][0]):
                continue
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

def get_M_z_boxes_fire(sample='main'):
    '''
    All halo masses calculated using the UM methods, always halo mass
    selection, boundaries from all non-bug ICs
    '''
    ## FIRE data
    if sample == 'main':
        simnames = sl.m12_hr_all2 + sl.m12_sr_all2 + sl.m12_f2md \
                   + sl.m13_hr_all2 + sl.m13_sr_all2 
    elif sample == 'm12_f3nobh_comp':
        simnames = sl.m12_nobh_clean2 + sl.m12_nobh_rest2 \
                   + sl.m12plus_f3nobh + sl.m12plus_f3nobh_lores
    elif sample == 'fire2':
        simnames = sl.m12_f2md
    elif sample == 'f3xtest':
        simnames = sl.m12_f2md + sl.m12_nobh_clean2 + sl.m12_nobh_rest2 \
                   + sl.m12_fire3x_tests
    for sn in sl.buglist2:
        if sn in simnames:
            simnames.remove(sn)
    meandef = 'BN98'
    firemasses, fireredshifts, fireradii, firesimnames = \
        readin_halodata(simnames, meandef=meandef,
                        zmin=0.45, zmax=1.05)
    firemasses = np.log10(firemasses)

    mass_minmax = {'m13': (np.inf, -np.inf),
                   'm12': (np.inf, -np.inf)}
    z_minmax = {'m13': (np.inf, -np.inf),
                'm12': (np.inf, -np.inf)}
    zmars = (0.05, 0.1)
    mmars = (0.2, 0.4)

    for simname, firemass, fireredshift \
            in zip(firesimnames, firemasses, fireredshifts):
        ic = sl.ic_from_simname(simname)
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

    mass_minmax = [{key: (mass_minmax[key][0] - mmar,
                            mass_minmax[key][1] + mmar)
                    for key in mass_minmax}
                    for mmar in mmars]
    z_minmax = [{key: (z_minmax[key][0] - zmar, z_minmax[key][1] + zmar)
                 for key in z_minmax}
                for zmar in zmars]
  
    return mass_minmax, z_minmax

def addobs_panel(ax, plotdata, obslabel, color, vlabel=True):
    _data = plotdata
    space = '\n' if vlabel else ' '
    noul_mainlabel = obslabel + space + '(det.)'
    ul_label = obslabel + space + '(UL)'   
    noul = _data['noul']
    isul = _data['isul']
    flag = _data['flag']
    noul_main = np.logical_and(noul, np.logical_not(flag))
    noul_flag = np.logical_and(noul, flag)
    if sum(noul_flag) > 0:
        noul_flaglabel = obslabel + space + '(det., !)'
    else:
        noul_flaglabel = None
    ax.errorbar(_data['z'][noul_main], _data['mh'][noul_main],
                yerr=_data['mherr'][:, noul_main], 
                linestyle='None', elinewidth=1.5, marker='o', 
                markersize=7, color=color, capsize=3,
                zorder=5, label=noul_mainlabel, alpha=1.)
    ax.errorbar(_data['z'][noul_flag], _data['mh'][noul_flag], 
                yerr=_data['mherr'][:, noul_flag], 
                color=color,
                linestyle='None', elinewidth=1.5, marker='o', 
                markersize=7, markeredgecolor=color, capsize=3,
                markerfacecolor='none', zorder=5, label=noul_flaglabel,
                alpha=1.)
    ax.errorbar(_data['z'][isul], _data['mh'][isul], 
                yerr=_data['mherr'][:, isul], 
                color=color,
                linestyle='None', elinewidth=1.5, marker='o', 
                markersize=7, markeredgecolor=color, capsize=3,
                markerfacecolor='none', zorder=5, label=ul_label,
                alpha=0.42)

def addfire_panel(ax, firesimnames, firemasses, fireredshifts,
                  mass_minmax, z_minmax, xcrop=None, xcropb=None):
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
            if xcropb is None:
                ax.plot(fireredshift, firemass, color='black',
                        linestyle=linestyle, linewidth=1.5,
                        marker=marker, markersize=ms,
                        label=label)
            else:
                transi = np.where(fireredshift > xcrop)[0][0]
                ax.plot(fireredshift[transi - 1:], firemass[transi - 1:],
                        color='black',
                        linestyle='dotted', linewidth=1.5,
                        marker=marker, markersize=ms,
                        label=None, alpha=0.7)
                ax.plot(fireredshift[:transi], firemass[:transi],
                        color='black',
                        linestyle=linestyle, linewidth=1.5,
                        marker=marker, markersize=ms,
                        label=label) 
            labeldone = True
        ics_used.add(ic)

    for mi, ls in enumerate(['solid']): # only plot one box
        for key in mass_minmax[mi]:
            if not np.isfinite(mass_minmax[mi][key][0]):
                continue
            if xcrop is None:
                line1 = [mass_minmax[mi][key][0], mass_minmax[mi][key][1], 
                         mass_minmax[mi][key][1], mass_minmax[mi][key][0], 
                         mass_minmax[mi][key][0]]
                line0 = [z_minmax[mi][key][0], z_minmax[mi][key][0], 
                         z_minmax[mi][key][1], z_minmax[mi][key][1], 
                         z_minmax[mi][key][0]]
            else:
                line1 = [mass_minmax[mi][key][1], mass_minmax[mi][key][1], 
                         mass_minmax[mi][key][0], mass_minmax[mi][key][0]]
                line0 = [min(z_minmax[mi][key][1], xcropb), 
                         z_minmax[mi][key][0], z_minmax[mi][key][0],
                         min(z_minmax[mi][key][1], xcropb)]
                for yi in (0, 1):
                    ax.plot([xcropb, z_minmax[mi][key][1]], 
                            [mass_minmax[mi][key][yi]] * 2,
                            linestyle='dotted',
                            color='gray', alpha=0.5, linewidth=2.)
                
            ax.plot(line0, line1, linestyle=ls, color='gray',
                    linewidth=2, alpha=0.5)
            
def plotMz_obs_fire_2panel(ricut_pkpc=450., sample='main'):
    '''
    All halo masses calculated using the UM methods, always halo mass
    selection, boundaries from all non-bug ICs, but only plot one 
    phys. model for each IC for legibility.
    '''
    obsdata = ('B+19', 'Q+23')
    plabels = ('CASBaH', 'CUBS')
    plotdata_obs = {}
    _colors = tc.tol_cset('high-contrast')
    obscolors = [_colors.blue, _colors.red]
    if 'B+19' in obsdata:
        data_bur = pd.read_csv(b19filen, sep='\t')
        z_bur = data_bur['zgal']
        m_bur = data_bur['logmvir_msun_bestest']
        m_bur_err = np.array([data_bur['logmvir_msun_bestest'] 
                                 - data_bur['logmvir_msun_lo'],
                              data_bur['logmvir_msun_hi']
                                 - data_bur['logmvir_msun_bestest']])
        isul_bur = data_bur['log_N_Ne8_isUL']
        noul_bur = np.logical_not(isul_bur)
        ri_bur = data_bur['impact_parameter_kpc']
        f1 = ri_bur <= ricut_pkpc
        flagged = data_bur['flagged_by_qu23']
        plotdata_obs['B+19'] = {'z': z_bur[f1],
                                'mh': m_bur[f1],
                                'mherr': m_bur_err[:, f1],
                                'isul': isul_bur[f1],
                                'noul': noul_bur[f1],
                                'flag': flagged[f1]}
    if 'Q+23' in obsdata:
        data_qu = pd.read_csv(q23filen, sep='\t')
        z_qu = data_qu['z_gal']
        m_qu = data_qu['logmvir_msun_bestest']
        m_qu_err = np.array([data_qu['logmvir_msun_bestest'] 
                                 - data_qu['logmvir_msun_lo'],
                              data_qu['logmvir_msun_hi']
                                 - data_qu['logmvir_msun_bestest']])
        isul_qu = data_qu['isul_ne8']
        noul_qu = np.logical_not(isul_qu)
        ri_qu = data_qu['impactpar_kpc']
        flagged = np.zeros(ri_qu.shape, dtype=bool)
        # can't compare to missing data
        f1 = np.logical_not(np.isnan(m_qu))
        f1 &= ri_qu <= ricut_pkpc
        plotdata_obs['Q+23'] = {'z': z_qu[f1],
                                'mh': m_qu[f1],
                                'mherr': m_qu_err[:, f1],
                                'isul': isul_qu[f1],
                                'noul': noul_qu[f1],
                                'flag': flagged[f1]}

    ## FIRE data
    if sample == 'main':
        simnames = sl.m12_hr_all2 + sl.m12_sr_all2 + sl.m12_f2md \
                   + sl.m13_hr_all2 + sl.m13_sr_all2     
    elif sample == 'm12_f3nobh_comp':
        simnames = sl.m12_nobh_clean2 + sl.m12_nobh_rest2 \
                   + sl.m12plus_f3nobh + sl.m12plus_f3nobh_lores 
    elif sample == 'fire2':
        simnames = sl.m12_f2md
    elif sample == 'f3xtest':
        simnames = sl.m12_f2md + sl.m12_nobh_clean2 + sl.m12_nobh_rest2 \
                   + sl.m12_fire3x_tests
    for sn in sl.buglist2:
        if sn in simnames:
            simnames.remove(sn)
    meandef = 'BN98'
    firemasses, fireredshifts, fireradii, firesimnames = \
        readin_halodata(simnames, meandef=meandef,
                        zmin=0.45, zmax=1.05)
    firemasses = np.log10(firemasses)

    fig = plt.figure(figsize=(5.5, 4.5))
    grid = gsp.GridSpec(nrows=2, ncols=2, hspace=0.3, wspace=0.0,
                        height_ratios=(3.5, 1.))
    axes = [fig.add_subplot(grid[0, i]) for i in range(len(obsdata))]
    lax = fig.add_subplot(grid[1, :])
    lax.axis('off')
    fontsize = 12
   
    for ax, obslabel, color, plabel in zip(axes, obsdata, obscolors, plabels):     
        addobs_panel(ax, plotdata_obs[obslabel], obslabel, color,
                     vlabel=False) #obslabel=='B+19'
        ax.set_title(plabel, fontsize=fontsize)
    
    # don't really need to replot, but keeps consistency with the
    # data plot selections
    mass_minmax, z_minmax = get_M_z_boxes_fire(sample=sample)
    addfire_panel(axes[0], firesimnames, firemasses, fireredshifts,
                  mass_minmax, z_minmax, xcrop=None, xcropb=None)
    addfire_panel(axes[1], firesimnames, firemasses, fireredshifts,
                  mass_minmax, z_minmax, xcrop=0.805, xcropb=0.78)

    axes[0].set_ylabel('$\\mathrm{M}_{\\mathrm{vir}} \\;'
                       ' [\\mathrm{M}_{\\odot}]$', fontsize=fontsize)
    axes[0].set_xlabel('redshift', fontsize=fontsize)
    axes[1].set_xlabel('redshift', fontsize=fontsize)
    axes[0].tick_params(labelsize=fontsize - 1., direction='in', 
                        which='both', top=True, right=True)
    axes[1].tick_params(labelsize=fontsize - 1., direction='in', 
                        which='both', top=True, right=True,
                        labelleft=False)
    
    axes[1].set_xlim((0.38, 0.85))
    ylims = [ax.get_ylim() for ax in axes]
    ymin = min([yl[0] for yl in ylims])
    ymax = max([yl[1] for yl in ylims])
    [ax.set_ylim((ymin, ymax)) for ax in axes]
    
    if np.isfinite(mass_minmax[0]['m13'][1]):
        axes[0].text(0.95, mass_minmax[0]['m13'][1],
                    'm13', fontsize=fontsize, color='gray',
                    horizontalalignment='center',
                    verticalalignment='bottom')
    axes[0].text(0.95, mass_minmax[0]['m12'][0] - 0.07,
                 'm12', fontsize=fontsize, color='gray',
                 horizontalalignment='center',
                 verticalalignment='top')
    if np.isfinite(mass_minmax[0]['m13'][0]):
        axes[1].text(0.78, mass_minmax[0]['m13'][1],
                    'm13', fontsize=fontsize, color='gray',
                    horizontalalignment='center',
                    verticalalignment='bottom')
    axes[1].text(0.78, mass_minmax[0]['m12'][0] - 0.07,
                 'm12', fontsize=fontsize, color='gray',
                 horizontalalignment='center',
                 verticalalignment='top')
        
    #l0 = axes[0].legend(fontsize=fontsize - 2, handlelength=1.0, ncol=1,
    #                    loc='upper right', handletextpad=0.4)
    handles0, _ = axes[0].get_legend_handles_labels()
    handles1, _ = axes[1].get_legend_handles_labels()
    #axes[1].legend(handles=handles[1:],
    #               fontsize=fontsize - 2, handlelength=0.7,
    #               ncol=1, loc='lower right', handletextpad=0.25,
    #               columnspacing=0.9, borderaxespad=0.18,
    #               labelspacing=0.28)
    #axes[0].add_artist(l0)
    handles = handles0[1:] + handles1[1:] + handles0[:1]
    lax.legend(handles=handles, 
               fontsize=fontsize, ncol=3, loc='upper center',
               handlelength=1.5, columnspacing=1.5)
    outdir = '/projects/b1026/nastasha/imgs/datacomp/'
    obss = '_'.join(obsdata)
    outfilen = outdir + ('recalc_halomass_z_selection_all2_simplified'
                         f'_{obss}_2panel_{sample}.pdf')
    plt.savefig(outfilen, bbox_inches='tight')
    return mass_minmax, z_minmax