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

def addcoldata(ax, data_bur, labelmeas=None, labelul=None):
    isul = data_bur['log_N_Ne8_isUL'].copy()
    notul = np.logical_not(isul)
    ax.errorbar(data_bur['impact_parameter_kpc'][notul], 
                data_bur['log_N_Ne8_pcm2'][notul],
                yerr=data_bur['log_N_Ne8_pcm2_err'][notul], 
                linestyle='None', elinewidth=1.5, marker='o', 
                markersize=3, color='black', capsize=3,
                label=labelmeas, zorder=5)
    ax.scatter(data_bur['impact_parameter_kpc'][isul], 
                data_bur['log_N_Ne8_pcm2'][isul],
                linestyle='None', marker='v', 
                s=10, facecolors='none', edgecolors='black', 
                label=labelul, zorder=5)

def addveldata(ax, data_bur, label=None, absvals=True):
    isul = data_bur['log_N_Ne8_isUL'].copy()
    notul = np.logical_not(isul)
    ydata = (data_bur['v_kmps'][notul]).astype(np.float)
    print(ydata)
    if absvals: 
        ydata = np.abs(ydata)
    ax.errorbar(data_bur['impact_parameter_kpc'][notul], 
                ydata,
                yerr=data_bur['v_kmps_err'][notul], 
                linestyle='None', elinewidth=1.5, marker='o', 
                markersize=3, color='black', capsize=3,
                label=label, zorder=5)

def cdprof_ne8_burchett19(filen_temp, simnames, rbins_pkpc,
                          showscatter='clean',
                          datafieldsels=None, outname=None,
                          plottype='coldens'):
    '''
    plottype: 'coldens' or 'abslosvel'
    '''
    oddir = '/projects/b1026/nastasha/extdata/'
    ofilen = oddir + 'data_burchett_etal_2019_table1.txt'

    fontsize = 12
    xlabel = '$\\mathrm{r}_{\\perp}$ [pkpc]'
    if plottype == 'coldens':
        ylabel = ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\, VIII}) \\;'
                  '[\\mathrm{cm}^{-2}]$')
    elif plottype == 'abslosvel':
        ylabel = ('$|v_{\\mathrm{los}}| \\; '
                  '[\\mathrm{km}\\, \\mathrm{s}^{-1}]$')
    axtitles = ['noBH', 'AGN-noCR', 'AGN-CR']
    phystab = {'noBH': lambda x: '_sdp1e10_' in x,
               'AGN-CR': lambda x: '_MHDCRspec1_' in x,
               'AGN-noCR': lambda x: ('_sdp1e10_' not in x 
                                      and '_MHDCRspec1_' not in x),
              }
    otherfills = [{'pax': 'x'}, {'pax': 'y'}, {'pax': 'z'}]
    sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2

    if showscatter == 'clean':
        ics = [simname.split('_')[0] for simname in simnames]
        icsset = np.array(list(set(ics)))
        cleansel = np.array([sum([_ic == ic for _ic in ics]) == 3
                            for ic in icsset])
        scatterics = icsset[cleansel]
    elif showscatter == 'all':
        scatterics = [simname.split('_')[0] for simname in simnames]
    else:
        scatterics = showscatter

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

    #title = 'Burchett et al. (2019) data vs. FIRE-3 z=0.5-1.0'
    #fig.suptitle(title, fontsize=fontsize)

    #print(modelfilelists)
    for simn in simnames:
        if simn in sims_sr:
            snapnums = sl.snaps_sr
        elif simn in sims_hr:
            snapnums = sl.snaps_hr
        else:
            msg = (f'No snap list for simname {simn}; options:'
                   f'{sims_hr},\n{sims_sr}')
            raise RuntimeError(msg)
        simlabel = simn.split('_')[0]
        filens = [filen_temp.format(simname=simn, snapnum=snap, **ofill)
                  for snap in snapnums for ofill in otherfills]
        if plottype == 'coldens':
            weightmap = True
            absvals = False
            unitconv = 1.
        elif plottype == 'abslosvel':
            weightmap = False
            absvals = True
            unitconv = 1e-5
        plo, pmed, phi = gpr.get_profile_massmap(filens, rbins_pkpc,
                                                 rbin_units='pkpc',
                                                 profiles=['perc-0.1', 
                                                           'perc-0.5', 
                                                           'perc-0.9'],
                                                 weightmap=weightmap,
                                                 absvals=absvals)
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
        ismain = simlabel in scatterics
        lw = lw_main if ismain else lw_sup
        ax = axes[np.where([plab == axt for axt in axtitles])[0][0]] 
        color = colors[simlabel]
        ls = 'solid'
        ax.plot(rcens, pmed * unitconv, color=color, linestyle=ls, linewidth=lw, 
                label=_label, path_effects=pu.getoutline(lw))
        if ismain:
            ax.fill_between(rcens, plo * unitconv, phi * unitconv, color=color, 
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
    if datafieldsels is not None:
        for seltuple in datafieldsels:
            field = seltuple[0]
            minv = seltuple[1]
            maxv = seltuple[2]
            data_bur = data_bur[data_bur[field] >= minv]
            data_bur = data_bur[data_bur[field] <= maxv]

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
        if plottype == 'coldens':
            addcoldata(ax, data_bur, labelmeas=_label, labelul=_ullabel)
        elif plottype == 'abslosvel':
            addveldata(ax, data_bur, label=_label, absvals=True)
    ylims = [ax.get_ylim() for ax in axes]
    if plottype == 'coldens':
        ymin = min([yl[0] for yl in ylims])
        ymin = max(ymin, 11.1)
    elif plottype == 'abslosvel':
        ymin = 0.
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
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotsets_ne8_burchett19(hsel='all'):
    mdir = '/projects/b1026/nastasha/maps/clean2_vlos/'
    filen_temp = ('vlos_by_coldens_Ne8_{simname}_snap{snapnum}'
                  '_shrink-sph-cen_BN98_depth_2.0rvir_{pax}-proj_v3.hdf5')
    datacomprange_m200m_msun = {'m13': (2437844520477.7627, 
                                        13425998015000.441),
                                'm12': (403416630932.0638, 
                                        1709940889606.5674)}
    datacomprange_z = {'m13': (0.4488065752755633, 1.0500000000106098),
                       'm12': (0.44880657526818074, 1.0500000000006244)}
    rbins_pkpc_m12 = np.linspace(0., 450., 50)
    rbins_pkpc_m13 = np.linspace(0., 600., 50)
    rbins = {'m12': rbins_pkpc_m12,
             'm13': rbins_pkpc_m13}
    simnames_all = {'m12': sl.m12_sr_all2 + sl.m12_hr_all2,
                    'm13': sl.m13_sr_all2 + sl.m13_hr_all2}
    simnames = dict()
    for key in simnames_all:
        _sns = simnames_all[key].copy()
        for sn in sl.buglist1:
            if sn in _sns:
                _sns.remove(sn)
        if hsel == 'all':
            simnames[key] = _sns
        elif hsel == 'clean':
            ics = [sn.split('_')[0] for sn in _sns]
            uniqueics = np.unique(ics)
            icsel = np.array([sum([_ic == ic for _ic in ics]) == 3 
                              for ic in uniqueics])
            ics_incl = uniqueics[icsel]
            simnames[key] = [sn for sn in _sns 
                             if sn.split('_')[0] in ics_incl]
    
    outdir = '/projects/b1026/nastasha/imgs/datacomp/'
    for mset in ['m12', 'm13']:
        datafieldsels = [('Mvir_Msun',) + datacomprange_m200m_msun[mset],
                         ('zgal',) + datacomprange_z[mset]]
        for plottype in ['coldens', 'abslosvel']:
            outname = outdir + f'{plottype}_Ne8comp_{mset}_{hsel}2.pdf'
    
            cdprof_ne8_burchett19(mdir + filen_temp, simnames[mset], 
                                  rbins[mset],
                                  showscatter='clean',
                                  datafieldsels=datafieldsels,
                                  outname=outname,
                                  plottype=plottype)
