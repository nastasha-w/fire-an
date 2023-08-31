import matplotlib.gridspec as gsp
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c

ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
totfilen = ddir + 'gas_Neon_Ne8_masses_rTcuts.dat'
avfilen = ddir + 'mean_ZNe_by_mass_volume_rcuts.dat'

mdir = '/projects/b1026/nastasha/imgs/cgmprop/'

physmodels = {'m12': ['FIRE-2', 'noBH', 'AGN-noCR', 'AGN-CR'],
              'm13': ['noBH', 'AGN-noCR', 'AGN-CR']}
solarmassfrac_Ne = 10**-2.9008431 #copied from PS20 tables

def readfgas(rrange_rvir=(0.1, 1.0), trange_logk=(-np.inf, np.inf),
             massset='m12'):
    data = pd.read_csv(totfilen, sep='\t')
    filter = np.isclose(data['rmin_rvir'], rrange_rvir[0])
    filter &= np.isclose(data['rmax_rvir'], rrange_rvir[1])
    filter &= np.isclose(data['tmin_logk'], trange_logk[0])
    filter &= np.isclose(data['tmax_logk'], trange_logk[1])
    filter &= data['weight'] == 'gasmass'
    data = data[filter]
    data['ic'] = np.array([sl.ic_from_simname(simname)
                           for simname in data['simname']])
    filter2 = np.array([ic.startswith(massset) for ic in data['ic']])
    data = data[filter2]

    data['physmodel'] = np.array([sl.physlabel_from_simname(simname)
                                  for simname in data['simname']])
    data['fgas'] = data['total [g or num. part.]'] \
                   / (data['Mvir_g'] * data['Omega_b'] / data['Omega_m'])
    return data

def readZNe_mass(rrange_rvir=(0.1, 1.0), trange_logk=(-np.inf, np.inf),
                 massset='m12'):
    data = pd.read_csv(totfilen, sep='\t')
    filter = np.isclose(data['rmin_rvir'], rrange_rvir[0])
    filter &= np.isclose(data['rmax_rvir'], rrange_rvir[1])
    filter &= np.isclose(data['tmin_logk'], trange_logk[0])
    filter &= np.isclose(data['tmax_logk'], trange_logk[1])
    filter &= data['weight'].isin(['gasmass', 'Neon'])
    data = data[filter]
    data['ic'] = np.array([sl.ic_from_simname(simname)
                           for simname in data['simname']])
    filter2 = np.array([ic.startswith(massset) for ic in data['ic']])
    data = data[filter2]
    data['physmodel'] = np.array([sl.physlabel_from_simname(simname)
                                  for simname in data['simname']])
    dfgas = data[data['weight'] == 'gasmass']
    dfne = data[data['weight'] == 'Neon']

    dfmodgas = dfgas.pivot_table(columns=['physmodel', 'ic', 'redshift'],
                                 index=('simname', 'snapnum'),
                                 values='total [g or num. part.]')
    dfmodne = dfne.pivot_table(columns=['physmodel', 'ic', 'redshift'],
                               index=('simname', 'snapnum'),
                               values='total [g or num. part.]')
    dfmodnefrac = dfmodne * c.atomw_Ne * c.u / (dfmodgas * solarmassfrac_Ne) 
    dfnefrac = dfmodnefrac.melt(value_name='ZNe_solar')
    # NaN values picked up from missing snapnums (different sets), 
    # missing physmodel/IC combinations
    dfnefrac = dfnefrac[np.logical_not(np.isnan(dfnefrac['ZNe_solar']))]
    return dfnefrac

def readNe8frac(rrange_rvir=(0.1, 1.0), trange_logk=(-np.inf, np.inf),
                massset='m12'):
    data = pd.read_csv(totfilen, sep='\t')
    filter = np.isclose(data['rmin_rvir'], rrange_rvir[0])
    filter &= np.isclose(data['rmax_rvir'], rrange_rvir[1])
    filter &= np.isclose(data['tmin_logk'], trange_logk[0])
    filter &= np.isclose(data['tmax_logk'], trange_logk[1])
    filter &= data['weight'].isin(['Ne8', 'Neon'])
    data = data[filter]
    data['ic'] = np.array([sl.ic_from_simname(simname)
                           for simname in data['simname']])
    filter2 = np.array([ic.startswith(massset) for ic in data['ic']])
    data = data[filter2]
    data['physmodel'] = np.array([sl.physlabel_from_simname(simname)
                                  for simname in data['simname']])
    dfne8 = data[data['weight'] == 'Ne8']
    dfne = data[data['weight'] == 'Neon']

    dfmodne8 = dfne8.pivot_table(columns=['physmodel', 'ic', 'redshift'],
                                 index=('simname', 'snapnum'),
                                 values='total [g or num. part.]')
    dfmodne = dfne.pivot_table(columns=['physmodel', 'ic', 'redshift'],
                               index=('simname', 'snapnum'),
                               values='total [g or num. part.]')
    dfmodne8frac = dfmodne8 / dfmodne
    dfne8frac = dfmodne8frac.melt(value_name='Ne8frac')
    # NaN values picked up from missing snapnums (different sets), 
    # missing physmodel/IC combinations
    dfne8frac = dfne8frac[np.logical_not(np.isnan(dfne8frac['Ne8frac']))]
    return dfne8frac

def plot_fgashists(rrange_rvir=(0.1, 1.0),
                   trange_logk=(-np.inf, np.inf),
                   massset='m12'):
    data = readfgas(rrange_rvir=rrange_rvir,
                    trange_logk=trange_logk,
                    massset=massset)
    vmin = data['fgas'].min()
    vmax = data['fgas'].max()
    bins = np.linspace(0.99 * vmin, 1.01 * vmax, 16)

    fig = plt.figure(figsize=(5.5, 5.))
    ax = fig.add_subplot(1, 1, 1)
    fontsize = 12
    hatches = ['\\', '/', '|', '-']

    for physmodel, hatch in zip(physmodels[massset], hatches):
        color = sl.physcolors[physmodel]
        label = sl.plotlabel_from_physlabel[physmodel]
        f1 = data['physmodel'] == physmodel
        ax.hist(data.loc[f1, 'fgas'], label=label, bins=bins, color=color,
                density=False, histtype='step', linewidth=2,
                linestyle='solid', hatch=hatch)
        #f2 = np.logical_and(f1, data['isclean'])
        #ax.hist(data.loc[f2, 'fgas'], label=None, bins=bins, color=color,
        #        density=False, histtype='stepfilled', alpha=0.5,
        #        linewidth=2, linestyle='dashed')
    
    handles0, _ = ax.get_legend_handles_labels()
    #handles1 = [mpatch.Patch((), (), label='full sample', color='gray',
    #                         linewidth=2, linestyle='dashed'),
    #            mpatch.Patch((), (), label='clean sample', color='gray',
    #                         alpha=0.5),
    #            ]
    ax.legend(handles=handles0, fontsize=fontsize - 1)

    ax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                   top=True, right=True)
    xlabel = ('$\\mathrm{M}_{\mathrm{gas}} \\,/\\, (\\Omega_{\\mathrm{b}} '
              ' \\mathrm{M}_{\mathrm{vir}} \\,/ \\, \\Omega_{\\mathrm{m}})$')
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel('number of snapshots', fontsize=fontsize)
    if trange_logk == (-np.inf, np.inf):
        title = massset + ', gas at ' \
                + (f'${rrange_rvir[0]} \\endash {rrange_rvir[1]}'
                    '\\, \\mathrm{R}_{\\mathrm{vir}}$')
    else:
        title = massset + ', gas at ' \
                + (f'${rrange_rvir[0]} \\endash {rrange_rvir[1]}'
                   '\\, \\mathrm{R}_{\\mathrm{vir}}, '
                   f' \\mathrm{{T}} > 10^{{{trange_logk[0]:.1f}}}'
                   '\\mathrm{{K}}$')
    fig.suptitle(title, fontsize=fontsize)

    outname = mdir + (f'fgascomp_{massset}_{rrange_rvir[0]}_to'
                      f'{rrange_rvir[1]}_Rvir_Tgas_ge_{trange_logk[0]:.1f}')
    outname = outname.replace('.', 'p') + '.pdf'
    plt.savefig(outname, bbox_inches='tight')

def plot_ZNehists(rrange_rvir=(0.1, 1.0),
                  trange_logk=(-np.inf, np.inf),
                  massset='m12'):
    data = readZNe_mass(rrange_rvir=rrange_rvir,
                        trange_logk=trange_logk,
                        massset=massset)
    vmin = np.log10(data['ZNe_solar'].min())
    vmax = np.log10(data['ZNe_solar'].max())
    bins = np.linspace(0.99 * vmin, 1.01 * vmax, 16)

    fig = plt.figure(figsize=(5.5, 5.))
    ax = fig.add_subplot(1, 1, 1)
    fontsize = 12
    hatches = ['\\', '/', '|', '-']

    for physmodel, hatch in zip(physmodels[massset], hatches):
        color = sl.physcolors[physmodel]
        label = sl.plotlabel_from_physlabel[physmodel]
        f1 = data['physmodel'] == physmodel
        ax.hist(np.log10(data.loc[f1, 'ZNe_solar']), label=label, 
                bins=bins, color=color,
                density=False, histtype='step', linewidth=2,
                linestyle='solid', hatch=hatch)
        #f2 = np.logical_and(f1, data['isclean'])
        #ax.hist(data.loc[f2, 'fgas'], label=None, bins=bins, color=color,
        #        density=False, histtype='stepfilled', alpha=0.5,
        #        linewidth=2, linestyle='dashed')
    
    handles0, _ = ax.get_legend_handles_labels()
    #handles1 = [mpatch.Patch((), (), label='full sample', color='gray',
    #                         linewidth=2, linestyle='dashed'),
    #            mpatch.Patch((), (), label='clean sample', color='gray',
    #                         alpha=0.5),
    #            ]
    ax.legend(handles=handles0, fontsize=fontsize - 1)

    ax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                   top=True, right=True)
    xlabel = ('$\\log_{10} \\, \\mathrm{Z}_{\\mathrm{Ne}}'
              ' \\; [\\mathrm{Z}_{\\mathrm{Ne}, \\odot}]$')
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel('number of snapshots', fontsize=fontsize)
    if trange_logk == (-np.inf, np.inf):
        title = massset + ', gas and Ne at ' \
                + (f'${rrange_rvir[0]} \\endash {rrange_rvir[1]}'
                    '\\, \\mathrm{R}_{\\mathrm{vir}}$')
    else:
        title = massset + ', gas and Ne at ' \
                + (f'${rrange_rvir[0]} \\endash {rrange_rvir[1]}'
                   '\\, \\mathrm{R}_{\\mathrm{vir}}, '
                   f' \\mathrm{{T}} > 10^{{{trange_logk[0]:.1f}}}'
                   '\\mathrm{{K}}$')
    fig.suptitle(title, fontsize=fontsize)

    outname = mdir + (f'ZNe_masswtd_comp_{massset}_{rrange_rvir[0]}_to'
                      f'{rrange_rvir[1]}_Rvir_Tgas_ge_{trange_logk[0]:.1f}')
    outname = outname.replace('.', 'p') + '.pdf'
    plt.savefig(outname, bbox_inches='tight')

def plot_Ne8frachists(rrange_rvir=(0.1, 1.0),
                      trange_logk=(-np.inf, np.inf),
                      massset='m12'):
    data = readNe8frac(rrange_rvir=rrange_rvir,
                       trange_logk=trange_logk,
                       massset=massset)
    vmin = np.log10(data['Ne8frac'].min())
    vmax = np.log10(data['Ne8frac'].max())
    bins = np.linspace(0.99 * vmin, 1.01 * vmax, 16)

    fig = plt.figure(figsize=(5.5, 5.))
    ax = fig.add_subplot(1, 1, 1)
    fontsize = 12
    hatches = ['\\', '/', '|', '-']

    for physmodel, hatch in zip(physmodels[massset], hatches):
        color = sl.physcolors[physmodel]
        label = sl.plotlabel_from_physlabel[physmodel]
        f1 = data['physmodel'] == physmodel
        ax.hist(np.log10(data.loc[f1, 'Ne8frac']), label=label, 
                bins=bins, color=color,
                density=False, histtype='step', linewidth=2,
                linestyle='solid', hatch=hatch)
        #f2 = np.logical_and(f1, data['isclean'])
        #ax.hist(data.loc[f2, 'fgas'], label=None, bins=bins, color=color,
        #        density=False, histtype='stepfilled', alpha=0.5,
        #        linewidth=2, linestyle='dashed')
    
    handles0, _ = ax.get_legend_handles_labels()
    #handles1 = [mpatch.Patch((), (), label='full sample', color='gray',
    #                         linewidth=2, linestyle='dashed'),
    #            mpatch.Patch((), (), label='clean sample', color='gray',
    #                         alpha=0.5),
    #            ]
    ax.legend(handles=handles0, fontsize=fontsize - 1)

    ax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                   top=True, right=True)
    xlabel = ('$\\log_{10} \\, \\mathrm{Ne\\,VIII} \\,/ \\,'
              ' \\mathrm{Ne}$')
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel('number of snapshots', fontsize=fontsize)
    if trange_logk == (-np.inf, np.inf):
        title = massset + ', Ne at ' \
                + (f'${rrange_rvir[0]} \\endash {rrange_rvir[1]}'
                    '\\, \\mathrm{R}_{\\mathrm{vir}}$')
    else:
        title = massset + ', Ne at ' \
                + (f'${rrange_rvir[0]} \\endash {rrange_rvir[1]}'
                   '\\, \\mathrm{R}_{\\mathrm{vir}}, '
                   f' \\mathrm{{T}} > 10^{{{trange_logk[0]:.1f}}}'
                   '\\mathrm{{K}}$')
    fig.suptitle(title, fontsize=fontsize)

    outname = mdir + (f'Ne8frac_comp_{massset}_{rrange_rvir[0]}_to'
                      f'{rrange_rvir[1]}_Rvir_Tgas_ge_{trange_logk[0]:.1f}')
    outname = outname.replace('.', 'p') + '.pdf'
    plt.savefig(outname, bbox_inches='tight')

def plot_prophists(qty='fgas',
                   rrange_rvir=(0.1, 1.0),
                   trange_logk=(-np.inf, np.inf),
                   massset='m12'):
    if qty == 'fgas':
        data = readfgas(rrange_rvir=rrange_rvir,
                        trange_logk=trange_logk,
                        massset=massset)
        datakey = 'fgas'
        xlabel = ('$\\mathrm{M}_{\mathrm{gas}} \\,/\\, '
                  '(\\Omega_{\\mathrm{b}} '
                  ' \\mathrm{M}_{\mathrm{vir}} \\,/ \\,'
                  '\\Omega_{\\mathrm{m}})$')
        outnamestr = 'fgas'
        titlestr = 'gas'
    

    vmin = np.log10(data[datakey].min())
    vmax = np.log10(data[datakey].max())
    bins = np.linspace(0.99 * vmin, 1.01 * vmax, 16)

    fig = plt.figure(figsize=(5.5, 5.))
    ax = fig.add_subplot(1, 1, 1)
    fontsize = 12
    hatches = ['\\', '/', '|', '-']

    for physmodel, hatch in zip(physmodels[massset], hatches):
        color = sl.physcolors[physmodel]
        label = sl.plotlabel_from_physlabel[physmodel]
        f1 = data['physmodel'] == physmodel
        ax.hist(data.loc[f1, datakey], label=label, bins=bins, color=color,
                density=False, histtype='step', linewidth=2,
                linestyle='solid', hatch=hatch)
        #f2 = np.logical_and(f1, data['isclean'])
        #ax.hist(data.loc[f2, 'fgas'], label=None, bins=bins, color=color,
        #        density=False, histtype='stepfilled', alpha=0.5,
        #        linewidth=2, linestyle='dashed')
    
    handles0, _ = ax.get_legend_handles_labels()
    #handles1 = [mpatch.Patch((), (), label='full sample', color='gray',
    #                         linewidth=2, linestyle='dashed'),
    #            mpatch.Patch((), (), label='clean sample', color='gray',
    #                         alpha=0.5),
    #            ]
    ax.legend(handles=handles0, fontsize=fontsize - 1)

    ax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                   top=True, right=True)
    
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel('number of snapshots', fontsize=fontsize)
    if trange_logk == (-np.inf, np.inf):
        title = massset + f', {titlestr} at ' \
                + (f'${rrange_rvir[0]} \\endash {rrange_rvir[1]}'
                    '\\, \\mathrm{R}_{\\mathrm{vir}}$')
    else:
        title = massset + f', {titlestr} at ' \
                + (f'${rrange_rvir[0]} \\endash {rrange_rvir[1]}'
                   '\\, \\mathrm{R}_{\\mathrm{vir}}, '
                   f' \\mathrm{{T}} > 10^{{{trange_logk[0]:.1f}}}'
                   '\\mathrm{{K}}$')
    fig.suptitle(title, fontsize=fontsize)

    outname = mdir + (f'{outnamestr}comp_{massset}_{rrange_rvir[0]}_to'
                      f'{rrange_rvir[1]}_Rvir_Tgas_ge_{trange_logk[0]:.1f}')
    outname = outname.replace('.', 'p') + '.pdf'
    plt.savefig(outname, bbox_inches='tight')

def compmodels_fgas(rrange_rvir=(0.1, 1.0),
                    trange_logk=(-np.inf, np.inf),
                    massset='m12'):
    data = readfgas(rrange_rvir=rrange_rvir,
                    trange_logk=trange_logk,
                    massset=massset)
    vmin = data['fgas'].min()
    vmax = data['fgas'].max()
    bins = np.linspace(0.99 * vmin, 1.01 * vmax, 16)
    
    physmodels_this = physmodels[massset]
    ncomp = len(physmodels_this) - 1
    panelsize = 2.
    figsize = (ncomp * panelsize, ) * 2 
    fig = plt.figure(figsize=figsize)
    grid = gsp.GridSpec(ncols=ncomp, nrows=ncomp, 
                        wspace=0.0, hspace=0.0)
    fontsize = 12
    xlabel = ('$\\mathrm{M}_{\mathrm{gas}} \\,/\\, (\\Omega_{\\mathrm{b}} '
              ' \\mathrm{M}_{\mathrm{vir}} \\,/ \\, \\Omega_{\\mathrm{m}})$')
    ylabel = 'number of snapshots'
    
    axes = []
    for pi, physmodel in enumerate(physmodels_this):
        if pi == 0:
            continue
        compmodels = physmodels_this[:pi]
        for ci, compmodel in enumerate(compmodels):
            ax = fig.add_subplot(grid[pi - 1, ci])
            axes.append(ax)
            ax.tick_params(which='both', direction='in', 
                           labelsize=fontsize - 1.,
                           top=True, right=True,
                           labelbottom=(pi == ncomp),
                           labelleft=(ci == 0))
            if ci == ncomp // 2 and pi == ncomp:
                ax.set_xlabel(xlabel, fontsize=fontsize)
            if pi == ncomp // 2 + 1 and ci == 0:
                ax.set_ylabel(ylabel, fontsize=fontsize)

            c1 = sl.physcolors[physmodel]
            c2 = sl.physcolors[compmodel]
            f1 = data['physmodel'] == physmodel
            f2 = data['physmodel'] == compmodel
            ax.hist(data.loc[f1, 'fgas'], bins=bins, color=c1,
                    density=False, histtype='step', linewidth=2,
                    linestyle='solid')
            ax.hist(data.loc[f2, 'fgas'], bins=bins, color=c2,
                    density=False, histtype='step', linewidth=2,
                    linestyle='dashed')
            ics1 = set(data.loc[f1, 'ic'])
            ics2 = set(data.loc[f2, 'ic'])
            commonics = ics1 & ics2
            fc = data['ic'].isin(commonics)
            f1c = np.logical_and(f1, fc)
            f2c = np.logical_and(f2, fc)
            ax.hist(data.loc[f1c, 'fgas'], label=None, bins=bins, color=c1,
                    density=False, histtype='stepfilled', alpha=0.5)
            ax.hist(data.loc[f2c, 'fgas'], label=None, bins=bins, color=c2,
                    density=False, histtype='stepfilled', alpha=0.5)
    ylims = [ax.get_ylim() for ax in axes]
    ymin = 0. #min([ylim[0] for ylim in ylims])
    ymax = max([ylim[1] for ylim in ylims])
    [ax.set_ylim((ymin, ymax)) for ax in axes]
    
    handles0  = [mpatch.Patch(label=sl.plotlabel_from_physlabel[pmodel],
                              ec=mcolors.to_rgb(sl.physcolors[pmodel]),
                              linewidth=2, linestyle='solid',
                              fc=mcolors.to_rgb(sl.physcolors[pmodel]) \
                                 + (0.5,))
                 for pmodel in physmodels_this]
    handles1 = [mpatch.Patch(label='all ICs', edgecolor='gray',
                             linewidth=2, linestyle='solid',
                             facecolor='none'),
                mpatch.Patch(label='shared ICs', 
                             fc=mcolors.to_rgb('gray') + (0.5,),
                             edgecolor='none'),
                ]
    ax = fig.add_subplot(grid[0, 1])
    ax.axis('off')
    ax.legend(handles=handles0 + handles1, fontsize=fontsize - 1,
              loc='upper left', ncol=1)

    if trange_logk == (-np.inf, np.inf):
        title = massset + ', gas at ' \
                + (f'${rrange_rvir[0]} \\endash {rrange_rvir[1]}'
                    '\\, \\mathrm{R}_{\\mathrm{vir}}$')
    else:
        title = massset + ', gas at ' \
                + (f'${rrange_rvir[0]} \\endash {rrange_rvir[1]}'
                   '\\, \\mathrm{R}_{\\mathrm{vir}}, '
                   f' \\mathrm{{T}} > 10^{{{trange_logk[0]:.1f}}}'
                   '\\mathrm{{K}}$')
    fig.suptitle(title, fontsize=fontsize)

    outname = mdir + (f'fgascomp_physpairs_{massset}_{rrange_rvir[0]}_to'
                      f'{rrange_rvir[1]}_Rvir_Tgas_ge_{trange_logk[0]:.1f}')
    outname = outname.replace('.', 'p') + '.pdf'
    plt.savefig(outname, bbox_inches='tight')