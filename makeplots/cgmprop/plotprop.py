import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fire_an.simlists as sl

ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
totfilen = ddir + 'gas_Neon_Ne8_masses_rTcuts.dat'
avfilen = ddir + 'mean_ZNe_by_mass_volume_rcuts.dat'

physmodels = {'m12': ['FIRE-2', 'noBH', 'AGN-noCR', 'AGN-CR'],
              'm13': ['noBH', 'AGN-noCR', 'AGN-CR']}

def plot_fgashists(rrange_rvir=(0.1, 1.0),
                   trange_logk=(-np.inf, np.inf),
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

    data['fgas'] = data['total [g or num. part.]'] \
                   / (data['Mvir_g'] * data['Omega_b'] / data['Omega_m'])
    data['physmodel'] = np.array([sl.physlabel_from_simname(simname)
                                  for simname in data['simname']])
    
    #data['isclean'] = data['ic'].isin(['m12f', 'm13h113', 'm13h206'])
    
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