import h5py
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
import fire_an.simlists as sl


ddir_clumpprof = '/projects/b1026/nastasha/clumps/'
outdir = '/projects/b1026/nastasha/imgs/clumpiness/'

def get_cprof_voltwtne8clumpfact(simname, snapnum):

    filen_temp = f'clumpines_measure_v1_{{mt}}_{simname}_{snapnum}'
    filen_temp = ddir_clumpprof + filen_temp + '.hdf5'

    with h5py.File(filen_temp.format(mt='Ne8num_Ne8dens')) as f:
        # num(Ne8) * n(Ne8) = V * n(ne8)**2
        sum_volne8sq = f['sums12/sums12'][:] 
        sum_volne8sq_tocgs =  f['sums12'].attrs['tocgs']
        rbins_rvir = f['rbins/rbins'][:]
    with h5py.File(filen_temp.format(mt='Vol_Ne8dens')) as f:
        # V * n(Ne8) = num(Ne8)
        sum_volne8 = f['sums12/sums12'][:] 
        sum_volne8_tocgs =  f['sums12'].attrs['tocgs']
        sum_vol =  f['sums1/sums1'][:]
        sum_vol_tocgs =  f['sums1'].attrs['tocgs']
    clumpiness = sum_volne8sq * sum_vol / sum_volne8**2
    clumpiness_tocgs = sum_volne8sq_tocgs * sum_vol_tocgs \
                       / sum_volne8_tocgs**2
    clumpiness = clumpiness * clumpiness_tocgs
    return rbins_rvir, clumpiness

def get_cprof_voltwtdensclumpfact(simname, snapnum):

    filen_temp = f'clumpines_measure_v1_{{mt}}_{simname}_{snapnum}'
    filen_temp = ddir_clumpprof + filen_temp + '.hdf5'

    with h5py.File(filen_temp.format(mt='Mass_dens')) as f:
        # num(Ne8) * n(Ne8) = V * n(ne8)**2
        sum_volrhosq = f['sums12/sums12'][:] 
        sum_volrhosq_tocgs =  f['sums12'].attrs['tocgs']
        rbins_rvir = f['rbins/rbins'][:]
    with h5py.File(filen_temp.format(mt='Vol_dens')) as f:
        # V * n(Ne8) = num(Ne8)
        sum_volrho = f['sums12/sums12'][:] 
        sum_volrho_tocgs =  f['sums12'].attrs['tocgs']
        sum_vol =  f['sums2/sums2'][:]
        sum_vol_tocgs =  f['sums2'].attrs['tocgs']
    clumpiness = sum_volrhosq * sum_vol / sum_volrho**2
    clumpiness_tocgs = sum_volrhosq_tocgs * sum_vol_tocgs \
                       / sum_volrho_tocgs**2
    clumpiness = clumpiness * clumpiness_tocgs
    return rbins_rvir, clumpiness


def plot_clumpiness_prof(target='Ne8'):
    '''
    target: {'Ne8', 'rho'}
        what to plot the clumpiness of: Ne VIII density or mass density
    '''
    showcolor_ics = ['m12f', 'm13h113', 'm13h206']
    icsnap_emph = [('m12f', sl.snaps_f2md[0]), ('m12f', sl.snaps_sr[0]),
                   ('m12f', sl.snaps_f2md[0]),
                   ('m13h113', sl.snaps_f2md[0]), ('m13h113', sl.snaps_sr[0]),
                   ]
    colors = sl.m12_iccolors.copy() 
    colors.update(sl.m13_iccolors)
    _colors = tc.tol_cset('muted')
    color_default = _colors.grey 

    axphys = [[('m12', 'FIRE-2'), ('m12', 'noBH'), ('m12', 'AGN-noCR'),
               ('m12', 'AGN-CR')],
              [None, ('m13', 'noBH'), ('m13', 'AGN-noCR'), 
               ('m13', 'AGN-CR')],
               ]
    _simnames = sl.m12_hr_all2 +  sl.m12_sr_all2 + sl.m12_f2md \
                + sl.m13_hr_all2 + sl.m13_sr_all2
    simnames = _simnames.copy()
    for sn in _simnames:
        if sn in sl.buglist1:
            simnames.remove(sn)
    sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    sims_f2md = sl.m12_f2md
    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    snaps_f2md = sl.snaps_f2md

    panelheight = 1.
    panelwidth = 1.3
    nrows = len(axphys)
    ncols = len(axphys[0])
    width_ratios = [panelwidth] * ncols
    height_ratios = [panelheight] * nrows
    wspace = 0.0
    hspace = 0.0
    width = sum(width_ratios)
    height = sum(height_ratios)
    fontsize = 12
    
    xlabel = '$r_{\\mathrm{3D}} \\; [\\mathrm{R}_{\\mathrm{vir}}]$'
    if target == 'Ne8':
        ylabel = ('$1 \\,/\\, f_{\\mathrm{c, V}}'
                  '(\mathrm{n}(\\mathrm{Ne\\,VIII}))$')
        cprof_func = get_cprof_voltwtne8clumpfact
    elif target == 'rho':
        ylabel = ('$1 \\,/\\, f_{\\mathrm{c, V}}'
                  '(\\rho)$')
        cprof_func = get_cprof_voltwtdensclumpfact

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows, wspace=wspace,
                        hspace=hspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    
    axes = []
    for ri, axphysrow in enumerate(axphys):
        for ci, axphys in enumerate(axphysrow):
            if axphys is None:
                continue
            ax = fig.add_subplot(grid[ri, ci])
            axes.append(ax)
            dobottom = (ri, ci) in [(0, 0), (1, 1), (1, 2), (1, 3)]
            doleft = (ri, ci) in [(0, 0), (1, 1)]
            ax.tick_params(which='both', direction='in',
                           labelsize=fontsize -1, top=True, right=True,
                           labelbottom=dobottom, labelleft=doleft)
            if dobottom:
                ax.set_xlabel(xlabel, fontsize=fontsize)
            if doleft:
                ax.set_ylabel(ylabel, fontsize=fontsize)
            axlabel = axphys[0] + ', ' \
                      + sl.plotlabel_from_physlabel(axphys[1])
            ax.text(0.05, 0.95, axlabel,
                    fontsize=fontsize, tranform=ax.transAxes,
                    horizontalalignment='left', verticalalignment='top')
            
            simnames_panel = [simname for simname in simnames 
                              if sl.ic_from_simname(simname)[:3] 
                                 == axphys[0]
                                 and sl.physlabel_from_simname(simname) 
                                     == axphys[1]]
            for sn in simnames_panel:
                snaps = (snaps_sr if sn in sims_sr
                         else snaps_hr if sn in sims_hr
                         else snaps_f2md if sn in sims_f2md
                         else None)
                ic = sl.ic_from_simname(sn)
                if ic in showcolor_ics:
                    color = colors[ic]
                    lw = 1.5
                    alpha = 1.
                    zo = 4
                else:
                    color = color_default
                    lw = 1.
                    alpha = 0.5
                    zo = 3
                for snap in snaps:
                    rbins, cf = cprof_func(sn, snap)
                    rcens = 0.5 * (rbins[:-1] + rbins[1:])
                    if (ic, snap) in icsnap_emph:
                        path_effects = pu.getoutline(1.5)
                        zo = 5
                    else:
                        path_effects = None
                    ax.plot(rcens, 1. / cf, linewidth=lw,
                            color=color, linestyle='solid',
                            alpha=alpha, path_effects=path_effects,
                            zorder = zo)
    [ax.set_ylim(0., 1.1) for ax in axes]
    [ax.set_xlim(0.1, 1.05) for ax in axes]

    lax = fig.add_subplot(grid[1, 0])
    lax.axis('off')
    handles = [mlines.Line2D((), (), linewidth=1.5, alpha=1., 
                             color=colors[ic], label=ic, linestyle='solid')
               for ic in showcolor_ics]
    handles = handles + [mlines.Line2D((), (), linewidth=1.0, alpha=0.5, 
                                       color=color_default, label='other ICs', 
                                       linestyle='solid')]
    lax.legend(handles=handles, loc='center', fontsize=fontsize)

    outname = outdir + f'volwtd_clumpfactor_prof_{target}.pdf'
    plt.savefig(outname, bbox_inches='tight')                    
              





    
    





