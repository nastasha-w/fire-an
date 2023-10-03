import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import fire_an.makeplots.get_2dprof as gpr
import fire_an.simlists as sl


def comp_crheatfix_crbug():
    outdir = '/projects/b1026/nastasha/imgs/datacomp/'
    outname = outdir + 'cr_heating_bug_fix_Ne8_coldensprof_comp.pdf'

    simnames = {'sn_bug_m12i': 'm12i_r7100',
                'sn_bug_m12f': 'm12f_r7100',
                'sn_fix_m12i': 'crheatfix_m12i_r7100',
                'sn_fix_m12f': 'crheatfix_m12f_r7100',
                }
    snapnums = sl.snaps_f2md
    mdir = '/projects/b1026/nastasha/maps/vdopmaps_all2/'
    filen_temp = ('vdoplos_by_coldens_Ne8_{simname}_snap{snapnum}'
                  '_shrink-sph-cen_BN98_depth_2.0rvir_{pax}-proj_v3.hdf5')
    
    rbins_pkpc = np.linspace(0., 450., 50)
    rcens = 0.5 * (rbins_pkpc[:-1] + rbins_pkpc[1:])
    otherfills = [{'pax': 'x'}, {'pax': 'y'}, {'pax': 'z'}]

    filens = {simkey: [mdir + filen_temp.format(simname=simnames[simkey], 
                                                snapnum=snap, **ofill)
                       for snap in snapnums for ofill in otherfills]
              for simkey in simnames}
    filens_z = {simkey: [[mdir + filen_temp.format(simname=simnames[simkey], 
                                                  snapnum=snap, **ofill)
                          for ofill in otherfills]
                         for snap in snapnums]
                for simkey in simnames}
    fig = plt.figure(figsize=(5.5, 3.))
    grid = gsp.GridSpec(ncols=2, nrows=1, hspace=0.0, wspace=0.0)
    axes = [fig.add_subplot(grid[0, i]) for i in range(2)]
    fontsize = 12
    colors = {'bug': 'gray',
              'fix': 'black'}

    for axi, ic in enumerate(['m12i', 'm12f']):
        ax = axes[axi]
        doleft = axi == 0
        ax.tick_params(which='both', direction='in', labelsize=fontsize - 1,
                       top=True, right=True, labelleft=doleft)
        if doleft:
            ax.set_ylabel('$\\log_{10}\\,\\mathrm{N}(\\mathrm{Ne\\,VIII})'
                          '\\; [\\mathrm{cm}^{-2}]$', fontsize=fontsize)
        ax.set_xlabel('$\\mathrm{r}_{\\perp} \\; [\\mathrm{pkpc}]$')
        ax.text(0.95, 0.95, ic, fontsize=fontsize,
                transform=ax.transAxes, horizontalalignment='right',
                verticalalignment='top')

        for opti, opt in enumerate(['bug', 'fix']):
            snkey = f'sn_{opt}_{ic}'
            _filens = filens[snkey]
            _filens_z = filens_z[snkey]
            plo, pmed, phi = gpr.get_profile_massmap(_filens, 
                                                     rbins_pkpc,
                                                     rbin_units='pkpc',
                                                     profiles=['perc-0.1', 
                                                               'perc-0.5', 
                                                               'perc-0.9'],
                                                     weightmap=True)
            meds_z = [gpr.get_profile_massmap(_fzs, 
                                              rbins_pkpc,
                                              rbin_units='pkpc',
                                              profiles=['perc-0.5'],
                                              weightmap=True)
                      for _fzs in _filens_z]
            meds_z = (np.array(meds_z)[:, 0, :])
            yerr = np.array([pmed - np.min(meds_z, axis=0), 
                             np.max(meds_z, axis=0) - pmed])
            ax.plot(rcens, plo, color=colors[opt], linestyle='dotted',
                    linewidth=1.)
            ax.plot(rcens, pmed, color=colors[opt], linestyle='solid',
                    linewidth=2.)
            ax.errorbar(rcens[opti::2], pmed[opti::2], yerr=yerr[:, opti::2],
                        color=colors[opt], linewidth=1.1, linestyle='none')
            ax.plot(rcens, phi, color=colors[opt], linestyle='dashed',
                    linewidth=1.5)
        if axi == 1:
            handles = [mlines.Line2D((), (), color='black',
                                     linewidth=1.5, linestyle='dashed',
                                     label='90%'),
                       mlines.Line2D((), (), color='black',
                                     linewidth=2., linestyle='solid',
                                     label='median'),
                       mlines.Line2D((), (), color='black',
                                     linewidth=1., linestyle='dotted',
                                     label='10%')]
        elif axi == 0:
            handles = [mlines.Line2D((), (), color=colors[opt],
                                     linewidth=1.5, linestyle='solid',
                                     label=f'CR heat. {opt}')
                       for opt in ['bug', 'fix']]
        ax.legend(handles=handles, fontsize=fontsize - 1, 
                  loc='lower left', handlelength=1.4,
                  frameon=False, bbox_to_anchor=(0.0, -0.02))
    ylims = [ax.get_ylim() for ax in axes]
    ymin = min([yl[0] for yl in ylims])
    ymax = max([yl[1] for yl in ylims])
    [ax.set_ylim((ymin, ymax)) for ax in axes]
    
    plt.savefig(outname, bbox_inches='tight')
            
