import h5py
import numpy as np

import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt

import fire_an.makeplots.tol_colors as tc
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c

def readdata_mprof(filen):
    with h5py.File(filen) as f:
        mass_g = 10**f['histogram/histogram'][:]
        rbins_rvir = f['axis_0/bins'][:]
        cosmopars = {key: val for key, val 
                     in f['Header/cosmopars'].attrs.items()}
    out = {'mass_g': mass_g, 
           'rbins_rvir': rbins_rvir, 
           'cosmopars': cosmopars}
    return out
        

def plot_cumul_mass(filen_temp, zfills, ptfills,
                    rmax_mstar_rvir=0.1, rmin_fb_rvir=0.1, rmax_fb_rvir=1.,
                    title=None, outname=None):
    '''
    cumulative mass vs. radius for one halo,
    different redshifts and particle types
    '''

    hdata = [[readdata_mprof(filen_temp.format(**(zfill | ptfill))) 
              for ptfill in ptfills]
             for zfill in zfills]
    nz = len(zfills)
    _colors = tc.tol_cmap('rainbow_discrete', nz)
    colors = _colors(np.linspace(1. - 1. / nz, 1. / nz, nz))
    styleargs = {0: {'linestyle': 'dashed', 'linewidth': 1.},
                 1: {'linestyle': 'solid', 'linewidth': 1.},
                 2: {'dashes': (8, 2), 'linewidth': 1.},
                 4: {'linestyle': 'dotted', 'linewidth': 1.},
                 5: {'linestyle': 'dashdot', 'linewidth': 1.},
                 't': {'linestyle': 'solid', 'linewidth': 2.},
                 'b': {'linestyle': 'dashed', 'linewidth': 2.},
                 }
    ptnames = {0: 'gas',
               1: 'hi-res DM',
               2: 'lo-res DM',
               4: 'stars',
               5: 'BHs',
               't': 'all zoom',
               'b': 'all baryons'}
    fig = plt.figure(figsize=(8.5, 5.5))
    fontsize = 12
    grid = gsp.GridSpec(ncols=2, nrows=2,
                        width_ratios=[5., 3.5], height_ratios=[4., 1.5],
                        hspace=0., wspace=0.)
    lax = fig.add_subplot(grid[1, :])
    lax.axis('off')
    tax = fig.add_subplot(grid[0, 1])
    tax.axis('off')
    ax = fig.add_subplot(grid[0, 0])
    ax.tick_params(direction='in', labelsize=fontsize - 1., top=True,
                   right=True, which='both')
    ax.set_yscale('log')
    ax.set_xscale('log')
    xlab = '$\\mathrm{r}_{\\mathrm{3D}} \\; [\\mathrm{R}_{\\mathrm{vir}}]$'
    ax.set_xlabel(xlab, fontsize=fontsize)
    ylab = 'cumulative mass $[\\mathrm{M}_{\\odot}]$'
    ax.set_ylabel(ylab, fontsize=fontsize)

    if title is not None:
        fig.suptitle(title, fontsize=fontsize)

    zi_zval = {}
    ymax = -np.inf
    texttop = 1.0
    textvspace = 0.155
    for zi in range(len(zfills)):
        profile_tot = None
        profile_b = None
        rx_all = None
        color = colors[zi]
        for pi, ptfill in enumerate(ptfills):
            _hdat = hdata[zi][pi]
            rx = _hdat['rbins_rvir'][1:]
            my = np.cumsum(_hdat['mass_g']) / c.solar_mass
            pt = ptfill['pt']
            _styleargs = styleargs[pt]
            cosmopars = _hdat['cosmopars']
            zv = cosmopars['z']
            if zi in zi_zval:
                if not np.isclose(zi_zval[zi], zv):
                    raise RuntimeError('redshift mismatch')
            else:
                zi_zval[zi] = zv
            if rx_all is None:
                rx_all = rx
            else:
                if not np.allclose(rx_all, rx):
                    raise RuntimeError(f'r bins mismatch pt{pt}')
            if pt == 2:
                loi = np.where(my > 0.)[0]
                if len(loi) > 0:
                    loi = loi[0] - 1
                    xp = rx[loi]
                    ax.axvline(xp, ymin=0., ymax=0.07, color=color)
            else:
                if profile_tot is None:
                    profile_tot = my.copy()
                else:
                    profile_tot += my
            if pt not in [1, 2]:
                if profile_b is None:
                    profile_b = my.copy()
                else:
                    profile_b += my
            if pt == 0:
                imin_fb = np.where(np.isclose(rx, rmin_fb_rvir))[0][0]
                imax_fb = np.where(np.isclose(rx, rmax_fb_rvir))[0][0] 
                mgas_cgm = my[imax_fb] - my[imin_fb]
            elif pt == 4:
                imax_ms = np.where(np.isclose(rx, rmax_mstar_rvir))[0][0]
                mstar = my[imax_ms]

            ax.plot(rx, my, color=color, **_styleargs)
        ax.plot(rx_all, profile_b, color=color, **(styleargs['b']))
        ax.plot(rx_all, profile_tot, color=color, **(styleargs['t']))
        mvir = profile_tot[imax_fb]
        mb = profile_b[imax_fb]
        info = (f'$z={zi_zval[zi]:.1f}:\\;'
                '\\mathrm{M}_{\\mathrm{vir}} = '
                f'10^{{{np.log10(mvir):.1f}}}'
                '\\mathrm{M}_{\\odot},\\,'
                '\\mathrm{M}_{*} = '
                f'10^{{{np.log10(mstar):.1f}}}'
                '\\mathrm{M}_{\\odot},$'
                '\n\t'
                '$\\mathrm{f}_{\\mathrm{b}} = '
                f'{mb / mvir:.3f},\\,'
                '\\mathrm{f}_{\\mathrm{CGM}} = '
                f'{mgas_cgm / mvir:.3f},\\,'
                '$')
        tax.text(0.05, texttop, info, fontsize=fontsize,
                 verticalalignment='top', horizontalalignment='left',
                 transform=tax.transAxes)
        texttop -= textvspace

        #print(profile_tot)
        #print(profile_b)
        ymax = max(ymax, profile_tot[-1])
    info = ('$\\Omega_{b} \\, / \\, \\Omega_{m}  = '
            f'{cosmopars["omegab"] / cosmopars["omegam"]:.3f}$')
    tax.text(0.05, texttop, info, fontsize=fontsize,
             verticalalignment='top', horizontalalignment='left',
             transform=tax.transAxes)
    ax.axvline(1., color='black')
    ax.set_ylim(1e-5 * ymax, 1.2 * ymax)

    pthandles = [mlines.Line2D((), (), color='black',
                               label=ptnames[key], 
                               **styleargs[key])
                 for key in ptnames]
    zhandles = [mlines.Line2D((), (), linestyle='solid', linewidth=1.5,
                              label=f'$z={zi_zval[zi]:.1f}$', 
                              color=colors[zi])
                for zi in range(len(zfills))]
    handles = pthandles + zhandles
    lax.legend(handles=handles, fontsize=fontsize, ncol=5,
               bbox_to_anchor=(0.5, 0.5), loc='upper center')

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')            


def plotset_cumul_mass():
    fdir = '/Users/nastasha/ciera/profiles/fire/ptmasses_all2/'
    _filen_temp = ('hist_r3D_by_mass_pt{{pt}}_{simname}_snap{{snap}}'
                   '_bins1_v1.hdf5')
    _filen_temp = fdir + _filen_temp
    outdir = ('/Users/nastasha/ciera/projects_lead/fire3_ionabs'
              '/massprof_model3_all2/')

    simnames_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    simnames_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    simnames = simnames_hr + simnames_sr
    zfills_hr = [{'snap': val} for val in sl.snaps_hr]
    zfills_sr = [{'snap': val} for val in sl.snaps_sr]
    zfills = [zfills_hr] * len(simnames_hr) + [zfills_sr] * len(simnames_sr)
    pts_nobh = (0, 1, 2, 4)
    pts_bh = (0, 1, 2, 4, 5)
    ptfills = [[{'pt': pt} for pt in pts_nobh] if '_sdp1e10_' in simname 
               else [{'pt': pt} for pt in pts_bh]
               for simname in simnames]
    for zf, ptf, simname in zip(zfills, ptfills, simnames):
        filen_temp = _filen_temp.format(simname=simname)
        title = simname.split('_')[0] + ', '
        title += ('noBH' if '_sdp1e10_' in simname
                  else 'AGN-CR' if '_MHDCRspec1_' in simname
                  else 'AGN-noCR')
        outname = outdir + f'massprof_cumul_{simname}.pdf'
        plot_cumul_mass(filen_temp, zf, ptf,
                        rmax_mstar_rvir=0.1, rmin_fb_rvir=0.1, 
                        rmax_fb_rvir=1.,
                        title=title, outname=outname)
        plt.close()