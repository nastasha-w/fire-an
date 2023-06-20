import h5py
import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np

import fire_an.analytic_halo.model_ionprof as mip
import fire_an.makeplots.plot_utils as pu
import fire_an.utils.constants_and_units as c
import fire_an.utils.opts_locs as ol

def readin_hm_data(filen, tomatch, ion):
    with h5py.File(mip.outdir_profiles + filen) as f:
        grpns_all = list(f.keys())
        grpns_use = [grpn for grpn in grpns_all if grpn.startswith(tomatch)]
        toplot = {}
        for grpn in grpns_use:
            grp = f[grpn]
            if bool(grp.attrs['failed']):
                continue
            lmvir = grp.attrs['logMvir_Msun_BN98']
            _sub = {}
            sgrp = grp[f'coldens_{ion}']
            _sub['b_kpc'] = sgrp['impactpar_cm'][:] / (c.cm_per_mpc * 1e-3)
            _sub['logN_cm2'] = np.log10(sgrp['coldens_cm2'][:])
            toplot[lmvir] = _sub
    return toplot


def plot_profiles_hm(filen, redshift, zsol, plind, mdot,
                     outname=None):
    tomatch = (f'z{redshift:.2f}_Zsolar{zsol:.2e}_vcplind{plind:.2f}'
               f'mdotMsunperYr{mdot:.2e}_logmvirMsun')
    data = readin_hm_data(filen, tomatch, 'Ne8')
    lmvs = np.sorted(list(data.keys()))
    vmin = lmvs[0]
    vmax = lmvs[-1]
    cmap = mcm.get_cmap('viridis')
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    fig = plt.figure(figsize=(5.5, 5.))
    grid = gsp.GridSpec(ncols=2, nrows=1, wspace=0.1, 
                        width_ratios=(10., 1.))
    ax = fig.add_subplot(grid[0])
    cax = fig.add_subplot(grid[1])
    fontsize = 12
    linewidth = 1.5
    path_effects = pu.getoutline(linewidth)

    for lmv in lmvs:
        xv = data[lmv]['b_kpc']
        yv = data[lmv]['logN_cm2']
        plt.plot(xv, yv, color=cmap((lmv - vmin) / (vmax - vmin)),
                 linewidth=linewidth, path_effects=path_effects)
    ax.set_xlabel('$\\mathrm{r}_{\\perp} \\; [\\mathrm{pkpc}]$',
                  fontsize=fontsize)
    ax.set_ylabel('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\,VIII})'
                  '\\; [\\mathrm{cm}^{-2}]$',
                  fontsize=fontsize)
    ax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                   top=True, right=True)

    plt.colorbar(mcm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                 orientation='vertical')
    cax.set_ylabel('$\\log_{10} \\, \\mathrm{M}_{\\mathrm{vir}}'
                   '\\; [\\mathrm{M}_{\\odot}]$',
                   fontsize=fontsize)
    cax.tick_params(labelsize=fontsize - 1.)
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

    