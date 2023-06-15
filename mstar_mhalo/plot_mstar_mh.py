import numpy as np
import matplotlib as mpl
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt

import fire_an.mstar_mhalo.loader_smdpl_sfr as ldsmdpl
import fire_an.utils.math_utils as mu
import fire_an.makeplots.plot_utils as pu

cosmopars_smdpl = {'omegalambda': 0.693,
                   'omegam': 0.307, 
                   'omegab': 0.048,
                   'h': 0.678}
def plot_mstar_mh(z_target):
    binsize = 0.1
    percvals = np.array([0.04, 0.2, 0.5, 0.8, 0.96])
    color_mhtoms = 'blue'
    color_mstomh = 'green'
    linestyles = ['dotted', 'dashed', 'solid', 'dashed', 'dotted']
    _cmap = mpl.cm.get_cmap('gist_yarg')
    cmap = pu.truncate_colormap(_cmap, minval=0., maxval=0.7)

    halos, a_used = ldsmdpl.loaddata(1. / (1. + z_target))
    z_used = 1. / a_used - 1.
    cosmopars =  cosmopars_smdpl.copy()
    cosmopars.update({'z': z_used, 'a': a_used})
    hsel = halos['upid'] == -1 # centrals
    mhalo_msun = halos['m'][hsel]  / cosmopars['h'] #BN98
    # for observed masses, would need to check the measurement method
    mstar_true_msun = halos['sm'][hsel] 

    ylabel = ('$\\log_10 \\, \\mathrm{M}_{\\star} \\;'
              ' [\\mathrm{M}_{\\mathrm{\\odot}}]$')
    xlabel = ('$\\log_10 \\, \\mathrm{M}_{\\mathrm{vir}} \\;'
              ' [\\mathrm{M}_{\\mathrm{\\odot}}]$')
    clabel = ('$\\log_10 \\, \\partial^2 \\mathrm{halo\\,frac.}'
              '\\, /\\, \\partial \\log_{10} \\mathrm{M}_{\\star}'
              '\\, /\\, \\partial \\log_{10} \\mathrm{M}_{\\mathrm{vir}}$')
    logmh = np.log10(mhalo_msun)
    logms = np.log10(mstar_true_msun)
    minmh = np.min(logmh)
    maxmh = np.max(logmh)
    b0 = np.floor(minmh / binsize) * binsize
    b1 = np.ceil(maxmh / binsize) * binsize
    mhbins = np.arange(b0, b1 + 0.5 * binsize, binsize)
    minms = np.min(logms)
    maxms = np.max(logms)
    b0 = np.floor(minms / binsize) * binsize
    b1 = np.ceil(maxms / binsize) * binsize
    msbins = np.arange(b0, b1 + 0.5 * binsize, binsize)

    hist, _bins = np.histogram([logmh, logms], bins=[mhbins, msbins])
    pvs = mu.percentiles_from_histogram(hist, msbins, axis=1,
                                        percentiles=percvals)
    mhcen = 0.5 * (mhbins[:-1] + mhbins[1:])
    
    fig = plt.figure((5.5, 5.))
    grid = gsp.GridSpec(ncols=2, nrows=1, wspace=0.1, 
                        width_ratios=[1., 0.1])
    ax = fig.add_subplot(grid[0])
    cax = fig.add_subplot(grid[1])
    fontsize = 12
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    
    phist = np.log10(hist / np.sum(hist) / binsize**2)
    pmax = np.max(phist)
    vmin = pmax - 7.
    img = ax.pcolormesh(mhbins, msbins, phist.T, cmap=cmap, vmin=vmin)
    plt.colorbar(img, cax=cax, extend='min')
    cax.set_ylabel(clabel, fontsize=fontsize)
    for pv, ls in zip(pvs, linestyles):
        ax.plot(mhcen, pv, linestyle=ls, color=color_mhtoms)
    
    pvs = mu.percentiles_from_histogram(hist, mhbins, axis=2,
                                        percentiles=percvals)
    mscen = 0.5 * (msbins[:-1] + msbins[1:])
    for pv, ls in zip(pvs, linestyles):
        ax.plot(pv, mscen, linestyle=ls, color=color_mstomh)

    

    




