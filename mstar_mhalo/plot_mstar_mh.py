import matplotlib as mpl
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sps

import fire_an.mstar_mhalo.analytical as an
import fire_an.mstar_mhalo.loader_smdpl_sfr as ldsmdpl
import fire_an.utils.math_utils as mu
import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc

imgdir = '/projects/b1026/nastasha/imgs/datacomp/shmh/'

cosmopars_smdpl = {'omegalambda': 0.693,
                   'omegam': 0.307, 
                   'omegab': 0.048,
                   'h': 0.678}

def gethist_mhalo_mstarcen(z_target, binsize=0.1):
    halos, a_used = ldsmdpl.loaddata(1. / (1. + z_target))
    z_used = 1. / a_used - 1.
    cosmopars =  cosmopars_smdpl.copy()
    cosmopars.update({'z': z_used, 'a': a_used})
    hsel = halos['upid'] == -1 # centrals
    mhalo_msun = halos['m'][hsel]  / cosmopars['h'] #BN98
    # for observed masses, would need to check the measurement method
    mstar_true_msun = halos['sm'][hsel] 
    logmh = np.log10(mhalo_msun)
    logms = np.log10(mstar_true_msun)
    minmh = np.min(logmh)
    maxmh = np.max(logmh)
    b0 = np.floor(minmh / binsize) * binsize
    b1 = np.ceil(maxmh / binsize) * binsize
    mhbins = np.arange(b0, b1 + 0.5 * binsize, binsize)
    minms = np.min(logms[np.isfinite(logms)])
    maxms = np.max(logms)
    b0 = np.floor(minms / binsize) * binsize
    b1 = np.ceil(maxms / binsize) * binsize
    msbins = np.arange(b0, b1 + 0.5 * binsize, binsize)
    msbins = np.append([-np.inf], msbins)
    _msbins = np.copy(msbins)
    _msbins[0] = _msbins[1] - binsize
    
    chunksize = 1_000_000
    arlen = len(logmh)
    nchunks = (arlen - 1) // chunksize + 1
    # bus error without feeding smaller bits of array into histogramdd
    for i in range(nchunks):
        sel = slice(chunksize * i, chunksize * (i + 1), None)
        _hist, _bins = np.histogramdd([logmh[sel], logms[sel]], 
                                      bins=[mhbins, msbins])
        if i == 0:
            hist = _hist
        else:
            hist += _hist
    return hist, msbins, _msbins, mhbins, cosmopars

def plot_mstar_mh(z_target):
    binsize = 0.1
    percvals = np.array([0.04, 0.2, 0.5, 0.8, 0.96])
    _colors = tc.tol_cset('vibrant')
    color_mhtoms = _colors[0]
    color_mstomh = _colors[1]
    color_moster13 = _colors[2]
    color_burchett19 = _colors[3]
    linestyles = ['dotted', 'dashed', 'solid', 'dashed', 'dotted']
    _cmap = mpl.cm.get_cmap('gist_yarg')
    cmap = pu.truncate_colormap(_cmap, minval=0., maxval=0.7)
    linewidth = 1.5
    path_effects = pu.getoutline(linewidth)

    hist, msbins, _msbins, mhbins, cosmopars = \
        gethist_mhalo_mstarcen(z_target, binsize=binsize)
    z_used = cosmopars['z']

    ylabel = ('$\\log_{10} \\, \\mathrm{M}_{\\star} \\;'
              ' [\\mathrm{M}_{\\mathrm{\\odot}}]$')
    xlabel = ('$\\log_{10} \\, \\mathrm{M}_{\\mathrm{vir}} \\;'
              ' [\\mathrm{M}_{\\mathrm{\\odot}}]$')
    clabel = ('$\\log_{10} \\, \\partial^2 \\mathrm{halo\\,frac.}'
              '\\, /\\, \\partial \\log_{10} \\mathrm{M}_{\\star}'
              '\\, /\\, \\partial \\log_{10} \\mathrm{M}_{\\mathrm{vir}}$')
    
    pvs = mu.percentiles_from_histogram(hist, msbins, axis=1,
                                        percentiles=percvals)
    mhcen = 0.5 * (mhbins[:-1] + mhbins[1:])
    
    fig = plt.figure(figsize=(5.5, 5.))
    grid = gsp.GridSpec(ncols=2, nrows=1, wspace=0.1, 
                        width_ratios=[1., 0.1])
    ax = fig.add_subplot(grid[0])
    cax = fig.add_subplot(grid[1])
    fontsize = 12
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(which='both', labelsize=fontsize - 1,
                   direction='in', top=True, right=True)
    
    phist = np.log10(hist / np.sum(hist) / binsize**2)
    pmax = np.max(phist)
    vmin = pmax - 6.
    img = ax.pcolormesh(mhbins, _msbins, phist.T, cmap=cmap, vmin=vmin)
    plt.colorbar(img, cax=cax, extend='min')
    cax.set_ylabel(clabel, fontsize=fontsize)
    for pv, ls in zip(pvs, linestyles):
        ax.plot(mhcen, pv, linestyle=ls, color=color_mhtoms,
                linewidth=linewidth, path_effects=path_effects)
    
    pvs = mu.percentiles_from_histogram(hist, mhbins, axis=0,
                                        percentiles=percvals)
    mscen = 0.5 * (msbins[:-1] + msbins[1:])
    mscen[0] = mscen[1] - binsize
    for pv, ls in zip(pvs, linestyles):
        ax.plot(pv, mscen, linestyle=ls, color=color_mstomh,
                linewidth=linewidth, path_effects=path_effects)
    
    xlims = ax.get_xlim()
    ax.set_xlim(10.5, xlims[1])
    ylims = ax.get_ylim()
    ax.set_ylim(7., ylims[1])

    yv_moster13 = np.log10(an.mstar_moster_etal_2013(10**mhcen, z_used))
    ax.plot(mhcen, yv_moster13, color=color_moster13, linewidth=linewidth,
            linestyle='dashdot', path_effects=path_effects)
    yv_burchett19 = np.log10(an.mstar_burchett_etal_2019(10**mhcen, z_used))
    ax.plot(mhcen, yv_burchett19, color=color_burchett19, linewidth=linewidth,
            linestyle='dashdot', path_effects=path_effects)
    
    label_mhtoms = '$\\mathrm{M}_{\\star}(\\mathrm{M}_{\\mathrm{vir}})$'
    label_mstomh = '$\\mathrm{M}_{\\mathrm{vir}}(\\mathrm{M}_{\\star})$'
    label_moster13 = ('M+13'
                      ' $\\mathrm{M}_{\\star}(\\mathrm{M}_{\\mathrm{200c}})$')
    label_burchett19 = ('B+19 '
                      ' $\\mathrm{M}_{\\star}(\\mathrm{M}_{\\mathrm{200c}})$')
    handles = [mlines.Line2D((), (), color='black',
                             linestyle=ls, label=f'{pv * 100.:.0f}%')
               for ls, pv in zip(linestyles, percvals)]
    handles = handles + \
              [mlines.Line2D((), (), color=color_mhtoms,
                             linestyle='solid', 
                             label=label_mhtoms,
                             linewidth=linewidth, 
                             path_effects=path_effects),
               mlines.Line2D((), (), color=color_mstomh,
                             linestyle='solid', 
                             label=label_mstomh,
                             linewidth=linewidth, 
                             path_effects=path_effects),
               mlines.Line2D((), (), color=color_moster13, 
                             linewidth=linewidth,
                             linestyle='dashdot', path_effects=path_effects,
                             label=label_moster13),
               mlines.Line2D((), (), color=color_burchett19, linewidth=linewidth,
                            linestyle='dashdot', path_effects=path_effects,
                            label=label_burchett19)]
    ax.legend(handles=handles, fontsize=fontsize -1, loc='lower right',
              ncol=2, bbox_to_anchor=(1.00, 0.00),
              handlelength=2., columnspacing=1.,
              handletextpad=0.4, framealpha=0.3)
    outname = f'mstar_mhalo_relation_universemachine_smdpl_z{z_used:.2f}'
    outname = imgdir + outname.replace('.', 'p') + '.pdf'
    plt.savefig(outname, bbox_inches='tight')

# kinda random from Burchett et al. (2019)
# (from the top of table 1, just got three values across 
# the stellar mass range)
# I'm going to assume it's 1 sigma for now...
# they don't seem to say in the paper.
_logmstar_examples = [(10.9, 0.1), (10.5, 0.1), (9.8, 0.2), (9.0, 0.2)]
    
def plot_mhdist_for_mstar(z_target):
    logmstar_examples = _logmstar_examples
    ls_examples = ['solid', 'dashed', 'dotted', 'dashdot']
    markers_examples = ['o', 's', 'P', 'h']
    _colors = tc.tol_cset('vibrant')
    color_mhtoms = _colors[0]
    color_mstomh = _colors[1]
    color_moster13 = _colors[2]
    color_burchett19 = _colors[3]
    color_mstomh_bayes = _colors[5]
    colors = [color_mhtoms, color_moster13, color_burchett19,
              color_mstomh, color_mstomh_bayes]
    colorlabels = ['med. $\\mathrm{M}_{\\star}(\\mathrm{M}_{\\mathrm{vir}})$',
                   ('M+13 $\\mathrm{M}_{\\star}'
                    '(\\mathrm{M}_{\\mathrm{200c}})$'),
                   ('B+19 $\\mathrm{M}_{\\star}'
                    '(\\mathrm{M}_{\\mathrm{200c}})$'),
                   'med. $\\mathrm{M}_{\\mathrm{vir}}(\\mathrm{M}_{\\star})$',
                   ('dist. $\\mathrm{M}_{\\mathrm{vir}}'
                    '(\\mathrm{M}_{\\star})$'),
                   ]

    hist, msbins, _msbins, mhbins, cosmopars = \
          gethist_mhalo_mstarcen(z_target, binsize=0.1)
    transmat_ms_to_mh = hist / np.sum(hist, axis=0)[np.newaxis, :]
    mhbins_fine = np.arange(mhbins[0], mhbins[-1] + 1e-4, 0.02)
    mscens = 0.5 * (_msbins[:-1] + _msbins[1:])
    mhcens = 0.5 * (mhbins[:-1] + mhbins[1:])
    # median Mh to Ms, solve for Mh
    medms_frommh = mu.percentiles_from_histogram(hist, msbins, axis=1,
        percentiles=np.array([0.5]))[0]
    medmh_fromms = mu.percentiles_from_histogram(hist, mhbins, axis=0,
        percentiles=np.array([0.5]))[0]
    fig = plt.figure(figsize=(5.5, 5.))
    ax = fig.add_subplot(1, 1, 1)
    fontsize = 12
    xlabel = ('$\\log_{10} \\, \\mathrm{M}_{\\mathrm{h}} \\;'
              ' [\\mathrm{M}_{\\odot}]$')
    ylabel = ('$\\partial \\mathrm{P} \\,/\\, '
              '\\partial \\log_{10} \\, \\mathrm{M}_{\\mathrm{h}}$')
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(which='both', labelsize=fontsize - 1,
                   direction='in', top=True, right=True)
    y_bestest = 0.8
    # last few halo mass bins get noisy, non-monotonic
    cutoff_mh = np.where(np.diff(medms_frommh) <= 0.)[0][0] + 1
    # for ms to mh, issues are at low and high masses
    cutoff_ms_lo = np.where(np.logical_and(np.diff(medmh_fromms) <= 0., 
                                           medmh_fromms[:-1] < 10))[0][-1]
    cutoff_ms_hi = np.where(np.logical_and(np.diff(medmh_fromms) <= 0., 
                                           medmh_fromms[:-1] > 10))[0] 
    if len(cutoff_ms_hi) >= 1:
        cutoff_ms_hi = cutoff_ms_hi[0] + 1
    else:
        cutoff_ms_hi = len(medmh_fromms)

    for lms_example, ls, mk in zip(logmstar_examples, ls_examples, 
                                   markers_examples):
        lms = lms_example[0]
        lms_err = lms_example[1]

        # Mh to Ms, propagate M* uncertainties through median Ms(Mh)
        _msi_max = np.where(_msbins > medms_frommh[cutoff_mh])[0][0]
        _msi_min = np.where(_msbins < medms_frommh[0])[0][-1] + 1

        mh_mhtoms = mu.linterpsolve(medms_frommh[:cutoff_mh], 
                                    mhcens[:cutoff_mh], lms)
        #print(medms_frommh[:cutoff_mh])
        #print(_msbins[_msi_min:_msi_max])
        mhvals = [mu.linterpsolve(medms_frommh[:cutoff_mh],
                                  mhcens[:cutoff_mh], _ms)
                  for _ms in _msbins[_msi_min:_msi_max]]
        mhvals = np.array(mhvals)
        pdist_ms = sps.erf((_msbins[_msi_min + 1:_msi_max] - lms) / lms_err) \
                   - sps.erf((_msbins[_msi_min :_msi_max-1] - lms) / lms_err)
        mh_pdist = pdist_ms / np.diff(mhvals)
        mh_pdistcens = 0.5 * (mhvals[:-1] + mhvals[1:])
        ax.plot(mh_pdistcens, mh_pdist, color=color_mhtoms,
                linestyle=ls)
        y_cenest = mu.linterpsolve(mh_pdistcens, mh_pdist, mh_mhtoms)
        ax.scatter([mh_mhtoms], [y_cenest],
                   color=color_mhtoms, marker=mk, s=20)
        
        # Ms to Mh, propagate M* uncertainties through median Mh(Ms)
        _msi_max = np.where(mhbins > medmh_fromms[cutoff_mh])[0][0]
        _msi_min = np.where(mhbins < medmh_fromms[0])[0][-1] + 1
        print(mscens[cutoff_ms_lo:cutoff_ms_hi])
        print(lms)
        mh_mstomh = mu.linterpsolve(mscens[cutoff_ms_lo:cutoff_ms_hi], 
                                    medmh_fromms[cutoff_ms_lo:cutoff_ms_hi],
                                    lms)
        mhvals = [mu.linterpsolve(mscens[cutoff_ms_lo:cutoff_ms_hi], 
                                  medmh_fromms[cutoff_ms_lo:cutoff_ms_hi],
                                  _ms)
                  for _ms in _msbins[cutoff_ms_lo + 1:cutoff_ms_hi - 1]]
        mhvals = np.array(mhvals)
        pdist_ms = sps.erf((_msbins[cutoff_ms_lo + 2: cutoff_ms_hi - 1] 
                            - lms) / lms_err) \
                   - sps.erf((_msbins[cutoff_ms_lo + 1: cutoff_ms_hi - 2] 
                              - lms) / lms_err)
        mh_pdist = pdist_ms / np.diff(mhvals)
        mh_pdistcens = 0.5 * (mhvals[:-1] + mhvals[1:])
        ax.plot(mh_pdistcens, mh_pdist, color=color_mstomh,
                linestyle=ls)
        y_cenest = mu.linterpsolve(mh_pdistcens, mh_pdist, mh_mstomh)
        ax.scatter([mh_mstomh], [y_cenest],
                   color=color_mstomh, marker=mk, s=20)
        
        # Ms to Mh, propagate M* uncertainties through Mh(Ms) distribution
        pdist_ms = sps.erf((msbins[1:] - lms) / lms_err) \
                   - sps.erf((_msbins[:-1] - lms) / lms_err)
        mh_pdist = np.sum(transmat_ms_to_mh * pdist_ms[np.newaxis, :], axis=1)
        mhvals = mhbins
        mh_pdist /= np.diff(mhvals)
        mh_pdistcens = 0.5 * (mhvals[:-1] + mhvals[1:])
        ax.plot(mh_pdistcens, mh_pdist, color=color_mstomh_bayes,
                linestyle=ls)
        
        # Mh to Ms (Moster et al. 2013), 
        # propagate M* uncertainties through median Ms(Mh)
        msfrommh = np.log10(an.mstar_moster_etal_2013(10**mhbins_fine, 
                                                      cosmopars['z']))
        mh_mhtoms = mu.linterpsolve(msfrommh, mhbins_fine, lms)
        mhvals = mhbins_fine
        mhvals = np.array(mhvals)
        pdist_ms = sps.erf((msfrommh[1:] - lms) / lms_err) \
                   - sps.erf((msfrommh[:-1] - lms) / lms_err)
        mh_pdist = pdist_ms / np.diff(mhvals)
        mh_pdistcens = 0.5 * (mhvals[:-1] + mhvals[1:])
        ax.plot(mh_pdistcens, mh_pdist, color=color_moster13,
                linestyle=ls)
        y_cenest = mu.linterpsolve(mh_pdistcens, mh_pdist, mh_mhtoms)
        ax.scatter([mh_mhtoms], [y_cenest],
                   color=color_moster13, marker=mk, s=20)
        
        # Mh to Ms (Burchett et al. 2019), 
        # propagate M* uncertainties through median Ms(Mh)
        msfrommh = np.log10(an.mstar_burchett_etal_2019(10**mhbins_fine,
                                                        cosmopars['z']))
        mh_mhtoms = mu.linterpsolve(msfrommh, mhbins_fine, lms)
        mhvals = mhbins_fine
        mhvals = np.array(mhvals)
        pdist_ms = sps.erf((msfrommh[1:] - lms) / lms_err) \
                   - sps.erf((msfrommh[:-1] - lms) / lms_err)
        mh_pdist = pdist_ms / np.diff(mhvals)
        mh_pdistcens = 0.5 * (mhvals[:-1] + mhvals[1:])
        ax.plot(mh_pdistcens, mh_pdist, color=color_burchett19,
                linestyle=ls)
        y_cenest = mu.linterpsolve(mh_pdistcens, mh_pdist, mh_mhtoms)
        ax.scatter([mh_mhtoms], [y_cenest],
                   color=color_burchett19, marker=mk, s=20)
    ax.set_xlim(10.7, 13.7)
    ylim = ax.get_ylim()
    yr = ylim[1] - ylim[0]
    ax.set_ylim(ylim[0], ylim[1] + 0.2 * yr)

    handles1 = [mlines.Line2D((), (), color=color, label=clabel)
                for color, clabel in zip(colors, colorlabels)]
    handles2 = [mlines.Line2D((), (), color='black', linestyle=ls,
                              label=(f'${lms[0]:.1f} \\pm {lms[1]:.1f}$'))
                for ls, lms in zip(ls_examples, logmstar_examples)]
    leg1 = ax.legend(handles=handles1, fontsize=fontsize - 2, ncol=3,
              loc='upper center', bbox_to_anchor=(0.5, 1.0),
              handlelength=1., columnspacing=0.7, handletextpad=0.4)
    leg2 = ax.legend(handles=handles2, fontsize=fontsize - 1, ncol=1,
              loc='upper right', bbox_to_anchor=(1.0, 0.8),
              handlelength=2., columnspacing=1., handletextpad=0.4,
              title=('$\\log_{10} \\, \\mathrm{M}_{\\star} \\;'
                    ' [\\mathrm{M}_{\\odot}]$'), 
              title_fontsize=fontsize - 1)
    ax.add_artist(leg1)
    outname = f'mhalo_from_example_mstar_z{cosmopars["z"]:.2f}'
    outname = imgdir + outname.replace('.', 'p') + '.pdf'
    plt.savefig(outname, bbox_inches='tight')
        

    
   




