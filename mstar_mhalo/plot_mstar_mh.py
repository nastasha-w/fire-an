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

## directory where the images go
imgdir = '/projects/b1026/nastasha/imgs/datacomp/smhm/'

def gethist_mhalo_mstarcen(x, binsize=3):
    return None

## copied from astropy.cosmology Planck15, 
## because that's what the Prochaska et al. (2019) CASBaH galaxy survey
## paper references.
cosmopars_p19 = {'omegam': 0.307, 'omegalambda': 1. - 0.307,
                 'omegab': 0.0486, 'h': 0.677}

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
        

def plot_mstar_mh(zs_target, variation='main'):
    '''
    zs_target: array-like of floats
        redshifts to compute the relation for. (Closest UM
        tabulated values are used.)
    variation: {'main', 'appendix'}
        show the UM relation + B+19 (main text) or the UM
        relation + median M* at Mh (appendix)
    '''
    binsize = 0.1
    # 1, 2 sigma
    sig1 = an.cumulgauss(1.) - an.cumulgauss(-1.)
    #sig2 = an.cumulgauss(2.) - an.cumulgauss(-2.)
    #percvals = np.array([0.5 - 0.5 * sig2, 0.5 - 0.5 * sig1, 0.5,
    #                     0.5 + 0.5 * sig1, 0.5 + 0.5 * sig2])
    percvals = np.array([0.5 - 0.5 * sig1, 0.5, 0.5 + 0.5 * sig1])
    _colors = tc.tol_cset('vibrant')
    color_mhtoms = _colors[1]
    color_mstomh = _colors[2]
    #color_moster13 = _colors[2]
    color_burchett19 = _colors[3]
    #linestyles = ['dotted', 'dashed', 'solid', 'dashed', 'dotted']
    linestyles = ['dashed', 'solid', 'dashed']
    _cmap = mpl.cm.get_cmap('gist_yarg')
    cmap = pu.truncate_colormap(_cmap, minval=0., maxval=0.7)
    linewidth = 1.5
    path_effects = pu.getoutline(linewidth)

    xlabel = ('$\\log_{10} \\, \\mathrm{M}_{\\star} \\;'
              ' [\\mathrm{M}_{\\mathrm{\\odot}}]$')
    ylabel = ('$\\log_{10} \\, \\mathrm{M}_{\\mathrm{vir}} \\;'
              ' [\\mathrm{M}_{\\mathrm{\\odot}}]$')
    #clabel = ('$\\log_{10} \\, \\partial^2 \\mathrm{halo\\,frac.}'
    #          '\\, /\\, \\partial \\log_{10} \\mathrm{M}_{\\star}'
    #          '\\, /\\, \\partial \\log_{10} \\mathrm{M}_{\\mathrm{vir}}$')
    clabel = ('$\\log_{10} \\, \\mathrm{halo\\,frac.}'
              '\\,/\\, \\mathrm{bin\\,size}$')
    
    npanels = len(zs_target)
    ncols = min(npanels, 3)
    nrows = (npanels - 1) // ncols + 1
    panelsize = 3.
    width_ratios = [panelsize] * ncols +  [0.1 * panelsize]
    height_ratios = [panelsize] * nrows
    fig = plt.figure(figsize=(sum(width_ratios), sum(height_ratios)))
    grid = gsp.GridSpec(ncols=ncols + 1, nrows=nrows, 
                        wspace=0.0, hspace=0.0,
                        width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols])
            for i in range(npanels)]
    cax = fig.add_subplot(grid[:, ncols])
    fontsize = 12
    vmax = 0.
    vmin = vmax - 6.
    
    zs_used = []
    for axi, (ax, z_target) in enumerate(zip(axes, zs_target)):
        dobottom = axi >= npanels - ncols
        doleft = axi % ncols == 0
        if dobottom:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        if doleft:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.tick_params(which='both', labelsize=fontsize - 1,
                       direction='in', top=True, right=True,
                       labelleft=doleft, labelbottom=dobottom)

        histobj = ldsmdpl.SMHMhists(np.array([z_target]), binsize=binsize)
        msbins_hist = histobj.getbins(z_target, mode='ms')
        mhbins_hist = histobj.getbins(z_target, mode='mh')
        hist_raw = histobj.gethist(z_target)
        z_used = list(histobj.hists.keys())[0]
        zs_used.append(z_used)
        ax.text(0.5, 0.97, f'$z = {z_used:.2f}$',
                fontsize=fontsize, transform=ax.transAxes,
                horizontalalignment='center', verticalalignment='top')

        phist = np.log10(hist_raw / np.sum(hist_raw) / binsize**2)
        img = ax.pcolormesh(msbins_hist, mhbins_hist, phist, cmap=cmap, 
                            vmin=vmin, vmax=vmax)
        
        msx, mhys = histobj.getperc_msmh(z_target, mode='mstomh', 
                                         percvals=percvals)
        for mhy, ls in zip(mhys, linestyles):
            ax.plot(msx, mhy, linestyle=ls, color=color_mstomh,
                    linewidth=linewidth, path_effects=path_effects)
        
        if variation == 'appendix':
            msxs, mhy = histobj.getperc_msmh(z_target, mode='mhtoms', 
                                            percvals=np.array([0.5]))
            ax.plot(msxs[0], mhy, linestyle='solid', color=color_mhtoms,
                    linewidth=linewidth, path_effects=path_effects)
        
        xlims = ax.get_xlim()
        ax.set_xlim(8., xlims[1])
        ylims = ax.get_ylim()
        ax.set_ylim(10.5, ylims[1])
        
        #xv_moster13 = np.log10(an.mstar_moster_etal_2013(10**mhbins_hist, z_used))
        #ax.plot(xv_moster13, mhbins_hist, color=color_moster13, linewidth=linewidth,
        #        linestyle='dashdot', path_effects=path_effects)
        if variation == 'main':
            logm200c = mhbins_hist
            cosmopars = cosmopars_p19.copy()
            cosmopars.update({'z': z_used, 'a': 1. / (1. + z_used)})
            c200_b19 = an.concentration_mass_m200c_dm14(10**logm200c, 
                                                        cosmopars)
            mh_b19 = np.array([an.m200c_to_mvir(m200c_msun, cosmopars, c200)
                               for m200c_msun, c200 in zip(10**logm200c,
                                                           c200_b19)])
            print(z_used, ':')
            print(logm200c)
            print(np.log10(mh_b19))
            ms_burchett19 = np.log10(an.mstar_burchett_etal_2019
                                     (10**logm200c, z_used))
            ax.plot(ms_burchett19, np.log10(mh_b19), color=color_burchett19, 
                    linewidth=linewidth,
                    linestyle='dashdot', path_effects=path_effects)
    xlims = [ax.get_xlim() for ax in axes]
    xmin = min([xlim[0] for xlim in xlims])
    xmax = max([xlim[1] for xlim in xlims])
    [ax.set_xlim((xmin, xmax)) for ax in axes]
    ylims = [ax.get_ylim() for ax in axes]
    ymin = min([ylim[0] for ylim in ylims])
    ymax = max([ylim[1] for ylim in ylims])
    [ax.set_ylim((ymin, ymax)) for ax in axes]

    plt.colorbar(img, cax=cax, extend='min')
    cax.set_ylabel(clabel, fontsize=fontsize)
    
    label_mhtoms = '$\\mathrm{M}_{\\star}(\\mathrm{M}_{\\mathrm{vir}})$'
    label_mstomh = '$\\mathrm{M}_{\\mathrm{vir}}(\\mathrm{M}_{\\star})$'
    #label_moster13 = ('M+13'
    #                  ' $\\mathrm{M}_{\\star}(\\mathrm{M}_{\\mathrm{200c}})$')
    label_burchett19 = ('B+19')
    #' $\\mathrm{M}_{\\star}(\\mathrm{M}_{\\mathrm{200c}})$')
    handles1 = [mlines.Line2D((), (), color='black',
                             linestyle=ls, label=f'{pv * 100.:.0f}%')
                for ls, pv in zip(linestyles, percvals)]
    handles2 = [mlines.Line2D((), (), color=color_mhtoms,
                              linestyle='solid', 
                              label=label_mhtoms,
                              linewidth=linewidth, 
                              path_effects=path_effects),
                mlines.Line2D((), (), color=color_mstomh,
                              linestyle='solid', 
                              label=label_mstomh,
                              linewidth=linewidth, 
                              path_effects=path_effects),
                mlines.Line2D((), (), color=color_burchett19, 
                              linewidth=linewidth,
                              linestyle='dashdot', 
                              path_effects=path_effects,
                              label=label_burchett19)]
    if variation == 'main':
        handles2 = handles2[1:]
    elif variation == 'appendix':
        handles2 = handles2[:2]
    if len(axes) == 1:
        axes[0].legend(handles=handles1 + handles2, 
                       fontsize=fontsize -1, 
                       loc='upper left',
                       ncol=2, bbox_to_anchor=(0.00, 0.9),
                       handlelength=2., columnspacing=1.,
                       handletextpad=0.4, framealpha=0.3)
    else:
        axes[0].legend(handles=handles1, 
                       fontsize=fontsize - 1, 
                       loc='upper left',
                       ncol=1, bbox_to_anchor=(0.00, 0.90),
                       handlelength=2., columnspacing=1.,
                       handletextpad=0.4, framealpha=0.3)
        axes[1].legend(handles=handles2, 
                       fontsize=fontsize - 1, 
                       loc='upper left',
                       ncol=1, bbox_to_anchor=(0.00, 0.90),
                       handlelength=2., columnspacing=1.,
                       handletextpad=0.4, framealpha=0.3)
    zstr = '_'.join([f'{z_used:.2f}' for z_used in zs_used])
    outname = ('mhalo_mstar_relation_universemachine_smdpl_z'
               f'_{zstr}_{variation}')
    outname = imgdir + outname.replace('.', 'p') + '.pdf'
    plt.savefig(outname, bbox_inches='tight')    
   

_logmstar_examples = [(10.9, 0.1), (10.5, 0.1), (9.8, 0.2), (9.0, 0.2)]
def plot_mhdist_for_mstar_example_scatterdecomp(z_target):
    binsize = 0.1
    logmstar_examples = _logmstar_examples
    ls_examples = ['solid', 'dashed', 'dotted', 'dashdot']
    markers_examples = ['o', 's', 'P', 'h']
    _colors = tc.tol_cset('vibrant')
    color_mhtoms = _colors[0]
    color_mstomh_obs_scat = _colors[1]
    #color_moster13 = _colors[2]
    #color_burchett19 = _colors[3]
    color_mstomh_obs_smhm_scat = _colors[5]
    color_mstomh_smhm_scat = _colors[2]

    colors = [color_mhtoms, 
              color_mstomh_obs_scat, 
              color_mstomh_obs_smhm_scat,
              color_mstomh_smhm_scat]
    colorlabels = [('med. $\\mathrm{M}_{\\star}(\\mathrm{M}_{\\mathrm{vir}}),'
                    ' \\sigma(\\mathrm{obs.})$'),
                   ('med. $\\mathrm{M}_{\\mathrm{vir}}(\\mathrm{M}_{\\star}),'
                    ' \\sigma(\\mathrm{obs.})$'),
                   #('M+13 $\\mathrm{M}_{\\star}'
                   # '(\\mathrm{M}_{\\mathrm{200c}})$'),
                   #('B+19 $\\mathrm{M}_{\\star}'
                   # '(\\mathrm{M}_{\\mathrm{200c}})$'),
                   ('$\\mathrm{M}_{\\mathrm{vir}}(\\mathrm{M}_{\\star}), '
                    '\\sigma(\\mathrm{SMHM}, \\mathrm{obs.})$'),
                   ('$\\mathrm{M}_{\\mathrm{vir}}(\\mathrm{M}_{\\star}), '
                    '\\sigma(\\mathrm{SMHM})$'),
                   ]
    histobj = ldsmdpl.SMHMhists(np.array([z_target]), binsize=binsize)
    msbins_hist = histobj.getbins(z_target, mode='ms')
    mhbins_hist = histobj.getbins(z_target, mode='mh')
    z_used = list(histobj.hists.keys())[0]

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
    for lms_example, ls, mk in zip(logmstar_examples, ls_examples, 
                                   markers_examples):
        lms = lms_example[0]
        lms_err = lms_example[1]
        msbins_p = msbins_hist
        msgaussdist = an.cumulgauss((msbins_p[1:] - lms) / lms_err) \
                      - an.cumulgauss((msbins_p[:-1] - lms) / lms_err)
        msbins_delta = np.array([lms - 1e-4, lms + 1e-4])
        msdeltadist = np.array([1.])

        # Mh to Ms, propagate M* uncertainties through median Ms(Mh)
        ms_mhtoms, mh_mhtoms = histobj.getperc_msmh(z_target, mode='mhtoms', 
                                                    percvals=np.array([0.5]))
        ms_mhtoms = ms_mhtoms[0]
        #print(ms_mhtoms, mh_mhtoms, msbins_p, msgaussdist)
        mhbins, mhpd = an.calcdist(ms_mhtoms, mh_mhtoms, msbins_p, 
                                   msgaussdist,
                                   filter_monorels=True,
                                   cutoff_xbins=True)
        mhdelta, _ = an.calcdist(ms_mhtoms, mh_mhtoms, msbins_delta, 
                                 msdeltadist,
                                 filter_monorels=True,
                                 cutoff_xbins=True)
        mhdelta = np.average(mhdelta)
        _mhcens = 0.5 * (mhbins[:-1] + mhbins[1:])

        yv_cenest = mu.linterpsolve(_mhcens, mhpd, mhdelta)
        ax.plot(_mhcens, mhpd, color=color_mhtoms, linestyle=ls)
        ax.scatter([mhdelta], [yv_cenest],
                   color=color_mhtoms, marker=mk, s=20)
        
        # Ms to Mh, propagate M* uncertainties through median Mh(Ms)
        ms_mstomh, mh_mstomh = histobj.getperc_msmh(z_target, mode='mstomh', 
                                                    percvals=np.array([0.5]))
        mh_mstomh = mh_mstomh[0]
        mhbins, mhpd = an.calcdist(ms_mstomh, mh_mstomh, msbins_p, 
                                   msgaussdist,
                                   filter_monorels=True,
                                   cutoff_xbins=True)
        mhdelta, _ = an.calcdist(ms_mstomh, mh_mstomh, msbins_delta, 
                                 msdeltadist,
                                 filter_monorels=True,
                                 cutoff_xbins=True)
        mhdelta = np.average(mhdelta)
        _mhcens = 0.5 * (mhbins[:-1] + mhbins[1:])

        yv_cenest = mu.linterpsolve(_mhcens, mhpd, mhdelta)
        ax.plot(_mhcens, mhpd, color=color_mstomh_obs_scat, linestyle=ls)
        ax.scatter([mhdelta], [yv_cenest],
                   color=color_mstomh_obs_scat, marker=mk, s=20)
        
        # Ms to Mh, propagate M* uncertainties through Mh(Ms) distribution
        mhp, _mhbins = histobj.matrixconv(msgaussdist, z_target, 
                                          mode='mstomh')
        _mhcens = 0.5 * (_mhbins[:-1] + _mhbins[1:])
        ax.plot(_mhcens, mhp / np.diff(_mhbins), 
                color=color_mstomh_obs_smhm_scat,
                linestyle=ls)
        # Ms to Mh, uncertainties only from Mh(Ms) distribution
        _msdist = np.zeros(len(msbins_hist) - 1, dtype=np.float)
        spikei = np.searchsorted(msbins_hist, lms) - 1
        _msdist[spikei] = 1.
        mhp, _mhbins = histobj.matrixconv(_msdist, z_target, 
                                          mode='mstomh')
        _mhcens = 0.5 * (_mhbins[:-1] + _mhbins[1:])
        ax.plot(_mhcens,  mhp / np.diff(_mhbins), 
                color=color_mstomh_smhm_scat, linestyle=ls)
    ax.set_xlim(10.7, 13.7)
    ylim = ax.get_ylim()
    yr = ylim[1] - ylim[0]
    ax.set_ylim(ylim[0], ylim[1] + 0.2 * yr)

    handles1 = [mlines.Line2D((), (), color=color, label=clabel)
                for color, clabel in zip(colors, colorlabels)]
    handles2 = [mlines.Line2D((), (), color='black', linestyle=ls,
                              label=(f'${lms[0]:.1f} \\pm {lms[1]:.1f}$'))
                for ls, lms in zip(ls_examples, logmstar_examples)]
    leg1 = ax.legend(handles=handles1, fontsize=fontsize - 2, ncol=2,
              loc='upper center', bbox_to_anchor=(0.5, 1.0),
              handlelength=1., columnspacing=0.7, handletextpad=0.4)
    leg2 = ax.legend(handles=handles2, fontsize=fontsize - 1, ncol=1,
              loc='upper right', bbox_to_anchor=(1.0, 0.8),
              handlelength=2., columnspacing=1., handletextpad=0.4,
              title=('$\\log_{10} \\, \\mathrm{M}_{\\star} \\;'
                    ' [\\mathrm{M}_{\\odot}]$'), 
              title_fontsize=fontsize - 1)
    ax.add_artist(leg1)
    outname = f'mhalo_from_example_mstar_z{z_used:.2f}'
    outname = imgdir + outname.replace('.', 'p') + '.pdf'
    plt.savefig(outname, bbox_inches='tight')


