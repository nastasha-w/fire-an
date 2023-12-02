'''
Rough comparison of the Burchett et al. (2019) data
to median column density curves from EAGLE (z=0.5)
in my 2020 paper with Joop and Ben.
'''

import h5py
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize as so

import fire_an.makeplots.litcomp.obsdataread as odr
import fire_an.utils.cosmo_utils as cu
import fire_an.utils.constants_and_units as c
import fire_an.utils.math_utils as mu

mdir = '/projects/b1026/nastasha/imgs/datacomp/eagle/'
eagledatadir = '/projects/b1026/nastasha/extdata/eaglepaper2/'
oddir = '/projects/b1026/nastasha/extdata/'
q23filen = oddir + 'plotdata_q23_nsigmas_1_2.dat'
b19filen = oddir + 'plotdata_b19_nsigmas_1_2.dat'


cosmopars_ea_23 = {'a': 0.6652884960735025,
                   'boxsize': 67.77,
                   'h': 0.6777,
                   'omegab': 0.0482519,
                   'omegalambda': 0.693,
                   'omegam': 0.307,
                   'z': 0.5031073074342141,
                   }

def readin_eagledata():
    '''
    R200c units, 
    0 - 2.5 R200c impact parameter, 0.1 R200c bin size
    0.5 dex mass bins (M200c)
    '''
    eagledatafilen = ('rdist_coldens_ne8_L0100N1504_23_test3.4_PtAb_C2Sm'
                      '_32000pix_6.25slice_zcen-all_z-projection_T4EOS'
                      '_1slice_to-100-pkpc-or-3-R200c_M200c-0p5dex-7000'
                      '_centrals_stored_profiles.hdf5')
    percdct = {}
    percentiles = [10., 50., 90.]
    pkeys = {_p: f'perc_{_p:.1f}' for _p in percentiles}
    subpath = '/R200c_bins/binset_0/'
    with h5py.File(eagledatadir + eagledatafilen, 'r') as f:
        galsets = list(f.keys())
        for galset in galsets:
            edges = f[galset + subpath + 'bin_edges'][:]
            pvals = {_p: f[galset + subpath + pkeys[_p]][:]
                     for _p in percentiles}
            mmin = f[galset].attrs['seltag'].decode()
            mmin = mmin.split('_')[0]
            mmin = float(mmin[3:])
            percdct[mmin] = {'edges': edges, 'pvals': pvals}
    return percdct

def mvir_to_m200c(mvir_msun):
    meandens_bn98 = cu.getmeandensity('BN98', cosmopars_ea_23)
    #meandens_200c = cu.getmeandensity('200c', cosmopars_ea_23)
    mvir_g = mvir_msun * c.solar_mass
    rbn98_cgs = (mvir_g / meandens_bn98 * 3. / (4. * np.pi))**(1./3.)
    
    rbins = np.linspace(0., 5. * rbn98_cgs, 5000.)
    rcens = 0.5 * (rbins[:-1] + rbins[1:])
    dr = np.average(np.diff(rbins))
    def minfunc(M200c):
        rhovals = cu.rho_NFW(rcens, M200c, delta=200, ref='rhocrit', 
                             z=cosmopars_ea_23['z'], 
                             cosmopars=cosmopars_ea_23, c='Schaller15')
        menc = np.cumsum(4. * np.pi * rhovals * rcens**2 * dr)
        rhomean = menc / (4. * np.pi / 3. * rbins[1:]**3)
        rhomean_rbn98 = mu.linterpsolve(rbins[1:], rhomean, rbn98_cgs)
        cost = np.abs(np.log10(rhomean_rbn98) - np.log10(meandens_bn98))
        return cost
    
    res = so.minimize(minfunc, mvir_g)
    if not res.succes:
        print(res)
        return res
    else:
        return res.x

def m200c_to_mvir(m200c_msun):
    meandens_bn98 = cu.getmeandensity('BN98', cosmopars_ea_23)
    meandens_200c = cu.getmeandensity('200c', cosmopars_ea_23)
    m200c_g = m200c_msun * c.solar_mass
    r200c_cgs = (m200c_g / meandens_200c * 3. / (4. * np.pi))**(1./3.)
    
    rbins = np.linspace(0., 5. * r200c_cgs, 5000)
    rcens = 0.5 * (rbins[:-1] + rbins[1:])
    dr = np.average(np.diff(rbins))
    rhovals = cu.rho_NFW(rcens, m200c_g, delta=200, ref='rhocrit', 
                         z=cosmopars_ea_23['z'],
                         cosmopars=cosmopars_ea_23, c='Schaller15')
    menc = np.cumsum(4. * np.pi  * rhovals * rcens**2 * dr)
    rhomean = menc / (4. * np.pi / 3. * rbins[1:]**3)
    rbn98_cgs = mu.linterpsolve(rhomean, rbins[1:], meandens_bn98)
    mbn98_g = 4. * np.pi / 3. * meandens_bn98 * rbn98_cgs**3
    mbn98_msun = mbn98_g / c.solar_mass
    return mbn98_msun 



## initial copy was from the paper 2 scripts
def plot_radprof_eagle_b19_comp():
    '''
    Very rough comparison! Do not put this in a paper. 
    '''
    fontsize = 12
    ylim = (12.0, 15.5)
    # for labeling, not passed to anything
    ypercs = [10., 50., 90.]
    nsigmas = [1, 2]

    imgname = 'ne8_b19_cols_vs_eaglez0p5_w20_roughcomp_v1.pdf'
    imgname = mdir + imgname

    kwa_pfills = {'color': (0.8, 0.8, 0.8)}
    kwa_med = {'color': 'black', 'linestyle': 'solid', 'linewidth': 2.}
    
    ylabel = ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\,VIII}) '
              '\\; [\\mathrm{cm}^{-2}]$')
    xlabel = '$r_{\\perp} \\; [\\mathrm{R}_{\\mathrm{200c}}]$'
    #clabel =( '$\\log_{10}\, \\mathrm{M}_{\\mathrm{200c}} '
    #         '\\; [\\mathrm{M}_{\\odot}]$')
    figtitle = ('rough match EAGLE data (Wijers et al. 2020)'
                ' to Burchett et al. (2019) obs.\n'
                'EAGLE masses: M200c categories, est. Mvir used'
                ' to compare to likely B+19 halo masses\n'
                'horiz. error bars are range of R200c fracs.'
                ' for each M200c range, given meas. pkpc')
    
    eagledat = readin_eagledata()
    data_bur = pd.read_csv(b19filen, sep='\t')
    # define used mass ranges
    deltaM200c = 0.5
    massbins_m200c_eagle = list(eagledat.keys())
    massbins_m200c_eagle.sort()
    massbins_m200c_eagle = massbins_m200c_eagle \
                           + [massbins_m200c_eagle[-1] + deltaM200c]
    print(massbins_m200c_eagle)
    massbins_mbn98_eagle = [np.log10(m200c_to_mvir(10**me)) 
                            for me in massbins_m200c_eagle]
    print(massbins_mbn98_eagle)
    npanels = len(massbins_m200c_eagle) - 1
    ncols = min(npanels, 4)
    nrows = (npanels - 1) // ncols + 1
    panelsize = 2.5
    height_ratios = [panelsize] * nrows
    width_ratios = [panelsize] * ncols
    width = sum(width_ratios)
    height = sum(height_ratios)

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(nrows=nrows, ncols=ncols,
                        hspace=0., wspace=0.,
                        width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols])
            for i in range(npanels)]
    fig.suptitle(figtitle, fontsize=fontsize)
    
    getlegax = None
    for axi in range(npanels):
        ax = axes[axi]
        dobottom = axi >= npanels - ncols
        doleft = axi % ncols == 0
        ax.tick_params(direction='in', which='both',
                       right=True, top=True, labelsize=fontsize - 1,
                       labelbottom=dobottom, labelleft=doleft)
        if dobottom:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        if doleft:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        
        mmin200c = massbins_m200c_eagle[axi]
        mmax200c = massbins_m200c_eagle[axi + 1]
        mminbn98 = massbins_mbn98_eagle[axi]
        mmaxbn98 = massbins_mbn98_eagle[axi + 1]
        r200cmin_cm = cu.Rhalo(10**mmin200c * c.solar_mass, 
                                 delta=200, ref='rhocrit', 
                                 z=cosmopars_ea_23['z'],
                                 cosmopars=cosmopars_ea_23)
        r200cmax_cm = cu.Rhalo(10**mmax200c * c.solar_mass, 
                                 delta=200, ref='rhocrit', 
                                 z=cosmopars_ea_23['z'],
                                 cosmopars=cosmopars_ea_23)
        
        axlabel = (f'$\\mathrm{{M}}_{{\\mathrm{{200c}}}}: {mmin200c:.1f}'
                   f'\\endash {mmax200c:.1f}$')
        ax.text(0.98, 0.98, axlabel, fontsize=fontsize,
                horizontalalignment='right', verticalalignment='top',
                transform=ax.transAxes)

        percdct = eagledat[mmin200c]
        ed = percdct['edges']
        cens = 0.5 * (ed[:-1] + ed[1:])
        med = percdct['pvals'][50.]
        plo = percdct['pvals'][10.]
        phi = percdct['pvals'][90.]
        
        ax.fill_between(cens, plo, phi, **kwa_pfills)
        ax.plot(cens, med, **kwa_med)
        
        detmsig1done = False
        detmsig0done = False
        detfsig1done = False
        detfsig0done = False
        ulsig1done = False
        ulsig0done = False
        for dbi in range(len(data_bur)):
            cloer = data_bur['logmvir_msun_loer'][dbi]
            chier = data_bur['logmvir_msun_hier'][dbi]
            #print(mminbn98, mmaxbn98)
            #print(cloer, chier)
            if cloer > mmaxbn98 or chier < mminbn98:
                continue

            _xv = data_bur['impact_parameter_kpc'][dbi]
            xlo = _xv * 1e-3 * c.cm_per_mpc / r200cmax_cm 
            xhi = _xv * 1e-3 * c.cm_per_mpc / r200cmin_cm 
            xmid = 0.5 * (xlo + xhi)
            xerr = xhi - xmid
            yv = data_bur['log_N_Ne8_pcm2'][dbi]
            isul = data_bur['log_N_Ne8_isUL'][dbi]
            yerr = data_bur['log_N_Ne8_pcm2_err'][dbi] if not isul else None
            #cbest = data_bur['logmvir_msun_bestest'][dbi]
            clo = data_bur['logmvir_msun_lo'][dbi]
            chi = data_bur['logmvir_msun_hi'][dbi]
            flag = data_bur['flagged_by_qu23'][dbi]
            #print(clo, chi)
            
            issig0 = not (clo > mmaxbn98 or chi < mminbn98)
            _label = None
            if issig0:
                _color = 'black'
                if isul and not ulsig0done:
                    _label = ('UL, $\\Delta\\mathrm{M}'
                              f' < {nsigmas[0]}\\sigma$')
                    ulsig0done = True
                elif (not isul and not flag) and not detmsig0done:
                    _label = ('det., $\\Delta\\mathrm{M}'
                              f' < {nsigmas[0]}\\sigma$')
                    detmsig0done = True
                elif (not isul and flag) and not detfsig0done:
                    _label = ('det., !, $\\Delta\\mathrm{M}'
                              f' < {nsigmas[0]}\\sigma$')
                    detfsig0done = True
            else:
                _color = 'gray'
                if isul and not ulsig1done:
                    _label = ('UL, $\\Delta\\mathrm{M}'
                              f' < {nsigmas[1]}\\sigma$')
                    ulsig1done = True
                elif (not isul and not flag) and not detmsig1done:
                    _label = ('det., $\\Delta\\mathrm{M}'
                              f' < {nsigmas[1]}\\sigma$')
                    detmsig1done = True
                elif (not isul and flag) and not detfsig1done:
                    _label = ('det., !, $\\Delta\\mathrm{M}'
                              f' < {nsigmas[1]}\\sigma$')
                    detfsig1done = True
            marker = 'v' if isul else 'o'
            markeredgecolor = _color if flag else 'black'
            markerfacecolor = 'none' if flag else _color
            markersize = 5
            zobase = 5. - 1. * isul

            ax.errorbar([xmid], [yv], yerr=yerr, xerr=xerr,
                        linestyle='none', elinewidth=1.5,
                        color=_color, capsize=3,
                        zorder=zobase,
                        marker=marker, markersize=markersize,
                        markeredgecolor=markeredgecolor, 
                        markerfacecolor=markerfacecolor,
                        markeredgewidth=1.0,
                        label=_label)
        if detmsig1done and detmsig0done and ulsig1done and ulsig0done \
                and detfsig1done and detfsig0done:
            getlegax = axi
    
    [ax.set_ylim(ylim) for ax in axes]
    # legend add
    handles, labels = axes[getlegax].get_legend_handles_labels()
    axes[-1].legend(handles=handles, fontsize=fontsize - 1,
                    handlelength=1., labelspacing=0.15, handletextpad=0.4,
                    loc='upper right', bbox_to_anchor=(1.0, 0.90))
    hndl1 = [mpatch.Patch(label=f'${ypercs[0]:.0f}\\endash{ypercs[-1]:.0f}$%',
                          **kwa_pfills),
             mlines.Line2D((), (), label='median', **kwa_med)]
    axes[-2].legend(handles=hndl1, fontsize=fontsize - 1,
                    handlelength=1., labelspacing=0.15, handletextpad=0.4,
                    loc='upper right', bbox_to_anchor=(1.0, 0.90))

    plt.savefig(imgname, format='pdf', bbox_inches='tight')


## initial copy was from the paper 2 scripts
def plot_radprof_eagle_q23_comp():
    '''
    Very rough comparison! Do not put this in a paper. 
    '''
    fontsize = 12
    ylim = (12.0, 15.5)
    # for labeling, not passed to anything
    ypercs = [10., 50., 90.]
    nsigmas = [1, 2]

    imgname = 'ne8_q23_cols_vs_eaglez0p5_w20_roughcomp_v1.pdf'
    imgname = mdir + imgname

    kwa_pfills = {'color': (0.8, 0.8, 0.8)}
    kwa_med = {'color': 'black', 'linestyle': 'solid', 'linewidth': 2.}
    
    ylabel = ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\,VIII}) '
              '\\; [\\mathrm{cm}^{-2}]$')
    xlabel = '$r_{\\perp} \\; [\\mathrm{R}_{\\mathrm{200c}}]$'
    #clabel =( '$\\log_{10}\, \\mathrm{M}_{\\mathrm{200c}} '
    #         '\\; [\\mathrm{M}_{\\odot}]$')
    figtitle = ('rough match EAGLE data (Wijers et al. 2020)'
                ' to Qu et al. (2023, in prep.) obs.\n'
                'EAGLE masses: M200c categories, est. Mvir used'
                ' to compare to likely Q+23 halo masses\n'
                'horiz. error bars are range of R200c fracs.'
                ' for each M200c range, given meas. pkpc')
    
    eagledat = readin_eagledata()
    data_cubs = pd.read_csv(q23filen, sep='\t')
    # define used mass ranges
    deltaM200c = 0.5
    massbins_m200c_eagle = list(eagledat.keys())
    massbins_m200c_eagle.sort()
    massbins_m200c_eagle = massbins_m200c_eagle \
                           + [massbins_m200c_eagle[-1] + deltaM200c]
    print(massbins_m200c_eagle)
    massbins_mbn98_eagle = [np.log10(m200c_to_mvir(10**me)) 
                            for me in massbins_m200c_eagle]
    print(massbins_mbn98_eagle)
    npanels = len(massbins_m200c_eagle) - 1
    ncols = min(npanels, 4)
    nrows = (npanels - 1) // ncols + 1
    panelsize = 2.5
    height_ratios = [panelsize] * nrows
    width_ratios = [panelsize] * ncols
    width = sum(width_ratios)
    height = sum(height_ratios)

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(nrows=nrows, ncols=ncols,
                        hspace=0., wspace=0.,
                        width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols])
            for i in range(npanels)]
    fig.suptitle(figtitle, fontsize=fontsize)
    
    getlegax = None
    print('data length CUBS: ', len(data_cubs))
    for axi in range(npanels):
        ax = axes[axi]
        dobottom = axi >= npanels - ncols
        doleft = axi % ncols == 0
        ax.tick_params(direction='in', which='both',
                       right=True, top=True, labelsize=fontsize - 1,
                       labelbottom=dobottom, labelleft=doleft)
        if dobottom:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        if doleft:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        
        mmin200c = massbins_m200c_eagle[axi]
        mmax200c = massbins_m200c_eagle[axi + 1]
        mminbn98 = massbins_mbn98_eagle[axi]
        mmaxbn98 = massbins_mbn98_eagle[axi + 1]
        r200cmin_cm = cu.Rhalo(10**mmin200c * c.solar_mass, 
                                 delta=200, ref='rhocrit', 
                                 z=cosmopars_ea_23['z'],
                                 cosmopars=cosmopars_ea_23)
        r200cmax_cm = cu.Rhalo(10**mmax200c * c.solar_mass, 
                                 delta=200, ref='rhocrit', 
                                 z=cosmopars_ea_23['z'],
                                 cosmopars=cosmopars_ea_23)
        
        axlabel = (f'$\\mathrm{{M}}_{{\\mathrm{{200c}}}}: {mmin200c:.1f}'
                   f'\\endash {mmax200c:.1f}$')
        ax.text(0.98, 0.98, axlabel, fontsize=fontsize,
                horizontalalignment='right', verticalalignment='top',
                transform=ax.transAxes)

        percdct = eagledat[mmin200c]
        ed = percdct['edges']
        cens = 0.5 * (ed[:-1] + ed[1:])
        med = percdct['pvals'][50.]
        plo = percdct['pvals'][10.]
        phi = percdct['pvals'][90.]
        
        ax.fill_between(cens, plo, phi, **kwa_pfills)
        ax.plot(cens, med, **kwa_med)
        
        detsig1done = False
        detsig0done = False
        ulsig1done = False
        ulsig0done = False
        for dbi in range(len(data_cubs['ne8col_logcm2'])):
            cloer = data_cubs['logmvir_msun_loer'][dbi]
            chier = data_cubs['logmvir_msun_hier'][dbi]
            #print(mminbn98, mmaxbn98)
            #print(cloer, chier)
            if cloer > mmaxbn98 or chier < mminbn98:
                continue

            _xv = data_cubs['impactpar_kpc'][dbi]
            xlo = _xv * 1e-3 * c.cm_per_mpc / r200cmax_cm 
            xhi = _xv * 1e-3 * c.cm_per_mpc / r200cmin_cm 
            xmid = 0.5 * (xlo + xhi)
            xerr = xhi - xmid
            yv = data_cubs['ne8col_logcm2'][dbi]
            isul = data_cubs['isul_ne8'][dbi]
            yerr = ([data_cubs['ne8col_2s_loerr_dex'][dbi]], 
                    [data_cubs['ne8col_2s_hierr_dex'][dbi]]) \
                    if not isul else None
            #cbest = data_bur['logmvir_msun_bestest'][dbi]
            clo = data_cubs['logmvir_msun_lo'][dbi]
            chi = data_cubs['logmvir_msun_hi'][dbi]
            #print(clo, chi)
            
            issig0 = not (clo > mmaxbn98 or chi < mminbn98)
            _label = None
            if issig0:
                _color = 'black'
                if isul and not ulsig0done:
                    _label = ('UL, $\\Delta\\mathrm{M}'
                              f' < {nsigmas[0]}\\sigma$')
                    ulsig0done = True
                elif not isul and not detsig0done:
                    _label = ('det., $\\Delta\\mathrm{M}'
                              f' < {nsigmas[0]}\\sigma$')
                    detsig0done = True
            else:
                _color = 'gray'
                if isul and not ulsig1done:
                    _label = ('UL, $\\Delta\\mathrm{M}'
                              f' < {nsigmas[1]}\\sigma$')
                    ulsig1done = True
                elif not isul and not detsig1done:
                    _label = ('det., $\\Delta\\mathrm{M}'
                              f' < {nsigmas[1]}\\sigma$')
                    detsig1done = True
            marker = 'v' if isul else 'o'
            markersize = 5
            zobase = 5. - 1. * isul
            ax.errorbar([xmid], [yv], yerr=yerr, xerr=xerr,
                        linestyle='none', elinewidth=1.5,
                        color=_color, capsize=3,
                        zorder=zobase,
                        marker=marker, markersize=markersize,
                        markeredgecolor='black', markeredgewidth=1.0,
                        label=_label, alpha=0.5)
        if detsig1done and detsig0done and ulsig1done and ulsig0done:
            getlegax = axi
    
    [ax.set_ylim(ylim) for ax in axes]
    # legend add
    handles, labels = axes[getlegax].get_legend_handles_labels()
    axes[-1].legend(handles=handles, fontsize=fontsize - 1,
                    handlelength=1., labelspacing=0.15, handletextpad=0.4,
                    loc='upper right', bbox_to_anchor=(1.0, 0.90))
    hndl1 = [mpatch.Patch(label=f'${ypercs[0]:.0f}\\endash{ypercs[-1]:.0f}$%',
                          **kwa_pfills),
             mlines.Line2D((), (), label='median', **kwa_med)]
    axes[-2].legend(handles=hndl1, fontsize=fontsize - 1,
                    handlelength=1., labelspacing=0.15, handletextpad=0.4,
                    loc='upper right', bbox_to_anchor=(1.0, 0.90))

    plt.savefig(imgname, format='pdf', bbox_inches='tight')

