import h5py
import matplotlib as mpl
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import fire_an.makeplots.plot_utils as pu
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.math_utils as mu

def get2dmap_r_vr(filen, minT=None):
    with h5py.File(filen, 'r') as f:
        hist = f['histogram/histogram'][:]
        islog = bool(f['histogram'].attrs['log'])
        if islog:
            hist = 10**hist
        if minT is not None:
            binsT =  f['axis_2/bins'][:]
            if bool(f['axis_2'].attrs['log']):
                binsT = 10**binsT
            iTmin = np.where(np.isclose(minT, binsT))[0][0]
        else:
            iTmin = 0
        hist = np.sum(hist[:, :, iTmin:], axis=2)
        rbins_rvir = f['axis_0/bins'][:]
        vradbins_cmps = f['axis_1/bins'][:]
        cosmopars = {key: val for key, val 
                     in f['Header/cosmopars'].attrs.items()}
        mvir_g = f['Header/inputpars/halodata'].attrs['Mvir_g']
        rvir_cm = f['Header/inputpars/halodata'].attrs['Rvir_cm']
    vradbins_kmps = vradbins_cmps * 1e-5
    normedhist = hist / np.sum(hist) / np.diff(rbins_rvir)[:, np.newaxis] \
                 / np.diff(vradbins_kmps)[np.newaxis, :]
    lognormedhist = np.log10(normedhist)
    out = {'rbins_rvir': rbins_rvir,
           'vradbins_kmps': vradbins_kmps,
           'cosmopars': cosmopars,
           'mvir_g': mvir_g,
           'rvir_cm': rvir_cm,
           'lognormedhist': lognormedhist,
           'hist': hist}
    return out

def plot_r_vr_weights(filen_temp, weightfills, weightlabels=None,
                      title=None, outname=None, minT=None):
    _cmap = mpl.cm.get_cmap('gist_yarg')
    cmap = pu.truncate_colormap(_cmap, minval=0., maxval=0.7)
    perccolor = 'orange'
    percvals = np.array([0.1, 0.5, 0.9])
    if minT is None:
        minT = [None] * len(weightfills)
    data = [get2dmap_r_vr(filen_temp.format(**fill), minT=_minT) 
            for _minT, fill in zip(minT, weightfills)]
    vmax = max([np.max(datum['lognormedhist']) for datum in data])
    vmin = min([np.min(datum['lognormedhist']) for datum in data])
    vmin = max(vmin, vmax - 5.)
    extend = 'min'
    clabel = ('$\\log_{10} \\; \\partial^2 \\mathrm{weight} \\,/\\,'
              '\\Sigma \\, \\mathrm{weight}'
              '\\,/\\, \\partial (r \\; [\\mathrm{R}_{\\mathrm{vir}}])'
              '\\,/\\, \\partial (v_{\\mathrm{r}}'
              '\\; [\\mathrm{km}\\, \\mathrm{s}^{-1}])$')
    xlabel = '$r \\; [\\mathrm{R}_{\\mathrm{vir}}]$'
    ylabel = '$v_{\\mathrm{r}}\\; [\\mathrm{km}\\, \\mathrm{s}^{-1}]$'
    fontsize = 12

    nw = len(weightfills)
    ncmax = 4
    panelsize = 2.5
    ncols = min(ncmax, nw)
    nrows = (nw - 1) // ncols + 1
    caxspace = 0.5
    width_ratios = [panelsize] * ncols + [caxspace]
    height_ratios = [panelsize] * nrows
    cbar_orientation = 'vertical'
    _ncols = ncols + 1
    width = sum(width_ratios)
    height = sum(height_ratios)
    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=_ncols, nrows=nrows,
                        height_ratios=height_ratios,
                        width_ratios=width_ratios,
                        hspace=0., wspace=0.)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols]) for i in range(nw)]
    cax = fig.add_subplot(grid[:, -1])

    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    # track y limits. excl. central galaxy
    ymin = np.inf
    ymax = -np.inf
    xr_ycheck_rvir = [0.2, 1.1]
    for i, (ax, datum) in enumerate(zip(axes, data)):
        below = i >= nw - ncols
        left = i % ncols == 0
        ax.tick_params(which='both', direction='in', top=True,
                       right=True, labelbottom=below,
                       labelleft=left, labelsize=fontsize - 1.)
        if left:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        if below:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        img = ax.pcolormesh(datum['rbins_rvir'], datum['vradbins_kmps'],
                            datum['lognormedhist'].T, cmap=cmap,
                            vmin=vmin, vmax=vmax, rasterized=True)
        ax.axvline(1., linestyle='dotted', color='black', ymax=0.8)
        vescvir = np.sqrt(c.gravity * datum['mvir_g'] / datum['rvir_cm']) \
                  * 1e-5
        ax.axhline(vescvir, linestyle='solid', color='black')
        ax.axhline(0., linestyle='solid', color='black')
        ax.axhline(-vescvir, linestyle='solid', color='black')

        ylo, ymed, yhi = mu.percentiles_from_histogram(
            datum['hist'], datum['vradbins_kmps'], axis=1, 
            percentiles=percvals)
        xc = 0.5 * (datum['rbins_rvir'][:-1] + datum['rbins_rvir'][1:])
        sclabel = f'{percvals[0] * 100:.0f}, {percvals[0-1] * 100:.0f}%'
        ax.plot(xc, ylo, color=perccolor, linestyle='dotted', label=sclabel)
        ax.plot(xc, yhi, color=perccolor, linestyle='dotted')
        ax.plot(xc, ymed, color=perccolor, linestyle='dashed', label='median')
        yc =  0.5 * (datum['vradbins_kmps'][:-1] + datum['vradbins_kmps'][1:])
        mean = np.sum(yc[np.newaxis, :] * datum['hist'], axis=1) \
               /  np.sum(datum['hist'], axis=1)
        ax.plot(xc, mean, color=perccolor, linestyle='solid', label='mean')  
        if i == 0:
            ax.legend(fontsize=fontsize - 2., loc='lower right', ncol=1,
                      handlelength=1., borderpad=0.3, handletextpad=0.6)
        if weightlabels is not None:
            ax.text(0.98, 0.98, weightlabels[i], color='black',
                    fontsize=fontsize, horizontalalignment='right',
                    verticalalignment='top', transform=ax.transAxes)
        xbins = datum['rbins_rvir']
        ybins = datum['vradbins_kmps']
        xilo = np.where(np.isclose(xr_ycheck_rvir[0], xbins))[0][0]
        xihi = np.where(np.isclose(xr_ycheck_rvir[-1], xbins))[0][0]
        checkpoints = np.sum(datum['hist'][xilo : xihi, :], axis=0)
        ysel = np.where(checkpoints > 0.)[0]
        _ymin = ybins[ysel[0]]
        _ymax = ybins[ysel[-1] + 1]
        ymin = min(ymin, _ymin)
        ymax = max(ymax, _ymax)
            
    plt.colorbar(img, cax=cax, extend=extend, orientation=cbar_orientation)
    cax.set_ylabel(clabel, fontsize=fontsize)
    
    xlims = [ax.get_xlim() for ax in axes]
    xmin = min([xlim[0] for xlim in xlims])
    xmax = max([xlim[1] for xlim in xlims])
    [ax.set_xlim((xmin, xmax)) for ax in axes]
    ymar = 0.05 * (ymax - ymin)    
    [ax.set_ylim((ymin - ymar, ymax + ymar)) for ax in axes]

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotset_r_vr_weights(hset='core'):
    if hset == 'core':
        ddir = '/projects/b1026/nastasha/hists/r_vr_clean2_nobug/'
        filetemp = ('hist_rcen_vcen_temperature_by_{{weight}}_{simname}'
                    '_snap{snapnum}_bins1_v1_hvcen.hdf5')
        weights = ['gasmass', 'gasvol', 'Oxygen', 'Neon',
                   'O6', 'Ne8', 'O7', 'Ne9', 'O8', 'Ne10']
        weightlabels = ['Mass', 'Volume', 'Oxygen > 1e5 K', 'Neon > 1e5 K',
                        'O VI', 'Ne VIII', 'O VII', 'Ne IX', 'O VIII', 'Ne X']
        minT = [None, None, 1e5, 1e5] + [None] * 6
        simnames = sl.m13_hr_clean2 + sl.m13_sr_clean2 +\
                   [('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp2e-4_gacc31_fa0.5'),
                    ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp2e-4_gacc31_fa0.5'),
                    ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5')]
        snaps_sr = [sl.snaps_sr[0], sl.snaps_sr[-1]]
        snaps_hr = [sl.snaps_hr[0], sl.snaps_hr[-1]]
        zs = [1.0, 0.5]
        zstrs = ['1p0', '0p5']
        sims_sr = sl.m13_sr_clean2 + sl.m12_sr_clean2
        sims_hr = sl.m13_hr_clean2 + sl.m12_hr_clean2

        outdir = '/projects/b1026/nastasha/imgs/r_vr_hists/'
        outname = 'weightdist_r_vr_ZTgeq5_{ic}_{phys}_{zstr}.pdf'
    weightfills = [{'weight': weight} for weight in weights]
    
    for simname in simnames:
        snaps = snaps_sr if simname in sims_sr \
                else snaps_hr if simname in sims_hr \
                else None 
        ic = simname.split('_')[0]
        phys = ('noBH' if '_sdp1e10_' in simname 
                else 'AGN-CR' if '_MHDCRspec1_' in simname 
                else 'AGN-noCR')
        for zi, snapnum in enumerate(snaps):
            filen_temp = ddir + filetemp.format(simname=simname, 
                                                snapnum=snapnum)
            zstr = zstrs[zi]
            zval = zs[zi]
            
            title = f'{ic}, {phys}, z={zval:.1f}'
            _outname = outdir + outname.format(ic=ic, phys=phys, zstr=zstr)

            plot_r_vr_weights(filen_temp, weightfills, 
                              weightlabels=weightlabels,
                              title=title, outname=_outname,
                              minT=minT)


def getmedmap_r_vr(filen):
    with h5py.File(filen, 'r') as f:
        hist = f['histogram/histogram'][:]
        islog = bool(f['histogram'].attrs['log'])
        if islog:
            hist = 10**hist
        bins2 =  f['axis_2/bins'][:]
        if bool(f['axis_2'].attrs['log']):
            bins2 = 10**bins2
        medianmap = mu.percentiles_from_histogram(hist, bins2, axis=2,
                                                  percentiles=np.array([0.5]))
        hist = np.sum(hist[:, :], axis=2)
        rbins_rvir = f['axis_0/bins'][:]
        vradbins_cmps = f['axis_1/bins'][:]
        cosmopars = {key: val for key, val 
                     in f['Header/cosmopars'].attrs.items()}
        mvir_g = f['Header/inputpars/halodata'].attrs['Mvir_g']
        rvir_cm = f['Header/inputpars/halodata'].attrs['Rvir_cm']
    vradbins_kmps = vradbins_cmps * 1e-5
    logmedianmap = np.log10(medianmap[0])
    out = {'rbins_rvir': rbins_rvir,
           'vradbins_kmps': vradbins_kmps,
           'cosmopars': cosmopars,
           'mvir_g': mvir_g,
           'rvir_cm': rvir_cm,
           'logmedianmap': logmedianmap,
           'hist': hist}
    return out

    
def plot_r_vr_medmap(filen_temp, weightfills, weightlabels=None,
                     title=None, outname=None, clabel=None):
    cmap = mpl.cm.get_cmap('viridis')
    cmap.set_bad('white', 1.)
    #cmap = pu.truncate_colormap(_cmap, minval=0., maxval=0.7)
    contourcolor = 'red'
    data = [getmedmap_r_vr(filen_temp.format(**fill)) 
            for fill in weightfills]
    #vmax = max([np.max(datum['logmedianmap'][
    #                np.isfinite(datum['logmedianmap'])])
    #            for datum in data])
    #vmin = min([np.min(datum['logmedianmap'][
    #                np.isfinite(datum['logmedianmap'])]) 
    #            for datum in data])
    if 'temperature' in filen_temp:
        vmin = 5.0
        vmax = 7.0
    elif 'density' in filen_temp:
        vmin = -31.
        vmax = -15.
    extend = 'both'
    xlabel = '$r \\; [\\mathrm{R}_{\\mathrm{vir}}]$'
    ylabel = '$v_{\\mathrm{r}}\\; [\\mathrm{km}\\, \\mathrm{s}^{-1}]$'
    fontsize = 12

    nw = len(weightfills)
    ncmax = 4
    panelsize = 2.5
    ncols = min(ncmax, nw)
    nrows = (nw - 1) // ncols + 1
    caxspace = 0.5
    width_ratios = [panelsize] * ncols + [caxspace]
    height_ratios = [panelsize] * nrows
    cbar_orientation = 'vertical'
    _ncols = ncols + 1
    width = sum(width_ratios)
    height = sum(height_ratios)
    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=_ncols, nrows=nrows,
                        height_ratios=height_ratios,
                        width_ratios=width_ratios,
                        hspace=0., wspace=0.)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols]) for i in range(nw)]
    cax = fig.add_subplot(grid[:, -1])

    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    # track y limits. excl. central galaxy
    ymin = np.inf
    ymax = -np.inf
    xr_ycheck_rvir = [0.2, 1.1]
    for i, (ax, datum) in enumerate(zip(axes, data)):
        below = i >= nw - ncols
        left = i % ncols == 0
        ax.tick_params(which='both', direction='in', top=True,
                       right=True, labelbottom=below,
                       labelleft=left, labelsize=fontsize - 1.)
        if left:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        if below:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        badsel = datum['hist'] <= 0.
        mapin = np.ma.masked_where(badsel, datum['logmedianmap'])
        img = ax.pcolormesh(datum['rbins_rvir'], 
                            datum['vradbins_kmps'],
                            mapin.T, 
                            cmap=cmap,
                            vmin=vmin, vmax=vmax, rasterized=True)
        ax.axvline(1., linestyle='dotted', color='black', ymax=0.8)
        vescvir = np.sqrt(c.gravity * datum['mvir_g'] / datum['rvir_cm']) \
                  * 1e-5
        ax.axhline(vescvir, linestyle='solid', color='black')
        ax.axhline(0., linestyle='solid', color='black')
        ax.axhline(-vescvir, linestyle='solid', color='black')
        
        edges = [datum['rbins_rvir'], datum['vradbins_kmps']]
        linestyles = ['solid', 'dashed', 'dotted']
        levels = [0.98, 0.9, 0.5]
        pu.add_2dhist_contours(ax, datum['hist'], edges, (0, 1),
                        mins=None, maxs=None, histlegend=False, 
                        fraclevels=True, levels=levels, 
                        legend=False, 
                        dimlabels=None, legendlabel=None,
                        legendlabel_pre=None, shiftx=0., shifty=0., 
                        dimshifts=None, colors=contourcolor,
                        linestyles=linestyles) 
        if i == 0:
            handles = [mlines.Line2D((), (), color=contourcolor,
                                    linestyle=ls, label=f'{lv*100:.1f}%')
                       for ls, lv in zip(linestyles, levels)]
            ax.legend(handles=handles, 
                      fontsize=fontsize - 2., loc='lower right', ncol=1,
                      handlelength=1., borderpad=0.3, handletextpad=0.6)
        if weightlabels is not None:
            ax.text(0.98, 0.98, weightlabels[i], color='black',
                    fontsize=fontsize, horizontalalignment='right',
                    verticalalignment='top', transform=ax.transAxes)
        xbins = datum['rbins_rvir']
        ybins = datum['vradbins_kmps']
        xilo = np.where(np.isclose(xr_ycheck_rvir[0], xbins))[0][0]
        xihi = np.where(np.isclose(xr_ycheck_rvir[-1], xbins))[0][0]
        checkpoints = np.sum(datum['hist'][xilo : xihi, :], axis=0)
        minv = np.min(np.max(datum['hist'], axis=1)) * 1e-3
        ysel = np.where(checkpoints > minv)[0]
        _ymin = ybins[ysel[0]]
        _ymax = ybins[ysel[-1] + 1]
        ymin = min(ymin, _ymin)
        ymax = max(ymax, _ymax)
            
    plt.colorbar(img, cax=cax, orientation=cbar_orientation, extend=extend)
    cax.set_ylabel(clabel, fontsize=fontsize)
    
    xlims = [ax.get_xlim() for ax in axes]
    xmin = min([xlim[0] for xlim in xlims])
    xmax = max([xlim[1] for xlim in xlims])
    [ax.set_xlim((xmin, xmax)) for ax in axes]
    ymar = 0.05 * (ymax - ymin)    
    [ax.set_ylim((ymin - ymar, ymax + ymar)) for ax in axes]

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotset_r_vr_weights(hset='core'):
    if hset == 'core':
        ddir = '/projects/b1026/nastasha/hists/r_vr_clean2_nobug/'
        filetemp = ('hist_rcen_vcen_{axis2}_by_{{weight}}_{simname}'
                    '_snap{snapnum}_bins1_v1_hvcen.hdf5')
        weights = ['gasmass', 'gasvol', 'Oxygen', 'Neon',
                   'O6', 'Ne8', 'O7', 'Ne9', 'O8', 'Ne10']
        weightlabels = ['Mass', 'Volume', 'Oxygen', 'Neon',
                        'O VI', 'Ne VIII', 'O VII', 'Ne IX', 'O VIII', 'Ne X']
        axes2 = ['temperature', 'density']
        clabels = ['median $\\log_{10} \\, \\mathrm{T} \\; [\\mathrm{K}]$',
                   ('median $\\log_{10} \\, \\rho \\; [\\mathrm{g} \\;'
                    ' \\mathrm{cm}^{-3}]$')]
        simnames = sl.m13_hr_clean2 + sl.m13_sr_clean2 +\
                   [('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp2e-4_gacc31_fa0.5'),
                    ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp2e-4_gacc31_fa0.5'),
                    ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5')]
        snaps_sr = [sl.snaps_sr[0], sl.snaps_sr[-1]]
        snaps_hr = [sl.snaps_hr[0], sl.snaps_hr[-1]]
        zs = [1.0, 0.5]
        zstrs = ['1p0', '0p5']
        sims_sr = sl.m13_sr_clean2 + sl.m12_sr_clean2
        sims_hr = sl.m13_hr_clean2 + sl.m12_hr_clean2

        outdir = '/projects/b1026/nastasha/imgs/r_vr_hists/'
        outname = 'propdist_r_vr_{axis2}_{ic}_{phys}_{zstr}.pdf'
    weightfills = [{'weight': weight} for weight in weights]
    
    for simname in simnames:
        snaps = snaps_sr if simname in sims_sr \
                else snaps_hr if simname in sims_hr \
                else None 
        ic = simname.split('_')[0]
        phys = ('noBH' if '_sdp1e10_' in simname 
                else 'AGN-CR' if '_MHDCRspec1_' in simname 
                else 'AGN-noCR')
        for zi, snapnum in enumerate(snaps):
            for axis2, clabel in zip(axes2, clabels):
                filen_temp = ddir + filetemp.format(simname=simname, 
                                                    snapnum=snapnum,
                                                    axis2=axis2)
                zstr = zstrs[zi]
                zval = zs[zi]
                
                title = f'{axis2}, {ic}, {phys}, z={zval:.1f}'
                _outname = outdir + outname.format(ic=ic, phys=phys,
                                                   zstr=zstr,
                                                   axis2=axis2)
                plot_r_vr_medmap(filen_temp, weightfills, 
                                 weightlabels=weightlabels,
                                 title=title, outname=_outname,
                                 clabel=clabel)
        break