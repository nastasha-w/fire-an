import h5py
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsp
import numpy as np

import fire_an.mainfunc.makehist as mh

def tryout_hist(index):
    outdir = '/projects/b1026/nastasha/tests/start_fire/hist_tests/'
    if index == 0:
        dirpath = '/projects/b1026/snapshots/fire3/m13h206_m3e5/' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000/'
        snapnum = 27
        simname = 'm13h206_m3e5__' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'
        outfilen = 'hist_Oxygen_by_Mass_0-1-2Rvir_{sc}_snap{sn}_shrink-sph-cen_BN98' + \
                   '_2rvir_v1.hdf5'

        axtypes = ['sim-direct']
        axtypes_args = [{'field': 'ElementAbundance/Oxygen'}]
        weighttype = 'Mass'
        weighttype_args = {}
        rbins = np.array([0., 1., 2.])
        runit = 'Rvir'

        outfilen = outfilen.format(sc=simname, sn=snapnum,)
    else:
        raise ValueError('invalid index: {}'.format(index))

    mh.histogram_radprof(dirpath, snapnum,
                         weighttype, weighttype_args, axtypes, axtypes_args,
                         particle_type=0, 
                         center='shrinksph', rbins=rbins, runit=runit,
                         logweights=True, logaxes=True, axbins=0.05,
                         outfilen=outdir + outfilen)
    
def tryout_coords_hist(ind):
    outdir = '/projects/b1026/nastasha/tests/hists_coords/'
    simname = ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
               '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000')
    simpath = 'm12f_m6e4/' + simname
    snapnum = 45

    weighttypes = ['Mass', 'ion', 'Metal', 'Volume']
    weighttype_argss = [dict(),
                        {'ion': 'Ne8',
                         'ps20depletion': False,
                         'density': False},
                        {'element': 'neon',
                         'ps20depletion': False,
                         'density': False},
                        dict(),
                        ]
    axtypess = ['coords', 'coords']
    axtypes_argss = [{'vel': 'vrad'},
                     {'vel': 'vtot'},
                     ]
    wtind = ind % len(weighttypes)
    atind = ind // len(weighttypes)
    weighttype = weighttypes[wtind]
    weighttype_args = weighttype_argss[wtind]
    axtypes = [axtypess[atind]]
    axtypes_args = [axtypes_argss[atind]]
    rbins = np.append(np.linspace(0., 0.095, 20), np.linspace(0.1, 2., 39))
    
    atstr = 'vel' + axtypes_args[0]['vel'] \
            if isinstance(axtypes_args[0]['vel'], int) else \
            axtypes_args[0]['vel']
    wtstr = 'gasmass' if weighttype == 'Mass' else\
            'gasvol' if weighttype == 'Volume' else\
            weighttype_args['ion'] if weighttype == 'ion' else \
            weighttype_args['element'] 
    outfilen = outdir +\
               f'hist_{atstr}_by_{wtstr}_{simname}_snap{snapnum}.hdf5'
    
    mh.histogram_radprof(simpath, snapnum,
                         weighttype, weighttype_args, axtypes, axtypes_args,
                         particle_type=0, 
                         center='shrinksph', rbins=rbins, runit='Rvir',
                         logweights=True, logaxes=False, axbins=[5e5],
                         outfilen=outfilen, overwrite=True)

def run_coords_hist():
    for ind in range(8):
        print('\n\n')
        tryout_coords_hist(ind)

def getdata_testhist(filen, yunit=1., xunit=1.):
    with h5py.File(filen, 'r') as f:
        hist = f['histogram/histogram'][:]
        islog = bool(f['histogram'].attrs['log'])
        if islog:
            hist = 10**hist
        rbins = f['axis_0/bins'][:] / xunit
        ybins = f['axis_1/bins'][:] / yunit   
    basehist = hist / np.diff(rbins)[:, np.newaxis] \
                    / np.diff(ybins)[np.newaxis, :] \
                    / np.sum(hist)
    # for contour levels
    norm1hist = hist / np.sum(hist, axis=1)[:, np.newaxis]
    addorder = np.argsort(basehist, axis=1)
    backorder = np.argsort(addorder, axis=1)
    csum = np.cumsum(np.take_along_axis(norm1hist, addorder, 1), axis=1)
    fracvals_y = np.take_along_axis(csum, backorder, 1)
    return basehist, fracvals_y, rbins, ybins


def plot_coords_hist():
    #ddir = '/projects/b1026/nastasha/tests/hists_coords/'
    ddir = '/Users/nastasha/ciera/tests/hists_coords/'
    simname = ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
               '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000')
    snapnum = 45
    
    dfilen = ddir +\
             f'hist_{{atstr}}_by_{{wtstr}}_{simname}_snap{snapnum}.hdf5'
    wtstrs = ['gasmass', 'gasvol', 'neon', 'Ne8']
    atstrs = ['vrad', 'vtot']
    vunit = 1e5
    fontsize = 12
    wttitles = {'gasmass': 'Gas mass-weighted', 
                'gasvol': 'Gas volume-weighted', 
                'neon': 'Neon-weighted', 
                'Ne8': 'Ne VIII-weighted'}
    attitles = {'vrad': ('$v_{\\mathrm{rad}} \\;'
                         '[\\mathrm{km} \\, \\mathrm{s}^{-1}]$'),
                'vtot': ('$v_{\\mathrm{tot}} \\;'
                         '[\\mathrm{km} \\, \\mathrm{s}^{-1}]$')}
    
    cmap = 'gist_yarg'
    lvls_contours = [0.04, 0.2]
    colors = 'fuchsia'
    linestyles = ['dashed', 'solid']
    
    lvlstr = ', '.join([f'{(1. - lvl) * 100:.1f}' for lvl in lvls_contours])
    title = ('histograms with coords test plot, m12f, AGN-CR, z=1.0\n'
             f'contour levels: {lvlstr}% of weight at fixed radius')
    outname = ddir + 'hists_coords_plot_sanity_check_m12f_AGN-CR_snap45.pdf'
    
    fig = plt.figure(figsize=(11., 5.))
    width_ratios = [2.5] * len(wtstrs) + [0.5]
    grid = gsp.GridSpec(ncols=len(wtstrs) + 1, nrows=len(atstrs),
                        wspace=0., hspace=0.,
                        width_ratios=width_ratios, top=0.85)
    
    cax = fig.add_subplot(grid[:, len(wtstrs)])
    clabel = ('$ \\log_{10} \\;'
              '\\partial^{2}\\, \\mathrm{weight} \\, / \\, '
              '\\left(\\Sigma(\\mathrm{weight}) \\;'
              ' \\partial \\, \\mathrm{r}'
              '\\; \\partial \\, v \\right) \\;'
              '[\\mathrm{R}_{\\mathrm{vir}}^{-1}'
              '\\, \\mathrm{km}^{-1}\\mathrm{s}]$')
    fig.suptitle(title, fontsize=fontsize)
    
    data = {(atstr, wtstr): 
            getdata_testhist(dfilen.format(atstr=atstr, wtstr=wtstr), 
                             yunit=vunit)
            for atstr in atstrs for wtstr in wtstrs}
    vmax = max([np.max(data[key][0]) for key in data])
    vmin = min([np.min(data[key][0]) for key in data])
    vmin = max(vmin, 1e-5 * vmax)
    vmin = np.log10(vmin)
    vmax = np.log10(vmax)

    axes = []
    for ati, atstr in enumerate(atstrs):
        for wti, wtstr in enumerate(wtstrs):
            
            hist, chist, rbins, ybins = data[(atstr, wtstr)]
            ax = fig.add_subplot(grid[ati, wti])
            axes.append(ax)
            dobottom = ati == len(atstrs) - 1
            doleft = wti == 0
            dotitle = ati == 0
            if dotitle:
                ax.set_title(wttitles[wtstr], fontsize=fontsize)
            if dobottom:
                ax.set_xlabel('$\\mathrm{r}_{\\mathrm{3D}} \\;'
                              ' [\\mathrm{R}_{\\mathrm{vir}}]$', 
                              fontsize=fontsize)
            if doleft:
                ax.set_ylabel(attitles[atstr], fontsize=fontsize)
            ax.tick_params(labelsize=fontsize - 1., which='both',
                           top=True, right=True, labelbottom=dobottom,
                           labelleft=doleft, direction='in')
            img = ax.pcolormesh(rbins, ybins, np.log10(hist).T, cmap=cmap,
                                vmin=vmin, vmax=vmax, rasterized=True)
            rc = 0.5 * (rbins[:-1] + rbins[1:])
            yc = 0.5 * (ybins[:-1] + ybins[1:])
            ax.contour(rc, yc, chist.T, levels=lvls_contours, colors=colors,
                       linestyles=linestyles)
    plt.colorbar(img, cax=cax, orientation='vertical', extend='min')
    cax.set_ylabel(clabel, fontsize=fontsize)
    cax.tick_params(labelsize=fontsize - 1)

    ylims1 = [ax.get_ylim() for ax in axes[:len(wtstrs)]]
    ymin = min([yl[0] for yl in ylims1])
    ymax = max([yl[1] for yl in ylims1])
    ymin = max(ymin, -350.)
    [ax.set_ylim((ymin, ymax)) for ax in axes[:len(wtstrs)]]
    ylims2 = [ax.get_ylim() for ax in axes[len(wtstrs):]]
    ymin = min([yl[0] for yl in ylims2])
    ymax = max([yl[1] for yl in ylims2])
    ymax = min(ymax, 510.)
    [ax.set_ylim((ymin, ymax)) for ax in axes[len(wtstrs):]]    

    plt.savefig(outname, bbox_inches='tight')    



    

    



