
import h5py
import matplotlib.cm as cm
import matplotlib.collections as mcol
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsp
import numpy as np

import fire_an.makeplots.plot_utils as pu
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c

def readmap(filen, weightmap=False):
    with h5py.File(filen, 'r') as f:
        if weightmap:
            mapkey = 'weightmap'
        else:
            mapkey = 'map'
        _map = f[mapkey][:]
        vmin = f[mapkey].attrs['minfinite']
        vmax = f[mapkey].attrs['max']

        box_cm = f['Header/inputpars'].attrs['diameter_used_cm']
        cosmopars = {key: val for key, val in \
                     f['Header/inputpars/cosmopars'].attrs.items()}
        #print(cosmopars)
        if 'Rvir_ckpcoverh' in f['Header/inputpars/halodata'].attrs:
            rvir_ckpcoverh = f['Header/inputpars/halodata'].attrs['Rvir_ckpcoverh']
            rvir_pkpc = rvir_ckpcoverh * cosmopars['a'] / cosmopars['h']
        elif 'Rvir_cm' in f['Header/inputpars/halodata'].attrs:
            rvir_cm = f['Header/inputpars/halodata'].attrs['Rvir_cm']
            rvir_pkpc = rvir_cm / (c.cm_per_mpc * 1e-3)
        xax = f['Header/inputpars'].attrs['Axis1']
        yax = f['Header/inputpars'].attrs['Axis2']
        box_pkpc = box_cm / (1e-3 * c.cm_per_mpc)
        extent = (-0.5 * box_pkpc[xax], 0.5 * box_pkpc[xax],
                  -0.5 * box_pkpc[yax], 0.5 * box_pkpc[yax])
    out = {'map': _map, 'vmin': vmin, 'vmax': vmax,
           'rvir_pkpc': rvir_pkpc, 'extent': extent}
    return out    

def plotmaps(filens, clabel, ctrans=13.3, sizeindic_pkpc=50., weightmap=False,
             axtitles=None, outname=None, vmin=None, vmax=None):

    mapdata = [readmap(filen, weightmap=weightmap) for filen in filens]
    if vmin is None:
        vmin = min([md['vmin'] for md in mapdata])
    if vmax is None:
        vmax = max([md['vmax'] for md in mapdata])
        vmin = max(vmin, vmax - 4.)
    if ctrans is None or ctrans < vmin or ctrans > vmax:
        cmap = cm.get_cmap('viridis')
    else:
        cmap = pu.paste_cmaps(['gist_yarg', 'viridis'], 
                              [vmin, ctrans, vmax])
    
    npanels = len(mapdata)
    ncols = min(npanels, 4)
    nrows = (npanels - 1) // ncols + 1
    panelsize = 2.5
    caxwidth = 0.5
    width_ratios = [panelsize] * ncols + [caxwidth]
    height_ratios = [panelsize] * nrows
    figsize = (sum(width_ratios), sum(height_ratios))
    fontsize = 12

    fig = plt.figure(figsize=figsize)
    grid = gsp.GridSpec(nrows=nrows, ncols=ncols + 1, hspace=0.0, 
                        wspace=0.0, width_ratios=width_ratios, 
                        height_ratios=height_ratios)
    axes = [fig.add_subplot(grid[i // ncols, i % ncols]) 
            for i in range(npanels)]
    cax = fig.add_subplot(grid[:2, ncols])
    
    xabsmin = min([md['extent'][1] for md in mapdata]) * 0.7
    yabsmin = min([md['extent'][3] for md in mapdata]) * 0.7
    for axi, (ax, mapd) in enumerate(zip(axes, mapdata)):
        img = ax.imshow(mapd['map'].T, extent=mapd['extent'], 
                        origin='lower', interpolation='nearest',
                        rasterized=True, vmin=vmin, vmax=vmax, 
                        cmap=cmap)
        ax.axis('off')
        ax.set_xlim((-xabsmin, xabsmin))
        ax.set_ylim((-yabsmin, yabsmin))

        rvir_pkpc = mapd['rvir_pkpc']
        patches = [mpatches.Circle((0., 0.), rvir_pkpc)]
        collection = mcol.PatchCollection(patches)
        collection.set(edgecolor=['red'], facecolor='none', linewidth=1.5)
        ax.add_collection(collection)
        if axi == 0:
            # label Rvir circle
            ax.text(1.05 * 2**-0.5 * rvir_pkpc, 1.05 * 2**-0.5 * rvir_pkpc, 
                    '$R_{\\mathrm{vir}}$',
                    color='red', fontsize=fontsize)
            
            # indicate size scale
            _xmin = -0.9 * xabsmin 
            _xmax = _xmin + sizeindic_pkpc
            _yv = -0.9 * yabsmin
            ax.plot([_xmin, _xmax], [_yv, _yv], color='red', linestyle='solid',
                    linewidth=2.5)
            ax.text(-0.9 * xabsmin, -0.87 * yabsmin, 
                    f'{sizeindic_pkpc:.0f} pkpc', fontsize=fontsize,
                    horizontalalignment='left', verticalalignment='bottom',
                    color='red')

        if axtitles is not None:
            ax.set_title(axtitles[axi], fontsize=fontsize)
    cbar = plt.colorbar(img, cax=cax, orientation='vertical', extend='both')
    cax.tick_params(which='both', labelsize=fontsize - 1.)
    cbar.set_label(clabel, fontsize=fontsize)        

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotoverview_ne8(ics, zis):
    clabel = ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\,VIII})'
              '\\; [\\mathrm{cm}^{-2}]$')
    ddir = '/projects/b1026/nastasha/maps/vdopmaps_all2/'
    filen_temp = ('vdoplos_by_coldens_Ne8_{simname}_snap{snapnum}'
                  '_shrink-sph-cen_BN98_depth_2.0rvir_{pax}-proj_v3.hdf5')
    outdir = '/projects/b1026/nastasha/imgs/summary_plots/ne8_thumbnails/'
    simnames_all = sl.m12_f2md + \
                   sl.m12_nobh_clean2 + sl.m12_nobh_rest2 +\
                   sl.m12_agnnocr_clean2 + sl.m12_agnnocr_rest2 +\
                   sl.m12_agncr_clean2 + sl.m12_agncr_rest2 +\
                   sl.m13_nobh_clean2 + sl.m13_nobh_rest2 +\
                   sl.m13_agnnocr_clean2 + sl.m13_agnnocr_rest2 +\
                   sl.m13_agncr_clean2 + sl.m13_agncr_rest2
    simnames_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    simnames_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    simnames_f2md = sl.m12_f2md
    zs_approx = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
    pax = 'z'
    for ic in ics:
        simnames = []
        for sn in simnames_all:
            if sl.ic_from_simname(sn) == ic and sn not in sl.buglist1:
                simnames.append(sn)
        if len(simnames) == 0:
            print(f'Found no (non-bug) simulations for IC {ic}')
            continue
        axtitles = [sl.plotlabel_from_physlabel[sl.physlabel_from_simname(sn)]
                    for sn in simnames]
        for zi in zis:
            z_approx = zs_approx[zi]
            snaps = [sl.snaps_hr[zi] if sn in simnames_hr 
                     else sl.snaps_sr[zi] if sn in simnames_sr
                     else sl.snaps_f2md[zi] if sn in simnames_f2md
                     else None
                     for sn in simnames]
            filens = [ddir + filen_temp.format(simname=simname, snapnum=snap, 
                                               pax=pax)
                      for simname, snap in zip(simnames, snaps)]
            if ic.startswith('m12'):
                sizeindic_pkpc = 50.
            else:
                sizeindic_pkpc = 100.
            outname = f'coldens_ne8_images_{ic}_z{z_approx:.1f}'
            outname = outname.replace('.', 'p')
            outname = outdir + outname + '.pdf'
            
            plotmaps(filens, clabel, ctrans=13.3, 
                     sizeindic_pkpc=sizeindic_pkpc, 
                     weightmap=True,
                     axtitles=axtitles, outname=outname,
                     vmin=11., vmax=15.)

def plotoverview_ne8_prop(ics, zis):
    '''
    simplified proposal version
    '''
    clabel = ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\,VIII})'
              '\\; [\\mathrm{cm}^{-2}]$')
    ddir = '/projects/b1026/nastasha/maps/vdopmaps_all2/'
    filen_temp = ('vdoplos_by_coldens_Ne8_{simname}_snap{snapnum}'
                  '_shrink-sph-cen_BN98_depth_2.0rvir_{pax}-proj_v3.hdf5')
    outdir = '/projects/b1026/nastasha/imgs/summary_plots/ne8_thumbnails/'
    simnames_all = sl.m12_f2md + \
                   sl.m12_nobh_clean2 + sl.m12_nobh_rest2 +\
                   sl.m12_agnnocr_clean2 + sl.m12_agnnocr_rest2 +\
                   sl.m12_agncr_clean2 + sl.m12_agncr_rest2 +\
                   sl.m13_nobh_clean2 + sl.m13_nobh_rest2 +\
                   sl.m13_agnnocr_clean2 + sl.m13_agnnocr_rest2 +\
                   sl.m13_agncr_clean2 + sl.m13_agncr_rest2
    simnames_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    simnames_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    simnames_f2md = sl.m12_f2md
    zs_approx = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
    pax = 'z'
    titlemap = {'noBH': 'no AGN',
                'AGN-noCR': 'AGN',
                'AGN-CR': 'AGN + cosmic rays',
                'FIRE-2': ''}
    for ic in ics:
        simnames = []
        for sn in simnames_all:
            if sl.ic_from_simname(sn) == ic and sn not in sl.buglist1 \
                    and sl.physlabel_from_simname(sn) not in ['FIRE-2']:
                simnames.append(sn)
        if len(simnames) == 0:
            print(f'Found no (non-bug) simulations for IC {ic}')
            continue
        axtitles = [titlemap[sl.physlabel_from_simname(sn)]
                    for sn in simnames]
        for zi in zis:
            z_approx = zs_approx[zi]
            snaps = [sl.snaps_hr[zi] if sn in simnames_hr 
                     else sl.snaps_sr[zi] if sn in simnames_sr
                     else sl.snaps_f2md[zi] if sn in simnames_f2md
                     else None
                     for sn in simnames]
            filens = [ddir + filen_temp.format(simname=simname, snapnum=snap, 
                                               pax=pax)
                      for simname, snap in zip(simnames, snaps)]
            if ic.startswith('m12'):
                sizeindic_pkpc = 50.
            else:
                sizeindic_pkpc = 100.
            outname = f'coldens_ne8_images_{ic}_z{z_approx:.1f}_F3only'
            outname = outname.replace('.', 'p')
            outname = outdir + outname + '.pdf'
            
            plotmaps(filens, clabel, ctrans=13.3, 
                     sizeindic_pkpc=sizeindic_pkpc, 
                     weightmap=True,
                     axtitles=axtitles, outname=outname,
                     vmin=11., vmax=14.5)

def plotoverview_ne8_fire2(ics, zis):
    '''
    talk version: just the FIRE-2 plot
    '''
    clabel = ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\,VIII})'
              '\\; [\\mathrm{cm}^{-2}]$')
    ddir = '/projects/b1026/nastasha/maps/vdopmaps_all2/'
    filen_temp = ('vdoplos_by_coldens_Ne8_{simname}_snap{snapnum}'
                  '_shrink-sph-cen_BN98_depth_2.0rvir_{pax}-proj_v3.hdf5')
    outdir = '/projects/b1026/nastasha/imgs/summary_plots/ne8_thumbnails/'
    simnames_all = sl.m12_f2md
    simnames_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    simnames_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    simnames_f2md = sl.m12_f2md
    zs_approx = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
    pax = 'z'
    for ic in ics:
        simnames = []
        for sn in simnames_all:
            if sl.ic_from_simname(sn) == ic and sn not in sl.buglist1:
                simnames.append(sn)
        if len(simnames) == 0:
            print(f'Found no (non-bug) simulations for IC {ic}')
            continue
        axtitles = [sl.plotlabel_from_physlabel[sl.physlabel_from_simname(sn)]
                    for sn in simnames]
        for zi in zis:
            z_approx = zs_approx[zi]
            snaps = [sl.snaps_hr[zi] if sn in simnames_hr 
                     else sl.snaps_sr[zi] if sn in simnames_sr
                     else sl.snaps_f2md[zi] if sn in simnames_f2md
                     else None
                     for sn in simnames]
            filens = [ddir + filen_temp.format(simname=simname, snapnum=snap, 
                                               pax=pax)
                      for simname, snap in zip(simnames, snaps)]
            if ic.startswith('m12'):
                sizeindic_pkpc = 50.
            else:
                sizeindic_pkpc = 100.
            outname = f'coldens_ne8_images_{ic}_z{z_approx:.1f}_fire2'
            outname = outname.replace('.', 'p')
            outname = outdir + outname + '.pdf'
            
            plotmaps(filens, clabel, ctrans=13.3, 
                     sizeindic_pkpc=sizeindic_pkpc, 
                     weightmap=True,
                     axtitles=axtitles, outname=outname,
                     vmin=11., vmax=15.)
    
