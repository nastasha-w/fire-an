
import h5py
import numpy as np

import matplotlib.collections as mcol
import matplotlib.gridspec as gsp
import matplotlib.patches as mpatch
import matplotlib.patheffects as mppe
import matplotlib.pyplot as plt

import fire_an.makeplots.plot_utils as pu
import fire_an.utils.constants_and_units as c


def test_ionsum_and_Z_maps():
    
    fdir = '/Users/nastasha/ciera/tests/fire_start/map_tests/'
    ionb = 'coldens_{ion}_m13h206_m3e5__m13h206_m3e5_MHDCRspec1_fire3' + \
           '_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5' +\
           '_fcr1e-3_vw3000_snap27_shrink-sph-cen_BN98_2rvir_v1.hdf5'
    ionfiles = [fdir + ionb.format(ion='O{}'.format(i)) for i in range(1, 10)]
    eltfile = fdir + ionb.format(ion='Oxygen')
    massfile = fdir + ionb.format(ion='gas-mass')
    fdir_h = '/Users/nastasha/ciera/tests/fire_start/hist_tests/'
    histZfile = fdir_h + 'hist_Oxygen_by_Mass_0-1-2Rvir_m13h206_m3e5__' + \
                'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021' +\
                '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000_snap27_' + \
                'shrink-sph-cen_BN98_2rvir_v1.hdf5'
    
    outfilen = fdir + 'O-sum-and-Z-frac-test_m13h206_m3e5_MHDCRspec1_fire3' +\
           '_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5' +\
           '_fcr1e-3_vw3000_snap27_shrink-sph-cen_BN98_2rvir_v1.pdf'

    ionmaps = []
    ion_mapext = []
    eltmass = c.atomw_O * c.u
    
    width_ratios = [1.] * 4 + [0.2]
    fig = plt.figure(figsize=(11., 11.))
    grid = gsp.GridSpec(nrows=4, ncols=5, hspace=0.05, wspace=0.05,
                        width_ratios=width_ratios)
    coordsax = fig.add_subplot(grid[:4, :4])
    ionaxes = [fig.add_subplot(grid[i // 4, i % 4]) for i in range(9)]
    ionsumax = fig.add_subplot(grid[2, 1])
    elttotax = fig.add_subplot(grid[2, 2])
    deltaionsumax = fig.add_subplot(grid[2, 3])
    massax = fig.add_subplot(grid[3, 0])
    metax = fig.add_subplot(grid[3, 1])
    
    # slightly smaller panel to add axes labels and stuff
    _histax = fig.add_subplot(grid[3, 2])
    _histax.set_axis_off()
    _histax.tick_params(which='both', labelleft=False, labelbottom=False,
                        left=False, bottom=False)
    hpos = _histax.get_position()
    fmarx = 0.3
    fmary = 0.2
    histax = fig.add_axes([hpos.x0 + fmarx * hpos.width, 
                           hpos.y0 + fmary * hpos.height,
                           hpos.width * (1. - fmarx), 
                           hpos.height * (1. - fmary)])

    _dhistax = fig.add_subplot(grid[3, 3])
    _dhistax.set_axis_off()
    _dhistax.tick_params(which='both', labelleft=False, labelbottom=False,
                         left=False, bottom=False)
    hpos = _dhistax.get_position()
    fmarx = 0.3
    fmary = 0.2
    dhistax = fig.add_axes([hpos.x0 + fmarx * hpos.width, 
                            hpos.y0 + fmary * hpos.height,
                            hpos.width * (1. - fmarx), 
                            hpos.height * (1. - fmary)])
    

    fontsize = 12
    

    cmap_cd = 'afmhot'
    cmap_gas = 'viridis'
    cmap_Z = 'plasma'
    cmap_delta = 'RdBu'

    cax_i = fig.add_subplot(grid[0, 4])
    cax_Z = fig.add_subplot(grid[1, 4])
    cax_delta = fig.add_subplot(grid[2, 4])
    cax_gas = fig.add_subplot(grid[3, 4])

    coordsax.spines['right'].set_visible(False)
    coordsax.spines['top'].set_visible(False)
    coordsax.spines['left'].set_visible(False)
    coordsax.spines['bottom'].set_visible(False)
    coordsax.tick_params(which='both', labelbottom=False, labelleft=False,
                         left=False, bottom=False)
    coordsax.set_xlabel('X [pkpc]', fontsize=fontsize, labelpad=20.)
    coordsax.set_ylabel('Y [pkpc]', fontsize=fontsize, labelpad=40.)

    patheff_text = [mppe.Stroke(linewidth=2.0, foreground="white"),
                    mppe.Stroke(linewidth=0.4, foreground="black"),
                    mppe.Normal()]  

    vmin_i = np.inf
    vmax_i = -np.inf
    for ionf in ionfiles:
        with h5py.File(ionf, 'r') as f:
            #print(ionf)
            _map = f['map'][:]
            ionmaps.append(_map)
            vmin = f['map'].attrs['minfinite']
            vmax = f['map'].attrs['max']
            vmin_i = min(vmin, vmin_i)
            vmax_i = max(vmax, vmax_i)

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
            ion_mapext.append(extent)

    with h5py.File(eltfile, 'r') as f:
        map_elt = f['map'][:]
        vmin = f['map'].attrs['minfinite']
        vmax = f['map'].attrs['max']
        vmin_i = min(vmin, vmin_i)
        vmax_i = max(vmax, vmax_i)

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
        extent_elt = (-0.5 * box_pkpc[xax], 0.5 * box_pkpc[xax],
                      -0.5 * box_pkpc[yax], 0.5 * box_pkpc[yax])
    
    with h5py.File(massfile, 'r') as f:
        map_mass = f['map'][:]
        vmin = f['map'].attrs['minfinite']
        vmax = f['map'].attrs['max']
        vmin_m = min(vmin, vmin_i)
        vmax_m = max(vmax, vmax_i)

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
        extent_mass = (-0.5 * box_pkpc[xax], 0.5 * box_pkpc[xax],
                       -0.5 * box_pkpc[yax], 0.5 * box_pkpc[yax])
        
    
    _vmin_i = max(vmin_i, vmax_i - 10.)
    extlow_i = 'neither' if _vmin_i >= vmin_i else 'min'
    vtrans = 12.5
    if vtrans > _vmin_i and vtrans < vmax_i:
        _cmap_cd = pu.paste_cmaps(['gist_yarg', cmap_cd], 
        [_vmin_i, vtrans, vmax_i])    
    else:
        _cmap_cd = cmap_cd
    patheff_circ = [mppe.Stroke(linewidth=2.0, foreground="white"),
                    mppe.Stroke(linewidth=1.5, foreground="black"),
                    mppe.Normal()]  

    isum = np.zeros(ionmaps[0].shape, dtype=ionmaps[0].dtype)
    for ii, (imap, iext) in enumerate(zip(ionmaps, ion_mapext)):
        ax = ionaxes[ii]
        ion = 'O{}'.format(ii + 1)
        cen = (0.5 * (iext[0] + iext[1]), 0.5 * (iext[2] + iext[3]))
        ynum = ii % 4 == 0
        ax.tick_params(axis='both', labelsize=fontsize-1, labelbottom=False,
                       labelleft=ynum, direction='out')

        img_i = ax.imshow(imap.T, origin='lower', interpolation='nearest',
                          extent=iext, vmin=_vmin_i, vmax=vmax_i, 
                          cmap=_cmap_cd)
        ax.text(0.05, 0.95, ion, fontsize=fontsize,
                horizontalalignment='left', verticalalignment='top',
                transform=ax.transAxes, color='blue', 
                path_effects=patheff_text)
        patches = [mpatch.Circle(cen, rvir_pkpc)]
        collection = mcol.PatchCollection(patches)
        collection.set(edgecolor=['blue'], facecolor='none', linewidth=1.5,
                       linestyle='dashed', path_effects=patheff_circ)
        ax.add_collection(collection)
        if ii == 0:
            ax.text(1.05 * 2**-0.5 * rvir_pkpc, 1.05 * 2**-0.5 * rvir_pkpc, 
                    '$R_{\\mathrm{vir}}$',
                    color='blue', fontsize=fontsize,
                    path_effects=patheff_text)
        isum += 10**imap
    
    plt.colorbar(img_i, cax=cax_i, extend=extlow_i, orientation='vertical')
    cax_i.set_ylabel('$\\log_{10} \\, \\mathrm{N} \\; [\\mathrm{cm}^{-2}]$',
                     fontsize=fontsize)

    isum = np.log10(isum)
    ax = ionsumax
    ion = 'ion sum'
    cen = (0.5 * (iext[0] + iext[1]), 0.5 * (iext[2] + iext[3]))
    ynum = False
    ax.tick_params(axis='both', labelsize=fontsize-1, labelbottom=False,
                   labelleft=ynum, direction='out')

    ax.imshow(isum.T, origin='lower', interpolation='nearest',
              extent=iext, vmin=_vmin_i, vmax=vmax_i, 
              cmap=_cmap_cd)
    ax.text(0.05, 0.95, ion, fontsize=fontsize,
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, color='blue', 
            path_effects=patheff_text)
    patches = [mpatch.Circle(cen, rvir_pkpc)]
    collection = mcol.PatchCollection(patches)
    collection.set(edgecolor=['blue'], facecolor='none', linewidth=1.5,
                    linestyle='dashed', path_effects=patheff_circ)
    ax.add_collection(collection)
    
    ax = elttotax
    ion = 'all O'
    cen = (0.5 * (iext[0] + iext[1]), 0.5 * (iext[2] + iext[3]))
    ynum = False
    ax.tick_params(axis='both', labelsize=fontsize-1, labelbottom=False,
                   labelleft=ynum, direction='out')
    #print(map_elt)
    ax.imshow(map_elt.T, origin='lower', interpolation='nearest',
              extent=extent_elt, vmin=_vmin_i, vmax=vmax_i, 
              cmap=_cmap_cd)
    ax.text(0.05, 0.95, ion, fontsize=fontsize,
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, color='blue', 
            path_effects=patheff_text)
    patches = [mpatch.Circle(cen, rvir_pkpc)]
    collection = mcol.PatchCollection(patches)
    collection.set(edgecolor=['blue'], facecolor='none', linewidth=1.5,
                    linestyle='dashed', path_effects=patheff_circ)
    ax.add_collection(collection)
     
    ax = deltaionsumax
    _map = isum - map_elt
    maxd = np.max(np.abs(_map[np.isfinite(_map)]))
    if np.any(np.abs(_map) > maxd):
        extend = 'both'
    else:
        extend = 'neither'
    ion = 'ion sum - all O'
    cen = (0.5 * (iext[0] + iext[1]), 0.5 * (iext[2] + iext[3]))
    ynum = False
    ax.tick_params(axis='both', labelsize=fontsize-1, labelbottom=False,
                   labelleft=ynum, direction='out')
    img_delta = ax.imshow(_map.T, origin='lower', interpolation='nearest',
                          extent=extent_elt, vmin=-maxd, vmax=maxd, 
                          cmap=cmap_delta)
    ax.text(0.05, 0.95, ion, fontsize=fontsize,
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, color='black', 
            path_effects=patheff_text)
    patches = [mpatch.Circle(cen, rvir_pkpc)]
    collection = mcol.PatchCollection(patches)
    collection.set(edgecolor=['black'], facecolor='none', linewidth=1.5,
                    linestyle='dashed', path_effects=patheff_circ)
    ax.add_collection(collection)
    
    plt.colorbar(img_delta, cax=cax_delta, extend=extend, orientation='vertical')
    cax_delta.set_ylabel('$\\Delta \\, \\log_{10} \\, \\mathrm{N}$',
                         fontsize=fontsize)

    ax = dhistax
    ax.set_xlabel('$\\Delta \\, \\log_{10} \\, \\mathrm{N}$',
                  fontsize=fontsize)
    ax.set_ylabel('# pixels', fontsize=fontsize)
    ax.hist(_map.flatten(), bins=100, log=True, histtype='stepfilled',
            color='blue')
    ax.text(0.05, 0.95, 'ion sum - all O', fontsize=fontsize,
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, color='black', 
            path_effects=patheff_text)
    ax.set_xlim(-1.05 * maxd, 1.05 * maxd)
    ax.tick_params(axis='both', labelsize=fontsize-1, labelbottom=True,
                   labelleft=True, direction='in')
    ax.tick_params(axis='x', which='both', rotation=45.)
    
    ax = massax
    _map = map_mass
    extend = 'neither'
    ion = 'gas'
    cen = (0.5 * (extent_mass[0] + extent_mass[1]), 
           0.5 * (extent_mass[2] + extent_mass[3]))
    ynum = True
    ax.tick_params(axis='both', labelsize=fontsize-1, labelbottom=True,
                   labelleft=ynum, direction='out')
    img_mass = ax.imshow(_map.T, origin='lower', interpolation='nearest',
                          extent=extent_mass, cmap=cmap_gas)
    ax.text(0.05, 0.95, ion, fontsize=fontsize,
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, color='red', 
            path_effects=patheff_text)
    patches = [mpatch.Circle(cen, rvir_pkpc)]
    collection = mcol.PatchCollection(patches)
    collection.set(edgecolor=['red'], facecolor='none', linewidth=1.5,
                    linestyle='dashed', path_effects=patheff_circ)
    ax.add_collection(collection)

    plt.colorbar(img_mass, cax=cax_gas, extend=extend, orientation='vertical')
    cax_gas.set_ylabel('$\\log_{10} \\, \\Sigma \\; ' + \
                         '[\\mathrm{g}\\,\\mathrm{cm}^{-2}]$',
                         fontsize=fontsize)
    
    ax = metax
    _map = map_elt + np.log10(eltmass) - map_mass
    extend = 'neither'
    ion = 'O mass frac.'
    cen = (0.5 * (extent_mass[0] + extent_mass[1]), 
           0.5 * (extent_mass[2] + extent_mass[3]))
    ynum = False
    ax.tick_params(axis='both', labelsize=fontsize-1, labelbottom=True,
                   labelleft=ynum, direction='out')
    img_Z = ax.imshow(_map.T, origin='lower', interpolation='nearest',
                      extent=extent_mass, cmap=cmap_Z)
    ax.text(0.05, 0.95, ion, fontsize=fontsize,
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, color='red', 
            path_effects=patheff_text)
    patches = [mpatch.Circle(cen, rvir_pkpc)]
    collection = mcol.PatchCollection(patches)
    collection.set(edgecolor=['red'], facecolor='none', linewidth=1.5,
                    linestyle='dashed', path_effects=patheff_circ)
    ax.add_collection(collection)

    plt.colorbar(img_Z, cax=cax_Z, extend=extend, orientation='vertical')
    cax_Z.set_ylabel('$\\log_{10} \\, \\mathrm{Z}$', fontsize=fontsize)
        
    with h5py.File(histZfile, 'r') as f:
        hist = f['histogram/histogram'][:]
        hist -= np.log10(c.solar_mass)
        xvals = f['axis_1/bins'][:]
    
    ax = histax
    ynum = True
    ax.tick_params(axis='both', labelsize=fontsize-1, labelbottom=True,
                   labelleft=ynum, direction='in')
    label = 'O/M hist.'
    ax.text(0.05, 0.95, label, fontsize=fontsize,
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, color='black', 
            path_effects=patheff_text)
    xlabel = '$\\log_{10} \\, \\mathrm{M}_{\\mathrm{O}} \\,/\\,' +\
             '\\mathrm{M}_{\\mathrm{gas}}$'
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel('$\\log_{10} \\, \\mathrm{M} \\; [\\mathrm{M}_{\\odot}]$',
                  fontsize=fontsize)
    
    _hist = np.empty((hist.shape[0], hist.shape[1] + 1), dtype=hist.dtype)
    _hist[:, :-1] = hist
    _hist[:, -1] = -np.inf
    maxy = np.max(_hist)
    miny = np.min(_hist[np.isfinite(_hist)])
    _hist[_hist == -np.inf] = -100.

    ax.step(xvals, _hist[0, :], color='black', linewidth=2., where='post')
    ax.step(xvals, _hist[1, :], color='blue', linewidth=1.5, 
            linestyle='dashed', where='post')
    ax.text(0.05, 0.70, '$< \\mathrm{R}_{\\mathrm{vir}}$', 
            fontsize=fontsize - 2,
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, color='black', 
            path_effects=patheff_text)
    ax.text(0.05, 0.80, '$1 \\endash 2 \\, \\mathrm{R}_{\\mathrm{vir}}$', 
            fontsize=fontsize - 2,
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, color='blue', 
            path_effects=patheff_text)
    ax.set_ylim(miny * 0.95, maxy * 1.1)
    
    plt.savefig(outfilen, bbox_inches='tight')

def check_h1maps():
    '''
    specific sanity test comparing total H, H I from the PS20 tables,
    and H I from the FIRE NeutralHydrogenAbundance field
    '''

    mapdir = '/Users/nastasha/ciera/tests/fire_start/h1_sim_test/'
    mapf_Htot = ('coldens_Hydrogen_m12m_m7e3_MHD_fire3_fireBH_Sep182021_hr'
                 '_crdiffc690_sdp1e10_gacc31_fa0.5_snap500_shrink-sph-cen'
                 '_BN98_2rvir_v2.hdf5')
    mapf_H1sim = ('coldens_H1-sim_m12m_m7e3_MHD_fire3_fireBH_Sep182021_hr'
                  '_crdiffc690_sdp1e10_gacc31_fa0.5_snap500_shrink-sph-cen'
                  '_BN98_2rvir_v2.hdf5')
    mapf_H1PS20 = ('coldens_H1-PS20_m12m_m7e3_MHD_fire3_fireBH_Sep182021_hr'
                   '_crdiffc690_sdp1e10_gacc31_fa0.5_snap500_shrink-sph-cen'
                   '_BN98_2rvir_v2.hdf5')

    maps = {}
    extents = {}
    rvirs = {}
    vmin = np.inf
    vmax = -np.inf
    mapkeys = ['htot', 'h1sim', 'h1ps20']
    maptitles = {'htot': '$\\mathrm{N}(\\mathrm{H})$',
                 'h1sim': '$\\mathrm{N}(\\mathrm{H I})$, FIRE',
                 'h1ps20': '$\\mathrm{N}(\\mathrm{H I})$, PS20 table',
                 }
    fractitles = {'h1sim': ('$\\mathrm{N}(\\mathrm{H I}) \\,/\\,'
                           ' \\mathrm{N}(\\mathrm{H})$, FIRE'), 
                  'h1ps20': ('$\\mathrm{N}(\\mathrm{H I}) \\,/\\,'
                             ' \\mathrm{N}(\\mathrm{H})$, PS20 tables'),
                  }

    for mapkey, mapf in zip(mapkeys, 
                            [mapf_Htot, mapf_H1sim, mapf_H1PS20]):
        with h5py.File(mapdir + mapf, 'r') as f:
            _map = f['map'][:]
            _vmin = f['map'].attrs['minfinite']
            _vmax = f['map'].attrs['max']
            print(mapf, _vmin, _vmax)

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
            
            maps[mapkey] = _map
            extents[mapkey] = extent
            rvirs[mapkey] = rvir_pkpc
            vmin = min(vmin, _vmin)
            vmax = max(vmax, _vmax)
    
    mincol = 13.6
    if mincol > vmin and mincol < vmax:
        cmap = pu.paste_cmaps(['gist_yarg', 'viridis'], [vmin, mincol, vmax])
    else:
        cmap = 'virdis'
    anyzero = np.min([np.min(maps[key]) for key in maps]) == vmin
    extend = 'neither' if anyzero else 'min'

    fig = plt.figure(figsize=(8.5, 5.))
    grid = gsp.GridSpec(nrows=2, ncols=5, hspace=0.2, wspace=0.2, 
                        width_ratios=[1., 5., 5., 5., 1.])
    mapaxes = [fig.add_subplot(grid[0, 1+ i]) for i in range(3)]
    mapcax = fig.add_subplot(grid[0, 4])
    fracaxes = [fig.add_subplot(grid[1, 2 + i]) for i in range(2)]
    fraccax = fig.add_subplot(grid[1, 4])
    diffax = fig.add_subplot(grid[1, 1])
    diffcax = fig.add_subplot(grid[1, 0])
    fontsize = 12
    
    # should be same for all three maps
    xlabel = ['X', 'Y', 'Z'][xax] + ' [pkpc]'
    ylabel = ['X', 'Y', 'Z'][yax] + ' [pkpc]'

    for mi, mapkey in enumerate(mapkeys):
        ax = mapaxes[mi]
        if mi == 0:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.set_title(maptitles[mapkey], fontsize=fontsize)

        img = ax.imshow(maps[mapkey].T, origin='lower', 
                       interpolation='nearest', 
                       vmin=vmin, vmax=vmax, cmap=cmap, 
                       extent=extents[mapkey]
                        )
        ax.tick_params(axis='both', labelsize=fontsize-1,
                       labelleft=mi == 0, labelbottom=False)

        patches = [mpatch.Circle((0., 0.), rvir_pkpc)]
        collection = mcol.PatchCollection(patches)
        collection.set(edgecolor=['red'], facecolor='none', linewidth=1.5)
        ax.add_collection(collection)
        if mi == 0:
            ax.text(1.05 * 2**-0.5 * rvir_pkpc, 1.05 * 2**-0.5 * rvir_pkpc, 
                   '$R_{\\mathrm{vir}}$',
                    color='red', fontsize=fontsize)
    plt.colorbar(img, cax=mapcax, extend=extend, orientation='vertical') 
    mapcax.set_ylabel('$\\log_{10} \\, \\mathrm{N} \\; [\\mathrm{cm}^{-2}]$',
                      fontsize=fontsize) 

    fmaps = {'h1sim': maps['h1sim'] - maps['htot'],
             'h1ps20': maps['h1ps20'] - maps['htot']
             }
    fvmin = np.min([np.min(fmaps[fkey][np.isfinite(fmaps[fkey])])\
                    for fkey in fmaps])
    fvmax = np.min([np.max(fmaps[fkey][np.isfinite(fmaps[fkey])])\
                    for fkey in fmaps])
    anyzero = np.min([np.min(fmaps[fkey]) for fkey in fmaps]) == fvmin
    extend = 'neither' if anyzero else 'min'
    fcmap = 'afmhot'
    
    for mi, mapkey in enumerate(['h1sim', 'h1ps20']):
        ax = fracaxes[mi]
        ax.set_xlabel(xlabel, fontsize=fontsize)
        #if mi == 0:
        #    ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.set_title(fractitles[mapkey], fontsize=fontsize)

        img = ax.imshow(fmaps[mapkey].T, origin='lower', 
                        interpolation='nearest', 
                        vmin=fvmin, vmax=fvmax, cmap=fcmap, 
                        extent=extents[mapkey]
                        )
        ax.tick_params(axis='both', labelsize=fontsize-1, 
                       labelleft=False, labelbottom=True)

        patches = [mpatch.Circle((0., 0.), rvir_pkpc)]
        collection = mcol.PatchCollection(patches)
        collection.set(edgecolor=['red'], facecolor='none', linewidth=1.5)
        ax.add_collection(collection)
    plt.colorbar(img, cax=fraccax, extend=extend, orientation='vertical')
    fraccax.set_ylabel('$\\log_{10} \\, \\mathrm{N}(\\mathrm{H I}) \\,/\\,'
                       ' \\mathrm{N}(\\mathrm{H})$', fontsize=fontsize) 
    
    ax = diffax
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_title('FIRE / PS20', fontsize=fontsize)
    diffmap = fmaps['h1sim'] - fmaps['h1ps20']
    vmax = np.max(np.abs(diffmap[np.isfinite(diffmap)]))
    vmin = -1. * vmax
    img = ax.imshow(diffmap.T, origin='lower', 
                    interpolation='nearest', 
                    vmin=vmin, vmax=vmax, cmap='RdBu', 
                    extent=extents['h1sim']
                    )
    ax.tick_params(axis='both', labelsize=fontsize-1,
                   labelleft=False)
    patches = [mpatch.Circle((0., 0.), rvir_pkpc)]
    collection = mcol.PatchCollection(patches)
    collection.set(edgecolor=['red'], facecolor='none', linewidth=1.5)
    ax.add_collection(collection)
    plt.colorbar(img, cax=diffcax, extend='neither', orientation='vertical')
    diffcax.set_ylabel('$\\Delta \\, \\log_{10} \\,'
                       ' \\mathrm{N}(\\mathrm{H I})$', fontsize=fontsize) 
    diffcax.yaxis.set_label_position('left')
    diffcax.yaxis.tick_left()

    outdir = '/Users/nastasha/ciera/tests/fire_start/h1_sim_test/'
    savename = outdir + 'mapcomp_H_H1-sim_H1-PS20.pdf'
    if savename is not None:
        plt.savefig(savename, bbox_inches='tight')