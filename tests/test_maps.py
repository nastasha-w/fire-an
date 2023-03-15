
import h5py
import numpy as np
import os

import matplotlib.collections as mcol
import matplotlib.gridspec as gsp
import matplotlib.patches as mpatch
import matplotlib.patheffects as mppe
import matplotlib.pyplot as plt

import mainfunc.makemap as mm 
import makeplots.plot_utils as pu
import utils.constants_and_units as c

# hard to do a true test, but check that projected masses and centering
# sort of make sense
def tryout_massmap(opt=1, center='AHFsmooth'):
    outdir = 'ls'
    _outfilen = 'mass_pt{pt}_{sc}_snap{sn}_ahf-cen_2rvir_v1.hdf5'
    if opt == 1:
        parttypes = [0, 1, 4]
        dirpath = '/projects/b1026/snapshots/metal_diffusion/m12i_res7100/'
        simcode = 'metal-diffusion-m12i-res7100'
        snapnum = 600
    elif opt == 2:
        parttypes = [0, 1, 4]
        dirpath = '/projects/b1026/snapshots/metal_diffusion/m12i_res7100/'
        simcode = 'metal-diffusion-m12i-res7100'
        snapnum = 399
    elif opt == 3:
        parttypes = [0, 1, 4]
        dirpath = '/projects/b1026/snapshots/metal_diffusion/m12i_res7100/'
        simcode = 'metal-diffusion-m12i-res7100'
        snapnum = 196

    for pt in parttypes:
        outfilen = outdir + _outfilen.format(pt=pt, sc=simcode, 
                                             sn=snapnum)
        mm.massmap(dirpath, snapnum, radius_rvir=2., particle_type=pt,
                   pixsize_pkpc=3., axis='z', outfilen=outfilen,
                   center=center)
        
def tryout_wholezoom(index):
    outdir = '/projects/b1026/nastasha/tests/start_fire/map_tests/'

    if index == 0:
        dirpath = '/projects/b1026/snapshots/fire3/m13h206_m3e5/' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000/' 
        simname = 'm13h206_m3e5__' + \
                  'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'                     
        snapnum = 27  
        outfilen_template = 'mass_pt{pt}_{sc}_snap{sn}_axis-{ax}_' + \
                            'wholezoom_v1.hdf5'
        _temp = outdir + outfilen_template 
        outfilens = {'outfilen_gas': _temp.format(pt=0, sc=simname, 
                                                 sn=snapnum, ax='{ax}'),
                     'outfilen_DM': _temp.format(pt=1, sc=simname, 
                                                 sn=snapnum, ax='{ax}'),
                     'outfilen_stars': _temp.format(pt=4, sc=simname, 
                                                 sn=snapnum, ax='{ax}'),
                     'outfilen_BH': _temp.format(pt=5, sc=simname, 
                                                 sn=snapnum, ax='{ax}'),                            
                    }

    mm.massmap_wholezoom(dirpath, snapnum, pixsize_pkpc=3.,
                         **outfilens)

def hasnan_map(filen):
    with h5py.File(filen, 'r') as f:
        if 'map' not in f:
            print('skipping {}'.format(filen))
            return False # don't flag anything in a file loop
        map_ = f['map'][:]
        if np.any(np.isnan(map_)):
            print('NaN values in map {}'.format(filen))
            return True
        else:
            return False

def checkdir_nanmap(dirn):
    filens = os.listdir(dirn)
    filens = [filen for filen in filens if filen.endswith('.hdf5')]
    anynan = False
    for filen in filens:
        anynan |= hasnan_map(dirn + filen)
    return anynan

def checkcenter_massmap(filen_template, savename=None, mincol=None,
                        center_simunits=None, Rvir_simunits=None):
    '''
    quick plot of the mass map in the file
    '''

    filens = {ax: filen_template.format(ax=ax) for ax in ['x', 'y', 'z']}
    
    fig = plt.figure(figsize=(8., 8.))
    grid = gsp.GridSpec(nrows=2, ncols=2, hspace=0.2, wspace=0.2, 
                        width_ratios=[1., 1.])
    axes = {}
    axes['z'] = fig.add_subplot(grid[0, 0]) 
    axes['y'] = fig.add_subplot(grid[0, 1])
    axes['x'] = fig.add_subplot(grid[1, 0])
    cax = fig.add_subplot(grid[1, 1])
    fontsize = 12

    axlabels = ['{} [sim. units: ckpc/h]'. format(_ax) for _ax in 'XYZ']
    
    massmaps = {}
    extents = {}
    xlabels = {}
    ylabels = {}
    xinds = {}
    yinds = {}
    vmin = np.inf
    vmax = -np.inf
    for ax in axes:
        filen = filens[ax]
        with h5py.File(filen, 'r') as f:
            _map = f['map'][:]
            vmin = min(vmin, f['map'].attrs['minfinite'])
            vmax = max(vmax, f['map'].attrs['max'])
            
            # error in file creation -- no actual conversion to cm
            region_simunits = f['Header/inputpars'].attrs['mapped_region_cm']
            #coords_to_CGS = f['Header/inputpars'].attrs['coords_toCGS']
            # region_simunits = region_cm / coords_to_CGS

            #box_cm = f['Header/inputpars'].attrs['diameter_used_cm']
            cosmopars = {key: val for key, val in \
                        f['Header/inputpars/cosmopars'].attrs.items()}
            _ax1 = f['Header/inputpars'].attrs['Axis1']
            _ax2 = f['Header/inputpars'].attrs['Axis2']
            _ax3 = f['Header/inputpars'].attrs['Axis3']
            if _ax3 == 2:
                xax = _ax1
                yax = _ax2
            elif _ax3 == 0:
                xax = _ax2
                yax = _ax1
                _map = _map.T
            elif _ax3 == 1:
                xax = _ax2
                yax = _ax1
                _map = _map.T
           
            extent = (region_simunits[xax][0], region_simunits[xax][1],
                      region_simunits[yax][0], region_simunits[yax][1])
            massmaps[ax] = _map
            extents[ax] = extent
            xlabels[ax] = axlabels[xax]
            ylabels[ax] = axlabels[yax]
            xinds[ax] = xax
            yinds[ax] = yax
    print('redshift: ', cosmopars['z'])

    if mincol is None:
        cmap = 'viridis'
    else:
        cmap = pu.paste_cmaps(['gist_yarg', 'viridis'], [vmin, mincol, vmax])
    extend = 'neither' if np.min(map) == vmin else 'min'
    
    for axn in axes:
        ax = axes[axn]
        ax.set_xlabel(xlabels[axn], fontsize=fontsize)
        ax.set_ylabel(ylabels[axn], fontsize=fontsize)

        img = ax.imshow(massmaps[axn].T, origin='lower', 
                        interpolation='nearest', vmin=vmin,
                        vmax=vmax, cmap=cmap, extent=extents[axn])
        ax.tick_params(axis='both', labelsize=fontsize-1)
        
        if center_simunits is not None:
            _cen = [center_simunits[xinds[axn]], center_simunits[yinds[axn]]]
            ax.scatter([_cen[0]], [_cen[1]], marker='.', color='red',
                        s=10)
            if Rvir_simunits is not None:
                patches = [mpatch.Circle(_cen, Rvir_simunits)]
                collection = mcol.PatchCollection(patches)
                collection.set(edgecolor=['red'], facecolor='none', 
                               linewidth=1.5)
                ax.add_collection(collection)
                ax.text(1.05 * 2**-0.5 * Rvir_simunits, 
                        1.05 * 2**-0.5 * Rvir_simunits, 
                        '$R_{\\mathrm{vir}}$',
                        color='red', fontsize=fontsize)
    
    cbar = plt.colorbar(img, cax=cax, extend=extend, orientation='horizontal',
                        aspect=10)
    clabel = 'surface density $[\\log_{10} \\mathrm{g}\\,\\mathrm{cm}^{-2}]$'
    cax.set_xlabel(clabel, fontsize=fontsize)
    cax.tick_params(labelsize=fontsize-1)
    
    if savename is not None:
        plt.savefig(savename, bbox_inches='tight')


def run_checkcenter_massmap(index, center=None, rvir=None,
                            masstype='gas'):
    outdir = '/projects/b1026/nastasha/tests/start_fire/map_tests/'
    cen = center
    mincols = {'gas': -5.,
               'DM': -5.,
               'stars': -7.,
               'BH': None}
    if index == 0:
        dirpath = '/projects/b1026/snapshots/fire3/m13h206_m3e5/' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000/' 
        simname = 'm13h206_m3e5__' + \
                  'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'                     
        snapnum = 27  
        outfilen_template = 'mass_pt{pt}_{sc}_snap{sn}_axis-{ax}_' + \
                            'wholezoom_v1.hdf5'
        _temp = outdir + outfilen_template 
        mapfilens = {'gas': _temp.format(pt=0, sc=simname, 
                                                 sn=snapnum, ax='{ax}'),
                     'DM': _temp.format(pt=1, sc=simname, 
                                                 sn=snapnum, ax='{ax}'),
                     'stars': _temp.format(pt=4, sc=simname, 
                                                 sn=snapnum, ax='{ax}'),
                     'BH': _temp.format(pt=5, sc=simname, 
                                                 sn=snapnum, ax='{ax}'),                            
                    }
        cen =  [48414.20743443, 49480.35333529, 48451.20700497]
    mapfile_template = mapfilens[masstype]
    
    checkcenter_massmap(mapfile_template, savename=None, 
                        mincol=mincols[masstype],
                        center_simunits=cen, Rvir_simunits=rvir)
    

def masstest_map(filens):
    '''
    files for all parttypes
    '''
    
    # mass in g will overflow 
    enclmass = np.float64(0.)
    for filen in filens:
        with h5py.File(filen, 'r') as f:
            map = 10**f['map'][:]

            box_cm = f['Header/inputpars'].attrs['diameter_used_cm']
            cosmopars = {key: val for key, val in \
                         f['Header/inputpars/cosmopars'].attrs.items()}
            #print(cosmopars)
            halopath = 'Header/inputpars/halodata'
            rvir_ckpcoverh = f[halopath].attrs['Rvir_ckpcoverh']
            mvir_msunoverh = np.float64(f[halopath].attrs['Mvir_Msunoverh'])
            pixsize_pkpc = f['Header/inputpars'].attrs['pixsize_pkpc']
            rvir_pkpc = rvir_ckpcoverh * cosmopars['a'] / cosmopars['h']
            xax = f['Header/inputpars'].attrs['Axis1']
            yax = f['Header/inputpars'].attrs['Axis2']
            box_pkpc = box_cm / (1e-3 * c.cm_per_mpc)
            xcminmax = (-0.5 * box_pkpc[xax] + 0.5 * pixsize_pkpc, 
                        0.5 * box_pkpc[xax] - 0.5 * pixsize_pkpc)
            ycminmax = (-0.5 * box_pkpc[yax] + 0.5 * pixsize_pkpc,
                        0.5 * box_pkpc[yax] - 0.5 * pixsize_pkpc)
            npix_x = map.shape[0]
            npix_y = map.shape[1]
            pixdist2_pkpc = np.linspace(xcminmax[0], xcminmax[1], npix_x)**2 +\
                            np.linspace(ycminmax[0], ycminmax[1], npix_y)**2
            mapsel = pixdist2_pkpc < rvir_pkpc**2
            pixarea_cm = (pixsize_pkpc * 1e-3 * c.cm_per_mpc)**2
            partmass = np.float64(np.sum(map[mapsel])) * pixarea_cm
            enclmass += partmass
    halomass_g = mvir_msunoverh * c.solar_mass / cosmopars['h']
    print('Found Mvir (AHF) = {:.3e} g'.format(halomass_g))
    print('Found enclosed mass in projection = {:.3e} g'.format(enclmass))

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
