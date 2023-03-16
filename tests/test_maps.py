
import h5py
import numpy as np
import os

import matplotlib.collections as mcol
import matplotlib.gridspec as gsp
import matplotlib.patches as mpatch
import matplotlib.patheffects as mppe
import matplotlib.pyplot as plt

import fire_an.mainfunc.makemap as mm 
import fire_an.makeplots.plot_utils as pu
import fire_an.utils.constants_and_units as c

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


