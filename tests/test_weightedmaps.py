import h5py
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np

import fire_an.mainfunc.makemap as mm
import fire_an.utils.constants_and_units as c

# quest
#tdir = '/projects/b1026/nastasha/tests/weightedmaps/'
# laptop
tdir = '/Users/nastasha/ciera/tests/weightedmaps/'

simname = ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
           '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000')
simpath = 'm12f_m6e4/' + simname
snapnum = 45

fn_stars_noweight = tdir + 'map_stars_m12f_m6e4_AGN-CR_snap45_noweight.hdf5'
fn_vlos_starwtd = tdir + \
                  'map_vlos_starwtd_starvcen_m12f_m6e4_AGN-CR_snap45.hdf5'
fn_gas_noweight = tdir + 'map_gas_m12f_m6e4_AGN-CR_snap45_noweight.hdf5'
fn_vlos_gaswtd = tdir + 'map_vlos_gaswtd_allvcen_m12f_m6e4_AGN-CR_snap45.hdf5'
fn_temp_gaswtd = tdir + \
                 'map_temperature_gaswtd_allvcen_m12f_m6e4_AGN-CR_snap45.hdf5'
fn_ne8_noweight = tdir + 'map_ne8_m12f_m6e4_AGN-CR_snap45_noweight.hdf5'
fn_vlos_ne8wtd = tdir + 'map_vlos_ne8wtd_allvcen_m12f_m6e4_AGN-CR_snap45.hdf5'
fn_temp_ne8wtd = tdir + \
                 'map_temperature_ne8wtd_allvcen_m12f_m6e4_AGN-CR_snap45.hdf5'

# debug until runs, get the data, 
# check outputs by eye (properties stored right)
def run_maps():
    print('\n\n')
    print('Trying to run unweighted star mass map')
    mm.massmap(simpath, snapnum, radius_rvir=2., particle_type=4,
               pixsize_pkpc=3., axis='x', outfilen=fn_stars_noweight,
               center='shrinksph', norm='pixsize_phys',
               maptype='Mass', maptype_args=None,
               weighttype=None, weighttype_args=None,
               save_weightmap=False, logmap=True,
               logweightmap=True)
    print('\n\n')
    print('Trying to run star-weighted vlos map')
    mm.massmap(simpath, snapnum, radius_rvir=2., particle_type=4,
               pixsize_pkpc=3., axis='x', outfilen=fn_vlos_starwtd,
               center='shrinksph', norm='pixsize_phys',
               maptype='coords', maptype_args={'vel': 'los', 
                    'vcen_cmps': {'parttypes': (4,), 'radius_rvir': 0.1}},
               weighttype='Mass', 
               weighttype_args=None,
               save_weightmap=True, logmap=False,
               logweightmap=True)
    print('\n\n')
    print('Trying to run unweighted gas mass map')
    mm.massmap(simpath, snapnum, radius_rvir=2., particle_type=0,
               pixsize_pkpc=3., axis='x', outfilen=fn_gas_noweight,
               center='shrinksph', norm='pixsize_phys',
               maptype='Mass', maptype_args=None,
               weighttype=None, weighttype_args=None,
               save_weightmap=False, logmap=True,
               logweightmap=False)
    print('\n\n')
    print('Trying to run gas-weighted vlos map')
    mm.massmap(simpath, snapnum, radius_rvir=2., particle_type=0,
               pixsize_pkpc=3., axis='x', outfilen=fn_vlos_gaswtd,
               center='shrinksph', norm='pixsize_phys',
               maptype='coords', maptype_args={'vel': 'los'},
               weighttype='Mass', 
               weighttype_args=None,
               save_weightmap=True, logmap=False,
               logweightmap=True)
    print('\n\n')
    print('Trying to run gas-weighted temperature map')
    mm.massmap(simpath, snapnum, radius_rvir=2., particle_type=0,
               pixsize_pkpc=3., axis='x', outfilen=fn_temp_gaswtd,
               center='shrinksph', norm='pixsize_phys',
               maptype='sim-direct', maptype_args={'field': 'Temperature'},
               weighttype='Mass', 
               weighttype_args=None,
               save_weightmap=False, logmap=True,
               logweightmap=True)
    print('\n\n')
    print('Trying to run unweighted ne8 map')
    mm.massmap(simpath, snapnum, radius_rvir=2., particle_type=0,
               pixsize_pkpc=3., axis='x', outfilen=fn_ne8_noweight,
               center='shrinksph', norm='pixsize_phys',
               maptype='ion', maptype_args={'ion': 'Ne8'},
               weighttype=None, weighttype_args=None,
               save_weightmap=False, logmap=True,
               logweightmap=False)
    print('\n\n')
    print('Trying to run ne8-weighted vlos map')
    mm.massmap(simpath, snapnum, radius_rvir=2., particle_type=0,
               pixsize_pkpc=3., axis='x', outfilen=fn_vlos_ne8wtd,
               center='shrinksph', norm='pixsize_phys',
               maptype='coords', maptype_args={'vel': 'los'},
               weighttype='ion', 
               weighttype_args={'ion': 'Ne8'},
               save_weightmap=True, logmap=False,
               logweightmap=True)
    print('\n\n')
    print('Trying to run ne8-weighted temperature map')
    mm.massmap(simpath, snapnum, radius_rvir=2., particle_type=0,
               pixsize_pkpc=3., axis='x', outfilen=fn_temp_ne8wtd,
               center='shrinksph', norm='pixsize_phys',
               maptype='sim-direct', maptype_args={'field': 'Temperature'},
               weighttype='ion',
               weighttype_args={'ion': 'Ne8'},
               save_weightmap=False, logmap=True,
               logweightmap=True)
    
def checkwsame():
    pairs = [(fn_stars_noweight, fn_vlos_starwtd),
             (fn_gas_noweight, fn_vlos_gaswtd),
             (fn_ne8_noweight, fn_vlos_ne8wtd),
             ]
    passed = True
    for now, wtd in pairs:
        with h5py.File(now, 'r') as f:
            map_now = f['map'][:]
        with h5py.File(wtd, 'r') as f:
            map_wtd = f['weightmap'][:]
        _passed = np.allclose(map_now, map_wtd)
        passed &= _passed
        if not _passed:
            print(f'Map mismatch:\n{now}\n{wtd}')
    return passed

def getmap(fn):
    with h5py.File(fn, 'r') as f:
        out = f['map'][:]
        pixs = f['Header/inputpars'].attrs['pixsize_pkpc']
        extent = (-0.5 * pixs * out.shape[0], 0.5 * pixs * out.shape[0],
                  -0.5 * pixs * out.shape[1], 0.5 * pixs * out.shape[1])
    return out, extent

def plot_maps():

    fontsize = 12 
    cmap_vel = 'RdBu'
    cmap_T = 'afmhot'
    cmap_gas = 'viridis'
    cmap_stars = 'gist_yarg'
    cmap_ne8 = 'plasma'
    
    title = 'maps and weighted maps test plot, m12f, AGN-CR, z=1.0, x-projection'
    outname = tdir + 'weightedmaps_plot_sanity_check_m12f_AGN-CR_snap45.pdf'

    fns = {'star': fn_stars_noweight,
           'starvel': fn_vlos_starwtd,
           'gas': fn_gas_noweight,
           'gasvel': fn_vlos_gaswtd,
           'gastemp': fn_temp_gaswtd,
           'ne8': fn_ne8_noweight,
           'ne8vel': fn_vlos_ne8wtd,
           'ne8temp': fn_temp_ne8wtd}
    
    maps = {key: getmap(fns[key]) for key in fns.keys()}
    vmin_vel = np.min([np.min(maps[key][0][np.isfinite(maps[key][0])]) 
                       for key in maps if key.endswith('vel')])
    vmax_vel = np.max([np.max(maps[key][0][np.isfinite(maps[key][0])]) 
                       for key in maps if key.endswith('vel')])
    vmax_vel = max(np.abs(vmax_vel), np.abs(vmin_vel))
    vmin_vel = -1. * vmax_vel
    vmin_t = np.min([np.min(maps[key][0][np.isfinite(maps[key][0])]) 
                     for key in maps if key.endswith('temp')])
    vmax_t = np.max([np.max(maps[key][0][np.isfinite(maps[key][0])]) 
                     for key in maps if key.endswith('temp')])
    
    panelsize = 2.5
    caxspace = 0.5     
    ncols = 3 
    nrows = 3
    wspace = 0.35
    hspace = 0.35
    width_ratios = [panelsize] * ncols + [caxspace]
    width = np.sum(width_ratios) * \
            (1. + wspace * ncols / np.sum(width_ratios))
    height_ratios = [panelsize] * nrows + [caxspace]
    height = np.sum(height_ratios) * \
            (1. + wspace * ncols / np.sum(height_ratios))
    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols + 1, nrows=nrows + 1,
                        width_ratios=width_ratios, 
                        height_ratios=height_ratios,
                        wspace=wspace, hspace=hspace)
    fig.suptitle(title, fontsize=fontsize)

    star_cax = fig.add_subplot(grid[0, ncols])
    gas_cax = fig.add_subplot(grid[1, ncols])
    ne8_cax = fig.add_subplot(grid[2, ncols])
    vel_cax = fig.add_subplot(grid[nrows, 1])
    temp_cax = fig.add_subplot(grid[nrows, 2])

    starax = fig.add_subplot(grid[0, 0])
    starunits = c.solar_mass / (c.cm_per_mpc * 1e-3)**2
    starimg = starax.imshow(maps['star'][0].T - np.log10(starunits), 
                            extent=maps['star'][1],
                            origin='lower', interpolation='nearest',
                            cmap=cmap_stars, vmin=0.)
    starax.set_ylabel('pkpc', fontsize=fontsize)
    plt.colorbar(starimg, cax=star_cax, orientation='vertical', extend='min')
    star_cax.set_ylabel('$\\log_{10} \\, \\Sigma_{*} \\, '
                        '[\\mathrm{M}_{\\odot} \\; \\mathrm{pkpc}^{-2}]$',
                        fontsize=fontsize)
    
    gasax = fig.add_subplot(grid[1, 0])
    gasunits = 0.752 / c.atomw_H * c.u
    gasimg = gasax.imshow(maps['gas'][0].T - np.log10(gasunits), 
                          extent=maps['gas'][1],
                          origin='lower', interpolation='nearest',
                          cmap=cmap_gas)
    gasax.set_ylabel('pkpc', fontsize=fontsize)
    plt.colorbar(gasimg, cax=gas_cax, orientation='vertical')
    gas_cax.set_ylabel('$\\log_{10} \\, \\Sigma_{\\mathrm{gas}} \\; '
                        '[\\mathrm{cm}^{-2}]$',
                        fontsize=fontsize)
    
    ne8ax = fig.add_subplot(grid[2, 0])
    ne8img = ne8ax.imshow(maps['ne8'][0].T, 
                          extent=maps['ne8'][1],
                          origin='lower', interpolation='nearest',
                          cmap=cmap_ne8)
    ne8ax.set_ylabel('pkpc', fontsize=fontsize)
    ne8ax.set_xlabel('pkpc', fontsize=fontsize)
    plt.colorbar(ne8img, cax=ne8_cax, orientation='vertical')
    ne8_cax.set_ylabel('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\,VIII}) \\;'
                        '[\\mathrm{cm}^{-2}]$',
                        fontsize=fontsize)

    svax = fig.add_subplot(grid[0, 1])
    vunits = 1e5
    vmin_vel /= vunits
    vmax_vel /= vunits
    vmin_vel = -150.
    vmax_vel = 150.
    svimg = svax.imshow(maps['starvel'][0].T / vunits, 
                        extent=maps['starvel'][1],
                        origin='lower', interpolation='nearest',
                        cmap=cmap_vel, vmin=vmin_vel, vmax=vmax_vel)
    plt.colorbar(svimg, cax=vel_cax, orientation='horizontal', extend='both')
    vel_cax.set_xlabel('$\\mathrm{v}_{\\mathrm{los}}\\, '
                       '[\\mathrm{km} \\; \\mathrm{s}^{-1}]$',
                       fontsize=fontsize)
    
    gvax = fig.add_subplot(grid[1, 1])
    gvimg = gvax.imshow(maps['gasvel'][0].T / vunits, 
                        extent=maps['gasvel'][1],
                        origin='lower', interpolation='nearest',
                        cmap=cmap_vel, vmin=vmin_vel, vmax=vmax_vel)
    
    nvax = fig.add_subplot(grid[2, 1])
    nvimg = nvax.imshow(maps['ne8vel'][0].T / vunits, 
                        extent=maps['ne8vel'][1],
                        origin='lower', interpolation='nearest',
                        cmap=cmap_vel, vmin=vmin_vel, vmax=vmax_vel)
    nvax.set_xlabel('pkpc', fontsize=fontsize)
    
    gtax = fig.add_subplot(grid[1, 2])
    gtimg = gtax.imshow(maps['gastemp'][0].T, 
                        extent=maps['gastemp'][1],
                        origin='lower', interpolation='nearest',
                        cmap=cmap_T, vmin=vmin_t, vmax=vmax_t)
    plt.colorbar(gtimg, cax=temp_cax, orientation='horizontal')
    temp_cax.set_xlabel('$\\log_{10} \\mathrm{T} \\; '
                        '[\\mathrm{K}]$',
                        fontsize=fontsize)
    
    ntax = fig.add_subplot(grid[2, 2])
    ntimg = ntax.imshow(maps['ne8temp'][0].T, 
                        extent=maps['ne8temp'][1],
                        origin='lower', interpolation='nearest',
                        cmap=cmap_T, vmin=vmin_t, vmax=vmax_t)
    ntax.set_xlabel('pkpc', fontsize=fontsize)
    
    plt.savefig(outname, bbox_inches='tight')
 