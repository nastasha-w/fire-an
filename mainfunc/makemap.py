import h5py
import numpy as np

import fire_an.mainfunc.get_qty as gq
import fire_an.mainfunc.haloprop as hp
import fire_an.readfire.readin_fire_data as rf
import fire_an.utils.constants_and_units as c
from fire_an.utils.projection import project

# AHF: sorta tested (enclosed 2D mass wasn't too far above Mvir)
# Rockstar: untested draft
# shrinking spheres: sort of tested (maps look right)
# mass maps: look ok
# ion/metal maps: tested sum of ions
def massmap(dirpath, snapnum, radius_rvir=2., particle_type=0,
            pixsize_pkpc=3., axis='z', outfilen=None,
            center='shrinksph', norm='pixsize_phys',
            maptype='Mass', maptype_args=None):
    '''
    Creates a mass map projected perpendicular to a line of sight axis
    by assuming the simulation resolution elements divide their mass 
    following a C2 SPH kernel.

    Parameters:
    -----------
    snapfile: str
        file (or example file, if split) containing the snapshot data
    dirpath: str
        path to the directory containing the 'output' directory with the
        snapshots
    snapnum: int
        snapshot number
    radius_rvir: float 
        radius of the cube to project in units of Rvir. Note that in the sky 
        plane, this will be (slightly) extended to get the exact pixel size.
    particle_type: int
        particle type to project (follows FIRE format)
    pixsize_pkpc: float
        size of the map pixels in proper kpc
    axis: str, 'x', 'y', or 'z'
        axis corresponding to the line of sight 
    outfilen: str or None. 
        if a string, the name of the file to save the output data to. The
        default is None, meaning the maps are returned as output
    center: str
        how to find the halo center.
        'AHFsmooth': use halo_00000_smooth.dat from AHF 
        'rockstar-maxmass': highest mass halo at snapshot from Rockstar
        'rockstar-mainprog': main progenitor of most massive halo at
                           final snapshot from Rockstar
        'rockstar-<int>': halo with snapshot halo catalogue index <int>
                          from Rockstar 
        'shrinksph': Imran's shrinking spheres method
    norm: {'pixsize_phys'}
        how to normalize the column values 
        'pixsize_phys': [quantity] / cm**2
    maptype: {'Mass', 'Metal', 'ion'}
        what sort of thing to map
        'Mass' -> g
        'Metal' -> number of nuclei of the selected element
        'ion' -> number of ions of the selected type
    maptype_args: dict or None
        see get_qty for parameters; options depend on maptype
    Output:
    -------
    massW: 2D array of floats
        projected mass image [log g/cm^-2]
    massQ: NaN array, for future work


    '''
    if axis == 'z':
        Axis1 = 0
        Axis2 = 1
        Axis3 = 2
    elif axis == 'x':
        Axis1 = 2
        Axis2 = 0
        Axis3 = 1
    elif axis == 'y':
        Axis1 = 1
        Axis2 = 2
        Axis3 = 0
    else:
        msg = 'axis should be "x", "y", or "z", not {}'
        raise ValueError(msg.format(axis))
    
    if center == 'AHFsmooth':
        halodat = hp.mainhalodata_AHFsmooth(dirpath, snapnum)
        snap = rf.get_Firesnap(dirpath, snapnum) 
        cen = np.array([halodat['Xc_ckpcoverh'], 
                        halodat['Yc_ckpcoverh'], 
                        halodat['Zc_ckpcoverh']])
        cen_cm = cen * snap.cosmopars.a * 1e-3 * c.cm_per_mpc \
                 / snap.cosmopars.h
        rvir_cm = halodat['Rvir_ckpcoverh'] * snap.cosmopars.a \
                  * 1e-3 * c.cm_per_mpc / snap.cosmopars.h
    elif center.startswith('rockstar'):
        select = center.split('-')[-1]
        if select not in ['maxmass', 'mainprog']:
            try:
                select = int(select)
            except ValueError:
                msg = 'invalid option for center: {}'.format(center)
                raise ValueError(msg)
        halodat, _csm_halo = hp.halodata_rockstar(dirpath, snapnum, 
                                                  select=select)
        snap = rf.get_Firesnap(dirpath, snapnum) 
        cen = np.array([halodat['Xc_ckpc'], 
                        halodat['Yc_ckpc'], 
                        halodat['Zc_ckpc']])
        cen_cm = cen * snap.cosmopars.a * 1e-3 * c.cm_per_mpc 
        rvir_cm = halodat['Rvir_cm'] 
    elif center ==  'shrinksph':
        halodat, _ = hp.gethalodata_shrinkingsphere(dirpath, snapnum,
                                                    meandef='BN98')
        cen_cm = np.array([halodat['Xc_cm'], 
                           halodat['Yc_cm'], 
                           halodat['Zc_cm']])
        rvir_cm = halodat['Rvir_cm']
        snap = rf.get_Firesnap(dirpath, snapnum) 
    else:
        raise ValueError('Invalid center option {}'.format(center))

    # calculate pixel numbers and projection region based
    # on target size and extended for integer pixel number
    target_size_cm = np.array([2. * radius_rvir * rvir_cm] * 3)
    pixel_cm = pixsize_pkpc * c.cm_per_mpc * 1e-3
    npix3 = (np.ceil(target_size_cm / pixel_cm)).astype(int)
    npix_x = npix3[Axis1]
    npix_y = npix3[Axis2]
    size_touse_cm = target_size_cm
    size_touse_cm[Axis1] = npix_x * pixel_cm
    size_touse_cm[Axis2] = npix_y * pixel_cm

    if norm == 'pixsize_phys':
        multipafter = 1. / pixel_cm**2
        norm_units = ' / (physical cm)**2'
    else:
        raise ValueError('Invalid norm option {}'.format(norm))

    basepath = 'PartType{}/'.format(particle_type)
    haslsmooth = particle_type == 0
    if haslsmooth: # gas
        lsmooth = snap.readarray_emulateEAGLE(basepath + 'SmoothingLength')
        lsmooth_toCGS = snap.toCGS

    coords = snap.readarray_emulateEAGLE(basepath + 'Coordinates')
    coords_toCGS = snap.toCGS
    # needed for projection step anyway
    coords -= cen_cm / coords_toCGS
    # select box region
    # zoom regions are generally centered -> don't worry
    # about edge overlap
    box_dims_coordunit = size_touse_cm / coords_toCGS

    if haslsmooth:
        # extreme values will occur at zoom region edges -> restrict
        filter_temp = np.all(np.abs((coords)) <= 0.5 * box_dims_coordunit, 
                             axis=1)
        lmax = np.max(lsmooth[filter_temp]) 
        conv = lsmooth_toCGS / coords_toCGS
        del filter_temp
        # might be lower-density stuff outside the region, but overlapping it
        lmargin = 2. * lmax * conv
        filter = np.all(np.abs((coords)) <= 0.5 * box_dims_coordunit \
                        + lmargin, axis=1)
        lsmooth = lsmooth[filter]
        if not np.isclose(conv, 1.):
            lsmooth *= conv
    
    else:
        filter = np.all(np.abs((coords)) <= 0.5 * box_dims_coordunit, axis=1)   
    
    coords = coords[filter]
    qW, toCGS, todoc = gq.get_qty(snap, particle_type, maptype, maptype_args,
                                  filterdct={'filter': filter})
    multipafter *= toCGS
    ## debugging: check for NaN values
    #naninds = np.where(np.isnan(qW))[0]
    #if len(naninds) > 0:
    #    print('Some qW values are NaN')
    #    print('Used {}, {}, {}'.format(particle_type, maptype, maptype_args))
    #    outfile_debug = outfilen.split('/')[:-1]
    #    outfile_debug.append('debug_qW_naninfo.hdf5')
    #    outfile_debug = '/'.join(outfile_debug)
    #
    #    numnan = len(naninds)
    #    minW = np.min(qW[np.isfinite(qW)])
    #    maxW = np.max(qW[np.isfinite(qW)])
    #    with h5py.File(outfile_debug, 'w') as f:
    #        hed = f.create_group('Header')
    #        hed.attrs.create('number of qW values', len(qW))
    #        hed.attrs.create('number of NaN values', numnan)
    #        hed.attrs.create('number of inf values', np.sum(qW == np.inf))
    #        hed.attrs.create('number of -inf values', np.sum(qW == -np.inf))
    #        hed.attrs.create('number of 0 values', np.sum(qW == 0))
    #        hed.attrs.create('number of values < 0', np.sum(qW < 0))
    #        hed.attrs.create('number of values > 0', np.sum(qW > 0))
    #        hed.attrs.create('qW_toCGS', toCGS)
    #        hed.attrs.create('multipafter', multipafter)
    #
    #        if minW > 0:
    #            bins = np.logspace(np.log10(minW), np.log10(maxW), 100)
    #        else:
    #            bins = np.linspace(minW, maxW, 100)
    #        # NaN, inf, -inf values just aren't counted
    #        hist, _ = np.histogram(qW, bins=bins)
    #
    #        f.create_dataset('qW_hist', data=hist)
    #        f.create_dataset('qW_hist_bins', data=bins)
    #
    #        # issues were for ion columns; check rho, T, Z of NaN values
    #        _filter = filter.copy()
    #        _filter[_filter] = np.isnan(qW)
    #
    #        _temp, _temp_toCGS, _temp_todoc = get_qty(snap, particle_type, 
    #                'sim-direct', {'field': 'Temperature'}, 
    #                filterdct={'filter': _filter})
    #        ds = f.create_dataset('Temperature_nanqW', data=_temp)
    #        ds.attrs.create('toCGS', _temp_toCGS)
    #        print('Temperature: ', _temp_todoc)
    #
    #        _dens, _dens_toCGS, _dens_todoc = get_qty(snap, particle_type, 
    #                'sim-direct', {'field': 'Density'}, 
    #                filterdct={'filter': _filter})
    #        ds = f.create_dataset('Density_nanqW', data=_dens)
    #        ds.attrs.create('toCGS', _dens_toCGS)
    #        print('Density: ', _dens_todoc)
    #
    #        _hden, _hden_toCGS, _emet_todoc = get_qty(snap, particle_type, 
    #                'Metal', {'element': 'Hydrogen', 'density': True}, 
    #                filterdct={'filter': _filter})
    #        ds = f.create_dataset('nH_nanqW', data=_hden)
    #        ds.attrs.create('toCGS', _hden_toCGS)
    #        print('Hydrogen number density: ', _emet_todoc)
    #
    #        if maptype == 'ion':
    #            ion = maptype_args['ion']
    #            dummytab = Linetable_PS20(ion, snap.cosmopars.z, 
    #                                      emission=False,
    #                                      vol=True, lintable=True)
    #            element = dummytab.element
    #            eltpath = basepath + 'ElementAbundance/' +\
    #                      string.capwords(element)
    #            _emet, _emet_toCGS, _emet_todoc = get_qty(snap, particle_type, 
    #                    'sim-direct', {'field': eltpath}, 
    #                    filterdct={'filter': _filter})
    #            ds = f.create_dataset('massfrac_{}_nanqW', data=_emet)
    #            ds.attrs.create('toCGS', _emet_toCGS)
    #            print(f'{element} mass fraction: ', _emet_todoc)
    #else:
    #    print('No NaN values in qW')
    # stars, black holes. DM: should do neighbour finding. Won't though.
    if not haslsmooth:
        # minimum smoothing length is set in the projection
        lsmooth = np.zeros(shape=(len(qW),), dtype=coords.dtype)
        lsmooth_toCGS = 1.
    
    tree = False
    periodic = False # zoom region
    NumPart = len(qW)
    dct = {'coords': coords, 'lsmooth': lsmooth, 
           'qW': qW, 
           'qQ': np.zeros(len(qW), dtype=np.float32)}
    Ls = box_dims_coordunit
    # cosmopars uses EAGLE-style cMpc/h units for the box
    box3 = [snap.cosmopars.boxsize * c.cm_per_mpc / snap.cosmopars.h \
            / coords_toCGS] * 3
    mapW, mapQ = project(NumPart, Ls, Axis1, Axis2, Axis3, box3,
                         periodic, npix_x, npix_y,
                         'C2', dct, tree, ompproj=True, 
                         projmin=None, projmax=None)
    lmapW = np.log10(mapW)
    ## debug NaN values in maps
    #if np.any(np.isnan(mapW)):
    #    print('NaN values in mapW after projection')
    #if np.any(mapW < 0.):
    #    print('values < 0 in mapW after projection')
    #if np.any(np.isnan(lmapW)):
    #    print('NaN values in log mapW before multipafter')

    lmapW += np.log10(multipafter)

    #if np.any(np.isnan(lmapW)):
    #    print('NaN values in log mapW after multipafter')
    #if outfilen is None:
    #    return lmapW, mapQ
    
    with h5py.File(outfilen, 'w') as f:
        # map (emulate make_maps format)
        f.create_dataset('map', data=lmapW)
        f['map'].attrs.create('log', True)
        minfinite = np.min(lmapW[np.isfinite(lmapW)])
        f['map'].attrs.create('minfinite', minfinite)
        f['map'].attrs.create('max', np.max(lmapW))
        
        # cosmopars (emulate make_maps format)
        hed = f.create_group('Header')
        cgrp = hed.create_group('inputpars/cosmopars')
        csm = snap.cosmopars.getdct()
        for key in csm:
            cgrp.attrs.create(key, csm[key])
        
        # direct input parameters
        igrp = hed['inputpars']
        igrp.attrs.create('snapfiles', np.array([np.string_(fn) for fn in snap.filens]))
        igrp.attrs.create('dirpath', np.string_(dirpath))
        igrp.attrs.create('radius_rvir', radius_rvir)
        igrp.attrs.create('particle_type', particle_type)
        igrp.attrs.create('pixsize_pkpc', pixsize_pkpc)
        igrp.attrs.create('axis', np.string_(axis))
        igrp.attrs.create('norm', np.string_(norm))
        igrp.attrs.create('outfilen', np.string_(outfilen))
        # useful derived/used stuff
        igrp.attrs.create('Axis1', Axis1)
        igrp.attrs.create('Axis2', Axis2)
        igrp.attrs.create('Axis3', Axis3)
        igrp.attrs.create('diameter_used_cm', np.array(size_touse_cm))
        if haslsmooth:
            igrp.attrs.create('margin_lsmooth_cm', lmargin * coords_toCGS)
        igrp.attrs.create('center', np.string_(center))
        _grp = igrp.create_group('halodata')
        for key in halodat:
            _grp.attrs.create(key, halodat[key])
        igrp.attrs.create('maptype', np.string_(maptype))
        if maptype_args is None:
            igrp.attrs.create('maptype_args', np.string_('None'))
        else:
            igrp.attrs.create('maptype_args', np.string_('dict'))
            _grp = igrp.create_group('maptype_args_dict')
            for key in maptype_args:
                val = maptype_args[key]
                if isinstance(val, type('')):
                    val = np.string_(val)
                _grp.attrs.create(key, val)
        for key in todoc:
            if key == 'units':
                val = todoc[key]
                val = val + norm_units
            else:
                val = todoc[key]
            if isinstance(val, type('')):
                val = np.string_(val)
            igrp.attrs.create(key, val)

def massmap_wholezoom(dirpath, snapnum, pixsize_pkpc=3.,
                      outfilen_DM='map_DM_{ax}-axis.hdf5',
                      outfilen_gas='map_gas_{ax}-axis.hdf5',
                      outfilen_stars='map_stars_{ax}-axis.hdf5',
                      outfilen_BH='map_BH_{ax}-axis.hdf5'):
    '''
    for debugging: make a mass map of basically the whole zoom region
    (for centering tests)
    '''
    parttype_outfilen = {0: outfilen_gas,
                         1: outfilen_DM,
                         4: outfilen_stars,
                         5: outfilen_BH}
    snap = rf.get_Firesnap(dirpath, snapnum)
    coords = {}
    mass = {}
    lsmooth = {}
    masspath = 'PartType{}/Masses'
    coordpath = 'PartType{}/Coordinates'
    lsmoothpath = 'PartType{}/SmoothingLength'
    coords_toCGS = None
    lsmooth_toCGS = None
    mass_toCGS = None
    maxl = -np.inf
    coordsbox = np.array([[np.inf, -np.inf]] * 3) 
    for pt in parttype_outfilen:
        try:
            coords[pt] = snap.readarray(coordpath.format(pt))
            _toCGS = snap.toCGS
            if coords_toCGS is None:
                coords_toCGS = _toCGS
            elif coords_toCGS != _toCGS:
                msg = ('Different particle types have different coordinate'
                       ' toCGS: {}, {}')
                raise RuntimeError(msg.format(coords_toCGS, _toCGS))
        except rf.FieldNotFoundError as err:
            print('PartType {} not found'.format(pt))
            print(err)
            continue
        mass[pt] = snap.readarray(masspath.format(pt))
        _toCGS = snap.toCGS
        if mass_toCGS is None:
            mass_toCGS = _toCGS
        elif mass_toCGS != _toCGS:
            msg = 'Different particle types have different mass'+\
                    ' toCGS: {}, {}'
            raise RuntimeError(msg.format(mass_toCGS, _toCGS))
        coordsbox[:, 0] = np.min([coordsbox[:, 0], 
                                  np.min(coords[pt], axis=0)], 
                                 axis=0)
        coordsbox[:, -1] = np.max([coordsbox[:, -1], 
                                  np.max(coords[pt], axis=0)], 
                                 axis=0)
        if pt == 0:
            lsmooth[pt] = snap.readarray(lsmoothpath.format(pt))
            _toCGS = snap.toCGS
            if lsmooth_toCGS is None:
                lsmooth_toCGS = _toCGS
            elif lsmooth_toCGS != _toCGS:
                msg = ('Different particle types have different smoothing '
                       'length toCGS: {}, {}')
                raise RuntimeError(msg.format(lsmooth_toCGS, _toCGS))
            maxl = max(maxl, np.max(lsmooth[pt]))
    coordsbox[:, 0] -= maxl * lsmooth_toCGS / coords_toCGS
    coordsbox[:, -1] += maxl * lsmooth_toCGS / coords_toCGS
    print('coordsbox before pixel adjustments: ', coordsbox)

    target_size_cm = (coordsbox[:, -1] - coordsbox[:, 0]) * coords_toCGS
    pixel_cm = pixsize_pkpc * c.cm_per_mpc * 1e-3
    npix3 = (np.ceil(target_size_cm / pixel_cm)).astype(int)
    center = 0.5 * np.sum(coordsbox, axis=1)
    Ls = npix3 * pixel_cm / coords_toCGS
    coordsbox = center[:, np.newaxis] \
                + Ls[:, np.newaxis] * np.array([-0.5, + 0.5])[np.newaxis, :]
    print('coordsbox: ', coordsbox)

    multipafter = mass_toCGS / pixel_cm**2
    units = 'g / (physical cm)**2'
    
    for pt in coords:
        print('Running particle type ', pt)
        qW = mass[pt]
        if pt in lsmooth:
            _lsmooth = lsmooth[pt]
            _lsmooth_toCGS = lsmooth_toCGS
        else:
            # minimum smoothing length is set in the projection
            _lsmooth = np.zeros(shape=(len(qW),), dtype=(coords[pt]).dtype)
            _lsmooth_toCGS = 1.
        tree = False
        periodic = False # zoom region
        NumPart = len(qW)
        dct = {'coords': coords[pt] - center, 'lsmooth': _lsmooth, 
               'qW': qW, 
               'qQ': np.zeros(len(qW), dtype=np.float32)}
        print('Extent of coordinates: ',
              np.min(coords[pt], axis=0), 
              ', ',
              np.max(coords[pt], axis=0))
        print('Extent of centered coordinates: ',
              np.min(dct['coords'], axis=0), 
              ', ',
              np.max(dct['coords'], axis=0))
        print('Ls: ', Ls)
        print('Coordinates in box: ', 
              np.sum(np.all(np.abs(dct['coords']) < Ls, axis=1)),
              ' / ', len(dct['coords']))

        # cosmopars uses EAGLE-style cMpc/h units for the box
        box3 = [snap.cosmopars.boxsize * c.cm_per_mpc / snap.cosmopars.h \
                / coords_toCGS] * 3
        for axis in ['x', 'y', 'z']:
            if axis == 'z':
                Axis1 = 0
                Axis2 = 1
                Axis3 = 2
            elif axis == 'x':
                Axis1 = 2
                Axis2 = 0
                Axis3 = 1
            elif axis == 'y':
                Axis1 = 1
                Axis2 = 2
                Axis3 = 0
            npix_x = npix3[Axis1]
            npix_y = npix3[Axis2]

            mapW, mapQ = project(NumPart, Ls, Axis1, Axis2, Axis3, box3,
                                 periodic, npix_x, npix_y,
                                 'C2', dct, tree, ompproj=True, 
                                 projmin=None, projmax=None)
            lmapW = np.log10(mapW)
            lmapW += np.log10(multipafter)
        
            outfilen = (parttype_outfilen[pt]).format(ax=axis)
    
            with h5py.File(outfilen, 'w') as f:
                # map (emulate make_maps format)
                f.create_dataset('map', data=lmapW)
                f['map'].attrs.create('log', True)
                minfinite = np.min(lmapW[np.isfinite(lmapW)])
                f['map'].attrs.create('minfinite', minfinite)
                f['map'].attrs.create('max', np.max(lmapW))
            
                # cosmopars (emulate make_maps format)
                hed = f.create_group('Header')
                cgrp = hed.create_group('inputpars/cosmopars')
                csm = snap.cosmopars.getdct()
                for key in csm:
                    cgrp.attrs.create(key, csm[key])
        
                # direct input parameters
                igrp = hed['inputpars']
                igrp.attrs.create('snapfiles', np.array([np.string_(fn) \
                                  for fn in snap.filens]))
                igrp.attrs.create('dirpath', np.string_(dirpath))
            
                igrp.attrs.create('particle_type', pt)
                igrp.attrs.create('pixsize_pkpc', pixsize_pkpc)
                igrp.attrs.create('axis', np.string_(axis))
                igrp.attrs.create('units', np.string_(units))
                igrp.attrs.create('outfilen', np.string_(outfilen))
                # useful derived/used stuff
                igrp.attrs.create('Axis1', Axis1)
                igrp.attrs.create('Axis2', Axis2)
                igrp.attrs.create('Axis3', Axis3)
                igrp.attrs.create('maptype', np.string_('Mass'))
                igrp.attrs.create('mapped_region_cm', coordsbox)
                igrp.attrs.create('maptype_args', np.string_('None'))
                # useful for plotting centers in sim units
                igrp.attrs.create('coords_toCGS', coords_toCGS)