import h5py
import numbers as num
import numpy as np

import fire_an.mainfunc.get_qty as gq
import fire_an.mainfunc.haloprop as hp
import fire_an.readfire.readin_fire_data as rf
import fire_an.utils.constants_and_units as c
from fire_an.utils.projection import project

def savedict_hdf5(grp, dct):
    for key in dct:
        val = dct[key]
        if isinstance(val, type('')):
            val = np.string_(val)
        elif val is None:
            val = np.string_('None')
        elif isinstance(val, dict):
            sgrp = grp.create_group(key + '_dict')
            _val = val.copy()
            savedict_hdf5(sgrp, _val)
            val = np.string_('dict')
        grp.attrs.create(key, val)

def process_typeargs_coords(simpath, snapnum, typeargs,
                            paxis=None):
    '''
    allows standard values or gethalodata_shrinkinssphere/get_vcom
    arguments to be used instead of specifiying vcen_cmps, rcen_cm
    explicitly.

    Parameters
    ----------
    typeargs: dict
        may match any maptype_args for maptype 'coords' in get_qty
        additional convenience options (by key):
        'rcen_cm': not included or None
            not included or None: use the BN98 overdensity definition,
            shrinking-sphere method with all particle types (except 2)
            and standard parameter values for centering.
            Otherwise, should match the get_qty options.
        'vcen_cmps': not included, None, or dict
            centering/virial radius always match the defaults under 
            'rcen_cm', and the virial radius is always the 'BN98' one.
            The defaults for 'radius_rvir' and 'parttypes' are 1. and
            'all', respectively. If a dict argument is given, these
            defaults may be overridden by including the 'radius_rvir'
            and/or 'parttypes' keywords, with options matching those
            arguments for haloprop.get_vcom 
            Otherwise, should match the get_qty options.
        'pos': 'los'
            get the line-of-sight position along whatever the 
            projection axis is. (requires the paxis argument to be set
            to 0, 1, or 2)
        'vel': 'los'
            get the line-of-sight velocity along whatever the 
            projection axis is. (requires the paxis argument to be set
            to 0, 1, or 2)
    paxis: {0, 1, 2, None}
        the projection axis. Required if a 'pos' or 'vel' coordinate
        is specified as 'los', otherwise ignored.
    '''
    needsv = 'vel' in typeargs or \
             ('multiple' in typeargs 
              and np.any(['vel' in dct for dct in typeargs['multiple']]))
    outdoc = {}
    typeargs_out = typeargs.copy()
    if ('center_cm' not in typeargs) or \
       ('center_cm' in typeargs and typeargs['center_cm'] is None):
        rhdat, rdoc = hp.gethalodata_shrinkingsphere(simpath, snapnum, 
                                                     meandef='BN98')
        outdoc.update({'coords_' + key: rdoc[key] for key in rdoc})
        rcen_cm = np.array([rhdat['Xc_cm'], rhdat['Yc_cm'], rhdat['Zc_cm']])
        typeargs_out.update({'center_cm': rcen_cm})
        outdoc.update({'coords_rcen_cm_in': 'default'})
        outdoc.update({'coords_center': 'shrinksph'})
    if needsv:
        if 'vcen_cmps' not in typeargs:
            typeargs['vcen_cmps'] = dict()
        elif typeargs['vcen_cmps'] is None:
            typeargs['vcen_cmps'] = dict()
        if (isinstance(typeargs['vcen_cmps'], dict)):
            tvdct = typeargs['vcen_cmps'].copy()
            if 'radius_rvir' in tvdct:
                radius_rvir = tvdct['radius_rvir']
            else:
                radius_rvir = 1.
            outdoc.update({'vcen_radius_rvir': radius_rvir})
            if 'parttypes' in tvdct:
                parttypes = tvdct['parttypes']
            else:
                parttypes = 'all'
            outdoc.update({'vcen_parttypes': parttypes})
            vdat, vdoc = hp.get_vcom(simpath, snapnum, 
                                     radius_rvir, meandef_rvir='BN98',
                                     parttypes=parttypes)
            outdoc.update({'vcen_' + key: vdoc[key] for key in vdoc})
            vcen_cmps = np.array([vdat['VXcom_cmps'],
                                  vdat['VYcom_cmps'],
                                  vdat['VZcom_cmps']])
            typeargs_out.update({'vcen_cmps': vcen_cmps})
        else:
            vcin = typeargs['vcen_cmps']
            if not (hasattr(vcin, '__len__') 
                    and len(vcin) == 3
                    and np.all([isinstance(vcin[i], num.Number) 
                                for i in range(3)])):
                raise ValueError('The "vcen_cmps" argument should be a'
                                 ' length 3 iterable of floats, '
                                 'a dictionary, or None')
    if 'vel' in typeargs and typeargs['vel'] == 'los':
        outdoc.update({'coords_vel_in': 'los'})
        typeargs_out['vel'] = paxis
    elif 'pos' in typeargs and typeargs['pos'] == 'los':
        outdoc.update({'coords_pos_in': 'los'})
        typeargs_out['pos'] = paxis
    elif 'multiple' in typeargs:
        cspec = typeargs['multiple'].copy()
        for di, dct in enumerate(cspec):
            if 'vel' in dct and dct['vel'] == 'los':
                outdoc.update({'coords_vel_in': 'los'})
                cspec[di].update({'vel': paxis})
            elif 'pos' in dct and dct['pos'] == 'los':
                outdoc.update({'coords_pos_in': 'los'})
                cspec[di].update({'pos': paxis})
        typeargs_out['multiple'] = cspec
    return typeargs_out, outdoc



# AHF: sorta tested (enclosed 2D mass wasn't too far above Mvir)
# Rockstar: untested draft
# shrinking spheres: sort of tested (maps look right)
# centers, masses agree ok with Imran's (DM only) code
# mass maps: look ok
# ion/metal maps: tested sum of ions
def massmap(dirpath, snapnum, radius_rvir=2., particle_type=0,
            pixsize_pkpc=3., axis='z', outfilen=None,
            center='shrinksph', norm='pixsize_phys',
            maptype='Mass', maptype_args=None,
            weighttype=None, weighttype_args=None,
            save_weightmap=False, logmap=True,
            logweightmap=True):
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
    maptype: {'Mass', 'Metal', 'ion', 'coords', 'sim-direct'}
        what sort of thing to map
        'Mass' -> g
        'Metal' -> number of nuclei of the selected element
        'ion' -> number of ions of the selected type
        'coords' -> positions and velocities
        'sim-direct' -> arrays stored in FIRE outputs
    maptype_args: dict or None
        see get_qty for parameters; options depend on maptype
    weighttype: same options as masstype, or None
        instead of projecting maptype directly, get the 
        weighttype-weighted average of maptype along the line of sight
    weighttype_args: dict or None
        same as maptype_args, but for the weighttype (if any)
    save_weightmap: bool
        if obtaining a weighted quantity map, also save the projected
        weights? (If true, both maps are stored in the same hdf5 file.)
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
        Axis1 = 1
        Axis2 = 2
        Axis3 = 0
    elif axis == 'y':
        Axis1 = 2
        Axis2 = 0
        Axis3 = 1
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
    elif center == 'shrinksph':
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
        multipafter_norm = 1. / pixel_cm**2
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
    if weighttype is None:
        if maptype == 'coords':
            maptype_args, todocW = process_typeargs_coords(dirpath, snapnum,
                                                           maptype_args,
                                                           paxis=Axis3)
        else:
            todocW = {}
        qW, toCGSW, _todocW = gq.get_qty(snap, particle_type, maptype,
                                      maptype_args,
                                      filterdct={'filter': filter})
        todocW.update(_todocW)
        multipafterW = toCGSW * multipafter_norm
        qQ = np.zeros(len(qW), dtype=np.float32)
    else:
        if maptype == 'coords':
            maptype_args, todocQ = process_typeargs_coords(dirpath, snapnum,
                                                           maptype_args,
                                                           paxis=Axis3)
        else:
            todocQ = {}
        qQ, toCGSQ, _todocQ = gq.get_qty(snap, particle_type, maptype,
                                        maptype_args,
                                        filterdct={'filter': filter})
        multipafterQ = toCGSQ
        todocQ.update(_todocQ)
        if weighttype == 'coords':
            maptype_args, todocW = process_typeargs_coords(dirpath, snapnum,
                                                           weighttype_args,
                                                           paxis=Axis3)
        else:
            todocW = {}
        qW, toCGSW, _todocW = gq.get_qty(snap, particle_type, weighttype,
                                         weighttype_args,
                                         filterdct={'filter': filter})
        todocW.update(_todocW)
        multipafterW = toCGSW * multipafter_norm
        
    if not haslsmooth:
        # minimum smoothing length is set in the projection
        lsmooth = np.zeros(shape=(len(qW),), dtype=coords.dtype)
        lsmooth_toCGS = 1.
    
    tree = False
    periodic = False # zoom region
    NumPart = len(qW)
    dct = {'coords': coords, 'lsmooth': lsmooth, 
           'qW': qW, 
           'qQ': qQ}
    Ls = box_dims_coordunit
    # cosmopars uses EAGLE-style cMpc/h units for the box
    box3 = [snap.cosmopars.boxsize * c.cm_per_mpc / snap.cosmopars.h \
            / coords_toCGS] * 3
    mapW, mapQ = project(NumPart, Ls, Axis1, Axis2, Axis3, box3,
                         periodic, npix_x, npix_y,
                         'C2', dct, tree, ompproj=True, 
                         projmin=None, projmax=None)
    if weighttype is None:
        if logmap:
            omapW = np.log10(mapW)
            omapW += np.log10(multipafterW)
        else:
            mapW *= multipafterW
            omapW = mapW
        mmap = omapW
        mdoc = todocW
        mlog = logmap
    else:
        if logmap:
            omapQ = np.log10(mapQ)
            omapQ += np.log10(multipafterQ)
        else:
            mapQ *= multipafterQ
            omapQ = mapQ
        if save_weightmap:
            if logweightmap:
                omapW = np.log10(mapW)
                omapW += np.log10(multipafterW)
            else:
                mapW *= multipafterW
                omapW = mapW
        mmap = omapQ
        mdoc = todocQ
        mlog = logmap

    with h5py.File(outfilen, 'w') as f:
        # map (emulate make_maps format)
        f.create_dataset('map', data=mmap)
        f['map'].attrs.create('log', mlog)
        if mlog:
            minfinite = np.min(mmap[np.isfinite(mmap)])
            f['map'].attrs.create('minfinite', minfinite)
        else:
            f['map'].attrs.create('min', np.min(mmap))
        f['map'].attrs.create('max', np.max(mmap))
        
        # cosmopars (emulate make_maps format)
        hed = f.create_group('Header')
        cgrp = hed.create_group('inputpars/cosmopars')
        csm = snap.cosmopars.getdct()
        savedict_hdf5(cgrp, csm)
        
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
        savedict_hdf5(_grp, halodat)
        igrp.attrs.create('maptype', np.string_(maptype))
        if maptype_args is None:
            igrp.attrs.create('maptype_args', np.string_('None'))
        else:
            igrp.attrs.create('maptype_args', np.string_('dict'))
            _grp = igrp.create_group('maptype_args_dict')
            savedict_hdf5(_grp, maptype_args)
        if weighttype is None and 'units' in mdoc:
            mdoc['units'] = mdoc['units'] + norm_units
        savedict_hdf5(igrp, mdoc)
        if weighttype is None:
            igrp.attrs.create('weighttype', np.string_('None'))
            igrp.attrs.create('weighttype_args', np.string_('None'))
        else:
            igrp.attrs.create('weighttype', np.string_(weighttype))
            if weighttype_args is None:
                igrp.attrs.create('weighttype_args', np.string_('None'))
            else:
                igrp.attrs.create('weighttype_args', np.string_('dict'))
                _grp = igrp.create_group('weighttype_args_dict')
                savedict_hdf5(_grp, weighttype_args)
                if 'units' in todocW:
                    todocW['units'] = todocW['units'] + norm_units
                savedict_hdf5(igrp, todocW)
            if save_weightmap:
                f.create_dataset('weightmap', data=omapW)
                f['weightmap'].attrs.create('log', logweightmap)
                if logweightmap:
                    minfinite = np.min(omapW[np.isfinite(omapW)])
                    f['weightmap'].attrs.create('minfinite', minfinite)
                else:
                    f['weightmap'].attrs.create('min', np.min(omapW))
                f['weightmap'].attrs.create('max', np.max(omapW))


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
                savedict_hdf5(cgrp, csm)
        
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