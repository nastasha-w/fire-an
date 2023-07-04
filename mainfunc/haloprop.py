
import glob
import h5py
import numpy as np
import os
import pandas as pd
import sys
import uuid

import matplotlib.pyplot as plt # debugging

import fire_an.simlists as sl


# Andrew Wetzel's Rockstar halo catalogue wrangler
try:
    import halo_analysis as ha
except ModuleNotFoundError:
    msg = ('Could not import module "halo_analysis";'
           ' Rockstar halo data read-in will fail.')
    print(msg)

import fire_an.readfire.readin_fire_data as rf
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu
import fire_an.utils.math_utils as mu
import fire_an.utils.opts_locs as ol



# setup from the internets
class NoStoredMatchError(Exception):
    def __init__(self, *args):
        if len(args) > 0:
            self.message = args[0]
        else:
            self.message = None
    def __str__(self):
        if self.message is not None:
            return 'NoStoredMatchError: {0} '.format(self.message)
        else:
            return 'NoStoredMatchError'


# seems to work for at least one halo 
# (m13 guinea pig at snapshot 27, comparing image to found center)
def calchalocen(coordsmassesdict, shrinkfrac=0.025, minparticles=1000, 
                initialradiusfactor=1.):
    '''
    from: https://github.com/isulta/massive-halos/blob/d2dc0dd3649f359c0cea7191bfefd11b3498eeda/scripts/halo_analysis_scripts.py#L164 
    Imran Sultan's method, citing Power et al. (2003):
    their parameter values: 
    shrinkpercent=2.5, minparticles=1000, initialradiusfactor=1

    '''
    coords = coordsmassesdict['coords']
    masses = coordsmassesdict['masses']
    totmass = np.sum(masses)
    com = np.sum(coords * masses[:, np.newaxis], axis=0) / totmass
    r2 = np.sum((coords - com[np.newaxis, :])**2, axis=1)
    searchrad2 = initialradiusfactor**2 * np.max(r2)
    Npart_conv = min(minparticles, len(masses) * 0.01)
 
    it = 0
    coords_it = coords.copy()
    masses_it = masses.copy()
    comlist = [com]
    radiuslist = [np.sqrt(searchrad2)]
    while len(masses_it) > Npart_conv:
        searchrad2 *= (1. - shrinkfrac)**2
        mask = r2 <= searchrad2
        coords_it = coords_it[mask]
        masses_it = masses_it[mask]
        com = np.sum(coords_it * masses_it[:, np.newaxis], axis=0) \
               / np.sum(masses_it)
        r2 = np.sum((coords_it - com[np.newaxis, :])**2, axis=1)

        it += 1
        comlist.append(com)
        radiuslist.append(np.sqrt(searchrad2))
    return com, comlist, radiuslist

# centering seems to work for at least one halo 
# (m13 guinea pig at snapshot 27, comparing image to found center)
def calchalodata_shrinkingsphere(path, snapshot, meandef=('200c', 'BN98')):
    '''
    Using Imran Sultan's shrinking spheres method, calculate the halo 
    center, then find the halo mass and radius for a given overdensity
    citerion

    Parameters:
    -----------
    path: str
        path containing the 'output' directory or the snapshot
        files/directories for the chosen simulation
    snapshot: int
        snapshot number
    meandef: str
        overdensity definition for the halo
        'BN98': Bryan & Norman 1998 fitting formula
        '<float>c': <float> times the critical density at the snapshot 
                    redshift
        '<float>m': <float> times the mean matter density at the 
                    snapshot redshift
        tuple of values -> return a list of Mvir and Rvir, in same order
    
    Returns:
    --------
    outdct: dict
        contains 
        'Xc_cm', 'Yc_cm', 'Zc_cm': floats 
            the coordinates of the halo center 
        'Rvir_cm': float or list of floats
            the virial radius (radii) according to the halo overdensity
            criterion. float or list matches string or iterable choice
            for the overdensity definition
        'Mvir_cm': float or list of floats
            the virial mass (masses) according to the halo overdensity
            criterion. float or list matches string or iterable choice
            for the overdensity definition
     todoc: dict
        contains information on parameter values and particle types used   
    
    '''
    minparticles = 1000
    minpart_halo = 1000
    snap = rf.get_Firesnap(path, snapshot)
    todoc = {}

    # get mass and coordinate data
    # use all zoom region particle types
    parttypes = [0, 1, 4, 5]
    dct_m = {}
    dct_c = {}
    toCGS_m = None
    toCGS_c = None
    for pt in parttypes:
        cpath = 'PartType{}/Coordinates'
        mpath = 'PartType{}/Mass'
        try:
            dct_c[pt] = snap.readarray_emulateEAGLE(cpath.format(pt))
            _toCGS_c = snap.toCGS
            dct_m[pt] = snap.readarray_emulateEAGLE(mpath.format(pt))
            _toCGS_m = snap.toCGS
        except (OSError, rf.FieldNotFoundError):
            msg = 'Skipping PartType {} in center calc: not present on file'
            print(msg.format(pt))
            continue
        if toCGS_m is None:
            toCGS_m = _toCGS_m
        elif not np.isclose(toCGS_m, _toCGS_m):
                msg = ('Different particle type masses have different'
                       ' CGS conversions in ' + snap.firstfilen)
                raise RuntimeError(msg)
        if toCGS_c is None:
            toCGS_c = _toCGS_c
        elif not np.isclose(toCGS_c, _toCGS_c):
                msg = ('Different particle type coordinates have different'
                       ' CGS conversions in ' + snap.firstfilen)
                raise RuntimeError(msg)
    pt_used = list(dct_m.keys())
    pt_used.sort()
    totlen = sum([len(dct_m[pt]) for pt in pt_used])
    masses = np.empty((totlen,), dtype=dct_m[pt_used[0]].dtype)
    coords = np.empty((totlen, dct_c[pt_used[0]].shape[1]), 
                      dtype=dct_c[pt_used[0]].dtype)
    todoc['parttypes_used'] = tuple(pt_used)
    start = 0
    for pt in pt_used:
        partlen = len(dct_m[pt])
        masses[start: start + partlen] = dct_m[pt]
        coords[start: start + partlen] = dct_c[pt]
        start += partlen

        del dct_m[pt]
        del dct_c[pt]
    coordsmassdict = {'masses': masses, 'coords': coords}
    com_simunits, comlist, radiuslist = \
        calchalocen(coordsmassdict, shrinkfrac=0.025, 
                    minparticles=minparticles, initialradiusfactor=1.)
    print('Found center of mass [sim units]: {}'.format(com_simunits))
    todoc.update({'shrinkfrac': 0.025, 
                  'minparticles': minparticles, 
                  'initialradiusfactor': 1.,
                  'minpart_halo': minpart_halo})
    # find Rvir/Mvir
    cosmopars = snap.cosmopars.getdct()
    todoc['cosmopars'] = cosmopars
    if isinstance(meandef, type('')):
        outputsingle = True
        dens_targets_cgs = [cu.getmeandensity(meandef, cosmopars)]
    else:
        outputsingle = False
        dens_targets_cgs = [cu.getmeandensity(md, cosmopars) 
                            for md in meandef]
        
    r2 = np.sum((coords - com_simunits[np.newaxis, :])**2, axis=1)
    del coords
    rorder = np.argsort(r2)
    r2_order = r2[rorder]
    # apparent truncation error issues in cumsum for some 
    # simulations/snapshots using float32. (enclosed mass plateaus)
    masses_order = np.asarray(masses[rorder], dtype=np.float64)
    del masses, r2, rorder
    dens_targets = [target / toCGS_m * toCGS_c**3 for target in \
                    dens_targets_cgs]
    dens2_order = (np.cumsum(masses_order)**2) \
                  / ((4. * np.pi / 3)**2 * r2_order**3)

    rsols_cgs = []
    msols_cgs = []
    xydct = {'x': r2_order, 'y': dens2_order}
    for dti, dens_target in enumerate(dens_targets):
        sols = mu.find_intercepts(None, None, dens_target**2, xydct=xydct)
        # no random low-density holes or anything
        sols = sols[sols >= r2_order[minpart_halo]]
        if len(sols) == 0:
            msg = 'No solutions found for density {}'.format(meandef[dti])
            print(msg)
        elif len(sols) == 1:
            rsol = np.sqrt(sols[0])
            msol = 4. * np.pi / 3. * rsol**3 * dens_target
            rsols_cgs.append(rsol * toCGS_c)
            msols_cgs.append(msol * toCGS_m)
            print(f'Found solution {dti}; r: {rsol}, m:{msol}')
        else:
            # technically a solution, but there will be some 
            # particle noise; smoothing?
            # on the other hand, those effects are proabably tiny
            sols_kpc = np.sqrt(sols) * toCGS_c / (1e-3 * c.cm_per_mpc)
            print('Found radius solution options [pkpc] {}'.format(sols_kpc))
            print('Selected first in list')
            rsol = np.sqrt(sols[0])
            msol = 4. * np.pi / 3. * rsol**3 * dens_target
            rsols_cgs.append(rsol * toCGS_c)
            msols_cgs.append(msol * toCGS_m)
    com_cgs = com_simunits * toCGS_c
    if outputsingle:
        rsols_cgs = rsols_cgs[0]
        msols_cgs = msols_cgs[0]
    outdct = {'Xc_cm': com_cgs[0], 'Yc_cm': com_cgs[1], 'Zc_cm': com_cgs[2],
              'Rvir_cm': rsols_cgs, 'Mvir_g': msols_cgs}
    return  outdct, todoc

def readhalodata_shrinkingsphere(path, snapshot, meandef=('200c', 'BN98')):
    '''
    returns the halo data, or a NoStoredMatchError if it is not in the
    main halo data file
    '''
    # this must contain *all* todoc entries from the previous function
    # except 'parttypes_used', which is assumed to be everything but
    # PartType2 (lo-res DM) present in the NumPart_Total table
    usedvals_calchalo = {'shrinkfrac': 0.025, 
                         'minparticles': 1000., 
                         'initialradiusfactor': 1.,
                         'minpart_halo': 1000.}
    snap = rf.get_Firesnap(path, snapshot, filetype='snap')
    with h5py.File(snap.firstfilen) as f:
        pts = list(f['Header'].attrs['NumPart_Total'])
    pts = [ind for ind in range(len(pts)) if pts[ind] > 0]
    pts.remove(2)
    pts.sort()
    usedvals_calchalo['parttypes_used'] = tuple(pts)
                                                
    filen_main = ol.filen_halocenrvir
    
    #pparts = path.split('/')
    #while '' in pparts:
    #    pparts.remove('')
    #if pparts[-1] == 'output':
    #    pparts = pparts[:-1]
    #simid = pparts[-1]
    simid = sl.simname_from_dirpath(path)
    
    with h5py.File(filen_main, 'r') as f:
        todoc = {}
        halodat = {}
        # check simulation run, snapshot
        if simid in f:
            smgrp = f[simid]
        else:
            raise NoStoredMatchError(f'Simulation {simid}')
        snn = f'snap_{snapshot}'
        if snn in smgrp:
            sngrp = smgrp[snn]
        else:
            raise NoStoredMatchError(f'Simulation {simid}, {snn}')
        cosmopars = {}
        for key, val in sngrp['cosmopars'].attrs.items():
            cosmopars[key] = val
        todoc['cosmopars'] = cosmopars
        # check center finding
        cengrpns = [grp for grp in sngrp.keys() if grp.startswith('cen')]
        for cengrpn in cengrpns:
            cgrp = sngrp[cengrpn]
            tomatch = usedvals_calchalo.keys()
            # using: 'parttypes_used' is a tuple, comparison to 
            # array gives boolean array, or False if different lengths
            if np.all([np.all(usedvals_calchalo[key] == cgrp.attrs[key])\
                        for key in tomatch]):
                halodat['Xc_cm'] = cgrp.attrs['Xc_cm']
                halodat['Yc_cm'] = cgrp.attrs['Yc_cm']
                halodat['Zc_cm'] = cgrp.attrs['Zc_cm']
                todoc.update(usedvals_calchalo)
                break
        if 'Xc_cm' not in halodat:
            msg = (f'Simulation {simid}, {snn}, '
                    f'center finding parameters {usedvals_calchalo}')
            raise NoStoredMatchError(msg)
        # check Mvir/Rvir def.
        outputsingle = False
        if isinstance(meandef, type('')):
            meandef = [meandef]
            outputsingle = True
        halodat['Rvir_cm'] = [] 
        halodat['Mvir_g'] = [] 
        for md in meandef:
            subgrpn = f'Rvir_{md}'
            if subgrpn in cgrp:
                sgrp = cgrp[subgrpn]
                halodat['Rvir_cm'].append(sgrp.attrs['Rvir_cm'])
                halodat['Mvir_g'].append(sgrp.attrs['Mvir_g'])
            else:
                msg = (f'Simulation {simid}, {snn}, '
                        f'center finding parameters {usedvals_calchalo}, '
                        f'overdensity definition {md}')
                raise NoStoredMatchError(msg)
        if outputsingle:
            halodat['Rvir_cm'] = halodat['Rvir_cm'][0]
            halodat['Mvir_g'] = halodat['Mvir_g'][0]
    print(f'Retrieved stored halo data from {filen_main}')
    return halodat, todoc  

def gethalodata_shrinkingsphere(path, snapshot, meandef=('200c', 'BN98')):
    '''
    same in/output as calchalodata_shrinkingsphere,
    but reads data from file if stored, and stores data to a temporary
    file if not.
    Run adddata_cenrvir() to add the temporary file data to the main 
    file. (Doing this during the main run could cause issues if multiple
    processes try to write to the same file at the same time.)
    '''
    # this must contain *all* todoc entries from the previous function
    # except 'parttypes_used', which is assumed to be everything but
    # PartType2 (lo-res DM) present in the NumPart_Total table
    usedvals_calchalo = {'shrinkfrac': 0.025, 
                         'minparticles': 1000., 
                         'initialradiusfactor': 1.,
                         'minpart_halo': 1000.}
    snap = rf.get_Firesnap(path, snapshot, filetype='snap')
    with h5py.File(snap.firstfilen) as f:
        pts = list(f['Header'].attrs['NumPart_Total'])
    pts = [ind for ind in range(len(pts)) if pts[ind] > 0]
    pts.remove(2)
    pts.sort()
    usedvals_calchalo['parttypes_used'] = tuple(pts)
                                                
    fdir = ol.dir_halodata
    
    newcalc = False
    #pparts = path.split('/')
    #while '' in pparts:
    #    pparts.remove('')
    #if pparts[-1] == 'output':
    #    pparts = pparts[:-1]
    #simid = pparts[-1]
    simid = sl.simname_from_dirpath(path)
    
    try:
        halodat, todoc = readhalodata_shrinkingsphere(
            path, snapshot, meandef=meandef)
        return halodat, todoc    
    except NoStoredMatchError as err:
        print(err)
        print('Center, Rvir were not stored')
        newcalc = True
    
    if newcalc:
        print('Calculating halo data...')
        halodat, todoc = calchalodata_shrinkingsphere(path, snapshot, 
                                                      meandef=meandef)
        print('Halo data calculated.')
        filen = fdir + f'temp_cen_rvir_{uuid.uuid1()}.hdf5'
        print(f'Saving data to file {filen}')
        if os.path.isfile(filen):
            msg = f'Temporary center/Rvir file {filen} already exists'
            raise RuntimeError(msg)
        with h5py.File(filen, 'w') as f:
            # sim, snap groups
            smgrp = f.create_group(simid)
            sngrp = smgrp.create_group(f'snap_{snapshot}')
            cmgrp = sngrp.create_group('cosmopars')
            for key in todoc['cosmopars']:
                cmgrp.attrs.create(key, todoc['cosmopars'][key])
            del todoc['cosmopars']
            # sim/snap subgroup for center pars.
            cengrp = sngrp.create_group('cen0')
            for cv in ['Xc_cm', 'Yc_cm', 'Zc_cm']:
                cengrp.attrs.create(cv, halodat[cv])
            for key in todoc:
                val = todoc[key]
                if isinstance(val, type('')):
                    val = np.string_(val)
                cengrp.attrs.create(key, val)
            # center pars. subgroups for mvir/rvir def.
            outputsingle = False
            if isinstance(meandef, type('')):
                meandef = [meandef]
                halodat['Rvir_cm'] = [halodat['Rvir_cm']] 
                halodat['Mvir_g'] = [halodat['Mvir_g']] 
                outputsingle = True
            for md, rv, mv in zip(meandef, halodat['Rvir_cm'], 
                                  halodat['Mvir_g']):
                gn = f'Rvir_{md}'
                vgrp = cengrp.create_group(gn)
                vgrp.attrs.create('Rvir_cm', rv)
                vgrp.attrs.create('Mvir_g', mv)
        print(f'Saved new halo data.')
        if outputsingle:
            halodat['Rvir_cm'] = halodat['Rvir_cm'][0]
            halodat['Mvir_g'] = halodat['Mvir_g'][0]
        return halodat, todoc
    
def adddata_cenrvir(rmtemp=False):
    '''
    put data in temporary cenrvir and vcom files into the main file

    Parameters:
    -----------
    rmtemp: bool
        remove temporary storage files after marking it as duplicate
        (This would be the second run of this function if a file was
        new -- this is done to ensure a file is not deleted if 
        something went wrong unexpectedly.)
    '''
    mainfilen =  ol.filen_halocenrvir 
    searchcrit = ol.dir_halodata + 'temp_cen_rvir_*.hdf5'
    tempfilens = glob.glob(searchcrit)
    searchcrit2 = ol.dir_halodata + 'temp_vcom_*.hdf5'
    tempfilens2 = glob.glob(searchcrit2)
    tempfilens = tempfilens + tempfilens2
    if len(tempfilens) == 0:
        print('No new data to add')
        return None
    with h5py.File(mainfilen, 'a') as fo:
        #print(mainfilen)
        for tfn in tempfilens:
            vcomfile = (tfn.split('/')[-1]).startswith('temp_vcom')
            with h5py.File(tfn, 'r') as fi:
                # should have one sim, snap, cen group
                # possibly multiple rvir definitions
                simid = next(iter(fi.keys()))
                #print('simid: ', simid)
                #print('main file keys: ', list(fo.keys()))
                if simid not in fo: #easy, copy whole thing
                    fi.copy(fi[simid], fo, name=simid)
                    print(f'Added file {tfn}:')
                    print(f'{simid}')
                    continue
                fo_smgrp = fo[simid]
                fi_smgrp = fi[simid]
                sngrpn = next(iter(fi_smgrp.keys()))
                #print('snap: ', sngrpn)
                #print('main file fo_smgrp keys: ', list(fo_smgrp.keys()))
                if sngrpn not in fo_smgrp: #easy, copy whole thing
                    fi.copy(fi_smgrp[sngrpn], fo_smgrp, name=sngrpn)
                    print(f'Added file {tfn}:')
                    print(f'{simid}, {sngrpn}')
                    continue
                # center matching/copy
                fo_sngrp = fo_smgrp[sngrpn]
                fi_sngrp = fi_smgrp[sngrpn]
                #print('main file cens: ', list(fo_sngrp.keys()))
                if 'cen0' not in fo_sngrp:
                    fi.copy(fi_sngrp['cen0'], fo_sngrp, name='cen0')
                    print(f'Added file {tfn}:')
                    print(f'{simid}, {sngrpn}, first center')
                    fo_cgrpn = 'cen0'
                    continue
                cens_fo = [grp for grp in fo_sngrp.keys() \
                           if grp.startswith('cen')]
                fi_cgrp = fi_sngrp['cen0']
                tocheck = ['Xc_cm', 'Yc_cm', 'Zc_cm']
                anymatch = False
                for cengrpn in cens_fo:
                    _fo_cgrp = fo_sngrp[cengrpn]
                    tomatch = set(fi_cgrp.attrs.keys()) - set(tocheck)
                    # using: 'parttypes_used' is a tuple, comparison to 
                    # array gives boolean array, or False if different lengths
                    if np.all([np.all(_fo_cgrp.attrs[key] \
                                      == fi_cgrp.attrs[key])\
                               for key in tomatch]):
                        fo_cgrp = _fo_cgrp
                        fo_cgrpn = cengrpn
                        anymatch = True
                        if not np.all([np.all(_fo_cgrp.attrs[key] \
                                       == fi_cgrp.attrs[key])\
                                       for key in tocheck]):
                            msg = (f'{mainfilen} and {tfn} have matching'
                                   f'simulation {simid}, {sngrpn}, '
                                   f'center finding, but different centers:\n'
                                   f'{fo_cgrp.attrs.items()},\n'
                                   f'{fi_cgrp.attrs.items()}')
                            raise RuntimeError(msg)
                if not anymatch:
                    fo_cgrpn = f'cen{len(cens_fo)}'
                    fi.copy(fi_cgrp, fo_sngrp, name=fo_cgrpn)
                    print(f'Added file {tfn}:')
                    print(f'{simid}, {sngrpn}, new center')
                    continue
                #print('main file densities: ', list(fo_cgrp.keys()))
                # mvir/rvir matching/copy
                fi_mrdefs = [grp for grp in fi_cgrp.keys() \
                             if grp.startswith('Rvir_')]
                #print('new densities: ', fi_mrdefs)
                for mdn in fi_mrdefs:
                    if mdn in fo_cgrp:
                        fi_dct = dict(fi_cgrp[mdn].attrs.items())
                        fo_dct = dict(fo_cgrp[mdn].attrs.items())
                        if fi_dct == fo_dct:
                            if not vcomfile:
                                break
                        else:
                            msg = (f'{mainfilen} and {tfn} have matching'
                                   f'simulation {simid}, {sngrpn}, '
                                   f'centers, but different Mvir or Rvir:\n'
                                   f'{fi_dct},\n'
                                   f'{fo_dct}')
                            raise RuntimeError(msg)
                    else:
                        fi.copy(fi_cgrp[mdn], fo_cgrp, name=mdn)
                        print(f'Added cen/Rvir file {tfn}:')
                        print(f'{simid}, {sngrpn}, {fi_mrdefs}')
                fi_rgrp = fi_cgrp[mdn]
                fo_rgrp = fo_cgrp[mdn]
                # Vcom copy, if any 
                if vcomfile:
                    if 'Vcom0' not in fo_rgrp:
                        fi.copy(fi_rgrp['Vcom0'], fo_rgrp, name='Vcom0')
                        print(f'Added Vcom file {tfn}:')
                        l2 = (f'{simid}, {sngrpn}, {fo_cgrpn}, {fi_mrdefs}'
                            ', first Vcom')
                        print(l2)
                        continue
                    vcom_fo = [grp for grp in fo_rgrp.keys() \
                               if grp.startswith('Vcom')]
                    fi_vgrp = fi_rgrp['Vcom0']
                    tocheck = ['VXcom_cmps', 'VYcom_cmps', 'VZcom_cmps']
                    anymatch = False
                    for vgrpn in vcom_fo:
                        _fo_vgrp = fo_rgrp[vgrpn]
                        tomatch = set(fi_vgrp.attrs.keys()) - set(tocheck)
                        # using: 'parttypes_used' is a tuple, 
                        # comparison to array gives boolean array, or
                        # False if different lengths
                        if np.all([np.all(_fo_vgrp.attrs[key] \
                                        == fi_vgrp.attrs[key])\
                                for key in tomatch]):
                            fo_vgrp = _fo_vgrp
                            anymatch = True
                            if not np.all([np.all(_fo_vgrp.attrs[key] \
                                        == fi_vgrp.attrs[key])\
                                        for key in tocheck]):
                                msg = (f'{mainfilen} and {tfn} have matching'
                                       f'simulation {simid}, {sngrpn}, '
                                       f'center finding, {mdn}'
                                       ' but different Vcom:\n'
                                       f'{fo_vgrp.attrs.items()},\n'
                                       f'{fi_vgrp.attrs.items()}')
                                raise RuntimeError(msg)
                    if not anymatch:
                        fo_vgrpn = f'Vcom{len(vcom_fo)}'
                        fi.copy(fi_vgrp, fo_rgrp, name=fo_vgrpn)
                        print(f'Added file {tfn}:')
                        print(f'{simid}, {sngrpn}, {mdn}, new Vcom')
                        continue
            print(f'skipped {tfn}; duplicate data')
            if rmtemp:
                print(f'deleting {tfn}')
                os.remove(tfn)

def calc_vcom(path, snapshot, radius_rvir, meandef_rvir='BN98',
              parttypes='all'):
    '''
    calculate center of mass velocity for specified particle 
    types in a specified fraction of the given virial radius

    Parameters:
    -----------
    path: str
        where to find the simulation data (include full path),
        directory should contain the snapshots or an 'ouput'
        subdirectory with the snapshots
    snapshot: int
        snapshot number
    radius_rvir: float
        radius of the sphere in which to include particles, in
        units of the virial radius
    meandef_rvir: ['BN98', '<overdensity>c', '<overdensity>m']
        which definition of the virial radius to use 
        (see calc_halodata)
    parttypes: tuple of ints (0, 1, 2, 4, and/or 5) or 'all'
        which particle types to include in the calculation.
        'all' uses all types except 2 (low-resolution dark matter),
        and 5 if the simulation does not contain black holes.
    '''
    halodat, todoc_cv = gethalodata_shrinkingsphere(path, snapshot, 
                                                    meandef=meandef_rvir)
    cen_cm = np.array([halodat['Xc_cm'], halodat['Yc_cm'], halodat['Zc_cm']])
    rvir_cm = halodat['Rvir_cm']
    snap = rf.get_Firesnap(path, snapshot)
    todoc = {'radius_rvir': radius_rvir}

    # get mass and coordinate data
    if parttypes == 'all':
        parttypes = (0, 1, 4, 5)
    else:
        parttypes = parttypes
    dct_m = {}
    dct_r = {}
    toCGS_m = None
    toCGS_c = None
    for pt in parttypes:
        cpath = 'PartType{}/Coordinates'
        mpath = 'PartType{}/Mass'
        try:
            ctemp = snap.readarray_emulateEAGLE(cpath.format(pt))
            _toCGS_c = snap.toCGS
            dct_m[pt] = snap.readarray_emulateEAGLE(mpath.format(pt))
            _toCGS_m = snap.toCGS
        except (OSError, rf.FieldNotFoundError):
            msg = 'Skipping PartType {} in COM vel. calc: not present on file'
            print(msg.format(pt))
            continue
        if toCGS_m is None:
            toCGS_m = _toCGS_m
        elif not np.isclose(toCGS_m, _toCGS_m):
                msg = ('Different particle type masses have different'
                       ' CGS conversions in ' + snap.firstfilen)
                raise RuntimeError(msg)
        if toCGS_c is None:
            toCGS_c = _toCGS_c
        elif not np.isclose(toCGS_c, _toCGS_c):
                msg = ('Different particle type coordinates have different'
                       ' CGS conversions in ' + snap.firstfilen)
                raise RuntimeError(msg)
        ctemp -= cen_cm / toCGS_c
        dct_r[pt] = np.sum(ctemp**2, axis=1)
        del ctemp

    pt_used = list(dct_m.keys())
    pt_used.sort()
    totlen = sum([len(dct_m[pt]) for pt in pt_used])
    masses = np.empty((totlen,), dtype=dct_m[pt_used[0]].dtype)
    rsq = np.empty((totlen,), dtype=dct_r[pt_used[0]].dtype)
    todoc['parttypes_used'] = tuple(pt_used)
    start = 0
    for pt in pt_used:
        partlen = len(dct_m[pt])
        masses[start: start + partlen] = dct_m[pt]
        rsq[start: start + partlen] = dct_r[pt]
        start += partlen

        del dct_m[pt]
        del dct_r[pt]
    rsq_max = (radius_rvir * rvir_cm / toCGS_c)**2
    psel = rsq <= rsq_max
    del rsq
    
    dct_v = {}
    toCGS_v = None
    # get velocities
    for pt in pt_used:
        vpath = 'PartType{}/Velocities'
        dct_v[pt] = snap.readarray_emulateEAGLE(vpath.format(pt))
        _toCGS_v = snap.toCGS
        if toCGS_v is None:
            toCGS_v = _toCGS_v
        elif not np.isclose(toCGS_v, _toCGS_v):
                msg = ('Different particle type velocities have different'
                       ' CGS conversions in ' + snap.firstfilen)
                raise RuntimeError(msg)
    totlen = sum([len(dct_v[pt]) for pt in pt_used])
    velocities = np.empty((totlen, 3), dtype=dct_v[pt_used[0]].dtype)
    start = 0
    for pt in pt_used:
        partlen = len(dct_v[pt])
        velocities[start: start + partlen] = dct_v[pt]
        start += partlen
        del dct_v[pt]
    masses = masses[psel]
    velocities = velocities[psel]
    vcom_codeunits = np.sum(masses[:, np.newaxis] * velocities, axis=0) \
                     / np.sum(masses)
    vcom_cmps = vcom_codeunits * toCGS_v
    todoc['units'] = 'cm * s**-1'
    out = {'VXcom_cmps': vcom_cmps[0], 
           'VYcom_cmps': vcom_cmps[1],
           'VZcom_cmps': vcom_cmps[2]}
    halodat.update(out)
    print('todoc output calc_vcom:')
    print(todoc)
    return halodat, todoc


def readdata_vcom(path, snapshot, radius_rvir, meandef_rvir='BN98',
                  parttypes='all', datafile=None):
    # raises NoStoredMatchError if data isn't present 
    # -> also no vcom data present
    halodat, todoc_cen = readhalodata_shrinkingsphere(path, snapshot,
                                                      meandef=meandef_rvir)
    
    todoc = {}
    usedvals_calchalo = {'shrinkfrac': 0.025, 
                         'minparticles': 1000., 
                         'initialradiusfactor': 1.,
                         'minpart_halo': 1000.}
    snap = rf.get_Firesnap(path, snapshot, filetype='snap')
    with h5py.File(snap.firstfilen) as f:
        pts = list(f['Header'].attrs['NumPart_Total'])
    pts = [ind for ind in range(len(pts)) if pts[ind] > 0]
    pts.remove(2)
    pts.sort()
    usedvals_calchalo['parttypes_used'] = tuple(pts)
    meandef = meandef_rvir

    if parttypes == 'all':
        pts_vcom = pts
    else:
        pts_vcom = list(parttypes)
        pts_vcom.sort()
        pts_vcom = tuple(pts_vcom)

    usedvalues_calcvcom = {'radius_rvir': radius_rvir,
                           'parttypes_used': pts_vcom}
    if datafile is None:                    
        filen_main = ol.filen_halocenrvir
    else:
        filen_main = datafile
    
    #pparts = path.split('/')
    #while '' in pparts:
    #    pparts.remove('')
    #if pparts[-1] == 'output':
    #    pparts = pparts[:-1]
    #simid = pparts[-1]
    simid = sl.simname_from_dirpath(path)

    with h5py.File(filen_main, 'r') as f:
        # find the halo center/rvir group
        smgrp = f[simid]
        snn = f'snap_{snapshot}'
        sngrp = smgrp[snn]
        cosmopars = {}
        for key, val in sngrp['cosmopars'].attrs.items():
            cosmopars[key] = val
        todoc['cosmopars'] = cosmopars
        cengrpns = [grp for grp in sngrp.keys() if grp.startswith('cen')]
        for cengrpn in cengrpns:
            cgrp = sngrp[cengrpn]
            tomatch = usedvals_calchalo.keys()
            # using: 'parttypes_used' is a tuple, comparison to 
            # array gives boolean array, or False if different lengths
            if np.all([np.all(usedvals_calchalo[key] == cgrp.attrs[key])\
                        for key in tomatch]):
                halodat['Xc_cm'] = cgrp.attrs['Xc_cm']
                halodat['Yc_cm'] = cgrp.attrs['Yc_cm']
                halodat['Zc_cm'] = cgrp.attrs['Zc_cm']
                todoc.update(usedvals_calchalo)
                break
        subgrpn = f'Rvir_{meandef}'
        sgrp = cgrp[subgrpn]
        halodat['Rvir_cm'] = sgrp.attrs['Rvir_cm']
        halodat['Mvir_g'] = sgrp.attrs['Mvir_g']
        # look for the center of mass
        for vcomgrpn in sgrp.keys():
            if not vcomgrpn.startswith('Vcom'):
                continue
            vgrp = sgrp[vcomgrpn]
            tomatch = usedvalues_calcvcom.keys()
            # using: 'parttypes_used' is a tuple, comparison to 
            # array gives boolean array, or False if different lengths
            if np.all([np.all(usedvalues_calcvcom[key] == vgrp.attrs[key])\
                        for key in tomatch]):
                halodat['VXcom_cmps'] = vgrp.attrs['VXcom_cmps']
                halodat['VYcom_cmps'] = vgrp.attrs['VYcom_cmps']
                halodat['VZcom_cmps'] = vgrp.attrs['VZcom_cmps']
                todoc.update(usedvalues_calcvcom)
                break
        if 'VXcom_cmps' not in halodat:
            msg = (f'No Vcom stored for {path}, snapshot {snapshot}, '
                    f'radius_rvir {radius_rvir}, meandef_rvir '
                    f'{meandef_rvir}, particle types {pts_vcom}')
            raise NoStoredMatchError(msg)
    print(f'Retrieved stored halo data from {filen_main}')
    return halodat, todoc 

def writedata_vcom(halodat, todoc,
                   path, snapshot, meandef_rvir='BN98',
                   datafile=None):
    usedvals_calchalo = {'shrinkfrac': 0.025, 
                         'minparticles': 1000., 
                         'initialradiusfactor': 1.,
                         'minpart_halo': 1000.}
    snap = rf.get_Firesnap(path, snapshot, filetype='snap')
    with h5py.File(snap.firstfilen) as f:
        pts = list(f['Header'].attrs['NumPart_Total'])
    pts = [ind for ind in range(len(pts)) if pts[ind] > 0]
    pts.remove(2)
    pts.sort()
    usedvals_calchalo['parttypes_used'] = tuple(pts)

    usedvalues_calcvcom = {'radius_rvir': todoc['radius_rvir'],
                           'parttypes_used': todoc['parttypes_used']}
    if datafile is None:                    
        filen = ol.filen_halocenrvir
    else:
        filen = datafile
    
    #pparts = path.split('/')
    #while '' in pparts:
    #    pparts.remove('')
    #if pparts[-1] == 'output':
    #    pparts = pparts[:-1]
    #simid = pparts[-1]
    simid = sl.simname_from_dirpath(path)

    with h5py.File(filen, 'a') as fo:
        # should have one sim, snap, cen group
        # possibly multiple rvir definitions
        #print('simid: ', simid)
        #print('main file keys: ', list(fo.keys()))
        if simid not in fo: #easy, copy whole thing
            sgrp = fo.create_group(simid)
        else:
            sgrp = fo[simid]
        sngrpn = f'snap_{snapshot}'
        if sngrpn in sgrp: #easy, copy whole thing
            sngrp = sgrp[sngrpn]
        else:
            sngrp = sgrp.create_group(sngrpn)
        #print('main file cens: ', list(fo_sngrp.keys()))
        cgrpn_opts = [grp for grp in sngrp.keys() 
                        if grp.startswith('cen')]
        tocheck = ['Xc_cm', 'Yc_cm', 'Zc_cm']
        anymatch = False
        for cengrpn in cgrpn_opts:
            _cgrp = sngrp[cengrpn]
            tomatch = set(_cgrp.attrs.keys()) - set(tocheck)
            # using: 'parttypes_used' is a tuple, comparison to 
            # array gives boolean array, or False if different lengths
            if np.all([np.all(_cgrp.attrs[key] \
                                == usedvals_calchalo[key])\
                        for key in tomatch]):
                cgrp = _cgrp
                anymatch = True
                if not np.all([np.all(_cgrp.attrs[key] \
                                == halodat[key])\
                                for key in tocheck]):
                    msg = (f'{filen} has matching'
                            f'simulation {simid}, {sngrpn}, '
                            f'center finding, but different centers:\n'
                            f'{_cgrp.attrs.items()},\n'
                            f'{halodat}')
                    raise RuntimeError(msg)
        if not anymatch:
            cgrpn = f'cen{len(cgrpn_opts)}'
            cgrp = sngrp.create_group(cgrpn)
            for key in usedvals_calchalo:
                cgrp.attrs.create(key, usedvals_calchalo[key])
            for key in ['Xc_cm', 'Yc_cm', 'Zc_cm']:
                cgrp.attrs.create(key, halodat[key])
        rgrpn = f'Rvir_{meandef_rvir}'
        #print('new densities: ', fi_mrdefs)
        if rgrpn in cgrp:
            rgrp = cgrp[rgrpn]
            if not np.all([np.all(rgrp.attrs[key] == halodat[key])\
                           for key in rgrp.attrs.keys()]):
                msg = (f'{filen} has matching'
                        f'simulation {simid}, {sngrpn}, '
                        f'center finding, but different Mv, Rv'
                        f' ({meandef_rvir}):\n'
                        f'{rgrp.attrs.items()},\n'
                        f'{halodat}')
                raise RuntimeError(msg)
            for key in ['Mvir_g', 'Rvir_cm']:
                rgrp.attrs.create(key, halodat[key])
        else:
            rgrp = cgrp.create_group(rgrpn)
            for key in ['Mvir_g', 'Rvir_cm']:
                rgrp.attrs.create(key, halodat[key])
        vmatch = False
        for vcomgrpn in rgrp.keys():
            if not vcomgrpn.startswith('Vcom'):
                continue
            vgrp = rgrp[vcomgrpn]
            tomatch = usedvalues_calcvcom.keys()
            # using: 'parttypes_used' is a tuple, comparison to 
            # array gives boolean array, or False if different lengths
            if np.all([np.all(usedvalues_calcvcom[key] == vgrp.attrs[key])\
                        for key in tomatch]):
                vmatch = True
                if not np.all([halodat[key] == vgrp.attrs[key]
                               for key in 
                               ['VXcom_cmps', 'VYcom_cmps', 'VZcom_cmps']]):
                    msg = (f'{filen} has matching'
                           f'simulation {simid}, {sngrpn}, '
                           f'center finding, Mv, Rv ({meandef_rvir})\n'
                           f'but different Vcom:\n'
                           f'{vgrp.attrs.items()},\n'
                           f'{halodat}')
                    raise RuntimeError(msg)(usedvalues_calcvcom)
                print(f'{filen} already contained data to store')
                break
        if not vmatch:
            grpnum = len([key for key in rgrp.keys() 
                          if key.startswith('Vcom')])
            vgrpn = f'Vcom{grpnum}'
            vgrp = rgrp.create_group(vgrpn)
            for key in usedvalues_calcvcom:
                vgrp.attrs.create(key, usedvalues_calcvcom[key])
            for key in ['VXcom_cmps', 'VYcom_cmps', 'VZcom_cmps']:
                vgrp.attrs.create(key, halodat[key])

def get_vcom(path, snapshot, radius_rvir, meandef_rvir='BN98',
             parttypes='all'):
    '''
    same in/output as calchalodata_shrinkingsphere,
    but reads data from file if stored, and stores data to a temporary
    file if not.
    Run adddata_cenrvir() to add the temporary file data to the main 
    file. (Doing this during the main run could cause issues if multiple
    processes try to write to the same file at the same time.)
    '''
    try:
        out = readdata_vcom(path, snapshot, radius_rvir, 
                            meandef_rvir=meandef_rvir,
                            parttypes=parttypes, datafile=None)
    except NoStoredMatchError as err:
        print(err)
        print('Calculating COM velocity')
        out = calc_vcom(path, snapshot, radius_rvir, meandef_rvir=meandef_rvir,
                        parttypes=parttypes)
        print('Vcom calculated.')
        filen = ol.dir_halodata + f'temp_vcom_{uuid.uuid1()}.hdf5'
        print(f'Saving data to file {filen}')
        if os.path.isfile(filen):
            msg = f'Temporary Vcom file {filen} already exists'
            raise RuntimeError(msg)
        writedata_vcom(out[0], out[1],
                       path, snapshot, meandef_rvir=meandef_rvir,
                       datafile=filen)
        print(f'Saved new halo data.')
    return out
        
def halodata_rockstar(path, snapnum, select='maxmass', 
                      masspath='mass.vir'):
    '''
    retrieve position, mass, and radius from rockstar halo data
    uses the Bryan and Norman overdensity mass (.vir in rockstar)
    
    Parameters:
    -----------
    path: str
        path to the directory containing the output and halo directories
        and the snapshot_times.txt file
    snapnum: int
        snapshot number
    select: {'maxmass', 'mainprog', int}
        how to select the halo to use
        'maxmass': highest mass in the snapshot
        'mainprog': main progenitor of the highest-mass halo at the
                    lowest redshift available
        int: index of the halo in the snapshot catalogue (Note: not the
             tree index that's unique across snapshots)
    masspath: path in hdf5 file to mass to use
        e.g., mass.mvir, mass.200c, mass.200m
    '''
    # options: 'BN98', '####c', '######m'
    if masspath == 'mass.vir':
        meandensdef = 'BN98'
    else:
        meandensdef = masspath.split('.')[-1]
    out = {}
    if select == 'maxmass' or isinstance(select, int):
        hal = ha.io.IO.read_catalogs('snapshot', snapnum, path)
        if select == 'maxmass':
            haloind = np.argmax(hal[masspath])
            print('Using snapshot halo index: {}'.format(haloind))
        else:
            haloind = select
        out['Mvir_Msun'] = hal[masspath][haloind]
        out['Xc_ckpc'], out['Yc_ckpc'], out['Zc_ckpc'] = \
            hal['position'][haloind]
        cosmopars = {}
        cosmopars['omegalambda'] = hal.Cosmology['omega_lambda']
        cosmopars['omegam'] = hal.Cosmology['omega_matter']
        cosmopars['omegab'] = hal.Cosmology['omega_baryon']
        cosmopars['h'] = hal.Cosmology['hubble']
        cosmopars['a'] = hal.snapshot['scalefactor']
        cosmopars['z'] = hal.snapshot['redshift']
    elif select == 'mainprog':
        halt = ha.io.IO.read_tree(simulation_directory=path, 
                                  species_snapshot_indices=[snapnum])
        # high-mass stuff isn't always run to z=0
        finalsnap = np.max(halt['snapshot'])
        wherefinalsnap = np.where(halt['snapshot'] == finalsnap)[0]
        whereind_maxmfinal = np.argmax(halt[masspath][wherefinalsnap])
        treeind_maxmfinal = wherefinalsnap[whereind_maxmfinal]
        prog_main_index = treeind_maxmfinal
        while prog_main_index >= 0:
            snap_current = halt['snapshot'][prog_main_index]
            if snap_current == snapnum:
                break
            if prog_main_index < 0:
                msg = 'No main progenitor at snapshot {} was found'
                raise RuntimeError(msg.format(snap_current + 1))
            prog_main_index = halt['progenitor.main.index'][prog_main_index]
        if bool(halt['am.phantom'][prog_main_index]):
            msg = 'This halo was not found by Rockstar,'+\
                  ' but interpolated'
            raise RuntimeError(msg)
        out['Mvir_Msun'] = halt['mass'][prog_main_index]
        out['Xc_ckpc'], out['Yc_ckpc'], out['Zc_ckpc'] = \
            halt['position'][prog_main_index]
        cosmopars = {}
        cosmopars['omegalambda'] = halt.Cosmology['omega_lambda']
        cosmopars['omegam'] = halt.Cosmology['omega_matter']
        cosmopars['omegab'] = halt.Cosmology['omega_baryon']
        cosmopars['h'] = halt.Cosmology['hubble']
        # get redshift from halo catalog
        hal = ha.io.IO.read_catalogs('snapshot', snapnum, path)
        cosmopars['a'] = hal.snapshot['scalefactor']
        cosmopars['z'] = hal.snapshot['redshift']
    
    meandens = cu.getmeandensity(meandensdef, cosmopars)
    #M = r_mean * 4/3 np.pi R63
    out['Rvir_cm'] = (3. / (4. * np.pi) * out['Mvir_Msun'] \
                      * c.solar_mass / meandens)**(1./3.)
    return out, cosmopars

def mainhalodata_AHFsmooth(path, snapnum):
    '''
    get properties of the main halo in the snapshot from halo_00000_smooth.dat
    assume units are intrinsic simulation units
    '''
    fn = path + '/halo/ahf/halo_00000_smooth.dat'
    df = pd.read_csv(fn, sep='\t')
    i = np.where(df['snum'] == snapnum)[0][0]
    out = {}
    # units from AHF docs: http://popia.ft.uam.es/AHF/files/AHF.pdf
    props = ['Mvir', 'Rvir', 'Xc', 'Yc', 'Zc']
    outprops = {'Mvir': 'Mvir_Msunoverh',
                'Rvir': 'Rvir_ckpcoverh',
                'Xc':   'Xc_ckpcoverh',
                'Yc':   'Yc_ckpcoverh',
                'Zc':   'Zc_ckpcoverh'}
    for prop in props:
        out[outprops[prop]] = df[prop][i]
    return out

if __name__ == '__main__':
    if len(sys.argv) > 1:
        mode = sys.argv[1]
        if mode == '--addstored':
            if len(sys.argv) > 2:
                rmtemp = bool(sys.argv[2])
            else:
                rmtemp = False
            adddata_cenrvir(rmtemp=rmtemp)
        else:
            raise ValueError(f'Invalid mode {mode}')
    else:
        print('No actions specified')