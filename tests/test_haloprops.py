

import h5py
import numpy as np
import os

import mainfunc.haloprop as hp
import readfire.readin_fire_data as rf
import utils.constants_and_units as c
import utils.cosmo_utils as cu


def test_mainhalodata_units_ahf(opt=1, dirpath=None, snapnum=None,
                            printfile=None):
    
    if opt == 1: # redshift 0 test
        dirpath = '/projects/b1026/snapshots/metal_diffusion/m12i_res7100/'
        snapfile = dirpath + 'output/snapdir_600/snapshot_600.0.hdf5'
        snapnum = 600
    elif opt == 2: # higher z test 
        dirpath = '/projects/b1026/snapshots/metal_diffusion/m12i_res7100/'
        snapfile = dirpath + 'output/snapdir_399/snapshot_399.0.hdf5'
        snapnum = 399
    elif opt == 3: # try other z
        dirpath = '/projects/b1026/snapshots/metal_diffusion/m12i_res7100/'
        snapfile = dirpath + 'output/snapdir_492/snapshot_492.0.hdf5'
        snapnum = 492
    elif opt is None:
        pathopts = ['output/snapdir_{sn:03d}/snapshot_{sn:03d}.0.hdf5',
                    'output/snapshot_{sn:03d}.hdf5']
        goodpath = False
        for pathopt in pathopts:
            snapfile = dirpath + pathopt.format(sn=snapnum)
            if os.path.isfile(snapfile):
                goodpath = True
                break
        if not goodpath:
            tried = [dirpath + pathopts.format()]
            msg = 'Could not find snapshot {} in {}. Tried:'.format(snapnum, dirpath)
            msg = msg + '\n' + '\n'.join(tried)
            raise RuntimeError(msg)
    else:
        msg = 'test_mainhalodata_units parameter opt = {} is invalid'
        raise ValueError(msg.format(opt))

    halodat = hp.mainhalodata_AHFsmooth(dirpath, snapnum)
    snap = rf.Firesnap(snapfile) 
    cen = np.array([halodat['Xc_ckpcoverh'], 
                    halodat['Yc_ckpcoverh'], 
                    halodat['Zc_ckpcoverh']])
    cen_cm = cen * snap.cosmopars.a * 1e-3 * c.cm_per_mpc / snap.cosmopars.h
    rvir_cm = halodat['Rvir_ckpcoverh'] * snap.cosmopars.a\
              * 1e-3 * c.cm_per_mpc / snap.cosmopars.h
    print('Cosmology:')
    print(snap.cosmopars.getdct())
    print('Center [AHF units]: {}'.format(cen))
    print('Rvir [AHF units]: {}'.format(halodat['Rvir_ckpcoverh']))
    print('Center [attempted cm]: {}'.format(cen_cm))
    print('Rvir [attempted cm]: {}'.format(rvir_cm))
    
    # gas
    coords_pt0 = snap.readarray_emulateEAGLE('PartType0/Coordinates')
    coords_pt0_toCGS = snap.toCGS
    masses_pt0 = snap.readarray_emulateEAGLE('PartType0/Masses')
    masses_pt0_toCGS = snap.toCGS
    # sanity check
    med_c = np.median(coords_pt0, axis=0)
    print('Median gas coords [sim units]: {}'.format(med_c))
    print('Median gas coordinates [cm]: {}'.format(med_c * coords_pt0_toCGS))

    d2 = np.sum((coords_pt0 - cen_cm / coords_pt0_toCGS)**2, axis=1)
    sel = d2 <= (rvir_cm / coords_pt0_toCGS) **2
    hm_pt0 = np.sum(masses_pt0[sel])
    print('Halo gas mass (sim units): ', hm_pt0)
    print('Selected {}/{} particles'.format(np.sum(sel), len(sel)))
    del coords_pt0
    del masses_pt0
    del d2
    del sel
    # dm (high-res)
    coords_pt1 = snap.readarray_emulateEAGLE('PartType1/Coordinates')
    coords_pt1_toCGS = snap.toCGS
    masses_pt1 = snap.readarray_emulateEAGLE('PartType1/Masses')
    masses_pt1_toCGS = snap.toCGS
    med_c = np.median(coords_pt1, axis=0)
    print('Median DM coords [sim units]: {}'.format(med_c))
    print('Median DM coordinates [cm]: {}'.format(med_c * coords_pt1_toCGS))
    d2 = np.sum((coords_pt1 - cen_cm / coords_pt1_toCGS)**2, axis=1)
    sel = d2 <= (rvir_cm / coords_pt1_toCGS) **2
    hm_pt1 = np.sum(masses_pt1[sel])
    print('Halo dm mass (sim units): ', hm_pt1)
    print('Selected {}/{} particles'.format(np.sum(sel), len(sel)))
    del coords_pt1
    del masses_pt1
    del d2
    del sel
    # stars
    coords_pt4 = snap.readarray_emulateEAGLE('PartType4/Coordinates')
    coords_pt4_toCGS = snap.toCGS
    masses_pt4 = snap.readarray_emulateEAGLE('PartType4/Masses')
    masses_pt4_toCGS = snap.toCGS
    med_c = np.median(coords_pt4, axis=0)
    print('Median star coords [sim units]: {}'.format(med_c))
    print('Median star coordinates [cm]: {}'.format(med_c * coords_pt4_toCGS))

    d2 = np.sum((coords_pt4 - cen_cm / coords_pt4_toCGS)**2, axis=1)
    sel = d2 <= (rvir_cm / coords_pt4_toCGS) **2
    hm_pt4 = np.sum(masses_pt4[sel])
    print('Halo stellar mass (sim units): ', hm_pt4)
    del coords_pt4
    del masses_pt4
    del d2
    del sel
    hm = hm_pt0 + hm_pt1 + hm_pt4

    msg = 'Got halo mass {hm}, listed Mvir is {Mvir}'
    hm_list_msun = halodat['Mvir_Msunoverh'] / snap.cosmopars.h
    hm_sum_msun = hm * (masses_pt0_toCGS / cu.c.solar_mass)
    print(msg.format(hm=hm_sum_msun, Mvir=hm_list_msun))
    hm_logmsun = np.log10(hm) + np.log10(masses_pt0_toCGS / cu.c.solar_mass)
    print('sum total is 10^{logm} Msun'.format(logm=hm_logmsun))

    if printfile is not None:
        new = not os.path.isfile(printfile)
        with open(printfile, 'a') as f:
            if new:
                columns = ['snapnum', 'redshift', 'Mvir_sum_Msun', 'Mvir_AHF_Msun']
                f.write('\t'.join(columns) + '\n')
            vals = [snapnum, snap.cosmopars.z, hm_sum_msun, hm_list_msun]
            f.write('\t'.join([str(val) for val in vals]) + '\n')

def test_mainhalodata_units_rockstar(opt=1, dirpath=None, snapnum=None,
                                     printfile=None, **kwargs):
    
    if opt == 1: # redshift 1 test
        dirpath = '/projects/b1026/snapshots/MassiveFIRE/h113_A4_res33000/'
        snapfile = dirpath + 'output/snapshot_277.hdf5'
        snapnum = 277
    elif opt == 2: # higher z test 
        dirpath = '/projects/b1026/snapshots/MassiveFIRE/h113_A4_res33000/'
        snapfile = dirpath + 'output/snapshot_200.hdf5'
        snapnum = 200
    elif opt == 3: # try other z
        dirpath = '/projects/b1026/snapshots/MassiveFIRE/h113_A4_res33000/'
        snapfile = dirpath + 'output/snapshot_100.hdf5'
        snapnum = 100
    elif opt is None:
        pathopts = ['output/snapdir_{sn:03d}/snapshot_{sn:03d}.0.hdf5',
                    'output/snapshot_{sn:03d}.hdf5']
        goodpath = False
        for pathopt in pathopts:
            snapfile = dirpath + pathopt.format(sn=snapnum)
            if os.path.isfile(snapfile):
                goodpath = True
                break
        if not goodpath:
            tried = [dirpath + pathopts.format()]
            msg = 'Could not find snapshot {} in {}. Tried:'.format(snapnum, dirpath)
            msg = msg + '\n' + '\n'.join(tried)
            raise RuntimeError(msg)
    else:
        msg = 'test_mainhalodata_units parameter opt = {} is invalid'
        raise ValueError(msg.format(opt))

    halodat, halo_cosmopars = hp.halodata_rockstar(dirpath, snapnum)
    snap = rf.get_Firesnap(dirpath, snapnum) 
    cen = np.array([halodat['Xc_ckpc'], 
                    halodat['Yc_ckpc'], 
                    halodat['Zc_ckpc']])
    cen_cm = cen * snap.cosmopars.a * 1e-3 * c.cm_per_mpc
    rvir_cm = halodat['Rvir_cm'] 
    print('Cosmology (snapshot):')
    print(snap.cosmopars.getdct())
    print('Cosmology (halo data):')
    print(halo_cosmopars)
    print('Center [rockstar units]: {}'.format(cen))
    print('Rvir [pkpc]: {}'.format(rvir_cm / (1e-3 * c.cm_per_mpc)))
    print('Center [attempted cm]: {}'.format(cen_cm))
    print('Rvir [attempted cm]: {}'.format(rvir_cm))
    
    # gas
    coords_pt0 = snap.readarray_emulateEAGLE('PartType0/Coordinates')
    coords_pt0_toCGS = snap.toCGS
    masses_pt0 = snap.readarray_emulateEAGLE('PartType0/Masses')
    masses_pt0_toCGS = snap.toCGS
    # sanity check
    med_c = np.median(coords_pt0, axis=0)
    print('Median gas coords [sim units]: {}'.format(med_c))
    print('Median gas coordinates [cm]: {}'.format(med_c * coords_pt0_toCGS))

    d2 = np.sum((coords_pt0 - cen_cm / coords_pt0_toCGS)**2, axis=1)
    sel = d2 <= (rvir_cm / coords_pt0_toCGS) **2
    hm_pt0 = np.sum(masses_pt0[sel])
    print('Halo gas mass (sim units): ', hm_pt0)
    print('Selected {}/{} particles'.format(np.sum(sel), len(sel)))
    del coords_pt0
    del masses_pt0
    del d2
    del sel
    # dm (high-res)
    coords_pt1 = snap.readarray_emulateEAGLE('PartType1/Coordinates')
    coords_pt1_toCGS = snap.toCGS
    masses_pt1 = snap.readarray_emulateEAGLE('PartType1/Masses')
    masses_pt1_toCGS = snap.toCGS
    med_c = np.median(coords_pt1, axis=0)
    print('Median DM coords [sim units]: {}'.format(med_c))
    print('Median DM coordinates [cm]: {}'.format(med_c * coords_pt1_toCGS))
    d2 = np.sum((coords_pt1 - cen_cm / coords_pt1_toCGS)**2, axis=1)
    sel = d2 <= (rvir_cm / coords_pt1_toCGS) **2
    hm_pt1 = np.sum(masses_pt1[sel])
    print('Halo dm mass (sim units): ', hm_pt1)
    print('Selected {}/{} particles'.format(np.sum(sel), len(sel)))
    del coords_pt1
    del masses_pt1
    del d2
    del sel
    # stars
    coords_pt4 = snap.readarray_emulateEAGLE('PartType4/Coordinates')
    coords_pt4_toCGS = snap.toCGS
    masses_pt4 = snap.readarray_emulateEAGLE('PartType4/Masses')
    masses_pt4_toCGS = snap.toCGS
    med_c = np.median(coords_pt4, axis=0)
    print('Median star coords [sim units]: {}'.format(med_c))
    print('Median star coordinates [cm]: {}'.format(med_c * coords_pt4_toCGS))

    d2 = np.sum((coords_pt4 - cen_cm / coords_pt4_toCGS)**2, axis=1)
    sel = d2 <= (rvir_cm / coords_pt4_toCGS) **2
    hm_pt4 = np.sum(masses_pt4[sel])
    print('Halo stellar mass (sim units): ', hm_pt4)
    del coords_pt4
    del masses_pt4
    del d2
    del sel
    hm = hm_pt0 + hm_pt1 + hm_pt4

    msg = 'Got halo mass {hm}, listed Mvir is {Mvir}'
    hm_list_msun = halodat['Mvir_Msun']
    hm_sum_msun = hm * (masses_pt0_toCGS / cu.c.solar_mass)
    print(msg.format(hm=hm_sum_msun, Mvir=hm_list_msun))
    hm_logmsun = np.log10(hm) + np.log10(masses_pt0_toCGS / cu.c.solar_mass)
    print('sum total is 10^{logm} Msun'.format(logm=hm_logmsun))

    if printfile is not None:
        new = not os.path.isfile(printfile)
        with open(printfile, 'a') as f:
            if new:
                columns = ['snapnum', 'redshift', 'Mvir_sum_Msun', 'Mvir_rockstar_Msun']
                f.write('\t'.join(columns) + '\n')
            vals = [snapnum, snap.cosmopars.z, hm_sum_msun, hm_list_msun]
            f.write('\t'.join([str(val) for val in vals]) + '\n')


# checkinh halo_0000_smooth.dat:
# Mvir is exactly flat over a large range of redshift values in that file
# might be an AHF issue?
def test_mainhalodata_units_multi(dirpath, printfile, version='ahf',
                                  **kwargs):
    print('running test_mainhalodata_units_multi')
    _snapdirs = os.listdir(dirpath + 'output/')
    snaps = []
    for _sd in _snapdirs:
        # looking for something like snapdir_196, extract 196
        if _sd.startswith('snapdir'):
            _snap = int(_sd.split('_')[-1])
            # special case, permissions error
            try: 
                os.listdir(dirpath + 'output/' + _sd)
                snaps.append(_snap)
            except PermissionError:
                # shows up seemingly randomly
                print('\nskipping snapshot {} due to permissions issues\n'.format(_snap))
                continue
        elif _sd.startswith('snapshot') and _sd.endswith('.hdf5'):
            # something like snapshot_164.hdf5
            print(_snap)
            _snap = int((_sd.split('_')[-1]).split('.')[0])
            try:
                f = h5py.File(dirpath + 'output/' + _sd, 'r')
                f.close()
            except Exception as err:
                print('\nSkipping snapshot {} due to h5py read issues:')
                print(err)
                print('\n')
                
    for snap in snaps:
        print('Snapshot ', snap)
        if version == 'ahf':
            test_mainhalodata_units_ahf(opt=None, dirpath=dirpath, 
                                        snapnum=snap,
                                        printfile=printfile)
        
        elif version == 'rockstar':
            test_mainhalodata_units_rockstar(opt=None, dirpath=dirpath, 
                                        snapnum=snap,
                                        printfile=printfile, **kwargs)
        else: 
            raise ValueError('invalid version option: {}'.format(version))
        print('\n')


def test_mainhalodata_units_multi_handler(opt=1):
    if opt == 1:
        dirpath = '/projects/b1026/snapshots/metal_diffusion/m12i_res7100/'
        printfile = '/projects/b1026/nastasha/tests/start_fire/AHF_unit_tests/'
        printfile += 'metal_diffusion__m12i_res7100.txt'
        version = 'ahf'
    elif opt == 2:
        dirpath = '/projects/b1026/snapshots/metal_diffusion/m11i_res7100/'
        printfile = '/projects/b1026/nastasha/tests/start_fire/AHF_unit_tests/'
        printfile += 'metal_diffusion__m11i_res7100.txt'
        version = 'ahf'
    else:
        raise ValueError('opt {} is not allowed'.format(opt))
    print('Running test_mainhalodata_units_multi(dirpath, printfile)')
    print('dirpath: ', dirpath)
    print('printfile: ', printfile)
    test_mainhalodata_units_multi(dirpath, printfile, version=version)
