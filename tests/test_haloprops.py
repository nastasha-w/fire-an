

import h5py
import numpy as np
import os
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gsp

import fire_an.mainfunc.haloprop as hp
import fire_an.makeplots.tol_colors as tc
import fire_an.readfire.readin_fire_data as rf
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu


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

def plot_halomasscheck(halofile, checkfile, imgname=None):
    '''
    compare the halo masses from the AHF/halo_0000_smooth.dat
    file and direct calculation of all mass within Rvir
    of the halo center Xc, Yc, Zc 

    Parameters:
    -----------
    halofile: str
        the halo_0000_smooth.dat file with full path
    checkfile: str
        file containing the calculated halo masses
    imgname: str
        name of the file to store the output image to
    
    Output:
    -------
    None
    '''

    checkdat = pd.read_csv(checkfile, sep='\t')
    ahfdat = pd.read_csv(halofile, sep='\t')
    # from log file of /projects/b1026/snapshots/metal_diffusion/m12i_res7100/
    # checks
    hconst =  0.702 

    fig = plt.figure(figsize=(11., 4.))
    grid = gsp.GridSpec(nrows=1, ncols=3, hspace=0.0, wspace=0.3, 
                        width_ratios=[1., 1., 1.])
    axes = [fig.add_subplot(grid[0, i]) for i in range(3)]
    fontsize = 12
    colors = tc.tol_cset('bright')
    size = 10.

    ax = axes[0]
    masslabel = 'Mvir [Msun]'
    xlabel = 'AHF ' + masslabel
    ylabel = masslabel + ' from AHF center and Rvir'
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(which='both', direction='in', labelsize=fontsize-1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    xv = np.array(checkdat['Mvir_AHF_Msun'])
    yv = np.array(checkdat['Mvir_sum_Msun'])
    minv = min(np.min(xv), np.min(yv))
    maxv = max(np.max(xv), np.max(yv))
    ax.plot([minv, maxv], [minv, maxv], color='gray', linestyle='dotted')
    ax.scatter(xv, yv, c=colors[0], s=size)

    ax = axes[1]
    xlabel = 'log (1 + redshift)'
    ylabel = masslabel
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(which='both', direction='in', labelsize=fontsize-1)
    ax.set_yscale('log')
    ahf_z = np.array(ahfdat['redshift'])
    _sort = np.argsort(ahf_z)
    ahf_z = ahf_z[_sort]
    ahf_mvir = np.array(ahfdat['Mvir'])[_sort] / hconst
    ahf_xv = np.log10(ahf_z + 1.)
    flat = np.diff(ahf_mvir) == 0.
    sectionbreaks = list(np.where(np.diff(flat) != 0)[0] + 1) 
    if flat[0]:
        sectionbreaks = [0] + sectionbreaks
    if flat[-1]:
        sectionbreaks = sectionbreaks + [len(flat)]
    sectionbreaks = np.array(sectionbreaks)
    sections = sectionbreaks.reshape((len(sectionbreaks) // 2, 2))
    flatx = [ahf_xv[sect[0]: sect[1] + 1] for sect in sections]
    flaty = [ahf_mvir[sect[0]: sect[1] + 1] for sect in sections]
    ax.plot(ahf_xv, ahf_mvir, color=colors[0])
    for _x, _y in zip(flatx, flaty):
        ax.plot(_x, _y, color='black')
    xv = np.log10(1. + checkdat['redshift'])
    ax.scatter(xv, checkdat['Mvir_AHF_Msun'], 
               color=colors[0], label='AHF mass', s=size)
    ax.scatter(xv, checkdat['Mvir_sum_Msun'], 
               color=colors[1], label='sum < AHF Rvir', s=size)
    ax.legend(fontsize=fontsize)
    ax.set_title('black: AHF halo mass is exactly flat', fontsize=fontsize)

    ax = axes[2]
    xlabel = 'log (1 + redshift)'
    ylabel = 'log abs ([M(< Rvir)] - [AHF]) / [AHF]'
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(which='both', direction='in', labelsize=fontsize-1)
    vahf = np.array(checkdat['Mvir_AHF_Msun'])
    vsum = np.array(checkdat['Mvir_sum_Msun'])
    yv = np.log10(np.abs((vsum - vahf) / vahf))
    xv = np.log10(1. + np.array(checkdat['redshift']))
    plt.scatter(xv, yv, c=colors[0], s=size)

    if imgname is not None:
        plt.savefig(imgname, bbox_inches='tight')

def runhalomasschecks(opt=1):
    checkdir = '/projects/b1026/nastasha/tests/start_fire/AHF_unit_tests/'
    if opt == 1:
        checkfile = checkdir +  'metal_diffusion__m12i_res7100.txt'
        halofile = '/projects/b1026/snapshots/metal_diffusion/m12i_res7100/halo/ahf/halo_00000_smooth.dat'
        imgname = checkdir + 'metal_diffusion__m11i_res7100_AHF-vs-sum.pdf'
        plot_halomasscheck(halofile, checkfile, imgname=imgname)
    else:
        raise ValueError('opt={} is not a valid option'.format(opt))