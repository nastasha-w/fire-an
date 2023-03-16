

import h5py
import numpy as np

import fire_an.readfire.readin_fire_data as rf


def checkfields_units(dirpath, snapnum, *args, numpart=100, 
                      outfilen='fields.hdf5'):
    '''
    Read in the data from the snapshot specified by dirpath
    and snap, read in the fields in args, convert to CGS, 
    save in file outfilen.
    '''
    snap = rf.get_Firesnap(dirpath, snapnum) 
    with h5py.File(outfilen, 'w') as f:
        hed = f.create_group('Header')
        cgrp = hed.create_group('cosmopars')
        cosmopars = snap.cosmopars.getdct()
        for key in cosmopars:
            cgrp.attrs.create(key, cosmopars[key])
        hed.attrs.create('snapnum', snapnum)
        hed.attrs.create('filepath_first', np.string_(snap.firstfilen))
        _info = 'datasets from FIRE hdf5 files stored in (physical) CGS units'
        hed.attrs.create('info', np.string_(_info))
        for arg in args:
            vals = snap.readarray_emulateEAGLE(arg)[:numpart]
            toCGS = snap.toCGS
            # e.g. masses overflow float32 in CGS
            # arrays are small anyway
            vals = vals.astype(np.float64) * toCGS 
            f.create_dataset(arg, data=vals)

def run_checkfields_units(index):
    dirpath1 = '/projects/b1026/snapshots/fire3/m13h206_m3e5/' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000/'
    simname1 = 'm13h206_m3e5__' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'
    snaps1 = [27, 45] # earliest and latest stored, for max a factor diffs
    fields1 = ('PartType0/Density', 
               'PartType0/Masses',
               'PartType0/Pressure',
               'PartType0/Temperature',
               'PartType0/Metallicity',
               'PartType0/ElementAbundance/Hydrogen',
               'PartType0/ElementAbundance/Oxygen',
               'PartType0/Coordinates',
               'PartType0/SmoothingLength')
    outdir1 = '/projects/b1026/nastasha/tests/start_fire/'
    outtemp = 'cgs_units_test_{simname}_snap{snap:03d}.hdf5'
    if index in [0, 1]:
        snap = snaps1[index]
        outfilen = outdir1 + outtemp.format(simname=simname1, snap=snap)
        checkfields_units(dirpath1, snap, *fields1, numpart=100, 
                          outfilen=outfilen)