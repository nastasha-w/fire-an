

import h5py

import readfire.readin_fire_data as rf


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