import numpy as np
import h5py

import fire_an.mainfunc.get_qty as gq
import fire_an.mainfunc.haloprop as hp
import fire_an.readfire.readin_fire_data as rfd
import fire_an.utils.constants_and_units as c
import fire_an.utils.h5utils as h5u

def wtdavnorm(dirpath, snapnum, parttype, 
              maptype1, maptype_args1, maptype2, maptype_args2,
              rbins=(0.1, 1.), rbins_unit='Rvir',
              savefile='example.hdf5'):
    snapobj = rfd.get_Firesnap(dirpath, snapnum)
    halodata, halo_todoc = hp.gethalodata_shrinkingsphere(dirpath, snapnum, 
                                                          meandef='BN98')
    center_cm = np.array([halodata['Xc_cm'], 
                          halodata['Yc_cm'], 
                          halodata['Zc_cm']])
    rvir_cm = halodata['Rvir_cm']

    rcen_simu, rcen_tocgs, rcen_todoc = gq.get_qty(
        snapobj, parttype, 'coords', {'pos': 'rcen', 'center_cm': center_cm}, 
        filterdct=None)
    if rbins_unit == 'Rvir':
        rbins_simu = np.array(rbins) * rvir_cm / rcen_tocgs
    elif rbins_unit == 'pkpc':
        rbins_simu = np.array(rbins) * (c.cm_per_mpc * 1e-3) / rcen_tocgs
    filter = rcen_simu <= rbins_simu[-1]
    filter &= rcen_simu >= rbins_simu[0]
    rcen_simu = rcen_simu[filter]
    filterdct = {'filter': filter}
    
    qty1, qty1_cgsconv, qty1_todoc = gq.get_qty(
        snapobj, parttype, maptype1, maptype_args1, filterdct=filterdct)
    qty2, qty2_cgsconv, qty2_todoc = gq.get_qty(
        snapobj, parttype, maptype2, maptype_args2, filterdct=filterdct)
    
    rinds = np.searchsorted(rbins_simu, rcen_simu) - 1
    sums1 = []
    sums2 = []
    sums12 = []
    navs = []
    counts = []
    for ri in range(len(rbins) - 1):
        rsel = rinds == ri
        _qty1 = qty1[rsel]
        _qty2 = qty2[rsel]
        count = np.sum(rsel)
        s1 = np.sum(_qty1)
        s2 = np.sum(_qty2)
        s12 = np.sum(_qty1 * _qty2)
        nav = s12 * float(count) / (s1 * s2)
        sums1.append(s1)
        sums2.append(s2)
        sums12.append(s12)
        navs.append(nav)
        counts.append(count)
    
    with h5py.File(savefile, 'w') as f:
        hed = f.create_group('Header')

        hgrp = hed.create_group('halodata')
        h5u.savedict_hdf5(hgrp, halodata)
        _hgrp = hgrp.create_group('halodata_todoc')
        h5u.savedict_hdf5(_hgrp, halo_todoc)

        rcgrp = hed.create_group('rcen_coordinates_doc')
        h5u.savedict_hdf5(rcgrp, rcen_todoc)
        
        q1grp = hed.create_group('qty1')
        q1grp.attrs.create('maptype', np.string_(maptype1))
        _q1grp = q1grp.create_group('maptype_args_dct')
        h5u.savedict_hdf5(_q1grp, maptype_args1)
        __q1grp = q1grp.create_group('todoc_dct')
        h5u.savedict_hdf5(__q1grp, qty1_todoc)

        q2grp = hed.create_group('qty2')
        q2grp.attrs.create('maptype', np.string_(maptype2))
        _q2grp = q2grp.create_group('maptype_args_dct')
        h5u.savedict_hdf5(_q2grp, maptype_args2)
        __q2grp = q2grp.create_group('todoc_dct')
        h5u.savedict_hdf5(__q2grp, qty2_todoc)

        s1grp = f.create_group('sums1')
        s1grp.attrs.create('tocgs', qty1_cgsconv)
        s1grp.create_dataset('sums1', data=sums1)

        s2grp = f.create_group('sums2')
        s2grp.attrs.create('tocgs', qty2_cgsconv)
        s2grp.create_dataset('sums2', data=sums2)

        s12grp = f.create_group('sums12')
        s12grp.attrs.create('tocgs', qty1_cgsconv * qty2_cgsconv)
        s12grp.create_dataset('sums12', data=sums12)

        ctgrp = f.create_group('particle_counts')
        ctgrp.create_dataset('particle_counts', data=counts)

        crgrp = f.create_group('normed_average')
        crgrp.create_dataset('normed_average', data=navs)

        rgrp = f.create_group('rbins')
        rgrp.attrs.create('units', np.string_(rbins_unit))
        rgrp.create_dataset('rbins', data=rbins)
        rgrp.create_dataset('rbins_cm', data=rbins_simu * rcen_tocgs)
    return navs







    


