import h5py
import numpy as np

import fire_an.mainfunc.get_qty as gq
import fire_an.mainfunc.haloprop as hp
import fire_an.readfire.readin_fire_data as rfd
import fire_an.utils.constants_and_units as c
import fire_an.utils.h5utils as h5u
import fire_an.utils.opts_locs as ol


def getpath(simname):
    dp2 = '_'.join(simname.split('_')[:2])
    if dp2.startswith('m13h02_'):
        dp2 = dp2.replace('m13h02', 'm13h002')
    simpath = '/'.join([ol.simdir_fire, dp2, simname]) 
    return simpath

def genfilter(simname: str, snapnum: int, 
              qtys: list, qtys_args: list,
              qtys_minmax_cgs: list,
              parttype=0, filterdct=None):
    '''
    selects >= min, < max each qty
    each element in qtys, qtys_args is passed to get_qty
    qtys: list[str], 
    qtys_args: list[dict],
    qtys_minmax_cgs: list[tuple]

    returns:
    --------
    The filter array (boolean, size equal to the number of True
    elements in the input filter array, otherwise equal to the 
    number of particles of the specified parttype)
    '''
    simpath = getpath(simname)
    snap = rfd.get_Firesnap(simpath, snapnum)
    filter = None
    for qty, qty_args, qty_minmax_cgs in zip(
            qtys, qtys_args, qtys_minmax_cgs):
        mincgs = qty_minmax_cgs[0]
        maxcgs = qty_minmax_cgs[1]
        if ((mincgs is None or mincgs == -np.inf)
             and (maxcgs is None or maxcgs == np.inf)):
            continue
        selvals, selvals_toCGS, selvals_todoc = gq.get_qty(
            snap, parttype, qty, qty_args, filterdct=filterdct)
        if mincgs is not None and mincgs != -np.inf:
            if filter is None:
                filter = (selvals >= mincgs / selvals_toCGS)
            else:
                filter &= (selvals >= mincgs / selvals_toCGS)
        if maxcgs is not None and maxcgs != np.inf:
            if filter is None:
                filter = (selvals >= maxcgs / selvals_toCGS)
            else:
                filter &= (selvals >= maxcgs / selvals_toCGS)
    return filter
    
def getweightfilter(simname, snapnum, selqty, selqty_args, 
                    filterdct=None, samplesize=1000,
                    parttype=0, 
                    strictselqtys=None, strictselqtys_args=None,
                    strictselqtys_minmax=None):
    '''
    select random resolution elements with probability equal
    to fraction of selqty in each particle 
    (so please don't use for temperatures or something)
    returns indices for selected elements
    '''
    simpath = getpath(simname)
    snap = rfd.get_Firesnap(simpath, snapnum)
    if strictselqtys is not None:
        filter = genfilter(simname, snapnum, strictselqtys, 
                           strictselqtys_args, strictselqtys_minmax,
                           parttype=parttype, filterdct=filterdct)
        setzero = np.logical_not(filter)
    else:
        setzero = slice(0, 0, 1) # select nothing
    selvals, selvals_toCGS, selvals_todoc = gq.get_qty(
        snap, parttype, selqty, selqty_args, filterdct=filterdct)
    selvals[setzero] = 0.
    normedvals = selvals / np.sum(selvals)
    out = np.random.choice(len(normedvals), size=(samplesize), 
                           replace=False, p=normedvals)
    return out

def getkininfo(simname, snapnum, filterdct=None, parttype=0, vr=False,
               vtot=False, rcen=False):
    '''
    returns positions in pkpc, velocities in km/s, and Rvir in pkpc
    '''
    simpath = getpath(simname)
    snap = rfd.get_Firesnap(simpath, snapnum)
    hdata = hp.get_vcom(simpath, snapnum, 1., meandef_rvir='BN98',
                        parttypes='all')
    vkeys = ['VXcom_cmps', 'VYcom_cmps', 'VZcom_cmps']
    vcen_cmps = np.array([hdata[0][key] for key in vkeys])
    pkeys = ['Xc_cm', 'Yc_cm', 'Zc_cm']
    cen_cm = np.array([hdata[0][key] for key in pkeys])
    rvir_cm = hdata[0]['Rvir_cm']
    rvir_pkpc = rvir_cm / (1e-3 * c.cm_per_mpc)
    
    qty = 'coords'
    coords_todo = [{'pos': 'allcart'}, {'vel': 'allcart'}]
    basei = 2
    if vr:
        coords_todo = coords_todo + [{'vel': 'vrad'}]
        vri = basei
        basei += 1
    if vtot:
        coords_todo = coords_todo + [{'vel': 'vtot'}]
        vti = basei
        basei += 1
    if rcen:
        coords_todo = coords_todo + [{'pos': 'rcen'}]
        rci = basei
    qty_args = {'center_cm': cen_cm, 'vcen_cmps': vcen_cmps,
                'multiple': coords_todo}
    valspv, toCGSpv, todoc_pv = gq.get_qty(
        snap, parttype, qty, qty_args, filterdct=filterdct)
    pos_pkpc = valspv[0] * (toCGSpv[0] / (c.cm_per_mpc * 1e-3))
    vel_kmps = valspv[1] * (toCGSpv[1] / (1e5))
    out = (pos_pkpc, vel_kmps)
    if vr:
        vr_kmps = valspv[vri] * (toCGSpv[vri] / (1e5))
        out = out + (vr_kmps,)
    if vtot:
        vt_kmps = valspv[vti] * (toCGSpv[vti] / (1e5))
        out = out + (vt_kmps,)
    if rcen:
        rc_pkpc = valspv[rci] * (toCGSpv[rci] /  (c.cm_per_mpc * 1e-3))
        out = out + (rc_pkpc,)
    return out +  (rvir_pkpc,)

def getspacefilter_box(simname, snapnum, boxradii_rvir, parttype=0):
    simpath = getpath(simname)
    snap = rfd.get_Firesnap(simpath, snapnum)
    hdata = hp.get_vcom(simpath, snapnum, 1., meandef_rvir='BN98',
                        parttypes='all')
    _hd = hdata[0]
    maxrad_cm = _hd['Rvir_cm'] * boxradii_rvir
    rcen_cm = np.array([_hd['Xc_cm'], _hd['Yc_cm'], _hd['Zc_cm']])

    pos = snap.readarray_emulateEAGLE(f'PartType{parttype}/Coordinates')
    pos_toCGS = snap.toCGS

    pos -= rcen_cm[np.newaxis, :] / pos_toCGS
    rsq = np.sum(pos**2, axis=1)
    filter = rsq <= ((maxrad_cm / pos_toCGS)**2)[np.newaxis, :]
    return filter, hdata

def save_selkin_box(simname, snapnum, boxradii_rvir, 
                    selqtys, selqtys_args, outname,
                    samplesize=100,
                    parttype=0, strictselqtyss=None, 
                    strictselqtyss_args=None,
                    strictselqtyss_minmax=None):
    sf, hdata = getspacefilter_box(simname, snapnum, boxradii_rvir, 
                                   parttype=parttype)
    sfd = {'filter': sf}

    pall_pkpc, vall_kmps, vrall_kmps, rvir_pkpc = getkininfo(
        simname, snapnum, filterdct=sfd, parttype=parttype, vr=True)
    pvs = []
    if strictselqtyss is None:
        strictselqtyss = [None] * len(selqtys)
        strictselqtyss_args = [None] * len(selqtys)
        strictselqtyss_minmax = [None] * len(selqtys)

    for selqty, selqty_args, sselqtys, sselqtys_args, sseltqys_minmax \
            in zip(selqtys, selqtys_args, strictselqtyss,
                   strictselqtyss_args, strictselqtyss_minmax):
        _filter = getweightfilter(simname, snapnum, selqty, selqty_args, 
                        filterdct=sfd, samplesize=samplesize,
                        parttype=parttype, strictselqtys=sselqtys, 
                        strictselqtys_args=sselqtys_args,
                        strictselqtys_minmax=sseltqys_minmax)
        pvs.append((pall_pkpc[_filter, :], vall_kmps[_filter, :],
                    vrall_kmps[_filter]))
    with h5py.File(outname, 'w') as f:
        hed = f.create_group('Header')
        hed.attrs.create('simname', np.string_(simname))
        hed.attrs.create('snapnum', snapnum)
        hed.attrs.create('boxradii_rvir', boxradii_rvir)
        hed.attrs.create('rvir_pkpc', rvir_pkpc)
        hed.attrs.create('samplesize', samplesize)
        hed.attrs.create('parttype', parttype)
        hgrp = hed.create_group('halodata')
        h5u.savedict_hdf5(hgrp, hdata)
        
        for si, (selqty, selqty_args, sselqtys, sselqtys_args, 
                 sseltqys_minmax) \
                in enumerate(zip(selqtys, selqtys_args, strictselqtyss,
                                 strictselqtyss_args, strictselqtyss_minmax)):
            sgrp = f.create_group(f'selection_{si}')
            sgrp.create_dataset('pos_norot_pkpc', data=pvs[si][0])
            sgrp.create_dataset('vel_norot_kmps', data=pvs[si][1])
            sgrp.create_dataset('vrad_kmps', data=pvs[si][2])

            sgrp.attrs.create('selqty', np.string_(selqty))
            if selqty_args is None:
                sgrp.attrs.create('selqty_args', np.string_('None'))
            else:
                sgrp.attrs.create('selqty_args', np.string_('dict'))
                _grp = sgrp.attrs.create('selqty_args_dict')
                h5u.savedict_hdf5(_grp, selqty_args)
            if sselqtys is None:
                sgrp.attrs.create('strictselqtys', np.string_('None'))
                sgrp.attrs.create('strictselqtys_args', 
                                  np.string_('None'))
                sgrp.attrs.create('strictselqtys_minmax', 
                                  np.string_('None'))
            else:
                sgrp.attrs.create('strictselqtys', 
                                  np.array([np.string_(qt) 
                                            for qt in sselqtys]))
                for ssi, (_ssel_args, _ssel_mm) in \
                        enumerate(zip(sselqtys_args, sseltqys_minmax)):
                    _grp = sgrp.create_group(f'strictsel_index_{ssi}')
                    if _ssel_args is None:
                        sgrp.attrs.create('strictselqty_args', 
                                          np.string_('None'))
                    else:
                        sgrp.attrs.create('strictselqty_args', 
                                          np.string_('dict'))
                        __grp = sgrp.attrs.create('strictselqty_args_dict')
                        h5u.savedict_hdf5(__grp, _ssel_args)
                        __grp.attrs.create('min_cgs', _ssel_mm[0])
                        __grp.attrs.create('max_cgs', _ssel_mm[1])