import h5py
import numpy as np

import fire_an.mainfunc.get_qty as gq
import fire_an.mainfunc.haloprop as hp
import fire_an.readfire.readin_fire_data as rfd
import fire_an.simlists as sl
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
    todocs_strict = []
    for qty, qty_args, qty_minmax_cgs in zip(
            qtys, qtys_args, qtys_minmax_cgs):
        mincgs = qty_minmax_cgs[0]
        maxcgs = qty_minmax_cgs[1]
        if ((mincgs is None or mincgs == -np.inf)
             and (maxcgs is None or maxcgs == np.inf)):
            continue
        selvals, selvals_toCGS, selvals_todoc = gq.get_qty(
            snap, parttype, qty, qty_args, filterdct=filterdct)
        todocs_strict.append(selvals_todoc)
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
    return filter, todocs_strict
    
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
        filter, strictsels_todoc = genfilter(simname, snapnum, strictselqtys,
            strictselqtys_args, strictselqtys_minmax, parttype=parttype, 
            filterdct=filterdct)
        setzero = np.logical_not(filter)
    else:
        setzero = slice(0, 0, 1) # select nothing
        strictsels_todoc = []
    selvals, selvals_toCGS, selvals_todoc = gq.get_qty(
        snap, parttype, selqty, selqty_args, filterdct=filterdct)
    selvals[setzero] = 0.
    normedvals = selvals / np.sum(selvals)
    out = np.random.choice(len(normedvals), size=(samplesize), 
                           replace=False, p=normedvals)
    return out, selvals_todoc, strictsels_todoc

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
    return out + (rvir_pkpc,)

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
    ra = np.abs(pos)
    filter = np.all(ra <= (maxrad_cm / pos_toCGS)[np.newaxis, :],
                    axis=1)
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
    selvalss_todoc = []
    strictselss_todoc = []
    if strictselqtyss is None:
        strictselqtyss = [None] * len(selqtys)
        strictselqtyss_args = [None] * len(selqtys)
        strictselqtyss_minmax = [None] * len(selqtys)
    
    for selqty, selqty_args, sselqtys, sselqtys_args, sseltqys_minmax \
            in zip(selqtys, selqtys_args, strictselqtyss,
                   strictselqtyss_args, strictselqtyss_minmax):
        _filter, selvals_todoc, strictsels_todoc = getweightfilter(
            simname, snapnum, selqty, selqty_args, filterdct=sfd, 
            samplesize=samplesize, parttype=parttype, 
            strictselqtys=sselqtys, strictselqtys_args=sselqtys_args,
            strictselqtys_minmax=sseltqys_minmax)
        selvalss_todoc.append(selvals_todoc)
        strictselss_todoc.append(strictsels_todoc)
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
        h5u.savedict_hdf5(hgrp, hdata[0])
        _hgrp = hgrp.create_group('doc')
        h5u.savedict_hdf5(_hgrp, hdata[1])
        
        for si, (selqty, selqty_args, selqty_todoc, 
                 sselqtys, sselqtys_args, sseltqys_minmax, sselqtys_todoc) \
                in enumerate(zip(selqtys, selqtys_args, selvalss_todoc,
                                 strictselqtyss, strictselqtyss_args, 
                                 strictselqtyss_minmax, strictselss_todoc)):
            sgrp = f.create_group(f'selection_{si}')
            sgrp.create_dataset('pos_norot_pkpc', data=pvs[si][0])
            sgrp.create_dataset('vel_norot_kmps', data=pvs[si][1])
            sgrp.create_dataset('vrad_kmps', data=pvs[si][2])

            sgrp.attrs.create('selqty', np.string_(selqty))
            if selqty_args is None:
                sgrp.attrs.create('selqty_args', np.string_('None'))
            else:
                sgrp.attrs.create('selqty_args', np.string_('dict'))
                _grp = sgrp.create_group('selqty_args_dict')
                h5u.savedict_hdf5(_grp, selqty_args)
            _grp = sgrp.create_group('selqty_doc')
            h5u.savedict_hdf5(_grp, selqty_todoc)
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
                for ssi, (_ssel_args, _ssel_mm, _ssel_doc) in \
                        enumerate(zip(sselqtys_args, sseltqys_minmax,
                                      sselqtys_todoc)):
                    _grp = sgrp.create_group(f'strictsel_index_{ssi}')
                    __grp = _grp.create_group('strictselqty_doc')
                    h5u.savedict_hdf5(__grp, _ssel_doc)
                    if _ssel_args is None:
                        sgrp.attrs.create('strictselqty_args', 
                                          np.string_('None'))
                    else:
                        sgrp.attrs.create('strictselqty_args', 
                                          np.string_('dict'))
                        __grp = sgrp.create_group('strictselqty_args_dict')
                        h5u.savedict_hdf5(__grp, _ssel_args)
                        __grp.attrs.create('min_cgs', _ssel_mm[0])
                        __grp.attrs.create('max_cgs', _ssel_mm[1])

def run_selkin_box_clean2(opt):
    # 3 indices per sim/snap
    boxradiis_rvir = [np.array([1.5, 1.5, 0.1]),
                      np.array([0.1, 1.5, 1.5]),
                      np.array([1.5, 0.1, 1.5]),
                      ]
    axes = ['z', 'x', 'y']
    selqtys = ['Mass', 'Volume', 'Metal', 'Metal'] + ['ion'] * 4
    selqtys_args = [None, None, 
                    {'element': 'Neon'}, {'element': 'Oxygen'},
                    {'ion': 'O6', 'density': False},
                    {'ion': 'O7', 'density': False},
                    {'ion': 'O8', 'density': False},
                    {'ion': 'Ne8', 'density': False}]
    samplesize = 150
    _outfilen = ('weightedsel_{samplesize}_{sc}_snap{sn}'
                 '_shrink-sph-cen_BN98'
                 '_depth_{llos:.1f}rvir_{los}-slice_v1.hdf5')
    #outdir = '/scratch1/08466/tg877653/output/slicepv_wtdsel_clean2/'
    outdir = './'
    if opt >= 0 and opt < 72:
        # m13sr: 72 indices
        ind = opt - 0
        simnames = sl.m13_sr_clean2 # len 4
        snaps = sl.snaplists['m13_sr'] # len 6
    elif opt >= 72 and opt < 108:
        # m13hr: 36 indices
        ind = opt - 72
        simnames = sl.m13_hr_clean2 # len 2
        snaps = sl.snaplists['m13_hr'] # len 6
    elif opt >= 108 and opt < 144:
        # m12sr: 36 indices
        ind = opt - 1728
        simnames = [('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000')
                    ]
        snaps = sl.snaplists['m12_sr'] # len 6
    elif opt >= 144 and opt < 216:
        # m12hr: 72 indices
        ind = opt - 144
        simnames = [('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp2e-4_gacc31_fa0.5'),
                    ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp2e-4_gacc31_fa0.5'),
                    ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ] # len 4
        snaps = sl.snaplists['m12_hr'] # len 6
        
    simi = ind // (len(snaps) * len(axes))
    snpi = (ind % (len(snaps) * len(axes))) \
           // ( len(axes))
    axi = ind % len(axes)
    simname = simnames[simi]
    snapnum = snaps[snpi]
    axis = axes[axi]
    boxradii_rvir = boxradiis_rvir[axi]
    
    outname = outdir + _outfilen.format(sc=simname, sn=snapnum,
                                        llos=0.1, los=axis,
                                        samplesize=samplesize)
    print('file to produce: ', outname)
    save_selkin_box(simname, snapnum, boxradii_rvir, 
                    selqtys, selqtys_args, outname,
                    samplesize=samplesize,
                    parttype=0, strictselqtyss=None, 
                    strictselqtyss_args=None,
                    strictselqtyss_minmax=None)