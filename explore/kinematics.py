'''
helper functions for explorations in matplotlib; 
possibly expandable into proper scripts
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import mpl_toolkits.mplot3d as m3d  

import fire_an.mainfunc.get_qty as gq
import fire_an.mainfunc.haloprop as hp
import fire_an.makeplots.plot_utils as pu
import fire_an.readfire.readin_fire_data as rfd
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.opts_locs as ol

# 10000 points is too many for 3D interactive plotting on quest
def getpath(simname):
    dp2 = '_'.join(simname.split('_')[:2])
    if dp2.startswith('m13h02_'):
        dp2 = dp2.replace('m13h02', 'm13h002')
    simpath = '/'.join([ol.simdir_fire, dp2, simname]) 
    return simpath

def getspacefilter(simname, snapnum, maxradius_rvir, parttype=0):
    simpath = getpath(simname)
    snap = rfd.get_Firesnap(simpath, snapnum)
    hdata = hp.get_vcom(simpath, snapnum, 1., meandef_rvir='BN98',
                        parttypes='all')
    _hd = hdata[0]
    maxrad_cm = _hd['Rvir_cm'] * maxradius_rvir
    rcen_cm = np.array([_hd['Xc_cm'], _hd['Yc_cm'], _hd['Zc_cm']])

    pos = snap.readarray_emulateEAGLE(f'PartType{parttype}/Coordinates')
    pos_toCGS = snap.toCGS

    pos -= rcen_cm[np.newaxis, :] / pos_toCGS
    rsq = np.sum(pos**2, axis=1)
    filter = rsq <= (maxrad_cm / pos_toCGS)**2
    return filter

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
               rcen=False):
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
    if rcen:
        rc_pkpc = valspv[rci] * (toCGSpv[rci] /  (c.cm_per_mpc * 1e-3))
        out = out + (rc_pkpc,)
    return out +  (rvir_pkpc,)

def get_selkin(simname, snapnum, maxradius_rvir, 
               selqtys, selqtys_args, samplesize=2000,
               parttype=0, strictselqtyss=None, 
               strictselqtyss_args=None,
               strictselqtyss_minmax=None):
    sf = getspacefilter(simname, snapnum, maxradius_rvir, parttype=parttype)
    sfd = {'filter': sf}

    pall_pkpc, vall_kmps, vrall_kmps, rvir_pkpc = getkininfo(
        simname, snapnum, filterdct=sfd, parttype=parttype, vr=True)
    pvs = []
    if strictselqtyss is None:
        strictselqtyss = [None] * len(selqtys)
        strictselqtyss_args = [None] * len(selqtys)
        strictselqtyss_minmax = [None] * len(selqtys)

    for selqty, selqty_args, sselqtys, sselqtys_args, sseltqys_minmax \
            in zip(selqtys, selqtys_args, strictselqtyss, strictselqtyss_args,
                   strictselqtyss_minmax ):
        _filter = getweightfilter(simname, snapnum, selqty, selqty_args, 
                        filterdct=sfd, samplesize=samplesize,
                        parttype=parttype, strictselqtys=sselqtys, 
                        strictselqtys_args=sselqtys_args,
                        strictselqtys_minmax=sseltqys_minmax)
        pvs.append((pall_pkpc[_filter, :], vall_kmps[_filter, :],
                    vrall_kmps[_filter]))
    return pvs, rvir_pkpc

def run_sel1(simname, snapnum):
    selqtys = ['Mass', 'Mass', 'Metal', 'Metal', 'Volume', 'ion']
    selqtys_args = [{}, {}, {'element': 'Neon'}, {'element': 'Neon'}, 
                    {}, {'ion': 'Ne8'}]
    sselqtys = [None, ['sim-direct'],
                 None, ['sim-direct'], 
                 None, None]
    sselqtys_args = [None, [{'field': 'Temperature'}], 
                     None, [{'field': 'Temperature'}], 
                     None, None]
    sselqtys_minmax = [None, [(1e5, np.inf)], 
                       None, [(1e5, np.inf)],
                       None, None]
    pv, rv = get_selkin(simname, snapnum, 1., 
                        selqtys, selqtys_args, samplesize=500,
                        parttype=0, strictselqtyss=sselqtys, 
                        strictselqtyss_args=sselqtys_args,
                        strictselqtyss_minmax=sselqtys_minmax)
    pvs, rvs = get_selkin(simname, snapnum, 1., 
                          ['Mass'], [{}], samplesize=300,
                          parttype=4)
    pvsc, rvsc = get_selkin(simname, snapnum, 0.05, 
                            ['Mass'], [{}], samplesize=100,
                            parttype=4)
    labels = ['Mass', 'Mass > 1e5 K', 'Neon', 'Neon > 1e5 K', 
              'Volume', 'Ne8', 'stars', 'stars < 0.05 Rvir']
    return pv + pvs + pvsc, labels

def run_sel2(simname, snapnum):
    selqtys = ['Volume', 'Mass', 'Mass', 'Metal'] + ['ion'] * 4
    selqtys_args = [{}, {}, {}, {'element': 'Neon'},
                    {'ion': 'Ne8'}, {'ion': 'O6'}, {'ion': 'O7'},
                    {'ion': 'O8'}]
    sselqtys = [None, None, ['sim-direct'], ['sim-direct']] + [None] * 4
    sselqtys_args = [None, None, [{'field': 'Temperature'}], 
                     [{'field': 'Temperature'}]] + [None] * 4
    sselqtys_minmax = [None, None, [(1e5, np.inf)], [(1e5, np.inf)]] \
                      + [None] * 4
    pv, rv = get_selkin(simname, snapnum, 1., 
                        selqtys, selqtys_args, samplesize=500,
                        parttype=0, strictselqtyss=sselqtys, 
                        strictselqtyss_args=sselqtys_args,
                        strictselqtyss_minmax=sselqtys_minmax)
    labels = ['Volume', 'Mass', 'Mass > 1e5 K', 'Neon > 1e5 K', 
              'Ne8', 'O6', 'O7', 'O8']
    return pv, labels

def quiverplot(pos, vel, alpha=0.2, vscale=0.1):
    fig = plt.figure()
    ax = m3d.Axes3D(fig)
    _vel = vscale * vel
    ax.quiver(pos[:, 0], pos[:, 1], pos[:, 2],
              _vel[:, 0], _vel[:, 1], _vel[:, 2],
              alpha=alpha)
    plt.show()

def quiverplots(posvels, alpha=0.2, vscales=0.1, axtitles=None,
                outname=None, title=None):
    ncmax = 4
    npanels = len(posvels)
    ncols = min(ncmax, npanels)
    nrows = (npanels - 1) // ncols + 1
    panelsize = 4.
    fig = plt.figure(figsize=(panelsize * ncols, panelsize * nrows))
    axes = [fig.add_subplot(nrows, ncols, i + 1, projection='3d') 
            for i in range(len(posvels))]
    if not hasattr(vscales, '__len__'):
        vscales = [vscales] * npanels
    vmax = max([np.max(np.abs(_pv[2])) for _pv in posvels])
    vmin = -1. * vmax 
    for i, (_p, _v, _vr) in enumerate(posvels):
        ax = axes[i]
        __v = vscales[i] * _v
        __p = np.copy(_p)
        cmap = mpl.cm.get_cmap('cool_r')
        cvals = cmap((_vr - vmin) / (vmax  - vmin))
        cvals[:, 3] = alpha
        out = ax.quiver(__p[:, 0], __p[:, 1], __p[:, 2],
                        __v[:, 0], __v[:, 1], __v[:, 2],
                        colors=cvals)
        # head color mismatch on quest + system anaconda3 
        # is a matplotlib bug.
        # np.repeat(cvals, 3, axis=0)
        # np.tile(cvals, (3, 1))
        if axtitles is not None:
            ax.set_title(axtitles[i])
    if title is not None:
        fig.suptitle(title)
    if outname is not None:
        outdir = '/projects/b1026/nastasha/imgs/vel3dcomp/3dplots_clean2/'
        plt.savefig(outdir + outname, bbox_inches='tight')
    plt.show()

def runplots_selset(hset='m12', selset=2):
    if hset == 'm12':
        sims = [#('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                # '_sdp1e10_gacc31_fa0.5'),
                #('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                # '_sdp2e-4_gacc31_fa0.5'),
                #('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                # '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                #('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                # '_sdp1e10_gacc31_fa0.5'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                 '_sdp2e-4_gacc31_fa0.5'),
                ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                 '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ]
    else:
        sims = sl.m13_nobh_clean2 + sl.m13_agnnocr_clean2 \
               + sl.m13_agncr_clean2
    snaps_hr = [sl.snaps_hr[0], sl.snaps_hr[-1]]
    snaps_sr = [sl.snaps_sr[0], sl.snaps_sr[-1]]
    hrset = sl.m12_hr_all2 + sl.m13_hr_all2
    srset = sl.m12_sr_all2 + sl.m13_sr_all2
    zs = [1.0, 0.5]
    zstrs = ['1p0', '0p5']

    #outdir = '/projects/b1026/nastasha/imgs/vel3dcomp/3dplots_clean2/'
    if selset == 1:
        outname_temp = 'pv3d_try2_{ic}_{phys}_z{zstr}_vscale_0p2_0p01.pdf'
        title_temp = ('{ic} {phys} z={z:.1f}, pos.: pkpc, '
                      'vel: km/s * (0.2, 0.01)')
        vscales = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.01]
    elif selset == 2:
        outname_temp = 'pv3d_try3_{ic}_{phys}_z{zstr}_vscale_0p2.pdf'
        title_temp = '{ic} {phys} z={z:.1f}, pos.: pkpc, vel: km/s * 0.2'
        vscales = [0.2] * 8
    for sim in sims:
        ic = sim.split('_')[0]
        phys = ('noBH' if '_sdp1e10_' in sim 
                else 'AGN-CR' if '_MHDCRspec1_' in sim 
                else 'AGN-noCR')
        for zi in range(len(zs)):
            zstr = zstrs[zi]
            zv = zs[zi]
            snap = (snaps_hr[zi] if sim in hrset 
                    else snaps_sr[zi] if sim in srset 
                    else None)
            outname = outname_temp.format(ic=ic, phys=phys, zstr=zstr)
            title = title_temp.format(ic=ic, phys=phys, z=zv)
            
            if selset == 1:
                pvs, labels = run_sel1(sim, snap)
            elif selset == 2:
                pvs, labels = run_sel2(sim, snap)
            quiverplots(pvs, alpha=0.2, vscales=vscales, axtitles=labels,
                        outname=outname, title=title)

def calcangmomprofile(simname, snapnum, rbins_rvir,
                      wtqtys, wtqtys_args, 
                      selqtys=None, selqtys_args=None, selqtys_minmax=None,
                      parttype=0):
    
    sf = getspacefilter(simname, snapnum, rbins_rvir[-1], parttype=parttype)
    sfd = {'filter': sf}
    pall_pkpc, vall_kmps, vrall_kmps, rcall_pkpc, rvir_pkpc = getkininfo(
        simname, snapnum, filterdct=sfd, parttype=parttype, vr=True, rcen=True)
    specL = np.cross(pall_pkpc, vall_kmps, axis=1)
    maxspecL = rcall_pkpc * vrall_kmps
    
    rbins_pkpc = np.array(rbins_rvir) * rvir_pkpc
    rinds = np.searchsorted(rbins_pkpc, rcall_pkpc) - 1 
    if selqtys is None:
        selqtys = [None] * len(wtqtys)
        selqtys_args = [None] * len(wtqtys)
        selqtys_minmax = [None] * len(wtqtys)
    wtd_specL = []
    wtd_maxspecL = []
    wts_tot = []
    for wtqty, wtqty_args, _selqtys, _selqtys_args, _selqtys_minmax in zip(
            wtqtys, wtqtys_args, selqtys, selqtys_args, selqtys_minmax):
        _specL_perbin = []
        _maxspecL_perbin = []
        _wtvals_perbin = []
        if _selqtys is not None:
            filter = genfilter(simname, snapnum, _selqtys, 
                               _selqtys_args, _selqtys_minmax, 
                               parttype=parttype, filterdct=sfd)
            setzero = np.logical_not(filter)
        else:
            setzero = slice(0, 0, 1) # select nothing
        wtvals, wtvals_toCGS, wtvals_todoc = gq.get_qty(
            snapnum, parttype, wtqty, wtqty_args, filterdct=sfd)
        wtvals[setzero] = 0.
        for ri in range(len(rbins_pkpc) - 1):
            rsel = rinds == ri
            _wt = wtvals[rsel]
            wtot = np.sum(_wt)
            _specL_perbin.append(np.sum(specL[rsel] * _wt) / wtot)
            _maxspecL_perbin.append(np.sum(maxspecL[rsel] * _wt) / wtot)
            _wtvals_perbin.append(wtot)
        wtd_specL.append(_specL_perbin)
        wtd_maxspecL.append(_maxspecL_perbin)
        wts_tot.append(_wtvals_perbin)
    return rbins_rvir, wtd_specL, wtd_maxspecL, wts_tot

    

    
    



    

    


    
    







