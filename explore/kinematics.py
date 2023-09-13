'''
helper functions for explorations in matplotlib; 
possibly expandable into proper scripts
'''
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import mpl_toolkits.mplot3d as m3d  
import numpy as np

import fire_an.explore.run_partsel_pv as rpv
import fire_an.mainfunc.get_qty as gq
import fire_an.mainfunc.haloprop as hp
import fire_an.makeplots.plot_utils as pu
import fire_an.readfire.readin_fire_data as rfd
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.opts_locs as ol

# 10000 points is too many for 3D interactive plotting on quest
def getpath(simname):
    return sl.dirpath_from_simname(simname)

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

def get_selkin(simname, snapnum, maxradius_rvir, 
               selqtys, selqtys_args, samplesize=2000,
               parttype=0, strictselqtyss=None, 
               strictselqtyss_args=None,
               strictselqtyss_minmax=None):
    sf = getspacefilter(simname, snapnum, maxradius_rvir, parttype=parttype)
    sfd = {'filter': sf}

    pall_pkpc, vall_kmps, vrall_kmps, rvir_pkpc = rpv.getkininfo(
        simname, snapnum, filterdct=sfd, parttype=parttype, vr=True)
    pvs = []
    if strictselqtyss is None:
        strictselqtyss = [None] * len(selqtys)
        strictselqtyss_args = [None] * len(selqtys)
        strictselqtyss_minmax = [None] * len(selqtys)

    for selqty, selqty_args, sselqtys, sselqtys_args, sseltqys_minmax \
            in zip(selqtys, selqtys_args, strictselqtyss, strictselqtyss_args,
                   strictselqtyss_minmax ):
        _filter = rpv.getweightfilter(simname, snapnum, selqty, selqty_args, 
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
        sims = sl.m13_agnnocr_clean2 \
               + sl.m13_agncr_clean2 #sl.m13_nobh_clean2 +
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



def plotcenter_check(pos, vel, cenpos, cenvel, vscale=0.005, alpha=0.2):
    fig = plt.figure()
    ax = m3d.Axes3D(fig)
    _vel = (vel - cenvel[np.newaxis, :]) * vscale
    _pos = pos - cenpos[np.newaxis, :]
    ax.quiver(_pos[:, 0], _pos[:, 1], _pos[:, 2],
              _vel[:, 0], _vel[:, 1], _vel[:, 2],
              alpha=alpha)
    plt.show()

def plotcenter_check2(pos, cenpos, alpha=0.2):
    fig = plt.figure()
    ax = m3d.Axes3D(fig)
    _pos = pos - cenpos[np.newaxis, :]
    ax.scatter(_pos[:, 0], _pos[:, 1], _pos[:, 2], color='gray',
               alpha=alpha)
    ax.scatter([0.], [0.], [0.], color='black')
    plt.show()

def checkcentering_stars(simname, snapnum, outname=None, title=None,
                         vscale=0.005, alpha=0.1):
    simpath = getpath(simname)
    pcen_cm, vcom_cmps, todoc = hp.getcengalcen(simpath, snapnum, 
                                                startrad_rvir=0.3,
                                                vcenrad_rvir=0.05)
    halodat = hp.gethalodata_shrinkingsphere(simpath, snapnum, meandef='BN98')
    rvir_cm = halodat[0]['Rvir_cm']

    snapobj = rfd.get_Firesnap(simpath, snapnum)
    spos_simu = snapobj.readarray('PartType4/Coordinates')
    spos_toCGS = snapobj.toCGS
    stard2 = np.sum((spos_simu - pcen_cm / spos_toCGS)**2, axis=1)  
    starsel = stard2 <= (0.3 * rvir_cm / spos_toCGS)
    spos_simu = spos_simu[starsel]
    spos_pkpc = spos_simu * (spos_toCGS / (1e-3 * c.cm_per_mpc))

    svel_simu = snapobj.readarray('PartType4/Velocities')[starsel]
    svel_toCGS = snapobj.toCGS
    svel_kmps = svel_simu * (svel_toCGS / 1e5)
    
    filterdct = {'filter': starsel}
    subfilter = rpv.getweightfilter(simname, snapnum, 'Mass', {}, 
                                    filterdct=filterdct, samplesize=200,
                                    parttype=4, 
                                    strictselqtys=None,
                                    strictselqtys_args=None,
                                    strictselqtys_minmax=None)
    
    fig = plt.figure(figsize=(9., 4.5))
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    if title is not None:
        fig.suptitle(title)

    _vel = (svel_kmps[subfilter, :] -  vcom_cmps[np.newaxis, :] * 1e-5) \
           * vscale
    _pos = spos_pkpc[subfilter, :] \
           - pcen_cm[np.newaxis, :] / (1e-3 * c.cm_per_mpc)
    ax1.quiver(_pos[:, 0], _pos[:, 1], _pos[:, 2],
               _vel[:, 0], _vel[:, 1], _vel[:, 2],
               alpha=alpha)
    ax1.scatter([0.], [0.], [0.], color='black')
    ax2.scatter(_pos[:, 0], _pos[:, 1], _pos[:, 2],
                alpha=alpha)
    ax2.scatter([0.], [0.], [0.], color='black')

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')
    print('\nShowing ', simname, snapnum)
    print()
    plt.show()

def runset_centeringcheck_stars(hset='m12'):
    if hset == 'm12':
        sims = [('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                 '_sdp1e10_gacc31_fa0.5'),
                ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                 '_sdp2e-4_gacc31_fa0.5'),
                ('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                 '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                 '_sdp1e10_gacc31_fa0.5'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                 '_sdp2e-4_gacc31_fa0.5'),
                ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                 '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ]
    else:
        sims = sl.m13_agnnocr_clean2 \
               + sl.m13_agncr_clean2 + sl.m13_nobh_clean2
    snaps_hr = [sl.snaps_hr[0], sl.snaps_hr[-1]]
    snaps_sr = [sl.snaps_sr[0], sl.snaps_sr[-1]]
    hrset = sl.m12_hr_all2 + sl.m13_hr_all2
    srset = sl.m12_sr_all2 + sl.m13_sr_all2
    zs = [1.0, 0.5]
    zstrs = ['1p0', '0p5']

    outdir = '/projects/b1026/nastasha/imgs/vel3dcomp/3dplots_clean2/'
    outname_temp = 'starrecen_check_try1_{ic}_{phys}_z{zstr}_vscale_0p005.pdf'
    title_temp = ('{ic} {phys} z={z:.1f}, pos.: pkpc, '
                    'vel: km/s * 0.005')
    vscale = 0.005
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
            outname = outdir + outname
            title = title_temp.format(ic=ic, phys=phys, z=zv)

            checkcentering_stars(sim, snap, outname=outname, title=title,
                                 vscale=vscale, alpha=0.1)

def calcangmomprofile_stars(simname, snapnum, rbins_rvir=None):
    if rbins_rvir is None:
        rbins_rvir = np.arange(0., 0.205, 0.01)
    simpath = getpath(simname)
    scen_cm, svcom_cmps, todoc = hp.getcengalcen(simpath, snapnum, 
                                                 startrad_rvir=0.3,
                                                 vcenrad_rvir=0.05)
    halodat = hp.gethalodata_shrinkingsphere(simpath, snapnum, meandef='BN98')
    rvir_cm = halodat[0]['Rvir_cm']
    snapobj = rfd.get_Firesnap(simpath, snapnum)
    spos_simu = snapobj.readarray('PartType4/Coordinates')
    spos_toCGS = snapobj.toCGS
    spos_simu -= scen_cm[np.newaxis, :] / spos_toCGS
    sel = np.sum(spos_simu**2, axis=1) \
          <= (rbins_rvir[-1] * rvir_cm / spos_toCGS)**2
    spos_simu = spos_simu[sel]
    svel_simu = snapobj.readarray('PartType4/Velocities')[sel]
    svel_toCGS = snapobj.toCGS
    svel_simu -= svcom_cmps[np.newaxis, :] / svel_toCGS

    specL = np.cross(spos_simu, svel_simu, axis=1)
    rcen_simu = np.sqrt(np.sum(spos_simu**2, axis=1))
    maxspecL = rcen_simu * np.sqrt(np.sum(svel_simu**2, axis=1))
    weight = snapobj.readarray('PartType4/Masses')[sel]
    wt_toCGS = snapobj.toCGS
    del spos_simu, svel_simu

    rvir_simu = rvir_cm / spos_toCGS
    
    print('Finding bin indices')
    rbins_pkpc = np.array(rbins_rvir) * rvir_simu
    rinds = np.searchsorted(rbins_pkpc, rcen_simu) - 1 
    print(rbins_pkpc)
    del rcen_simu

    specL_perbin = []
    maxspecL_perbin = []
    wtvals_perbin = []
    for ri in range(len(rbins_pkpc) - 1):
        rsel = rinds == ri
        _wt = weight[rsel]
        wtot = np.sum(_wt)
        specL_perbin.append(np.sum(specL[rsel] * _wt[:, np.newaxis], axis=0) 
                                / wtot)
        maxspecL_perbin.append(np.sum(maxspecL[rsel] * _wt) / wtot)
        wtvals_perbin.append(wtot)
    specL_perbin = np.array(specL_perbin)
    maxspecL_perbin = np.array(maxspecL_perbin)
    wtvals_perbin = np.array(wtvals_perbin)
    specL_perbin *= (spos_toCGS * svel_toCGS) / (1e-3 * c.cm_per_mpc * 1e5)
    maxspecL_perbin *= (spos_toCGS * svel_toCGS) / (1e-3 * c.cm_per_mpc * 1e5)
    wtvals_perbin *= (wt_toCGS / c.solar_mass)
    return specL_perbin, maxspecL_perbin, wtvals_perbin, rbins_rvir


def calcangmomprofile(simname, snapnum, rbins_rvir,
                      wtqtys, wtqtys_args, 
                      selqtys=None, selqtys_args=None, selqtys_minmax=None,
                      parttype=0):
    simpath = getpath(simname)
    snap = rfd.get_Firesnap(simpath, snapnum)
    sf = getspacefilter(simname, snapnum, rbins_rvir[-1], parttype=parttype)
    sfd = {'filter': sf}
    pall_pkpc, vall_kmps, vtotall_kmps, rcall_pkpc, rvir_pkpc =\
        rpv.getkininfo(simname, snapnum, filterdct=sfd, parttype=parttype, 
                       vr=False, vtot=True, rcen=True)
    print('Calculating L, max L')
    specL = np.cross(pall_pkpc, vall_kmps, axis=1)
    maxspecL = rcall_pkpc * vtotall_kmps
    del pall_pkpc, vall_kmps, vtotall_kmps
    
    print('Finding bin indices')
    rbins_pkpc = np.array(rbins_rvir) * rvir_pkpc
    rinds = np.searchsorted(rbins_pkpc, rcall_pkpc) - 1 
    print(rbins_pkpc)
    del rcall_pkpc

    print('Starting wt loop')
    if selqtys is None:
        selqtys = [None] * len(wtqtys)
        selqtys_args = [None] * len(wtqtys)
        selqtys_minmax = [None] * len(wtqtys)
    wtd_specL = []
    wtd_maxspecL = []
    wts_tot = []
    for wtqty, wtqty_args, _selqtys, _selqtys_args, _selqtys_minmax in zip(
            wtqtys, wtqtys_args, selqtys, selqtys_args, selqtys_minmax):
        print(f'Loop: {wtqty}, {wtqty_args}, subsample: {_selqtys}')
        _specL_perbin = []
        _maxspecL_perbin = []
        _wtvals_perbin = []
        if _selqtys is not None:
            filter = rpv.genfilter(simname, snapnum, _selqtys, 
                                   _selqtys_args, _selqtys_minmax, 
                                   parttype=parttype, filterdct=sfd)
            setzero = np.logical_not(filter)
        else:
            setzero = slice(0, 0, 1) # select nothing
        wtvals, wtvals_toCGS, wtvals_todoc = gq.get_qty(
            snap, parttype, wtqty, wtqty_args, filterdct=sfd)
        wtvals[setzero] = 0.
        for ri in range(len(rbins_pkpc) - 1):
            rsel = rinds == ri
            _wt = wtvals[rsel]
            wtot = np.sum(_wt)
            _specL_perbin.append(np.sum(specL[rsel] * _wt[:, np.newaxis], axis=0) 
                                 / wtot)
            _maxspecL_perbin.append(np.sum(maxspecL[rsel] * _wt) / wtot)
            _wtvals_perbin.append(wtot)
        wtd_specL.append(_specL_perbin)
        wtd_maxspecL.append(_maxspecL_perbin)
        wts_tot.append(_wtvals_perbin)
    return rbins_rvir, wtd_specL, wtd_maxspecL, wts_tot

def getangmomprofile_gvs(simname, snapnum):

    wtqtys_g = ['Mass', 'Mass', 'Volume']
    wtqtys_args_g = [{}, {}, {}]
    selqtys_g = [None, ['sim-direct'], None]
    selqtys_args_g = [None, [{'field': 'Temperature'}], None]
    selqtys_minmax_g = [None, [(1e5, np.inf)], None]
    rbins = np.array([0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 
                      0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])
    
    rbins_rvir, wtd_specL_g, wtd_maxspecL_g, wts_tot_g = \
        calcangmomprofile(simname, snapnum, rbins,
                          wtqtys_g, wtqtys_args_g, 
                          selqtys=selqtys_g, selqtys_args=selqtys_args_g, 
                          selqtys_minmax=selqtys_minmax_g,
                          parttype=0)
    rbins_rvir, wtd_specL_s, wtd_maxspecL_s, wts_tot_s = \
        calcangmomprofile(simname, snapnum, rbins,
                          ['Mass'], [{}], 
                          selqtys=None, selqtys_args=None, 
                          selqtys_minmax=None,
                          parttype=4)
    labels = ['Mass', 'Mass > 1e5 K', 'Volume', 'Stars']
    out = (rbins_rvir, 
           np.array(wtd_specL_g + wtd_specL_s), 
           np.array(wtd_maxspecL_g + wtd_maxspecL_s),
           np.array(wts_tot_g +  wts_tot_s),
           labels)
    return out

def plot_angmomprofile(rbins_rvir, wtd_specL, wtd_maxspecL, wts_tot,
                       specL_ref,
                       wtlabels, wtcolors, wtstyles, title=None):
    # plot delta(wt, ref wt at cen)
    # plot delta(wt, prev. weight)
    # plot planar fractions
    # plot radial dist.
    # different weights in different colors
    # different plots for different haloes
    # m12f AGN-CR: looks like a dang mess.
    fontsize = 12

    fig = plt.figure(figsize=(7., 6.))
    grid = gsp.GridSpec(ncols=2, nrows=3, hspace=0.2,
                        wspace=0.3,
                        height_ratios=[3., 3., 1.])
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    rc = 0.5 * (rbins_rvir[:-1] + rbins_rvir[1:])

    ax = fig.add_subplot(grid[0, 0])
    for i, _sL in enumerate(wtd_specL):
        print(_sL)
        print(specL_ref)
        normin = np.einsum('ij,j->i', _sL, specL_ref) \
                 / np.sqrt(np.sum(_sL**2, axis=1) * np.sum(specL_ref**2))
        angle = np.arccos(normin)
        ax.plot(rc, angle, label=wtlabels[i], color=wtcolors[i], 
                **wtstyles[i])
    ax.set_ylabel('angle with ref. L [rad]')
    ax.tick_params(which='both', top=True, right=True, 
                   labelsize=fontsize - 1., direction='in')
    
    ax = fig.add_subplot(grid[0, 1])
    for i, _sL in enumerate(wtd_specL):
        normin = np.einsum('ij,ij->i', _sL[1:], _sL[:-1]) \
                 / np.sqrt(np.sum(_sL[1:]**2, axis=1) \
                           * np.sum(_sL[:-1]**2, axis=1))
        angle = np.arccos(normin)
        print(normin)
        ax.plot(rc[1:], angle, label=wtlabels[i], color=wtcolors[i], 
                **wtstyles[i])
    ax.set_ylabel('L angle with prev. radius [rad]')
    ax.tick_params(which='both', top=True, right=True, 
                   labelsize=fontsize - 1., direction='in')
    
    ax = fig.add_subplot(grid[1, 0])
    for i, (_sL, _sLm) in enumerate(zip(wtd_specL, wtd_maxspecL)):
        frac = np.sqrt(np.sum(_sL**2, axis=1)) / _sLm
        ax.plot(rc, frac, label=wtlabels[i], color=wtcolors[i], 
                **wtstyles[i])
    ax.set_ylabel('$\\Sigma_{i} \\, w_i \\, \\vec{r} \\times '
                  '\\vec{v} \\,/\\,'
                  '\\Sigma_{i} \\, w_i \\, |\\vec{r}| \\, |\\vec{v}|$',
                  fontsize=fontsize)
    ax.tick_params(which='both', top=True, right=True, 
                   labelsize=fontsize - 1., direction='in')
    ax.set_xlabel('$\\mathrm{r}_{\\mathrm{3D}} \\,'
                  ' [\\mathrm{R}_{\\mathrm{vir}}]$',
                  fontsize=fontsize)
    
    ax = fig.add_subplot(grid[1, 1])
    for i, wt in enumerate(wts_tot):
        frac = wt / np.sum(wt) / np.diff(rbins_rvir)
        ax.plot(rc, frac, label=wtlabels[i], color=wtcolors[i], 
                **wtstyles[i])
    ax.set_ylabel('$\\Delta \\, \\mathrm{weight}\\,\\mathrm{frac.} \\,//\\,'
                  '\\Delta \\, \\mathrm{r}$',
                  fontsize=fontsize)
    ax.tick_params(which='both', top=True, right=True, 
                   labelsize=fontsize - 1., direction='in')
    ax.set_xlabel('$\\mathrm{r}_{\\mathrm{3D}} \\,'
                  ' [\\mathrm{R}_{\\mathrm{vir}}]$',
                  fontsize=fontsize)
    ax.set_yscale('log')

    lax = fig.add_subplot(grid[2, :])
    lax.axis('off')
    handles = [mlines.Line2D((), (), color=c, label=lab, **st)
               for lab, c, st in zip(wtlabels, wtcolors, wtstyles)]
    lax.legend(handles=handles, fontsize=fontsize, ncol=3, 
               loc='upper center', bbox_to_anchor=(0.5, 0.65))

def runset_angmomprof_starrecen(hset='m12'):
    if hset == 'm12':
        sims = [('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                 '_sdp1e10_gacc31_fa0.5'),
                ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                 '_sdp2e-4_gacc31_fa0.5'),
                ('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                 '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                 '_sdp1e10_gacc31_fa0.5'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                 '_sdp2e-4_gacc31_fa0.5'),
                ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                 '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ]
    else:
        sims = sl.m13_agnnocr_clean2 \
               + sl.m13_agncr_clean2 + sl.m13_nobh_clean2
    snaps_hr = [sl.snaps_hr[0], sl.snaps_hr[-1]]
    snaps_sr = [sl.snaps_sr[0], sl.snaps_sr[-1]]
    hrset = sl.m12_hr_all2 + sl.m13_hr_all2
    srset = sl.m12_sr_all2 + sl.m13_sr_all2
    zs = [1.0, 0.5]
    zstrs = ['1p0', '0p5']

    outdir = '/projects/b1026/nastasha/imgs/vel3dcomp/3dplots_clean2/'
    outname_temp = 'angmomprof_starrecen_try1_{ic}_{phys}_z{zstr}.pdf'
    title_temp = ('{ic} {phys} z={z:.1f}; ref: weighted. av. < 0.1 Rvir')
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
            outname = outdir + outname
            title = title_temp.format(ic=ic, phys=phys, z=zv)
            
            specL_perbin, maxspecL_perbin, wtvals_perbin, rbins_rvir = \
                calcangmomprofile_stars(sim, snap, rbins_rvir=None)
            ref = np.sum(specL_perbin[:10, :] 
                         * wtvals_perbin[:10, np.newaxis], axis=0)\
                  / np.sum(wtvals_perbin[:10])
            plot_angmomprofile(rbins_rvir, [specL_perbin], [maxspecL_perbin], 
                               [wtvals_perbin], ref,
                               ['stars'], ['blue'], [{}], title=title)
            plt.savefig(outname, bbox_inches='tight')
            plt.show()

def setup_vizcheck_outflow_Lstarrecen(simname, snapnum, galrad_rvir=0.1):
    '''
    Just do Volume and Mass > 1e5 K, those seemed to show
    outflows best in terms of basic properties.
    '''
    specL_perbin, maxspecL_perbin, wtvals_perbin, rbins_rvir = \
        calcangmomprofile_stars(simname, snapnum, 
                                rbins_rvir=np.array([0., galrad_rvir]))
    Ldir = specL_perbin[0] / np.sqrt(np.sum(specL_perbin**2))
    
    maxradius_rvir = 1.1
    selqtys = ['Mass', 'Volume']
    selqtys_args = [{}, {}]
    strictselqtyss = [['sim-direct'], None]
    strictselqtyss_args = [[{'field': 'Temperature'}], None]
    strictselqtyss_minmax = [[(1e5, np.inf)], None]
    pvs, rvir_pkpc = get_selkin(simname, snapnum, maxradius_rvir, 
                                selqtys, selqtys_args, samplesize=300,
                                parttype=0, strictselqtyss=strictselqtyss, 
                                strictselqtyss_args=strictselqtyss_args,
                                strictselqtyss_minmax=strictselqtyss_minmax)
    return pvs, Ldir, rvir_pkpc

def plot_vizcheck_outflow_Lstarrecen(posvels, Ldir, rvir_pkpc,
                                     title=None, axtitles=None, outname=None,
                                     alpha=0.2, vscales=0.1):
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
        Lvec = rvir_pkpc * Ldir
        Lvec2 = - rvir_pkpc * Ldir
        ax.quiver([0., 0.], [0., 0.], [0., 0.],
                  [Lvec[0], Lvec2[0]], 
                  [Lvec[1], Lvec2[1]], 
                  [Lvec[2], Lvec2[2]], color='black')
        if axtitles is not None:
            ax.set_title(axtitles[i])
    if title is not None:
        fig.suptitle(title)
    if outname is not None:
        outdir = '/projects/b1026/nastasha/imgs/vel3dcomp/3dplots_clean2/'
        plt.savefig(outdir + outname, bbox_inches='tight')
    plt.show()   

def run_vizcheck_outflow_Lstarrecen(hset='m12'):
    if hset == 'm12':
        sims = [('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                 '_sdp1e10_gacc31_fa0.5'),
                ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                 '_sdp2e-4_gacc31_fa0.5'),
                ('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                 '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                 '_sdp1e10_gacc31_fa0.5'),
                ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                 '_sdp2e-4_gacc31_fa0.5'),
                ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                 '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                ]
    else:
        sims = sl.m13_agnnocr_clean2 \
               + sl.m13_agncr_clean2 + sl.m13_nobh_clean2
    snaps_hr = [sl.snaps_hr[0], sl.snaps_hr[-1]]
    snaps_sr = [sl.snaps_sr[0], sl.snaps_sr[-1]]
    hrset = sl.m12_hr_all2 + sl.m13_hr_all2
    srset = sl.m12_sr_all2 + sl.m13_sr_all2
    zs = [1.0, 0.5]
    zstrs = ['1p0', '0p5']
    axtitles = ['Mass < 1e5 K', 'Volume']

    outdir = '/projects/b1026/nastasha/imgs/vel3dcomp/3dplots_clean2/'
    outname_temp = 'dircheck_outflow_Lstarrecen_try3_{ic}_{phys}_z{zstr}.pdf'
    title_temp = ('{ic} {phys} z={z:.1f}; pos: pkpc, vel: km/s * 0.2;'
                  ' black: star L (< 0.025 Rvir)')
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
            #outname = outdir + outname
            title = title_temp.format(ic=ic, phys=phys, z=zv)
            
            pvs, Ldir, rvir_pkpc = setup_vizcheck_outflow_Lstarrecen(
                sim, snap, galrad_rvir=0.025)
            plot_vizcheck_outflow_Lstarrecen(pvs, Ldir, rvir_pkpc,
                                     title=title, axtitles=axtitles, 
                                     outname=outname,
                                     alpha=0.2, vscales=0.2)
                        

            


     


            


    

    
    



    

    


    
    







