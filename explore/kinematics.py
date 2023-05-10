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

def getweightfilter(simname, snapnum, selqty, selqty_args, 
                    filterdct=None, samplesize=1000,
                    parttype=0):
    '''
    select random resolution elements with probability equal
    to fraction of selqty in each particle 
    (so please don't use for temperatures or something)
    returns indices for selected elements
    '''
    simpath = getpath(simname)
    snap = rfd.get_Firesnap(simpath, snapnum)

    selvals, selvals_toCGS, selvals_todoc = gq.get_qty(
        snap, parttype, selqty, selqty_args, filterdct=filterdct)
    normedvals = selvals / np.sum(selvals)
    out = np.random.choice(len(normedvals), size=(samplesize), 
                           replace=False, p=normedvals)
    return out

def getkininfo(simname, snapnum, filterdct=None, parttype=0, vr=False):
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
    if vr:
        coords_todo = coords_todo + [{'vel': 'vrad'}]
    qty_args = {'center_cm': cen_cm, 'vcen_cmps': vcen_cmps,
                'multiple': coords_todo}
    valspv, toCGSpv, todoc_pv = gq.get_qty(
        snap, parttype, qty, qty_args, filterdct=filterdct)
    pos_pkpc = valspv[0] * (toCGSpv[0] / (c.cm_per_mpc * 1e-3))
    vel_kmps = valspv[1] * (toCGSpv[1] / (1e5))
    if vr:
        vr_kmps = valspv[2] * (toCGSpv[2] / (1e5))
        return pos_pkpc, vel_kmps, vr_kmps, rvir_pkpc
    else:
        return pos_pkpc, vel_kmps, rvir_pkpc

def get_selkin(simname, snapnum, maxradius_rvir, 
               selqtys, selqtys_args, samplesize=2000,
               parttype=0):
    sf = getspacefilter(simname, snapnum, maxradius_rvir, parttype=parttype)
    sfd = {'filter': sf}

    pall_pkpc, vall_kmps, vrall_kmps, rvir_pkpc = getkininfo(
        simname, snapnum, filterdct=sfd, parttype=parttype, vr=True)
    pvs = []
    for selqty, selqty_args in zip(selqtys, selqtys_args):
        _filter = getweightfilter(simname, snapnum, selqty, selqty_args, 
                        filterdct=sfd, samplesize=samplesize,
                        parttype=parttype)
        pvs.append((pall_pkpc[_filter, :], vall_kmps[_filter, :],
                    vrall_kmps[_filter]))
    return pvs, rvir_pkpc

def run_sel1(simname, snapnum):
    selqtys = ['Mass', 'Volume', 'ion']
    selqtys_args = [{}, {}, {'ion': 'Ne8'}]
    pv, rv = get_selkin(simname, snapnum, 1., 
                        selqtys, selqtys_args, samplesize=500,
                        parttype=0)
    pvs, rvs = get_selkin(simname, snapnum, 1., 
                          ['Mass'], [{}], samplesize=500,
                          parttype=4)
    pvsc, rvsc = get_selkin(simname, snapnum, 0.05, 
                            ['Mass'], [{}], samplesize=100,
                            parttype=4)
    labels = ['Mass', 'Volume', 'Ne8', 'stars', 'stars < 0.05 Rvir']
    return pv + pvs + pvsc, labels

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
    ncmax = 3
    npanels = len(posvels)
    ncols = min(ncmax, npanels)
    nrows = (npanels - 1) // ncols + 1
    panelsize = 4.
    fig = plt.figure(figsize=(panelsize * ncols, panelsize * nrows))
    axes = [fig.add_subplot(nrows, ncols, i + 1, projection='3d') 
            for i in range(len(posvels))]
    if not hasattr(vscales, '__len__'):
        vscales = [vscales] * npanels
    for i, (_p, _v, _vr) in enumerate(posvels):
        ax = axes[i]
        __v = vscales[i] * _v
        __p = np.copy(_p)
        vmax = np.max(np.abs(_vr))
        vmin = -1. * vmax 
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
        plt.savefig(outdir + outname)
    plt.show()


    
    



    

    


    
    







