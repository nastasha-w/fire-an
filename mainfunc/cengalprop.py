'''
starting from haloprop halo center and rvir, 
get central galaxy center, velocity, stellar mass
'''

import glob
import h5py
import numpy as np
import os
import sys
import uuid

import fire_an.mainfunc.haloprop as hp
import fire_an.readfire.readin_fire_data as rf
import fire_an.simlists as sl
import fire_an.utils.h5utils as h5u
import fire_an.utils.opts_locs as ol

def calccengalcen(simpath, snapnum, startrad_rvir=0.3,
                  vcenrad_rvir=0.05, mstarrad_rvir=0.1):
    '''
    starting from the all-types halo center of mass, find the 
    central galaxy center of mass and center of velocity
    (halo centering gets close, but ~1 kpc off is too much for 
    angular momentum calculations)

    This runs pretty fast.
    '''
    todoc = {}

    halodat = hp.gethalodata_shrinkingsphere(simpath, snapnum, meandef='BN98')
    posstart_cm = np.array([halodat[0]['Xc_cm'], halodat[0]['Yc_cm'], 
                            halodat[0]['Zc_cm']])
    rvir_cm = halodat[0]['Rvir_cm']
    todoc['halodata'] = halodat[0]
    todoc['halodata_doc'] = halodat[1]
    snapobj = rf.get_Firesnap(simpath, snapnum)
    spos_simu = snapobj.readarray('PartType4/Coordinates')
    spos_toCGS = snapobj.toCGS
    stard2 = np.sum((spos_simu - posstart_cm / spos_toCGS)**2, axis=1)  
    starsel = stard2 <= (startrad_rvir * rvir_cm / spos_toCGS)
    del stard2
    spos_simu = spos_simu[starsel]
    smass_simu = snapobj.readarray('PartType4/Masses')[starsel]
    smass_toCGS = snapobj.toCGS
    coordsmassesdict = {'coords': spos_simu,
                        'masses': smass_simu}

    kwargs_calccen = {'shrinkfrac': 0.025, 
                      'minparticles': 1000, 
                      'initialradiusfactor': 1.}
    todoc['kwargs_calchalocen_stars'] = kwargs_calccen
    todoc['startrad_rvir'] = startrad_rvir
    todoc['vcenrad_rvir'] = vcenrad_rvir

    scen_simu, _, _ = hp.calchalocen(coordsmassesdict, **kwargs_calccen)
    stard2 = np.sum((spos_simu - scen_simu)**2, axis=1)
    starsel2 = stard2 <= (vcenrad_rvir * rvir_cm / spos_toCGS)
    starsel[starsel] = starsel2
    smass_simu = smass_simu[starsel2]
    svel_simu = snapobj.readarray('PartType4/Velocities')[starsel]
    svel_toCGS = snapobj.toCGS
    vcom_simu = np.sum(svel_simu * smass_simu[:, np.newaxis], axis=0) \
                / np.sum(smass_simu)
    vcom_cmps = vcom_simu * svel_toCGS
    pcen_cm = scen_simu * spos_toCGS
    todoc['starcen_cm'] = pcen_cm
    todoc['starvcom_cmps'] = vcom_cmps
    
    spos_simu = snapobj.readarray('PartType4/Coordinates')
    stard2 = np.sum((spos_simu - scen_simu)**2, axis=1)
    starsel = stard2 <= (mstarrad_rvir * rvir_cm / spos_toCGS)
    mass_simu = snapobj.readarray('PartType4/Masses')[starsel]
    mass_toCGS = snapobj.toCGS
    starmass_g = np.sum(mass_simu) * mass_toCGS

    todoc['mstarrad_rvir'] = mstarrad_rvir
    todoc['mstar_gal_g'] = starmass_g
    return pcen_cm, vcom_cmps, todoc

def savedata_cengalcen(simpath, snapnum, pcen_cm, vcom_cmps, todoc):
    filen = ol.dir_halodata + f'temp_pvcengal_{uuid.uuid1()}.hdf5'
    if os.path.isfile(filen):
        raise RuntimeError(f'temp. file {filen} already exists')
    #if simpath.endswith('output'):
    #    simname = simpath.split('/')[-2]
    #elif simpath.endswith('output/'):
    #    simname = simpath.split('/')[-3]
    #elif simpath.endswith('/'):
    #    simname = simpath.split('/')[-2]
    #else:
    #    simname = simpath.split('/')[-1]
    simname = sl.simname_from_dirpath(simpath)

    with h5py.File(filen, 'w') as f:
        g1 = f.create_group(simname)
        g2 = g1.create_group(f'snap_{snapnum}')
        g3 = g2.create_group('pv0')
        gd = g3.create_group('doc')
        h5u.savedict_hdf5(gd, todoc)
        g3.attrs.create('pcen_cm', pcen_cm)
        g3.attrs.create('vcom_cmps', vcom_cmps)

def readdata_cengalcen(simpath, snapnum, startrad_rvir=0.3,
                       vcenrad_rvir=0.05, mstarrad_rvir=0.1):
    filen = ol.dir_halodata + 'pvcengal.hdf5'
    #simname = simpath.split('/')[-1]
    simname = sl.simname_from_dirpath(simpath)
    with h5py.File(filen, 'r') as f:
        if simname not in f:
            raise hp.NoStoredMatchError(f'simname {simname} data not stored')
        g1 = f[simname]
        if f'snap_{snapnum}' not in g1:
            raise hp.NoStoredMatchError(f'snapshot {snapnum} data not'
                                     f' stored for {simname}')
        g2 = g1[f'snap_{snapnum}']
        g3nopts = g2.keys()
        tomatch = {'startrad_rvir': startrad_rvir,
                   'vcenrad_rvir': vcenrad_rvir,
                   'mstarrad_rvir': mstarrad_rvir}
        mkeys = list(tomatch.keys())
        g3n = None
        for g3nopt in g3nopts:
            if np.all([np.isclose(tomatch[key], g2[g3nopt]['doc'].attrs[key]) 
                       for key in mkeys]):
                g3n = g3nopt
                break
        if g3n is None:
            raise hp.NoStoredMatchError(f'{simname}, snapshot {snapnum}'
                                        f' data with parameters {tomatch}'
                                        ' not found')
        g3 = g2[g3n]
        todoc = h5u.readgrp_todict(g3['doc'], subgroups=True)
        pcen_cm = g3.attrs['pcen_cm']
        vcom_cmps = g3.attrs['vcom_cmps']
    return pcen_cm, vcom_cmps, todoc
        
def getcengalcen(simpath, snapnum, startrad_rvir=0.3,
                 vcenrad_rvir=0.05, mstarrad_rvir=0.1):
    try:
        out = readdata_cengalcen(simpath, snapnum, 
                                 startrad_rvir=startrad_rvir,
                                 vcenrad_rvir=vcenrad_rvir, 
                                 mstarrad_rvir=mstarrad_rvir)
    except hp.NoStoredMatchError:
        out = calccengalcen(simpath, snapnum, 
                            startrad_rvir=startrad_rvir,
                            vcenrad_rvir=vcenrad_rvir, 
                            mstarrad_rvir=mstarrad_rvir)
        savedata_cengalcen(simpath, snapnum, *out)
    return out

def adddata_cengalcen(rmtemp=False):
    mainfilen =  ol.dir_halodata + 'pvcengal.hdf5'
    searchcrit = ol.dir_halodata + 'temp_pvcengal_*.hdf5'
    tempfilens = glob.glob(searchcrit)
    if len(tempfilens) == 0:
        print('No new data to add')
        return None
    with h5py.File(mainfilen, 'a') as fo:
        for tfn in tempfilens:
            with h5py.File(tfn, 'r') as fi:
                simid = next(iter(fi.keys()))
                #print('simid: ', simid)
                #print('main file keys: ', list(fo.keys()))
                if simid not in fo: #easy, copy whole thing
                    fi.copy(fi[simid], fo, name=simid)
                    print(f'Added file {tfn}:')
                    print(f'{simid}')
                    continue
                fo_smgrp = fo[simid]
                fi_smgrp = fi[simid]
                sngrpn = next(iter(fi_smgrp.keys()))
                #print('snap: ', sngrpn)
                #print('main file fo_smgrp keys: ', list(fo_smgrp.keys()))
                if sngrpn not in fo_smgrp: #easy, copy whole thing
                    fi.copy(fi_smgrp[sngrpn], fo_smgrp, name=sngrpn)
                    print(f'Added file {tfn}:')
                    print(f'{simid}, {sngrpn}')
                    continue
                # center matching/copy
                fo_sngrp = fo_smgrp[sngrpn]
                fi_sngrp = fi_smgrp[sngrpn]
                #print('main file cens: ', list(fo_sngrp.keys()))
                if 'pv0' not in fo_sngrp:
                    fi.copy(fi_sngrp['pv0'], fo_sngrp, name='pv0')
                    print(f'Added file {tfn}:')
                    print(f'{simid}, {sngrpn}, first center')
                    continue
                cens_fo = [grp for grp in fo_sngrp.keys() \
                           if grp.startswith('pv')]
                fi_cgrp = fi_sngrp['pv0']
                tocheck = ['vcom_cmps', 'pcen_cm']
                anymatch = False
                for cengrpn in cens_fo:
                    _fo_cgrp = fo_sngrp[cengrpn]
                    ismatch = h5u.checkattrsmatch(_fo_cgrp['doc'],
                                                  fi_cgrp['doc'],
                                                  verbose=True)
                    if ismatch:
                        fo_cgrp = _fo_cgrp
                        anymatch = True
                        print(fo_cgrp)
                        print(fi_cgrp)
                        if not np.all([np.allclose(fo_cgrp.attrs[key],
                                                   fi_cgrp.attrs[key])
                                       for key in tocheck]):
                            msg = (f'{mainfilen} and {tfn} have matching'
                                   f'simulation {simid}, {sngrpn}, '
                                   f'center finding, but different centers:\n'
                                   f'{fo_cgrp.attrs.items()},\n'
                                   f'{fi_cgrp.attrs.items()}')
                            raise RuntimeError(msg)
                        # no continue if everything just matches -> goes to
                        # delete file check
                if not anymatch:
                    fo_cgrpn = f'pv{len(cens_fo)}'
                    fi.copy(fi_cgrp, fo_sngrp, name=fo_cgrpn)
                    print(f'Added file {tfn}:')
                    print(f'{simid}, {sngrpn}, new center')
                    continue
            print(f'skipped {tfn}; duplicate data')
            if rmtemp:
                print(f'deleting {tfn}')
                os.remove(tfn)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        mode = sys.argv[1]
        if mode == '--addstored':
            if len(sys.argv) > 2:
                rmtemp = bool(sys.argv[2])
            else:
                rmtemp = False
            adddata_cengalcen(rmtemp=rmtemp)
        else:
            raise ValueError(f'Invalid mode {mode}')
    else:
        print('No actions specified')