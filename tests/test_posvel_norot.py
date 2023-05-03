'''
test get_qty / CoordinateWrangler position and velocity calculations
testing direct values (from read_fire) against get_qty values
and total/radial velocity

Uses a specific simulation and snapshot; change simname, simpath,
and snapshot to test on a different one.
(This one was chosen to be a relatively small FIRE-3 dataset within
the sample I'm analysing.)
'''

import numpy as np

import fire_an.mainfunc.get_qty as gq
import fire_an.mainfunc.haloprop as hp
import fire_an.readfire.readin_fire_data as rfd
import fire_an.utils.constants_and_units as c

simname = ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
           '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000')
simpath = 'm12f_m6e4/' + simname
snapnum = 45

snapobj = rfd.get_Firesnap(simpath, snapnum)

vcen_all = hp.get_vcom(simpath, snapnum, 1., meandef_rvir='BN98',
                       parttypes='all')
vcen_all_cmps = (vcen_all[0]['VXcom_cmps'], vcen_all[0]['VYcom_cmps'],
                 vcen_all[0]['VZcom_cmps'])
pcen_cm = (vcen_all[0]['Xc_cm'], vcen_all[0]['Yc_cm'],
           vcen_all[0]['Zc_cm'])
#vcen_gal = hp.get_vcom(simpath, snapnum, 1., meandef_rvir='BN98',
#                       parttypes='all')
#vcen_gal_cmps = (vcen_gal[0]['VXcom_cmps'], vcen_gal[0]['VYcom_cmps'],
#                 vcen_gal[0]['VZcom_cmps'])
maptype_args_all = {'vcen_cmps': vcen_all_cmps, 'center_cm': pcen_cm}
rvir_cm = vcen_all[0]['Rvir_cm']
pcen_cm = np.array(pcen_cm)
vcen_all_cmps = np.array(vcen_all_cmps)

# baseline data for value/calculation tests
posdirect_simu = \
    snapobj.readarray_emulateEAGLE('PartType0/Coordinates')
posdirect_tocgs = snapobj.toCGS
veldirect_simu = \
    snapobj.readarray_emulateEAGLE('PartType0/Velocities')
veldirect_tocgs = snapobj.toCGS


def test_cgsconv():
    passed = True
    ckpch = 1e-3 * c.cm_per_mpc * snapobj.cosmopars.a / snapobj.cosmopars.h
    if not np.isclose(posdirect_tocgs, ckpch):
        print(f'Unexpected position units: {posdirect_tocgs}')
        passed = False
    if not np.isclose(veldirect_tocgs, 1e5 * np.sqrt(snapobj.cosmopars.a)):
        print(f'Unexpected velocity units: {veldirect_tocgs}')
        passed = False
    return passed
    
def test_cartesian_values():
    print('Testing whether direct and get_qty positions and velocities match')
    passed = True
    specargs = {'multiple': [{'pos': 'allcart'}, {'vel': 'allcart'},
                             {'pos': 0}, {'pos': 1}, {'pos': 2},
                             {'vel': 0}, {'vel': 1}, {'vel': 2}]}
    specargs.update(maptype_args_all)
    posvel_gq_simu, posvel_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                                    specargs,
                                                    filterdct=None)
    if not np.isclose(posvel_gq_tocgs[0], posdirect_tocgs):
        print('Position unit conversion mismatch between get_qty and direct')
        passed = False
    if not np.isclose(posvel_gq_tocgs[1], veldirect_tocgs):
        print('Velocity unit conversion mismatch between get_qty and direct')
        passed = False
    if not np.allclose(posdirect_simu - pcen_cm / posdirect_tocgs, 
                       posvel_gq_simu[0]):
        print('Mismatch between get_qty and direct positions')
        passed = False
    if not np.allclose(veldirect_simu - vcen_all_cmps / veldirect_tocgs, 
                       posvel_gq_simu[1]):
        print(veldirect_simu - vcen_all_cmps / veldirect_tocgs)
        print(posvel_gq_simu[1])
        printfilter = np.where(np.logical_not(np.isclose(veldirect_simu - vcen_all_cmps / veldirect_tocgs, 
                       posvel_gq_simu[1], rtol=1e-4)))
        print()
        print((veldirect_simu- vcen_all_cmps)[printfilter] / veldirect_tocgs)
        print(posvel_gq_simu[1][printfilter])
        print('Mismatch between get_qty and direct velocities')
        passed = False
    for ci in range(3):
        if not np.isclose(posvel_gq_tocgs[2 + ci], posdirect_tocgs):
            print(f'Position unit conversion (index {ci}) mismatch between'
                  'get_qty and direct')
            passed = False
        if not np.isclose(posvel_gq_tocgs[5 + ci], veldirect_tocgs):
            print(f'Velocity unit conversion (index {ci}) mismatch between'
                  'get_qty and direct')
            passed = False
        if not np.all(posvel_gq_simu[0][:, ci] == posvel_gq_simu[2 + ci]):
            print(f'Mismatch between get_qty allcart and index {ci} '
                  'positions')
            passed = False
        if not np.all(posvel_gq_simu[1][:, ci] == posvel_gq_simu[5 + ci]):
            print(f'Mismatch between get_qty allcart and index {ci} '
                  'velocities')
            passed = False
    print()
    return passed

def test_multiple_vs_single():
    passed = True
    print('Testing multiples vs. vel, pos calls')
    specargs = {'multiple': [{'pos': 'allcart'}, {'pos': 0}, {'pos': 1},
                             {'pos': 2}]}
    specargs.update(maptype_args_all)
    pall_gq_simu, pall_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                                specargs)
    if not np.allclose(pall_gq_simu[0], 
                       posdirect_simu - pcen_cm / posdirect_tocgs):
        print('direct vs. get_qty mismatch (pos, mulitple)')
        passed = False
    if not np.allclose(pall_gq_simu[1], 
                       (posdirect_simu - pcen_cm / posdirect_tocgs)[:, 0]):
        print('direct vs. get_qty mismatch (pos0, mulitple)')
        passed = False
    if not np.allclose(pall_gq_simu[2], 
                       (posdirect_simu - pcen_cm / posdirect_tocgs)[:, 1]):
        print('direct vs. get_qty mismatch (pos1, mulitple)')
        passed = False
    if not np.allclose(pall_gq_simu[3], 
                       (posdirect_simu - pcen_cm / posdirect_tocgs)[:, 2]):
        print('direct vs. get_qty mismatch (pos2, mulitple)')
        passed = False
    specargs = {'multiple': [{'vel': 'allcart'}, {'vel': 0}, {'vel': 1},
                             {'vel': 2}]}
    specargs.update(maptype_args_all)
    vall_gq_simu, vall_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                                specargs)
    if not np.allclose(vall_gq_simu[0], 
                       veldirect_simu - vcen_all_cmps / veldirect_tocgs):
        print('direct vs. get_qty mismatch (vel, mulitple)')
        passed = False
    if not np.allclose(vall_gq_simu[1], 
                       (veldirect_simu 
                        - vcen_all_cmps / veldirect_tocgs)[:, 0]):
        print('direct vs. get_qty mismatch (vel0, mulitple)')
        passed = False
    if not np.allclose(vall_gq_simu[2], 
                       (veldirect_simu 
                        - vcen_all_cmps / veldirect_tocgs)[:, 1]):
        print('direct vs. get_qty mismatch (vel1, mulitple)')
        passed = False
    if not np.allclose(vall_gq_simu[3], 
                       (veldirect_simu 
                        - vcen_all_cmps / veldirect_tocgs)[:, 2]):
        print('direct vs. get_qty mismatch (vel2, mulitple)')
        passed = False
    specargs = {'multiple': [{'pos': 1}, {'pos': 2}]}
    specargs.update(maptype_args_all)
    p12_gq_simu, p12_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                                specargs)
    if not np.allclose(np.array(p12_gq_simu).T, 
                       (posdirect_simu - pcen_cm / posdirect_tocgs)[:, 1:]):
        print('direct vs. get_qty mismatch (pos, 1, 2, multiple)')
        passed = False
    specargs = {'multiple': [{'vel': 1}, {'vel': 2}]}
    specargs.update(maptype_args_all)
    v12_gq_simu, v12_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                              specargs)
    if not np.allclose(np.array(v12_gq_simu).T, 
                       (veldirect_simu 
                        - vcen_all_cmps / veldirect_tocgs)[:, 1:]):
        print('direct vs. get_qty mismatch (vel, 1, 2, multiple)')
        passed = False
    specargs = {'pos': 0}
    specargs.update(maptype_args_all)
    p0_gq_simu, p0_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                            specargs)
    if not np.allclose(p0_gq_simu, 
                       (posdirect_simu - pcen_cm / posdirect_tocgs)[:, 0]):
        print('direct vs. get_qty mismatch (pos, 0, single)')
        passed = False
    specargs = {'vel': 0}
    specargs.update(maptype_args_all)
    v0_gq_simu, v0_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                            specargs)
    if not np.allclose(v0_gq_simu, 
                       (veldirect_simu 
                        - vcen_all_cmps / veldirect_tocgs)[:, 0]):
        print('direct vs. get_qty mismatch (vel, 0, single)')
        passed = False
    return passed

def test_calc_values():
    passed = True
    print('Testing whether direct and get_qty calculated quantities match')
    specargs = {'multiple': [{'pos': 'rcen'}, {'vel': 'vrad'},
                             {'vel': 'vtot'}]}
    specargs.update(maptype_args_all)
    gq_simu, gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                      specargs, filterdct=None)
    if not np.isclose(gq_tocgs[0], posdirect_tocgs):
        print('Rcen unit conversion mismatch between get_qty and direct')
        passed = False
    if not np.isclose(gq_tocgs[1], veldirect_tocgs):
        print('Vrad unit conversion mismatch between get_qty and direct')
        passed = False
    if not np.isclose(gq_tocgs[2], veldirect_tocgs):
        print('Vtot unit conversion mismatch between get_qty and direct')
        passed = False
    rcen_direct = np.sum((posdirect_simu - pcen_cm / posdirect_tocgs)**2, 
                         axis=1)
    rcen_direct = np.sqrt(rcen_direct)
    if not np.allclose(gq_simu[0], rcen_direct):
        print('Rcen mismatch between get_qty and direct')
        passed = False
    rdir = (posdirect_simu - pcen_cm / posdirect_tocgs) \
           / rcen_direct[:, np.newaxis]
    del rcen_direct
    vrad_direct = veldirect_simu - vcen_all_cmps / veldirect_tocgs
    vrad_direct = np.sum(rdir * vrad_direct, axis=1)
    if not np.allclose(gq_simu[1], vrad_direct):
        print('Vrad mismatch between get_qty and direct')
        passed = False
    del rdir, vrad_direct
    vtot_direct = np.sum((veldirect_simu 
                          - vcen_all_cmps / veldirect_tocgs)**2,
                         axis=1)
    vtot_direct = np.sqrt(vtot_direct)
    if not np.allclose(gq_simu[2], vtot_direct):
        print('Rcen mismatch between get_qty and direct')
        passed = False
    print()
    return passed

def test_filterdct():
    print('Testing filterdct')
    passed = True
    filter = np.array([1, 6, 131])
    filterdct = {'filter': filter}
    specargs = {'pos': 1}
    specargs.update(maptype_args_all)
    pos1_gq_simu, pos1_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                                specargs,filterdct=filterdct)
    if not np.allclose(pos1_gq_simu, 
                       posdirect_simu[filter, 1] 
                        - pcen_cm[1] / posdirect_tocgs):
        print(pos1_gq_simu)
        print(posdirect_simu[filter, 1] - pcen_cm[1] / posdirect_tocgs)
        print('Filterdct selection direct vs. get_qty mismatch (pos1)')
        passed = False
    specargs = {'multiple': [{'pos': 0}, {'pos': 2}]}
    specargs.update(maptype_args_all)
    pos02_gq_simu, pos02_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                                  specargs,
                                                  filterdct=filterdct)
    if not np.allclose(np.array(pos02_gq_simu).T, 
                       posdirect_simu[filter, :][:, np.array([0, 2])]
                        - pcen_cm[np.array([0, 2])] / posdirect_tocgs):
        print(pos02_gq_simu)
        print(posdirect_simu[filter, :][:, np.array([0, 2])])
        print('Filterdct selection direct vs. get_qty mismatch (pos0, 2)')
        passed = False
    specargs = {'vel': 2}
    specargs.update(maptype_args_all)
    vel2_gq_simu, vel2_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                                specargs,
                                                filterdct=filterdct)
    if not np.allclose(vel2_gq_simu, 
                       veldirect_simu[filter, 2]
                        - vcen_all_cmps[2] / veldirect_tocgs):
        print('Filterdct selection direct vs. get_qty mismatch (vel2)')
        passed = False
    specargs = {'multiple': [{'vel': 1}, {'vel': 2}]}
    specargs.update(maptype_args_all)
    vel12_gq_simu, vel12_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                                  specargs,
                                                  filterdct=filterdct)
    if not np.allclose(np.array(vel12_gq_simu).T, 
                       veldirect_simu[filter, :][:, np.array([1, 2])]
                        - vcen_all_cmps[np.array([1, 2])] 
                          / veldirect_tocgs):
        print(vel12_gq_simu)
        print(veldirect_simu[filter, :][:, np.array([1, 2])])
        print('Filterdct selection direct vs. get_qty mismatch (vel1, 2)')
        passed = False
    specargs = {'pos': 'allcart'}
    specargs.update(maptype_args_all)
    pos_gq_simu, pos_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                              specargs,
                                              filterdct=filterdct)
    if not np.allclose(pos_gq_simu, 
                       posdirect_simu[filter] - pcen_cm / posdirect_tocgs):
        print(pos_gq_simu)
        print(posdirect_simu[filter, :])
        print('Filterdct selection direct vs. get_qty mismatch (pos)')
        passed = False
    specargs = {'vel': 'allcart'}
    specargs.update(maptype_args_all)
    vel_gq_simu, vel_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                              specargs, filterdct=filterdct)
    if not np.allclose(vel_gq_simu, 
                       veldirect_simu[filter, :] 
                        - vcen_all_cmps / veldirect_tocgs):
        print(vel_gq_simu)
        print(veldirect_simu[filter, :])
        print('Filterdct selection direct vs. get_qty mismatch (vel)')
        passed = False
    specargs = {'vel': 'vrad'}
    specargs.update(maptype_args_all)
    vrd_gq_simu, vrd_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                              specargs, filterdct=filterdct)
    rcen_direct = np.sum((posdirect_simu[filter, :] 
                          - pcen_cm / posdirect_tocgs)**2, 
                         axis=1)
    rcen_direct = np.sqrt(rcen_direct)
    rdir = (posdirect_simu[filter, :] - pcen_cm / posdirect_tocgs) \
           / rcen_direct[:, np.newaxis]
    del rcen_direct
    vrad_direct = veldirect_simu[filter, :] - vcen_all_cmps / veldirect_tocgs
    vrad_direct = np.sum(rdir * vrad_direct, axis=1)
    del rdir
    if not np.allclose(vrd_gq_simu, vrad_direct):
        print('Filterdct selection direct vs. get_qty mismatch (vrad)')
        passed = False
    return passed

def test_all():
    allpassed = True
    allpassed &= test_cgsconv()
    allpassed &= test_cartesian_values()
    allpassed &= test_calc_values()
    allpassed &= test_multiple_vs_single()
    allpassed &= test_filterdct()

    if allpassed:
        print('All tests passed')
    else:
        print('Some tests failed')
    return allpassed

if __name__ == '__main__':
    test_all()
