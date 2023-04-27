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
vcen_gal = hp.get_vcom(simpath, snapnum, 1., meandef_rvir='BN98',
                       parttypes='all')
vcen_gal_cmps = (vcen_gal[0]['VXcom_cmps'], vcen_gal[0]['VYcom_cmps'],
                 vcen_gal[0]['VZcom_cmps'])
maptype_args_all = {'vcen_cmps': vcen_all_cmps, 'cen_cm': pcen_cm}
rvir_cm = vcen_all[0]['Rvir_cm']

# baseline data for value/calculation tests
posdirect_simu, posdirect_tocgs = \
    snapobj.readarray_emulateEAGLE('PartType0/Coordinates')
veldirect_simu, veldirect_tocgs = \
    snapobj.readarray_emulateEAGLE('PartType0/Velocities')

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
                             {'pos': 0}, {'pos', 1}, {'pos': 2},
                             {'vel', 0}, {'vel', 1}, {'vel': 2}]}
    specargs.update(maptype_args_all)
    posvel_gq_simu, posvel_gq_tocgs, _ = gq.get_qty(snapobj, 0, 'coords',
                                                    specargs,
                                                    filterdct=None)
    if not np.isclose(posvel_gq_tocgs[0], posdirect_tocgs):
        print('Position unit conversion mismatch between get_qty and direct')
        passed = False
    if not np.isclose(posvel_gq_tocgs[1], veldirect_tocgs):
        print('Position unit conversion mismatch between get_qty and direct')
        passed = False
    if not np.allclose(posdirect_simu - pcen_cm / posdirect_tocgs, 
                       posvel_gq_simu[0]):
        print('Mismatch between get_qty and direct positions')
        passed = False
    if not np.allclose(veldirect_simu - vcen_all_cmps / veldirect_tocgs, 
                       posvel_gq_simu[1]):
        print('Mismatch between get_qty and direct positions')
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
                         axis=0)
    rcen_direct = np.sqrt(rcen_direct)
    if not np.isclose(gq_simu[0], rcen_direct):
        print('Rcen mismatch between get_qty and direct')
        passed = False
    rdir = (posdirect_simu - pcen_cm / posdirect_tocgs) \
           / rcen_direct[:, np.newaxis]
    del rcen_direct
    vrad_direct = veldirect_simu - vcen_all / veldirect_tocgs
    vrad_direct = np.sum(rdir * vrad_direct, axis=1)
    if not np.isclose(gq_simu[1], vrad_direct):
        print('Vrad mismatch between get_qty and direct')
        passed = False
    del rdir, vrad_direct
    vcen_direct = np.sum((veldirect_simu 
                          - vcen_all_cmps / veldirect_tocgs)**2,
                         axis=0)
    vcen_direct = np.sqrt(vcen_direct)
    if not np.isclose(gq_simu[1], vcen_direct):
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
    if not np.allclose(pos1_gq_simu, posdirect_simu[filter, :]):
        print('Filterdct selection direct vs. get_qty mismatch')
        passed = False
    return passed

if __name__ == '__main__':
    allpassed = True
    allpassed &= test_cgsconv()
    allpassed &= test_cartesian_values()
    allpassed &= test_calc_values()
    allpassed &= test_filterdct()

    if allpassed:
        print('All tests passed')
    else:
        print('Some tests failed')
