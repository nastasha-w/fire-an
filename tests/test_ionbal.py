
import h5py
import numpy as np

from ionrad.ion_utils import Linetable_PS20
import mainfunc.get_qty as gq
import readfire.readin_fire_data as rf
import utils.constants_and_units as c



def test_ionbal_calc(dirpath, snapnum, ion, target_Z=0.01, delta_Z=0.001,
                     ps20depletion=False, outfilen='ionbal_test.hdf5',
                     lintable=False):
    snap = rf.get_Firesnap(dirpath, snapnum)
    cosmopars = snap.cosmopars.getdct()
    
    # filter sim. particles and calculate ion balances, rho, T, Z
    metallicity = snap.readarray_emulateEAGLE('PartType0/Metallicity')
    zfilter = metallicity >= target_Z - delta_Z
    zfilter &= metallicity <= target_Z + delta_Z
    metallicity = metallicity[zfilter]
    indct = {'filter': zfilter}
    ionbals = gq.get_ionfrac(snap, ion, indct=indct, table='PS20', 
                             simtype='fire',
                             ps20depletion=ps20depletion, lintable=lintable)
    temperature = snap.readarray_emulateEAGLE('PartType0/Temperature')[zfilter]
    temperature *= snap.toCGS
    hdens = snap.readarray_emulateEAGLE('PartType0/Density')[zfilter]
    hconv = snap.toCGS
    hdens *= snap.readarray_emulateEAGLE('PartType0/ElementAbundance/Hydrogen')[zfilter]
    hconv *= snap.toCGS
    hconv /= (c.atomw_H * c.u)
    hdens *= hconv
    
    # get corresponding ion balance table
    # for table read-in only, lin/log shouldn't matter; problems weren't there
    iontab = Linetable_PS20(ion, cosmopars['z'], emission=False, vol=True,
                            lintable=False)
    iontab.findiontable()
    tab_logT = iontab.logTK
    tab_lognH = iontab.lognHcm3
    tab_logZ = iontab.logZsol + np.log10(iontab.solarZ)
    tab_ionbal_T_Z_nH = iontab.iontable_T_Z_nH.copy()
    if ps20depletion:
        tab_ionbal_T_Z_nH = 10**tab_ionbal_T_Z_nH
        iontab.finddepletiontable()
        tab_depletion_T_Z_nH = iontab.depletiontable_T_Z_nH.copy()
        tab_ionbal_T_Z_nH *= (1. - 10**tab_depletion_T_Z_nH)
        tab_ionbal_T_Z_nH = np.log10(tab_ionbal_T_Z_nH)
        tab_ionbal_T_Z_nH[tab_ionbal_T_Z_nH < -50.] = -50.

    interpvalZ = np.log10(target_Z)
    iZhi = np.where(tab_logZ >= interpvalZ)[0][0]
    iZlo = np.where(tab_logZ <= interpvalZ)[0][-1]
    if iZlo == iZhi:
        tab_ionbal_T_nH = tab_ionbal_T_Z_nH[:, iZlo, :] 
    else:
        hiZ = tab_logZ[iZhi]
        loZ = tab_logZ[iZlo]
        tab_ionbal_T_nH = (hiZ - interpvalZ) / (hiZ - loZ) * tab_ionbal_T_Z_nH[:, iZlo, :] +\
                          (interpvalZ - loZ) / (hiZ - loZ) * tab_ionbal_T_Z_nH[:, iZhi, :]
    tab_ionbal_T_nH = 10**tab_ionbal_T_nH

    # save data
    with h5py.File(outfilen, 'w') as f:
        hed = f.create_group('Header')
        cgrp = hed.create_group('cosmopars')
        cosmopars = snap.cosmopars.getdct()
        for key in cosmopars:
            cgrp.attrs.create(key, cosmopars[key])
        hed.attrs.create('snapnum', snapnum)
        hed.attrs.create('filepath_first', np.string_(snap.firstfilen))
        _info = 'FIRE calculated ion balances and the underlying ion balance table'
        hed.attrs.create('info', np.string_(_info))
        hed.attrs.create('target_Z', target_Z)
        hed.attrs.create('delta_Z', delta_Z)
        hed.attrs.create('ion', np.string_(ion))
        hed.attrs.create('ps20depletion', ps20depletion)
        hed.attrs.create('lintable', lintable)
        
        gsim = f.create_group('simulation_data')
        gsim.create_dataset('ionbal', data=ionbals)
        gsim.create_dataset('T_K', data=temperature)
        gsim.create_dataset('nH_cm**-3', data=hdens)
        gsim.create_dataset('metallicity_abs_mass_frac', data=metallicity)
        
        gtab = f.create_group('iontab_data')
        print('About to save tab_ionbal_T_nH')
        print('{} / {} NaN'.format(np.sum(np.isnan(tab_ionbal_T_nH)), 
                                   np.prod(tab_ionbal_T_nH.shape)))
        print(tab_ionbal_T_nH)
        gtab.create_dataset('ionbal_T_nH', data=tab_ionbal_T_nH)
        gtab.create_dataset('logT_K', data=tab_logT)
        gtab.create_dataset('lognH_cm**-3', data=tab_lognH)

def run_ionbal_test(opt=0):
    dirpath1 = '/projects/b1026/snapshots/fire3/m13h206_m3e5/' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000/'
    simname1 = 'm13h206_m3e5__' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'
    snaps1 = [27, 45]
    ions1 = ['O{}'.format(i) for i in range(1, 10)]

    outdir =  '/projects/b1026/nastasha/tests/start_fire/ionbal_tests/'
    outtemplate1 = outdir + 'ionbal_test_PS20_{ion}_depletion-{dp}_Z-{Z}' + \
                            '_snap{snap:03d}_{sim}.hdf5'
    outtemplate2 = outdir + 'ionbal_test_PS20_{ion}_depletion-{dp}_Z-{Z}' + \
                            '_snap{snap:03d}_lintable-{lintable}_{sim}.hdf5'
    
    if opt >= 0 and opt < 6:
        dirpath = dirpath1
        simname = simname1
        ions = ions1
        ps20depletion = bool(opt % 2)
        snapnum = snaps1[opt // 4]
        target_Z = [0.01, 0.0001][(opt  // 2) % 2]
        delta_Z = 0.1 * target_Z
        outtemplate = outtemplate1
        dolintable = False
    if opt >= 6 and opt < 18:
        _opt = opt - 6
        __opt = _opt % 6
        lintable = bool(_opt // 6)
        dirpath = dirpath1
        simname = simname1
        ions = ions1
        ps20depletion = bool(__opt % 2)
        snapnum = snaps1[__opt // 4]
        target_Z = [0.01, 0.0001][(__opt  // 2) % 2]
        delta_Z = 0.1 * target_Z
        outtemplate = outtemplate2
        dolintable = True
    else:
        raise ValueError('Invalid opt {}'.format(opt))
    for ion in ions:
        if dolintable:
            outfilen = outtemplate.format(ion=ion, dp=ps20depletion, 
                                          Z=target_Z, sim=simname, 
                                          snap=snapnum, lintable=lintable)
            test_ionbal_calc(dirpath, snapnum, ion, target_Z=target_Z, 
                             delta_Z=delta_Z, ps20depletion=ps20depletion, 
                             outfilen=outfilen, lintable=lintable)
        else:
            outfilen = outtemplate.format(ion=ion, dp=ps20depletion, 
                                          Z=target_Z, sim=simname, 
                                          snap=snapnum)
            test_ionbal_calc(dirpath, snapnum, ion, target_Z=target_Z, 
                             delta_Z=delta_Z, ps20depletion=ps20depletion, 
                             outfilen=outfilen)
