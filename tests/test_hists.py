
import numpy as np

import mainfunc.makehist as mh

def tryout_hist(index):
    outdir = '/projects/b1026/nastasha/tests/start_fire/hist_tests/'
    if index == 0:
        dirpath = '/projects/b1026/snapshots/fire3/m13h206_m3e5/' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000/'
        snapnum = 27
        simname = 'm13h206_m3e5__' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'
        outfilen = 'hist_Oxygen_by_Mass_0-1-2Rvir_{sc}_snap{sn}_shrink-sph-cen_BN98' + \
                   '_2rvir_v1.hdf5'

        axtypes = ['sim-direct']
        axtypes_args = [{'field': 'ElementAbundance/Oxygen'}]
        weighttype = 'Mass'
        weighttype_args = {}
        rbins = np.array([0., 1., 2.])
        runit = 'Rvir'

        outfilen = outfilen.format(sc=simname, sn=snapnum,)
    else:
        raise ValueError('invalid index: {}'.format(index))

    mh.histogram_radprof(dirpath, snapnum,
                         weighttype, weighttype_args, axtypes, axtypes_args,
                         particle_type=0, 
                         center='shrinksph', rbins=rbins, runit=runit,
                         logweights=True, logaxes=True, axbins=0.05,
                         outfilen=outdir + outfilen)