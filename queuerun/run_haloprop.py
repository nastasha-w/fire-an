
import numpy as np

import mainfunc.haloprop as hp

def run_halodata(opt):
    # test cases
    if opt >= 0 and opt < 6:
        ind = opt - 0
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = [('m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ]
        snaps = [45, 50]
        meandef = ('BN98', '200c', '200m', '500c')
    # clean samples 1 and 2
    elif opt >= 6 and opt < 30:
        ind = opt - 6
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                   ]
        snaps = [45, 46, 47, 48, 49, 50]
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
    elif opt >= 30 and opt < 42:
        ind = opt - 30
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31_fa0.5',
                   ]
        snaps = [186, 197, 210, 224, 240, 258]
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
    elif opt >= 42 and opt < 54:
        ind = opt - 42
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
                    ]
        snaps = [186, 197, 210, 224, 240, 258]
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
    elif opt >= 54 and opt < 60:
        ind = opt - 54
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    ]
        snaps = [45, 46, 47, 48, 49, 50]
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
    simi = ind // len(snaps)
    snapi = ind % len(snaps)
    simname = simnames[simi]
    snapshot = snaps[snapi]

    dp2 = '_'.join(simname.split('_')[:2])
    if dp2.startswith('m13h02_'):
        dp2 = dp2.replace('m13h02', 'm13h002')
    dirpath = '/'.join([_dirpath, dp2, simname]) 
    
    hp.gethalodata_shrinkingsphere(dirpath, snapshot, meandef=meandef)