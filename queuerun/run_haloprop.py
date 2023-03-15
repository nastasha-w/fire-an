
import numpy as np

import mainfunc.haloprop as hp
import simlists as sl

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
    # redo runs after truncation/saturation cumsum bug
    elif opt >= 60 and opt < 126:
        ind = opt - 60
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = sl.m13_sr_all1 # len 11
        snaps = sl.snaplists['m13_sr'] # len 6
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
    elif opt >= 126 and opt < 174:
        ind = opt - 126
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
        simnames = sl.m13_hr_all1 # len 8
        snaps = sl.snaplists['m13_hr'] # len 6
    elif opt >= 174 and opt < 198:
        ind = opt - 174
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
        simnames = sl.m12_sr_all1 # len 4
        snaps = sl.snaplists['m12_sr'] # len 6
    elif opt >= 198 and opt < 294:
        ind = opt - 198
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
        simnames = sl.m12_hr_all1 # len 16
        snaps = sl.snaplists['m12_hr'] # len 6

    simi = ind // len(snaps)
    snapi = ind % len(snaps)
    simname = simnames[simi]
    snapshot = snaps[snapi]

    dp2 = '_'.join(simname.split('_')[:2])
    if dp2.startswith('m13h02_'):
        dp2 = dp2.replace('m13h02', 'm13h002')
    dirpath = '/'.join([_dirpath, dp2, simname]) 
    
    hp.gethalodata_shrinkingsphere(dirpath, snapshot, meandef=meandef)
        