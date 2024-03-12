
import fire_an.mainfunc.cengalprop as cgp
import fire_an.mainfunc.haloprop as hp
import fire_an.simlists as sl
import fire_an.utils.opts_locs as ol

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
    # and add new ICs
    elif opt >= 60 and opt < 150:
        # m13-sr, 90 indices
        ind = opt - 60
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = sl.m13_sr_all1 # len 15
        snaps = sl.snaplists['m13_sr'] # len 6
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
    elif opt >= 150 and opt < 162:
        # m13-hr, 12 indices
        ind = opt - 150
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
        simnames = sl.m13_hr_all1 # len 2
        snaps = sl.snaplists['m13_hr'] # len 6
    elif opt >= 162 and opt < 186:
        # m12-sr, 24 indices
        ind = opt - 162
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
        simnames = sl.m12_sr_all1 # len 4
        snaps = sl.snaplists['m12_sr'] # len 6
    elif opt >= 186 and opt < 282:
        # m12-hr, 96 indices
        ind = opt - 186
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
        simnames = sl.m12_hr_all1 # len 16
        snaps = sl.snaplists['m12_hr'] # len 6
    elif opt >= 282 and opt < 294:
        # m12-hr, 12 indices (all1 -> all2 supplement)
        ind = opt - 282
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
        simnames = [('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp2e-4_gacc31_fa0.5')
                   ]
        snaps = sl.snaplists['m12_hr'] # len 6
    elif opt >= 294 and opt < 304:
        # m13-sr z=0 
        # for C-series profiles z=0., 0.5, 1., supplement clean sample
        ind = opt - 294
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
        simnames = sl.m13_sr_all2_z0 # len 10
        snaps = [60]
    elif opt >= 304 and opt < 308:
        # m12-sr z=0 
        # for C-series profiles z=0., 0.5, 1., supplement clean sample
        ind = opt - 304
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
        simnames = sl.m12_sr_all2_z0 # len 4
        snaps = [60]
    # m13_hr_all2_z0 : length 0
    elif opt >= 308 and opt < 325:
        # m12-hr z=0 
        # for C-series profiles z=0., 0.5, 1., supplement clean sample
        ind = opt - 308
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
        simnames = sl.m12_hr_all2_z0 # len 17
        snaps = [500]
    elif opt >= 325 and opt < 352:
        # m11-hr
        ind = opt - 325
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
        simnames = sl.m11_hr_set1 # len 9
        snaps = sl.snaps_hr_051 # len 3
    elif opt >= 352 and opt < 421:
        # m11-sr
        ind = opt - 352
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # might as well; extra overdensities are cheap
        meandef = ('BN98', '200c', '200m', '500c', '500m', 
                   '2500c', '2500m', '178c', '178m', '100c', '100m')
        simnames = sl.m11_sr_set1 # len 23
        snaps = sl.snaps_sr_051 # len 3

    simi = ind // len(snaps)
    snapi = ind % len(snaps)
    simname = simnames[simi]
    snapshot = snaps[snapi]

    dp2 = '_'.join(simname.split('_')[:2])
    if dp2.startswith('m13h02_'):
        dp2 = dp2.replace('m13h02', 'm13h002')
    dirpath = '/'.join([_dirpath, dp2, simname]) 
    
    hp.gethalodata_shrinkingsphere(dirpath, snapshot, meandef=meandef)

def run_vcen1(opt):
    '''
    1 index runs the Halo (all particles but type 2, 1 Rvir) and
    the galaxy (parttype 4, 0.1 Rvir) velocity calculation.
    '''
    _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
    meandef = 'BN98'
    if opt >= 0 and opt < 90:
        ind = opt - 0
        simnames = sl.m13_sr_all2 # 15
        snaps = sl.snaps_sr # 6
    elif opt >= 90 and opt < 102:
        ind = opt - 90
        simnames = sl.m13_hr_all2 # 2
        snaps = sl.snaps_hr # 6
    elif opt >= 102 and opt < 126:
        ind = opt - 102
        simnames = sl.m12_sr_all2 # 4
        snaps = sl.snaps_sr # 6
    elif opt >= 126 and opt < 234:
        ind = opt - 126
        simnames = sl.m12_hr_all2 # 18
        snaps = sl.snaps_hr # 6
    else:
        raise ValueError(f'Nothing specified for index {ind}')
    simi = ind // len(snaps)
    snapi = ind % len(snaps)
    simname = simnames[simi]
    snapnum = snaps[snapi]

    dp2 = '_'.join(simname.split('_')[:2])
    if dp2.startswith('m13h02_'):
        dp2 = dp2.replace('m13h02', 'm13h002')
    dirpath = '/'.join([_dirpath, dp2, simname]) 

    print(f'Whole halo, {simname}, snap {snapnum}')
    hp.get_vcom(dirpath, snapnum, 1., meandef_rvir='BN98',
                       parttypes='all')
    print('\n\n')
    print(f'Galaxy, {simname}, snap {snapnum}')
    hp.get_vcom(dirpath, snapnum, 0.1, meandef_rvir=meandef,
                parttypes=(4,))
        

def runcengal1(opt=0):
    _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
    if opt >= 0 and opt < 90:
        ind = opt - 0
        simnames = sl.m13_sr_all2 # 15
        snaps = sl.snaps_sr # 6
    elif opt >= 90 and opt < 102:
        ind = opt - 90
        simnames = sl.m13_hr_all2 # 2
        snaps = sl.snaps_hr # 6
    elif opt >= 102 and opt < 126:
        ind = opt - 102
        simnames = sl.m12_sr_all2 # 4
        snaps = sl.snaps_sr # 6
    elif opt >= 126 and opt < 234:
        ind = opt - 126
        simnames = sl.m12_hr_all2 # 18
        snaps = sl.snaps_hr # 6
    else:
        raise ValueError(f'Nothing specified for index {ind}')
    simi = ind // len(snaps)
    snapi = ind % len(snaps)
    simname = simnames[simi]
    snapnum = snaps[snapi]

    dp2 = '_'.join(simname.split('_')[:2])
    if dp2.startswith('m13h02_'):
        dp2 = dp2.replace('m13h02', 'm13h002')
    dirpath = '/'.join([_dirpath, dp2, simname]) 

    cgp.getcengalcen(dirpath, snapnum, startrad_rvir=0.3,
                     vcenrad_rvir=0.05, mstarrad_rvir=0.1)
    
def run_haloprop_f2md(opt):
    '''
    1 index runs the Halo (all particles but type 2, 1 Rvir) and
    the galaxy (parttype 4, 0.1 Rvir) velocity calculation,
    and the central galaxy re-centering (parttype 4 in 0.3 Rvir)
    '''
    _dirpath = '/scratch/projects/xsede/GalaxiesOnFIRE/metal_diffusion/'
    meandef = 'BN98'

    # 48 indices (1st run)
    # + 12 indices, 2nd run with crheatfix added
    ind = opt - 0
    simnames = sl.m12_f2md # 8
    snaps = sl.snaps_f2md # 6

    simi = ind // len(snaps)
    snapi = ind % len(snaps)
    simname = simnames[simi]
    snapnum = snaps[snapi]

    #dirpath = '/'.join([_dirpath, simname]) 
    dirpath = sl.dirpath_from_simname(simname)

    #print(f'Whole halo center, {simname}, snap {snapnum}')
    #hp.gethalodata_shrinkingsphere(dirpath, snapnum, 
    #                               meandef=('BN98', '200c', '200m'))

    print(f'Whole halo Vcom, {simname}, snap {snapnum}')
    hp.get_vcom(dirpath, snapnum, 1., meandef_rvir=meandef,
                parttypes='all')
    print('\n\n')
    print(f'Galaxy Vcom, {simname}, snap {snapnum}')
    hp.get_vcom(dirpath, snapnum, 0.1, meandef_rvir=meandef,
                parttypes=(4,)) 
    print('\n\n')
    print(f'Galaxy re-centering, {simname}, snap {snapnum}')
    cgp.getcengalcen(dirpath, snapnum, startrad_rvir=0.3,
                     vcenrad_rvir=0.05, mstarrad_rvir=0.1)

def run_haloprop_fire3_m12new(opt):
    '''
    1 index runs the Halo (all particles but type 2, 1 Rvir) and
    the galaxy (parttype 4, 0.1 Rvir) velocity calculation,
    and the central galaxy re-centering (parttype 4 in 0.3 Rvir)
    '''
    # 1st run does not include m12g_r7100: 
    # Pratik is re-running after files were corrupted
    _dirpath = ol.simdir_fire3_m12plus
    meandef = 'BN98'

    if opt >= 0 and opt < 18:
        ind = opt - 0
        simnames = sl.m12plus_f3nobh # 3
        snaps = sl.snaps_hr # 6
    elif opt >= 18 and opt < 72:
        # 54 indices
        ind = opt - 18
        simnames = sl.m12plus_f3nobh_lores # 9
        snaps = sl.snaps_hr # 6

    simi = ind // len(snaps)
    snapi = ind % len(snaps)
    simname = simnames[simi]
    snapnum = snaps[snapi]

    #dirpath = '/'.join([_dirpath, simname]) 
    dirpath = sl.dirpath_from_simname(simname)

    #print(f'Whole halo center, {simname}, snap {snapnum}')
    #hp.gethalodata_shrinkingsphere(dirpath, snapnum, 
    #                               meandef=('BN98', '200c', '200m'))

    print(f'Whole halo Vcom, {simname}, snap {snapnum}')
    hp.get_vcom(dirpath, snapnum, 1., meandef_rvir=meandef,
                parttypes='all')
    print('\n\n')
    print(f'Galaxy Vcom, {simname}, snap {snapnum}')
    hp.get_vcom(dirpath, snapnum, 0.1, meandef_rvir=meandef,
                parttypes=(4,)) 
    print('\n\n')
    print(f'Galaxy re-centering, {simname}, snap {snapnum}')
    cgp.getcengalcen(dirpath, snapnum, startrad_rvir=0.3,
                     vcenrad_rvir=0.05, mstarrad_rvir=0.1)
    

def run_haloprop_fire3_m12new(opt):
    '''
    1 index runs the Halo (all particles but type 2, 1 Rvir) and
    the galaxy (parttype 4, 0.1 Rvir) velocity calculation,
    and the central galaxy re-centering (parttype 4 in 0.3 Rvir)
    '''
    # 1st run does not include m12g_r7100: 
    # Pratik is re-running after files were corrupted
    _dirpath = ol.simdir_fire3x_tests
    meandef = 'BN98'

    if opt >= 0 and opt < 24:
        ind = opt - 0
        simnames = sl.m12_fire3x_tests #4
        snaps = sl.snaps_sr # 6

    simi = ind // len(snaps)
    snapi = ind % len(snaps)
    simname = simnames[simi]
    snapnum = snaps[snapi]

    #dirpath = '/'.join([_dirpath, simname]) 
    dirpath = sl.dirpath_from_simname(simname)

    #print(f'Whole halo center, {simname}, snap {snapnum}')
    #hp.gethalodata_shrinkingsphere(dirpath, snapnum, 
    #                               meandef=('BN98', '200c', '200m'))

    print(f'Whole halo Vcom, {simname}, snap {snapnum}')
    hp.get_vcom(dirpath, snapnum, 1., meandef_rvir=meandef,
                parttypes='all')
    print('\n\n')
    print(f'Galaxy Vcom, {simname}, snap {snapnum}')
    hp.get_vcom(dirpath, snapnum, 0.1, meandef_rvir=meandef,
                parttypes=(4,)) 
    print('\n\n')
    print(f'Galaxy re-centering, {simname}, snap {snapnum}')
    cgp.getcengalcen(dirpath, snapnum, startrad_rvir=0.3,
                     vcenrad_rvir=0.05, mstarrad_rvir=0.1)
