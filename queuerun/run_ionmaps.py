
import os

import fire_an.mainfunc.makemap as mm
import fire_an.simlists as sl
import fire_an.utils.opts_locs as ol

def tryout_ionmap(opt=1):
    outdir = '/projects/b1026/nastasha/tests/start_fire/map_tests/'
    dirpath1 = '/projects/b1026/snapshots/fire3/m13h206_m3e5/' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000/'
    simname1 = 'm13h206_m3e5__' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'
    # version number depends on code edits; some indices might have
    # been run with previous versions
    _outfilen = 'coldens_{qt}_{sc}_snap{sn}_shrink-sph-cen_BN98' + \
                '_2rvir{depl}_v2.hdf5'
    checkfileflag = False
    if opt == 1:
        simname = simname1
        dirpath = dirpath1
        ions = ['o6']
        maptype = 'ion'
        _maptype_args = {'ps20depletion': False}
        snapnum = 27
    
        maptype_argss = [{key: _maptype_args[key] for key in _maptype_args} \
                         for ion in ions]
        [maptype_args.update({'ion': ion}) \
         for maptype_args, ion in zip(maptype_argss, ions)]
    
    elif opt > 1 and opt <= 9:
        dirpath = '/projects/b1026/isultan/m12i_noAGNfb/'
        simname = 'm12i_noAGNfb_CR-diff-coeff-690_FIRE-2'
        _ions = ['si4', 'n5', 'o6', 'ne8']
        ions = [_ions[(opt - 2) % 4]]
        maptype = 'ion'
        _maptype_args = {'ps20depletion': False}
        snapshots = [277, 600]
        snapnum = snapshots[(opt - 2) // 4]
        
        maptype_argss = [{key: _maptype_args[key] for key in _maptype_args} \
                         for ion in ions]
        [maptype_args.update({'ion': ion}) \
         for maptype_args, ion in zip(maptype_argss, ions)]
    
    elif opt >= 10 and opt < 21:
        # check ion sum, los metallicity
        simname = simname1
        dirpath = dirpath1
        snapnum = 27
        if opt < 19:
            ions = ['O{}'.format(opt - 9)]
            maptype = 'ion'
            _maptype_args = {'ps20depletion': False}
    
            maptype_argss = [{key: _maptype_args[key] for key in _maptype_args} \
                              for ion in ions]
            [maptype_args.update({'ion': ion}) \
             for maptype_args, ion in zip(maptype_argss, ions)]
        elif opt == 19:
            maptype = 'Metal'
            _maptype_args = {'element': 'Oxygen'}
            
            maptype_argss = [_maptype_args]
        elif opt == 20:
            maptype = 'Mass'
            _maptype_args = {}
            
            maptype_argss = [_maptype_args]
    elif opt >= 21 and opt < 57:
        # 36 indices; frontera paths
        outdir = '/scratch1/08466/tg877653/output/maps/set1_BH_noBH/'
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['Mass', 'O6', 'Ne8', 'N5', 'C2', 'Si2', 'Fe2', 'Mg2', 'Mg10']
        # z=0.4, 0.5, 0.6, redshifts match exactly at least for same ICs
        snaps = [50] # 49, 51
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # the m13 sdp1e-4 run (with BH) was only run down to z=1
        simnames = ['m12i_m6e4_MHD_fire3_fireBH_Sep052021_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm12i_m6e4_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h206_m3e5_MHD_fire3_fireBH_Sep052021_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                   ]
        ind = opt - 21
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()]

    elif opt >= 57 and opt < 93:
        # 36 indices; frontera paths
        outdir = '/scratch1/08466/tg877653/output/maps/set2_BH_noBH/'
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['Mass', 'O6', 'Ne8', 'N5', 'C2', 'Si2', 'Fe2', 'Mg2', 'Mg10']
        # z=0.4, 0.5, 0.6, redshifts match exactly at least for same ICs
        snaps = [50] # 49, 51
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # the m13s noBH run down to z=0 are all crdiff690
        # m13s with BH: m13h29 has crdiffc690, but noBH counterpart is not down to z=0
        #               same for m13h113, m13h236
        # m13h206 has m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_sdp1e-4_gacc31_fa0.5
        #             m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR0_vcr1000_Oct252021_crdiffc690_sdp1e-4_gacc31_fa0.5_fcr3e-4_vw3000
        # to z=0, with BH
        # noBH h206:  m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5
        #             m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_sdp1e10_gacc31_fa0.5  
        # -> ONE m13 option down to z=0 for same physics model BH/noBH comp.             
        # also only has one m12 match, for m12i
        # checked all have 60 snaps 
        simnames = ['m12i_m6e4_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm12i_m6e4_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_sdp1e10_gacc31_fa0.5',
                   ]
        ind = opt - 57
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()]

    elif opt >= 93 and opt < 189:
        # 96 indices; frontera paths
        ind = opt - 93
        outdir = '/scratch1/08466/tg877653/output/maps/set3_BH_noBH/'
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['Mass', 'O6', 'Ne8', 'Mg10', 'N5', 'Mg2']
        # z=0.4, 0.5, 0.6, redshifts match exactly at least for same ICs
        snaps = [49, 51, 45, 60] # 49, 51
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # the m13s noBH run down to z=0 are all crdiff690
        # m13s with BH: m13h29 has crdiffc690, but noBH counterpart is not down to z=0
        #               same for m13h113, m13h236
        # m13h206 has m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_sdp1e-4_gacc31_fa0.5
        #             m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR0_vcr1000_Oct252021_crdiffc690_sdp1e-4_gacc31_fa0.5_fcr3e-4_vw3000
        # to z=0, with BH
        # noBH h206:  m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5
        #             m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_sdp1e10_gacc31_fa0.5  
        # -> ONE m13 option down to z=0 for same physics model BH/noBH comp.             
        # also only has one m12 match, for m12i
        # checked all have 60 snaps 
        simnames = ['m12i_m6e4_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm12i_m6e4_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_sdp1e10_gacc31_fa0.5',
                   ]
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()]
    
    elif opt >= 189 and opt < 261:
        # split into m12 and m13 as two sets: different snapshot ranges
        ind = opt - 189
        outdir = '/scratch1/08466/tg877653/output/maps/set4_BH_noBH/'
        checkfileflag = True
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['Mass', 'O6', 'Ne8', 'Mg10', 'N5', 'Mg2']
        # z=0.0, 0.5, 1.0
        snaps = [500, 258, 186] #z=0.0, 0.50, 1.0
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # from Lindsey's selection, hr.
        # m12f since m12i noBH counterpart is crossed out
        simnames = ['m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
                    'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm12m_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
                    'm12m_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
                   ]
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()]
    
    elif opt >= 261 and opt < 333:
        # split into m12 and m13 as two sets: different snapshot ranges
        ind = opt - 261
        outdir = '/scratch1/08466/tg877653/output/maps/set5_BH_noBH/'
        checkfileflag = True
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['Mass', 'O6', 'Ne8', 'Mg10', 'N5', 'Mg2']
        # z=0.0, 0.5, 1.0
        snaps = [60, 50, 45] #z=0.0, 0.50, 1.0
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # from Lindsey's selection, sr.
        # these m13s selected because noBH ran down to z=0
        simnames = ['m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h007_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h007_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                   ]
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()]

    elif opt >= 333 and opt < 357:
        # split into m12 and m13 as two sets: different snapshot ranges
        ind = opt - 357
        outdir = '/scratch1/08466/tg877653/output/maps/set6_BH_noBH/'
        checkfileflag = True
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['Mass', 'O6', 'Ne8', 'Mg10', 'N5', 'Mg2']
        # z=0.0, 0.5, 1.0
        snaps = [50, 45] #z=0.50, 1.0
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # from Lindsey's selection, sr.
        # for comparison to the to z=0.0 noBH m13s: did run to z=0.5,
        #     but second-fewest snapshots of the noBH m13s
        simnames = ['m13h002_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h002_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                   ]
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()]
    elif opt == 357:
        # example of map with NaN values for debugging
        outdir = '/scratch1/08466/tg877653/output/maps/debug_mapnan/'
        checkfileflag = False
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simname = 'm13h002_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5'
        snapnum = 45
        ion = 'Mg10'
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)
        maptype = 'ion'
        _maptype_args = {'ps20depletion': False}
        _maptype_args.update({'ion': ion})
        maptype_argss = [_maptype_args.copy()]
    elif opt >= 358 and opt < 366:
        # split into m12 and m13 as two sets: different snapshot ranges
        ind = opt - 358
        outdir = '/scratch1/08466/tg877653/output/maps/clean_set1/'
        checkfileflag = True
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['Mass', 'O6', 'Ne8', 'Mg10']
        # z=0.0, 0.5, 1.0
        snaps = [50] #z=0.50,
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # from Lindsey's selection, sr.
        # for comparison to the to z=0.0 noBH m13s: did run to z=0.5,
        #     but second-fewest snapshots of the noBH m13s
        simnames = ['m12m_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                   ]
        # already have the AGN with CR and no BH versions of this in set4
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()]
    elif opt >= 366 and opt < 378:
        # split into m12 and m13 as two sets: different snapshot ranges
        ind = opt - 366
        outdir = '/scratch1/08466/tg877653/output/maps/clean_set1/'
        checkfileflag = True
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['Mass', 'O6', 'Ne8', 'Mg10']
        # z=0.0, 0.5, 1.0
        snaps = [258] #z=0.50,
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # from Lindsey's selection, sr.
        # for comparison to the to z=0.0 noBH m13s: did run to z=0.5,
        #     but second-fewest snapshots of the noBH m13s
        simnames = ['m13h029_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
                    'm13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31_fa0.5',
                   ]
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()]    
    elif opt >= 378 and opt < 394:
        # split into m12 and m13 as two sets: different snapshot ranges
        ind = opt - 378
        outdir = '/scratch1/08466/tg877653/output/maps/clean_set1/'
        checkfileflag = True
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['Mass', 'O6', 'Ne8', 'Mg10']
        # z=0.0, 0.5, 1.0
        snaps = [50] #z=0.50,
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # from Lindsey's selection, sr.
        # for comparison to the to z=0.0 noBH m13s: did run to z=0.5,
        #     but second-fewest snapshots of the noBH m13s
        simnames = ['m13h029_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h029_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                   ]
        # m13h206: already have AGN-no CR and no BH at z=0.5 (set 5)
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()]  
    elif opt >= 395 and opt < 398:
        # quest
        outdir = '/projects/b1026/nastasha/tests/start_fire/h1sim_tests/'
        checkfileflag = True
        # standard res M12, M13 w and w/o BH
        _dirpath = '/projects/b1026/snapshots/fire3/'
        # from Lindsey's selection, sr.
        # for comparison to the to z=0.0 noBH m13s: did run to z=0.5,
        #     but second-fewest snapshots of the noBH m13s
        simname = 'm12m_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5'
        # m13h206: already have AGN-no CR and no BH at z=0.5 (set 5)
        snapnum = 500
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)
        if opt == 395:
            maptype = 'Metal'
            _maptype_args = {'element': 'Hydrogen'}
            _outfilen = 'coldens_{qt}_{sc}_snap{sn}_shrink-sph-cen_BN98' + \
                '_2rvir{depl}_v2.hdf5'
        elif opt == 396:
            maptype = 'ion'
            _maptype_args = {'ion': 'H1', 'ionfrac-method': 'PS20',
                             'ps20depletion': False}
            _outfilen = 'coldens_{qt}-PS20_{sc}_snap{sn}_shrink-sph-cen_BN98' + \
                '_2rvir{depl}_v2.hdf5'
        elif opt == 397:
            maptype = 'ion'
            _maptype_args = {'ion': 'H1', 'ionfrac-method': 'sim',
                             'ps20depletion': False}
            _outfilen = 'coldens_{qt}-sim_{sc}_snap{sn}_shrink-sph-cen_BN98' + \
                '_2rvir{depl}_v2.hdf5'
        maptype_argss = [_maptype_args.copy()]  

    elif opt >= 398 and opt < 402:
        # add H 1 to clean sample (excluding m12m -- has a bug)
        ind = opt - 398
        outdir = '/scratch1/08466/tg877653/output/maps/clean_set1/'
        checkfileflag = True
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['H1']
        # z=0.0, 0.5, 1.0
        snaps = [50] #z=0.50,
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # from Lindsey's selection, sr.
        # for comparison to the to z=0.0 noBH m13s: did run to z=0.5,
        #     but second-fewest snapshots of the noBH m13s
        simnames = ['m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                   ]
        # m13h206: already have AGN-no CR and no BH at z=0.5 (set 5)
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            if ion == 'H1':
                _maptype_args = {'ps20depletion': False, 
                                 'ionfrac-method': 'sim'}
            else:
                _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()] 
    elif opt >= 402 and opt < 404:
        # add H 1 to clean sample (excluding m12m -- has a bug)
        ind = opt - 402
        outdir = '/scratch1/08466/tg877653/output/maps/clean_set1/'
        checkfileflag = True
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['H1']
        # z=0.0, 0.5, 1.0
        snaps = [258] #z=0.50,
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31_fa0.5',
                   ]
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            if ion == 'H1':
                _maptype_args = {'ps20depletion': False, 
                                 'ionfrac-method': 'sim'}
            else:
                _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()] 
    elif opt >= 404 and opt < 406:
        # add H 1 to clean sample (excluding m12m -- has a bug)
        ind = opt - 404
        outdir = '/scratch1/08466/tg877653/output/maps/clean_set1/'
        checkfileflag = True
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['H1']
        # z=0.0, 0.5, 1.0
        snaps = [258] #z=0.50,
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # from Lindsey's selection, sr.
        # for comparison to the to z=0.0 noBH m13s: did run to z=0.5,
        #     but second-fewest snapshots of the noBH m13s
        simnames = ['m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
                    ]
        # m13h206: already have AGN-no CR and no BH at z=0.5 (set 5)
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            if ion == 'H1':
                _maptype_args = {'ps20depletion': False, 
                                 'ionfrac-method': 'sim'}
            else:
                _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()] 
    elif opt >= 406 and opt < 407:
        # add H 1 to clean sample (excluding m12m -- has a bug)
        ind = opt - 406
        outdir = '/scratch1/08466/tg877653/output/maps/clean_set1/'
        checkfileflag = True
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['H1']
        # z=0.0, 0.5, 1.0
        snaps = [50] #z=0.50,
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # from Lindsey's selection, sr.
        # for comparison to the to z=0.0 noBH m13s: did run to z=0.5,
        #     but second-fewest snapshots of the noBH m13s
        simnames = ['m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    ]
        # m13h206: already have AGN-no CR and no BH at z=0.5 (set 5)
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            if ion == 'H1':
                _maptype_args = {'ps20depletion': False, 
                                 'ionfrac-method': 'sim'}
            else:
                _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()] 
    
    elif opt >= 407 and opt < 507:
        # add H 1 to clean sample (excluding m12m -- has a bug)
        ind = opt - 407
        outdir = '/scratch1/08466/tg877653/output/maps/clean_set2/'
        checkfileflag = True
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['Mass', 'O6', 'Ne8', 'Mg10', 'H1']
        # z=0.0, 0.5, 1.0
        snaps = [49, 48, 47, 46, 45] #z=0.6, 0.7, 0.8, 0.9, 1.0
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # from Lindsey's selection, sr.
        # for comparison to the to z=0.0 noBH m13s: did run to z=0.5,
        #     but second-fewest snapshots of the noBH m13s
        simnames = ['m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                   ]
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            if ion == 'H1':
                _maptype_args = {'ps20depletion': False, 
                                 'ionfrac-method': 'sim'}
            else:
                _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()] 
    elif opt >= 507 and opt < 557:
        # add z=0.6 - 1.0, 0.1 steps, to clean sample (set2)
        ind = opt - 507
        outdir = '/scratch1/08466/tg877653/output/maps/clean_set2/'
        checkfileflag = True
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['Mass', 'O6', 'Ne8', 'Mg10', 'H1']
        snaps = [240, 224, 210, 197, 186] #z=0.6, 0.7, 0.8, 0.9, 1.0
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31_fa0.5',
                   ]
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            if ion == 'H1':
                _maptype_args = {'ps20depletion': False, 
                                 'ionfrac-method': 'sim'}
            else:
                _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()] 
    elif opt >= 557 and opt < 607:
        # add H 1 to clean sample (excluding m12m -- has a bug)
        ind = opt - 557
        outdir = '/scratch1/08466/tg877653/output/maps/clean_set2/'
        checkfileflag = True
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['Mass', 'O6', 'Ne8', 'Mg10', 'H1']
        snaps = [240, 224, 210, 197, 186] #z=0.6, 0.7, 0.8, 0.9, 1.0
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # from Lindsey's selection, sr.
        # for comparison to the to z=0.0 noBH m13s: did run to z=0.5,
        #     but second-fewest snapshots of the noBH m13s
        simnames = ['m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
                    ]
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            if ion == 'H1':
                _maptype_args = {'ps20depletion': False, 
                                 'ionfrac-method': 'sim'}
            else:
                _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()] 
    elif opt >= 607 and opt < 632:
        # add H 1 to clean sample (excluding m12m -- has a bug)
        ind = opt - 607
        outdir = '/scratch1/08466/tg877653/output/maps/clean_set2/'
        checkfileflag = True
        # CUBS https://arxiv.org/pdf/2209.01228.pdf: 
        # At 𝑧≈1, HST/COS FUV spectra cover a wide
        # range of ions, including 
        # H i, He i, Cii, N ii to N iv, O i to O v, S ii to
        # S v, Ne iv to Ne vi, Ne viii, and Mg x
        # kinda random subset of those, H I not yet FIRE-consistent
        ions = ['Mass', 'O6', 'Ne8', 'Mg10', 'H1']
        snaps = [49, 48, 47, 46, 45] #z=0.6, 0.7, 0.8, 0.9, 1.0 #z=0.50,
        # standard res M12, M13 w and w/o BH
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        # from Lindsey's selection, sr.
        # for comparison to the to z=0.0 noBH m13s: did run to z=0.5,
        #     but second-fewest snapshots of the noBH m13s
        simnames = ['m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    ]
        # m13h206: already have AGN-no CR and no BH at z=0.5 (set 5)
        simi = ind // (len(snaps) * len(ions))
        snpi = (ind % (len(snaps) * len(ions))) // len(ions)
        ioni = ind % len(ions)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        ion = ions[ioni]
        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])
        #print(dirpath)

        if ion == 'Mass':
            maptype = 'Mass'
            maptype_argss = [{}]
        else:
            maptype = 'ion'
            if ion == 'H1':
                _maptype_args = {'ps20depletion': False, 
                                 'ionfrac-method': 'sim'}
            else:
                _maptype_args = {'ps20depletion': False}
            _maptype_args.update({'ion': ion})
            maptype_argss = [_maptype_args.copy()] 

    for maptype_args in maptype_argss:
        depl = ''
        if maptype == 'ion':
            qt = maptype_args['ion']
            if 'ionfrac-method' in maptype_args:
                if maptype_args['ionfrac-method'] == 'sim':
                    depl = '_ionfrac-fromsim'
                else:
                    _depl = maptype_args['ps20depletion']
                    if _depl:
                        depl = '_ps20-depl'
            else:
                _depl = maptype_args['ps20depletion']
                if _depl:
                    depl = '_ps20-depl'
        elif maptype == 'Metal':
            qt = maptype_args['element']
        elif maptype == 'Mass':
            qt = 'gas-mass'

        outfilen = outdir + _outfilen.format(sc=simname, sn=snapnum, 
                                             depl=depl, qt=qt)
        if checkfileflag:
            if os.path.isfile(outfilen):
                msg = 'For opt {}, output file already exists:\n{}'
                print(msg.format(opt, outfilen))
                print('Not running this map again')
                return None

        mm.massmap(dirpath, snapnum, radius_rvir=2., particle_type=0,
                   pixsize_pkpc=3., axis='z', outfilen=outfilen,
                   center='shrinksph', norm='pixsize_phys',
                   maptype=maptype, maptype_args=maptype_args)

def run_ionmap_xyz(opt=0):
    # version number depends on code edits; some indices might have
    # been run with previous versions
    _outfilen = ('coldens_{qt}_{sc}_snap{sn}_shrink-sph-cen_BN98'
                 '_2rvir{depl}_{los}-proj_v3.hdf5')
    _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
    checkfileflag = False
    if opt >= 0 and opt < 1350:
        # m13sr: 1350 indices
        ind = opt - 0
        outdir = '/scratch1/08466/tg877653/output/maps/set3_model3/'
        checkfileflag = True
        ions = ['Mass', 'O6', 'Ne8', 'Mg10', 'Neon']
        loss = ['x', 'y', 'z']
        simnames = sl.m13_sr_all1 # len 15
        snaps = sl.snaplists['m13_sr'] # len 6
    elif opt >= 1350 and opt < 1530:
        # m13hr: 180 indices
        ind = opt - 1350
        outdir = '/scratch1/08466/tg877653/output/maps/set3_model3/'
        checkfileflag = True
        ions = ['Mass', 'O6', 'Ne8', 'Mg10', 'Neon']
        loss = ['x', 'y', 'z']
        simnames = sl.m13_hr_all1 # len 2
        snaps = sl.snaplists['m13_hr'] # len 6
        checkfileflag = True
    elif opt >= 1530 and opt < 1890:
        # m12sr: 360 indices
        ind = opt - 1530
        outdir = '/scratch1/08466/tg877653/output/maps/set3_model3/'
        checkfileflag = True
        ions = ['Mass', 'O6', 'Ne8', 'Mg10', 'Neon']
        loss = ['x', 'y', 'z']
        simnames = sl.m12_sr_all1 # len 4
        snaps = sl.snaplists['m12_sr'] # len 6
    elif opt >= 1890 and opt < 3330:
        # m12hr: 1440 indices
        ind = opt - 1890
        outdir = '/scratch1/08466/tg877653/output/maps/set3_model3/'
        checkfileflag = True
        ions = ['Mass', 'O6', 'Ne8', 'Mg10', 'Neon']
        loss = ['x', 'y', 'z']
        simnames = sl.m12_hr_all1 # len 16
        snaps = sl.snaplists['m12_hr'] # len 6
    # m12q runs: all1 -> all2
    elif opt >= 3330 and opt < 3510:
        # m12q additional: 180 indices
        ind = opt - 3330
        outdir = '/scratch1/08466/tg877653/output/maps/set3_model3/'
        checkfileflag = True
        ions = ['Mass', 'O6', 'Ne8', 'Mg10', 'Neon']
        loss = ['x', 'y', 'z']
        simnames = [('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp2e-4_gacc31_fa0.5')
                   ] # len 2
        snaps = sl.snaplists['m12_hr'] # len 6
        
    simi = ind // (len(snaps) * len(ions) * len(loss))
    snpi = (ind % (len(snaps) * len(ions) * len(loss))) \
           // (len(ions) * len(loss))
    ioni = (ind % (len(ions) * len(loss))) // len(loss)
    losi = ind % len(loss)
    simname = simnames[simi]
    snapnum = snaps[snpi]
    ion = ions[ioni]
    los = loss[losi]

    # directory is halo name + resolution 
    dp2 = '_'.join(simname.split('_')[:2])
    if dp2.startswith('m13h02_'):
        dp2 = dp2.replace('m13h02', 'm13h002')
    dirpath = '/'.join([_dirpath, dp2, simname])
    #print(dirpath)

    if ion == 'Mass':
        maptype = 'Mass'
        maptype_argss = [{}]
    elif ion in ['Hydrogen', 'Helium', 'Carbon', 'Nitrogen', 'Oxygen',
                 'Neon', 'Magnesium', 'Silicon', 'Iron']:
        maptype = 'Metal'
        maptype_argss = [{'element': ion, 'density': False}]
    else:
        maptype = 'ion'
        if ion == 'H1':
            _maptype_args = {'ps20depletion': False, 
                             'ionfrac-method': 'sim'}
        else:
            _maptype_args = {'ps20depletion': False}
        _maptype_args.update({'ion': ion})
        maptype_argss = [_maptype_args.copy()] 

    for maptype_args in maptype_argss:
        depl = ''
        if maptype == 'ion':
            qt = maptype_args['ion']
            if 'ionfrac-method' in maptype_args:
                if maptype_args['ionfrac-method'] == 'sim':
                    depl = '_ionfrac-fromsim'
                else:
                    _depl = maptype_args['ps20depletion']
                    if _depl:
                        depl = '_ps20-depl'
            else:
                _depl = maptype_args['ps20depletion']
                if _depl:
                    depl = '_ps20-depl'
        elif maptype == 'Metal':
            qt = maptype_args['element']
        elif maptype == 'Mass':
            qt = 'gas-mass'

        outfilen = outdir + _outfilen.format(sc=simname, sn=snapnum, 
                                             depl=depl, qt=qt, los=los)
        if checkfileflag:
            if os.path.isfile(outfilen):
                msg = 'For opt {}, output file already exists:\n{}'
                print(msg.format(opt, outfilen))
                print('Not running this map again')
                return None

        mm.massmap(dirpath, snapnum, radius_rvir=2., particle_type=0,
                   pixsize_pkpc=3., axis=los, outfilen=outfilen,
                   center='shrinksph', norm='pixsize_phys',
                   maptype=maptype, maptype_args=maptype_args)

def run_vlosmaps(opt=0):
    # version number depends on code edits; some indices might have
    # been run with previous versions
    _outfilen = ('vlos_by_coldens_{qt}_{sc}_snap{sn}_shrink-sph-cen_BN98'
                 '_depth_{llos:.1f}rvir{depl}_{los}-proj_v3.hdf5')
    _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
    outdir = '/scratch1/08466/tg877653/output/maps/velmaps_clean2/'
    checkfileflag = True
    # 48 maps per halo/snap
    ions = ['Mass', 'Volume', 'Oxygen', 'Neon',
            'O6', 'O7', 'O8', 'Ne8']
    loss = ['x', 'y', 'z']
    rloss_rvir = [0.1, 2.] 
    wtdtype = 'coords'
    wtdtype_args = {'vel': 'los'}
    if opt >= 0 and opt < 1152:
        # m13sr: 1152 indices
        ind = opt - 0
        simnames = sl.m13_sr_clean2 # len 4
        snaps = sl.snaplists['m13_sr'] # len 6
    elif opt >= 1152 and opt < 1728:
        # m13hr: 576 indices
        ind = opt - 1152
        simnames = sl.m13_hr_clean2 # len 2
        snaps = sl.snaplists['m13_hr'] # len 6
    elif opt >= 1728 and opt < 2304:
        # m12sr: 576 indices
        ind = opt - 1728
        simnames = [('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
                    ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
                     '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000')
                    ]
        snaps = sl.snaplists['m12_sr'] # len 6
    elif opt >= 2304 and opt < 3456:
        # m12hr: 1152 indices
        ind = opt - 2304
        simnames = [('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp2e-4_gacc31_fa0.5'),
                    ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp2e-4_gacc31_fa0.5'),
                    ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
                     '_sdp1e10_gacc31_fa0.5'),
                    ] # len 4
        snaps = sl.snaplists['m12_hr'] # len 6
    ## m12, Ne8 full sample
    elif opt >= 3456 and opt < 3492:
        # 3 axes, 6 snaps, 2 haloes -> 36 indices
        ind = opt - 3456
        ions = ['Ne8']
        rloss_rvir = [2.]
        simnames_done = [
            ('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
             '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
            ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
             '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000')] # len 2
        simnames = sl.m12_sr_all2
        for _sn in simnames_done:
            simnames.remove(_sn)
        # len 2
        snaps = sl.snaps_sr
    elif opt >= 3492 and opt < 3744:
        # 3 axes, 6 snaps, 14 haloes -> 252 indices
        ind = opt - 3492
        ions = ['Ne8']
        rloss_rvir = [2.]
        simnames_done = [
            ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
             '_sdp2e-4_gacc31_fa0.5'),
            ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
             '_sdp1e10_gacc31_fa0.5'),
            ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
             '_sdp2e-4_gacc31_fa0.5'),
            ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
             '_sdp1e10_gacc31_fa0.5')]
        simnames = sl.m12_hr_all2
        for _sn in simnames_done:
            simnames.remove(_sn)
        # len 14
        snaps = sl.snaps_hr
    # m13, Ne8 full sample
    # m13-hr: all in the clean sample
    elif opt >= 3744 and opt < 3942:
        # 3 axes, 6 snaps, 11 haloes -> 198 indices
        ind = opt - 3744
        ions = ['Ne8']
        rloss_rvir = [2.]
        # z = 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
        # m13-sr
        # 1320 indices
        simnames = sl.m13_sr_all2
        for _sn in sl.m13_sr_clean2:
            simnames.remove(_sn)
        # len 11
        snaps = sl.snaps_sr # len 6
        
    simi = ind // (len(snaps) * len(ions) * len(loss) * len(rloss_rvir))
    snpi = (ind % (len(snaps) * len(ions) * len(loss) * len(rloss_rvir))) \
           // (len(ions) * len(loss) * len(rloss_rvir))
    ioni = (ind % (len(ions) * len(loss) * len(rloss_rvir))) \
            // (len(loss) * len(rloss_rvir))
    losi = (ind % len(loss) * len(rloss_rvir)) // len(rloss_rvir)
    rlosi = ind % len(rloss_rvir)
    simname = simnames[simi]
    snapnum = snaps[snpi]
    ion = ions[ioni]
    los = loss[losi]
    rlos_rvir = rloss_rvir[rlosi]

    # directory is halo name + resolution 
    dp2 = '_'.join(simname.split('_')[:2])
    if dp2.startswith('m13h02_'):
        dp2 = dp2.replace('m13h02', 'm13h002')
    dirpath = '/'.join([_dirpath, dp2, simname])
    #print(dirpath)

    if ion == 'Mass':
        maptype = 'Mass'
        maptype_argss = [{}]
    elif ion == 'Volume':
        maptype = 'Volume'
        maptype_argss = [{}]
    elif ion in ['Hydrogen', 'Helium', 'Carbon', 'Nitrogen', 'Oxygen',
                 'Neon', 'Magnesium', 'Silicon', 'Iron']:
        maptype = 'Metal'
        maptype_argss = [{'element': ion, 'density': False}]
    else:
        maptype = 'ion'
        if ion == 'H1':
            _maptype_args = {'ps20depletion': False, 
                             'ionfrac-method': 'sim'}
        else:
            _maptype_args = {'ps20depletion': False}
        _maptype_args.update({'ion': ion})
        maptype_argss = [_maptype_args.copy()] 

    for maptype_args in maptype_argss:
        depl = ''
        if maptype == 'ion':
            qt = maptype_args['ion']
            if 'ionfrac-method' in maptype_args:
                if maptype_args['ionfrac-method'] == 'sim':
                    depl = '_ionfrac-fromsim'
                else:
                    _depl = maptype_args['ps20depletion']
                    if _depl:
                        depl = '_ps20-depl'
            else:
                _depl = maptype_args['ps20depletion']
                if _depl:
                    depl = '_ps20-depl'
        elif maptype == 'Metal':
            qt = maptype_args['element']
        elif maptype == 'Mass':
            qt = 'gas-mass'
        elif maptype == 'Volume':
            qt = 'gas-vol'

        outfilen = outdir + _outfilen.format(sc=simname, sn=snapnum, 
                                             depl=depl, qt=qt, los=los,
                                             llos=rlos_rvir)
        if checkfileflag:
            if os.path.isfile(outfilen):
                msg = 'For opt {}, output file already exists:\n{}'
                print(msg.format(opt, outfilen))
                print('Not running this map again')
                return None

        mm.massmap(dirpath, snapnum, radius_rvir=2., particle_type=0,
                   pixsize_pkpc=3., axis=los, outfilen=outfilen,
                   center='shrinksph', norm='pixsize_phys',
                   weighttype=maptype, weighttype_args=maptype_args,
                   maptype=wtdtype, maptype_args=wtdtype_args,
                   save_weightmap=True, logmap=False,
                   logweightmap=True, losradius_rvir=rlos_rvir)
        

def run_vdoplosmaps(opt=0):
    # version number depends on code edits; some indices might have
    # been run with previous versions
    _outfilen = ('vdoplos_by_coldens_{qt}_{sc}_snap{sn}_shrink-sph-cen_BN98'
                 '_depth_{llos:.1f}rvir{depl}_{los}-proj_v3.hdf5')
    _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
    outdir = '/scratch1/08466/tg877653/output/maps/vdopmaps_clean2/'
    checkfileflag = True
    # 24 maps per halo/snap
    ions = ['Mass', 'Volume', 'Neon', 'Ne8']
    loss = ['x', 'y', 'z']
    rloss_rvir = [0.1, 2.]
    wtdtype = 'coords'
    wtdtype_args = {'vel': 'doplos'}
    if opt >= 0 and opt < 2160:
        # m13sr: 2160 indices
        ind = opt - 0
        simnames = sl.m13_sr_all2 # len 15
        snaps = sl.snaplists['m13_sr'] # len 6
    elif opt >= 2160 and opt < 2448:
        # m13hr: 288 indices
        ind = opt - 2160
        simnames = sl.m13_hr_all2 # len 2
        snaps = sl.snaplists['m13_hr'] # len 6
    elif opt >= 2448 and opt < 3024:
        # m12sr: 576 indices
        ind = opt - 2448
        simnames = sl.m12_sr_all2 # len 4
        snaps = sl.snaplists['m12_sr'] # len 6
    elif opt >= 3024 and opt < 5616:
        # m12hr: 2592 indices
        ind = opt - 3024
        simnames = sl.m12_hr_all2 # len 18
        snaps = sl.snaplists['m12_hr'] # len 6
        
    simi = ind // (len(snaps) * len(ions) * len(loss) * len(rloss_rvir))
    snpi = (ind % (len(snaps) * len(ions) * len(loss) * len(rloss_rvir))) \
           // (len(ions) * len(loss) * len(rloss_rvir))
    ioni = (ind % (len(ions) * len(loss) * len(rloss_rvir))) \
            // (len(loss) * len(rloss_rvir))
    losi = (ind % len(loss) * len(rloss_rvir)) // len(rloss_rvir)
    rlosi = ind % len(rloss_rvir)
    simname = simnames[simi]
    snapnum = snaps[snpi]
    ion = ions[ioni]
    los = loss[losi]
    rlos_rvir = rloss_rvir[rlosi]

    # directory is halo name + resolution 
    dp2 = '_'.join(simname.split('_')[:2])
    if dp2.startswith('m13h02_'):
        dp2 = dp2.replace('m13h02', 'm13h002')
    dirpath = '/'.join([_dirpath, dp2, simname])
    #print(dirpath)

    if ion == 'Mass':
        maptype = 'Mass'
        maptype_argss = [{}]
    elif ion == 'Volume':
        maptype = 'Volume'
        maptype_argss = [{}]
    elif ion in ['Hydrogen', 'Helium', 'Carbon', 'Nitrogen', 'Oxygen',
                 'Neon', 'Magnesium', 'Silicon', 'Iron']:
        maptype = 'Metal'
        maptype_argss = [{'element': ion, 'density': False}]
    else:
        maptype = 'ion'
        if ion == 'H1':
            _maptype_args = {'ps20depletion': False, 
                             'ionfrac-method': 'sim'}
        else:
            _maptype_args = {'ps20depletion': False}
        _maptype_args.update({'ion': ion})
        maptype_argss = [_maptype_args.copy()] 

    for maptype_args in maptype_argss:
        depl = ''
        if maptype == 'ion':
            qt = maptype_args['ion']
            if 'ionfrac-method' in maptype_args:
                if maptype_args['ionfrac-method'] == 'sim':
                    depl = '_ionfrac-fromsim'
                else:
                    _depl = maptype_args['ps20depletion']
                    if _depl:
                        depl = '_ps20-depl'
            else:
                _depl = maptype_args['ps20depletion']
                if _depl:
                    depl = '_ps20-depl'
        elif maptype == 'Metal':
            qt = maptype_args['element']
        elif maptype == 'Mass':
            qt = 'gas-mass'
        elif maptype == 'Volume':
            qt = 'gas-vol'

        outfilen = outdir + _outfilen.format(sc=simname, sn=snapnum, 
                                             depl=depl, qt=qt, los=los,
                                             llos=rlos_rvir)
        if checkfileflag:
            if os.path.isfile(outfilen):
                msg = 'For opt {}, output file already exists:\n{}'
                print(msg.format(opt, outfilen))
                print('Not running this map again')
                return None

        mm.massmap(dirpath, snapnum, radius_rvir=2., particle_type=0,
                   pixsize_pkpc=3., axis=los, outfilen=outfilen,
                   center='shrinksph', norm='pixsize_phys',
                   weighttype=maptype, weighttype_args=maptype_args,
                   maptype=wtdtype, maptype_args=wtdtype_args,
                   save_weightmap=True, logmap=False,
                   logweightmap=True, losradius_rvir=rlos_rvir)

def run_vdoplosmaps_f2md(opt=0):
    # version number depends on code edits; some indices might have
    # been run with previous versions
    _outfilen = ('vdoplos_by_coldens_{qt}_{sc}_snap{sn}_shrink-sph-cen_BN98'
                 '_depth_{llos:.1f}rvir{depl}_{los}-proj_v3.hdf5')
    _dirpath = '/scratch/projects/xsede/GalaxiesOnFIRE/metal_diffusion/'
    outdir = '/scratch/08466/tg877653/output/maps/vdopmaps/'
    checkfileflag = True
    # 24 maps per halo/snap
    ions = ['Mass', 'Volume', 'Neon', 'Ne8']
    loss = ['x', 'y', 'z']
    rloss_rvir = [0.1, 2.]
    wtdtype = 'coords'
    wtdtype_args = {'vel': 'doplos'}
    
    # 48 haloes (0 - 1151)
    # + 12 haloes (crheatingfix; 1152 - 1440, 288 indices)
    ind = opt - 0
    simnames = sl.m12_f2md # len 8
    snaps = sl.snaps_f2md # len 6
        
    simi = ind // (len(snaps) * len(ions) * len(loss) * len(rloss_rvir))
    snpi = (ind % (len(snaps) * len(ions) * len(loss) * len(rloss_rvir))) \
           // (len(ions) * len(loss) * len(rloss_rvir))
    ioni = (ind % (len(ions) * len(loss) * len(rloss_rvir))) \
            // (len(loss) * len(rloss_rvir))
    losi = (ind % len(loss) * len(rloss_rvir)) // len(rloss_rvir)
    rlosi = ind % len(rloss_rvir)
    simname = simnames[simi]
    snapnum = snaps[snpi]
    ion = ions[ioni]
    los = loss[losi]
    rlos_rvir = rloss_rvir[rlosi]

    dirpath = sl.dirpath_from_simname(simname)

    if ion == 'Mass':
        maptype = 'Mass'
        maptype_argss = [{}]
    elif ion == 'Volume':
        maptype = 'Volume'
        maptype_argss = [{}]
    elif ion in ['Hydrogen', 'Helium', 'Carbon', 'Nitrogen', 'Oxygen',
                 'Neon', 'Magnesium', 'Silicon', 'Iron']:
        maptype = 'Metal'
        maptype_argss = [{'element': ion, 'density': False}]
    else:
        maptype = 'ion'
        if ion == 'H1':
            _maptype_args = {'ps20depletion': False, 
                             'ionfrac-method': 'sim'}
        else:
            _maptype_args = {'ps20depletion': False}
        _maptype_args.update({'ion': ion})
        maptype_argss = [_maptype_args.copy()] 

    for maptype_args in maptype_argss:
        depl = ''
        if maptype == 'ion':
            qt = maptype_args['ion']
            if 'ionfrac-method' in maptype_args:
                if maptype_args['ionfrac-method'] == 'sim':
                    depl = '_ionfrac-fromsim'
                else:
                    _depl = maptype_args['ps20depletion']
                    if _depl:
                        depl = '_ps20-depl'
            else:
                _depl = maptype_args['ps20depletion']
                if _depl:
                    depl = '_ps20-depl'
        elif maptype == 'Metal':
            qt = maptype_args['element']
        elif maptype == 'Mass':
            qt = 'gas-mass'
        elif maptype == 'Volume':
            qt = 'gas-vol'

        outfilen = outdir + _outfilen.format(sc=simname, sn=snapnum, 
                                             depl=depl, qt=qt, los=los,
                                             llos=rlos_rvir)
        if checkfileflag:
            if os.path.isfile(outfilen):
                msg = 'For opt {}, output file already exists:\n{}'
                print(msg.format(opt, outfilen))
                print('Not running this map again')
                return None

        mm.massmap(dirpath, snapnum, radius_rvir=2., particle_type=0,
                   pixsize_pkpc=3., axis=los, outfilen=outfilen,
                   center='shrinksph', norm='pixsize_phys',
                   weighttype=maptype, weighttype_args=maptype_args,
                   maptype=wtdtype, maptype_args=wtdtype_args,
                   save_weightmap=True, logmap=False,
                   logweightmap=True, losradius_rvir=rlos_rvir)

def run_vdoplosmaps_m12new(opt=0):
    # version number depends on code edits; some indices might have
    # been run with previous versions
    _outfilen = ('vdoplos_by_coldens_{qt}_{sc}_snap{sn}_shrink-sph-cen_BN98'
                 '_depth_{llos:.1f}rvir{depl}_{los}-proj_v3.hdf5')
    _dirpath = ol.simdir_fire3_m12plus
    outdir = '/projects/b1026/nastasha/maps/vdopmaps_all2/'
    checkfileflag = True
    # 6 maps per halo/snap
    ions = ['Ne8']
    loss = ['x', 'y', 'z']
    rloss_rvir = [0.1, 2.]
    wtdtype = 'coords'
    wtdtype_args = {'vel': 'doplos'}
    
    if opt >=0 and opt < 108:
        # 18 haloes 
        # m12g not included in the first run
        ind = opt - 0
        simnames = sl.m12plus_f3nobh # len 3
        snaps = sl.snaps_hr # len 6
    elif opt >= 108 and opt < 594:
        # 54 haloes, 324 indices
        ind = opt - 108
        simnames = sl.m12plus_f3nobh_lores # len 9
        snaps = sl.snaps_hr # len 6
    elif opt >= 594 and opt < 738:
        # frontera
        outdir = '/scratch1/08466/tg877653/output/maps/vdopmaps_f3x/'
        # 24 halos, 144 indices
        ind = opt - 594
        simnames = sl.m12_fire3x_tests
        snaps = sl.snaps_sr

        
    simi = ind // (len(snaps) * len(ions) * len(loss) * len(rloss_rvir))
    snpi = (ind % (len(snaps) * len(ions) * len(loss) * len(rloss_rvir))) \
           // (len(ions) * len(loss) * len(rloss_rvir))
    ioni = (ind % (len(ions) * len(loss) * len(rloss_rvir))) \
            // (len(loss) * len(rloss_rvir))
    losi = (ind % len(loss) * len(rloss_rvir)) // len(rloss_rvir)
    rlosi = ind % len(rloss_rvir)
    simname = simnames[simi]
    snapnum = snaps[snpi]
    ion = ions[ioni]
    los = loss[losi]
    rlos_rvir = rloss_rvir[rlosi]

    dirpath = sl.dirpath_from_simname(simname)

    if ion == 'Mass':
        maptype = 'Mass'
        maptype_argss = [{}]
    elif ion == 'Volume':
        maptype = 'Volume'
        maptype_argss = [{}]
    elif ion in ['Hydrogen', 'Helium', 'Carbon', 'Nitrogen', 'Oxygen',
                 'Neon', 'Magnesium', 'Silicon', 'Iron']:
        maptype = 'Metal'
        maptype_argss = [{'element': ion, 'density': False}]
    else:
        maptype = 'ion'
        if ion == 'H1':
            _maptype_args = {'ps20depletion': False, 
                             'ionfrac-method': 'sim'}
        else:
            _maptype_args = {'ps20depletion': False}
        _maptype_args.update({'ion': ion})
        maptype_argss = [_maptype_args.copy()] 

    for maptype_args in maptype_argss:
        depl = ''
        if maptype == 'ion':
            qt = maptype_args['ion']
            if 'ionfrac-method' in maptype_args:
                if maptype_args['ionfrac-method'] == 'sim':
                    depl = '_ionfrac-fromsim'
                else:
                    _depl = maptype_args['ps20depletion']
                    if _depl:
                        depl = '_ps20-depl'
            else:
                _depl = maptype_args['ps20depletion']
                if _depl:
                    depl = '_ps20-depl'
        elif maptype == 'Metal':
            qt = maptype_args['element']
        elif maptype == 'Mass':
            qt = 'gas-mass'
        elif maptype == 'Volume':
            qt = 'gas-vol'

        outfilen = outdir + _outfilen.format(sc=simname, sn=snapnum, 
                                             depl=depl, qt=qt, los=los,
                                             llos=rlos_rvir)
        if checkfileflag:
            if os.path.isfile(outfilen):
                msg = 'For opt {}, output file already exists:\n{}'
                print(msg.format(opt, outfilen))
                print('Not running this map again')
                return None

        mm.massmap(dirpath, snapnum, radius_rvir=2., particle_type=0,
                   pixsize_pkpc=3., axis=los, outfilen=outfilen,
                   center='shrinksph', norm='pixsize_phys',
                   weighttype=maptype, weighttype_args=maptype_args,
                   maptype=wtdtype, maptype_args=wtdtype_args,
                   save_weightmap=True, logmap=False,
                   logweightmap=True, losradius_rvir=rlos_rvir)