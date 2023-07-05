
import sys

import fire_an.explore.run_partsel_pv as spv
import fire_an.mainfunc.haloprop as hp
import fire_an.tests.test_haloprops as th
import fire_an.tests.test_maps as tm
import fire_an.tests.test_readfire as trf
import fire_an.tests.test_ionbal as tib
import fire_an.queuerun.run_ionmaps as rim
import fire_an.queuerun.run_hists as rhs
import fire_an.queuerun.run_haloprop as rhp


def fromcommandline(index):
    '''
    This mapping is just based on the order in which I (first) ran things,
    and will not generally follow any kind of logic
    '''
    print('Running fire_maps.py process {}'.format(index))
    if index > 0 and index < 4:
        th.test_mainhalodata_units_ahf(opt=index)
    elif index == 4:
        # test a whole lot of snapshots in one go
        th.test_mainhalodata_units_multi_handler(opt=1)
    elif index == 5:
        th.test_mainhalodata_units_multi_handler(opt=2)
    elif index == 6:
        tm.tryout_massmap(opt=1)
    elif index == 7:
        tm.tryout_massmap(opt=2)
    elif index == 8:
        tm.tryout_massmap(opt=3)
    elif index > 8 and index <= 11: # opt starts at 1
        opt = index - 8
        msg = 'Calling test_mainhalodata_units_rockstar(opt={})'
        print(msg.format(opt))
        th.test_mainhalodata_units_rockstar(opt=opt)
    elif index in [12, 13]:
        opt = index - 12
        trf.run_checkfields_units(opt)
    elif index >= 14 and index < 20:
        tib.run_ionbal_test(opt=index - 14)
    elif index == 20:
        tm.tryout_ionmap(opt=1)
    elif index >= 21 and index < 33:
        # opt in [6, 18)
        tib.run_ionbal_test(opt=index - 15)
    elif index >= 33 and index < 41:
        tm.tryout_ionmap(opt=index - 31)
    elif index >= 41 and index < 52:
        tm.tryout_ionmap(opt=index - 31)
    elif index == 52:
        th.tryout_hist(0)
    elif index >= 53 and index < 58:
        # launcher + script loading test
        print('Hello from index {}'.format(index))
    elif index >= 58 and index < 94:
        # set 1 maps -- 1 sim failed
        rim.tryout_ionmap(opt=index - 58 + 21)
    elif index >= 94 and index < 130:
        # set 2 maps 
        rim.tryout_ionmap(opt=index - 94 + 57) 
    elif index >= 130 and index < 226:
        # set 3 maps
        rim.tryout_ionmap(opt=index - 130 + 93)
    elif index >= 226 and index < 394:
        # sets 4, 5, 6
        # set 4 (m12, high-res): 226 - 297 (72 inds)
        # sets 5,6 (m13, standard-res): 298 - 393 (96 inds)
        rim.tryout_ionmap(opt=index - 226 + 189)
    elif index == 394:
        # NaN values in maps debug: single map example
        rim.tryout_ionmap(opt=357)
    elif index >= 395 and index < 431:
        rim.tryout_ionmap(opt=index - 395 + 358)
        # clean sample set 1: 
        # the parts that weren't already in sets 4-6
        # 395 - 402: m12 lower-res
        # 403 - 414: m13 higher-res
        # 415 - 430: m13 standard-res
    elif index >= 431 and index < 434:
        rim.tryout_ionmap(opt=index - 431 + 395)
        # two H I methods and H total: sanity check for H1-sim impl.
    elif index >= 434 and index < 443:
        rim.tryout_ionmap(opt=index - 434 + 398)
        # H I maps for clean sample z=0.5
        # 434 - 437: m13 standard-res
        # 438 - 439: m13 hi-res
        # 440 - 441: m12 hi-res
        # 442:       m12 standard-res
    elif index >= 443 and index < 668:
        rim.tryout_ionmap(opt=index - 443 + 407)
        # 4 ions + mass for 5 redshifts 0.6 - 1.0, clean sample 
        # (no m12m)
        # 443 - 542: m13-SR (4 IC/phys)
        # 543 - 592: m13-HR (2 phys)
        # 593 - 642: m12-HR (2 phys)
        # 643 - 667: m12-SR (1 IC/phys)
        # opts [407, 632) (9 x 25 indices)
    elif index >= 668 and index < 884:
        rhs.run_hist(index - 668 + 0)
        # z=0.5 only
        # Mass, Volume, H I: (T, rho, O, Ne, Mg) profiles
        # O6, Ne8, Mg10: (T, rho, parent element) profiles
        # 668 - 727: m13-SR (4 IC/phys), Mass, Volume, HI
        # 728 - 763: m13-SR (4 IC/phys), O6, Ne8, Mg10
        # 764 - 793: m12-HR (2 IC/phys), Mass, Volume, HI
        # 794 - 811: m12-HR (2 IC/phys), O6, Ne8, Mg10
        # 812 - 841: m13-HR (2 IC/phys), Mass, Volume, HI
        # 842 - 859: m13-HR (2 IC/phys), O6, Ne8, Mg10
        # 860 - 874: m12-SR (1 IC/phys), Mass, Volume, HI
        # 875 - 883: m12-SR (1 IC/phys), O6, Ne8, Mg10
    elif index >= 884 and index < 1964:
        rhs.run_hist(index - 884 + 216)
        # z=0.6, 0.7, 0.8, 0.9, 1.0
        # Mass, Volume, H I: (T, rho, O, Ne, Mg) profiles
        # O6, Ne8, Mg10: (T, rho, parent element) profiles
        #  884 - 1183: m13-SR (4 IC/phys), Mass, Volume, HI
        # 1184 - 1363: m13-SR (4 IC/phys), O6, Ne8, Mg10
        # (480 inds)
        # 1364 - 1513: m12-HR (2 IC/phys), Mass, Volume, HI
        # 1514 - 1603: m12-HR (2 IC/phys), O6, Ne8, Mg10
        # (240 inds)
        # 1604 - 1753: m13-HR (2 IC/phys), Mass, Volume, HI
        # 1754 - 1843: m13-HR (2 IC/phys), O6, Ne8, Mg10
        # (240 inds)
        # 1844 - 1918: m12-SR (1 IC/phys), Mass, Volume, HI
        # 1919 - 1963: m12-SR (1 IC/phys), O6, Ne8, Mg10
        # (120 inds)
    elif index >= 1964 and index < 1970:
        # test halo centering script
        rhp.run_halodata(index - 1964)
    elif index >= 1970 and index < 2024:
        # clean samples 1/2 
        # 1970 - 1993: m13-SR (24 inds)
        # 1994 - 2005: m13-HR (12 inds)
        # 2006 - 2017: m12-HR (12 inds)
        # 2018 - 2023: m12-SR (6 inds)
        rhp.run_halodata(index - 1970 + 6)
    elif index >= 2024 and index < 2030:
        # debugging Mvir/Rvir finder: fp same Mvir values
        snaps = [186, 197, 210, 224, 240, 258]
        snapshot = snaps[index - 2024]
        path = ('/scratch3/01799/phopkins/fire3_suite_done/m12f_m7e3/'
                'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr'
                '_crdiffc690_sdp1e10_gacc31_fa0.5')
        res = hp.calchalodata_shrinkingsphere(path, snapshot, 
                                              meandef=('200c', 'BN98'))
        print(path)
        print(snapshot)
        print(res)
    elif index >= 2030 and index < 2252:
        rhp.run_halodata(index - 2030 + 60)
        # all 3model m12/m13 runs that got to z=0.5
        # from Lindsey's spreadsheet
        # 2030 - 2119: m13-SR (90 inds)
        # 2120 - 2131: m13-HR (12 inds) 
        # 2132 - 2155: m12-SR (24 inds)
        # 2156 - 2251: m12-HR (96 inds)
    elif index >= 2252 and index < 5582:
        rim.run_ionmap_xyz(index - 2252)
        # all 3model m12/m13 runs that got to z=0.5
        # from Lindsey's spreadsheet
        # z = 1.0 - 0.5, Mass, Ne8, Neon, O6, Mg10
        # 2252 - 3601: m13-SR (1350 inds)
        # 3602 - 3781: m13-HR ( 180 inds) 
        # 3782 - 4141: m12-SR ( 360 inds) 
        # 4142 - 5581: m12-HR (1440 inds)
    elif index >= 5582 and index < 5594:
        rhp.run_halodata(index - 5582 + 282)
        # m12q noBH and AGN-noCR
    elif index >= 5594 and index < 5721:
        rhp.run_halodata(index - 5582 + 282)
        # 5594 - 5603: m13-sr z=0 (10 indices)
        # 5604 - 5607: m12-sr z=0 (4 indices)
        # 5608 - 5624: m12-hr z=0 (17 indices)
        # no z=0 m13-hr runs
        # 5625 - 5720: # m11 z=0, 0.5, 1 (96 indices)
    elif index >= 5721 and index < 7156:
        rhs.run_hist_carbonions_z051(index - 5721)
        # 5721 - 6392: m11 z=0, 0.5, 1 (672 indices)
        # 6393 - 6511: m12-hr z=0      (119 indices)
        # 6512 - 6763: m12-hr z=0.5, 1 (252 indices)
        # 6764 - 6791: m12-sr z=0      (28 indices)
        # 6792 - 6847: m12-sr z=0.5, 1 (56 indices)
        # m13hr z=0: no sims
        # 6848 - 6875: m13-hr z=0.5, 1 (28 indices)
        # 6876 - 6945: m13-sr z=0      (70 indices)
        # 6946 - 7155: m13-sr z=0.5, 1 (210 indices)
    elif index >= 7156 and index < 8326:
        rhs.run_hist_ptmasses_all2(index - 7156)
        # 7156 - 7605: m13-sr (450 inds)
        # 7606 - 7665: m13-hr ( 60 inds)
        # 7666 - 7785: m12-sr (120 inds)
        # 7786 - 8325: m12-hr (540 inds)
    elif index >= 8326 and index < 8506:
        rim.run_ionmap_xyz(index - 8326 + 3330)
        # added m12q haloes (all1 -> all2)
        # z = 1.0 - 0.5, Mass, Ne8, Neon, O6, Mg10
        # (180 inds)
    elif index >= 8506 and index < 11926:
        rhs.run_hist_o6ne8mg10_tosample2(index - 8506)
        # add most of all2 sample 3D profiles
        #  8506 - 10329: m12-hr (1824 inds)
        # 10330 - 10671: m12-sr (342 inds)
        # 10672 - 11925: m13-sr (1254 inds)
    elif index >= 11926 and index < 12160:
        rhp.run_vcen1(index - 11926)
        # get all2 centers of velocity (halo and gal for each index)
        # 11926 - 12015: m13-sr (90 inds)
        # 12016 - 12027: m13-hr (12 inds)
        # 12028 - 12051: m12-sr (24 inds)
        # 12052 - 12159: m12-hr (108 inds) 
    elif index >= 12160 and index < 14968:
        # get all2 vrad/vtot for different weights
        # 12160 - 13239: m13-sr (1080 indices)
        # 13240 - 13383: m13-hr (144 indices)
        # 13384 - 13671: m12-sr (288 indices)
        # 13672 - 14967: m12-hr (1296 indices)
        rhs.run_hist_vtotrad(index - 12160)
    elif index >= 14968 and index < 16408:
        # clean2_nobug r, vr, nH/T for different weights
        # 14968 - 15447: m13-sr (480 inds)
        # 15448 - 15687: m13-hr (240 inds)
        # 15688 - 15927: m12-sr (240 inds)
        # 15928 - 16407: m12-hr (480 inds)
        rhs.run_hist_rad_vrad_weighted(index - 14968)
    elif index >= 16408 and index < 19864:
        # clean2_nobug, vlos maps with different weights
        # 16408 - 17559: m13-sr (1152 inds)
        # 17560 - 18135: m13-hr (576 inds)
        # 18136 - 18711: m12-sr (576 inds)
        # 18712 - 19863: m12-hr (1152 inds)
        rim.run_vlosmaps(opt=index - 16408)
    elif index >= 19864 and index < 20080:
        spv.run_selkin_box_clean2(index - 19864)
    elif index >= 20080 and index < 28000:
        ## supplement r_vr_T/nH clean2 -> all2
        # 20080 - 21399: m13-sr (1320 inds)
        # 21400 - 21639: m12-sr ( 240 inds)
        # 21640 - 23319: m12-hr (1680 inds)
        ## add r_vr_O/Ne all2
        # 23320 - 25119: m13-sr (1800 inds)
        # 25120 - 25359: m13-hr ( 240 inds)
        # 25360 - 25839: m12-sr ( 480 inds)
        # 25840 - 27999: m12-hr (2160 inds)
        rhs.run_hist_rad_vrad_weighted(index - 20080 + 1440)
    elif index >= 28000 and index < 28288:
        rim.run_vlosmaps(opt=index - 28000 + 3456)
        # Ne8 2 Rvir depth los maps for clean -> all m12 haloes
        # indices 28000 - 28287
    elif index >= 28288 and index < 28522:
        rhp.runcengal1(opt=index - 28288)
        # get all2 stellar center/vcom
    elif index >= 28522 and index < 28720:
        # Ne8 vlos maps m13 clean2 -> all2
        rim.run_vlosmaps(opt=index - 28522 + 3744)
    elif index >= 28720 and index < 34336:
        rim.run_vdoplosmaps(opt=index - 28720)
        # Ne8, Ne, M, V vdop maps, 3 axis proj, thin + thick slices 
        # 28720 - 30879: m13sr (2160 indices)
        # 30880 - 31167: m13hr ( 288 indices)
        # 31168 - 31743: m12sr ( 576 indices)
        # 31744 - 34335: m12hr (2592 indices)
    elif index >= 34336 and index < 37144:
        rhs.run_hist_ppv(index - 34336)
        # Ne8 ppv 3 axes
        # 34336 - 35415: m13sr (1080 indices)
        # 35416 - 35559: m13hr ( 144 indices)
        # 35560 - 35847: m12sr ( 288 indices)
        # 35848 - 37143: m12hr (1296 indices)
    elif index >= 37144 and index < 39952:
        rhs.run_hist_vdoplos_vrad(index - 37144)
        # Ne8 vdoplos - vrad 3 axes
        # 37144 - 38223: m13sr (1080 indices)
        # 38224 - 38367: m13hr ( 144 indices)
        # 38368 - 38655: m12sr ( 288 indices)
        # 38656 - 39951: m12hr (1296 indices)
    elif index >= 39952 and index < 40000:
        rhp.run_haloprop_f2md(index - 39952)
        # 48 indices
    elif index >= 40000 and index < 41152:
        rim.run_vdoplosmaps_f2md(opt=index - 40000)
        # 1152 indices
    elif index >= 41152 and index < 41728:
        # 576 indices
        rhs.run_hist_rad_vrad_weighted(index - 41152)
    elif index >= 41728 and index < 41740:
        rhp.run_haloprop_f2md(index - 41728 + 48)
        # 12 indices; adding crheatfix snaps
    elif index >= 41740 and index < 42028:
        # + 12 haloes (crheatingfix; 1152 - 1440, 288 indices)
        rim.run_vdoplosmaps_f2md(opt=index - 41740 + 1152)
    elif index >= 42028 and index < 42172:
        # + 12 haloes (crheatingfix; 144 indices, starting at 576)
        rhs.run_hist_rad_vrad_weighted(opt=index - 42048 + 576)
    else:
        raise ValueError('Nothing specified for index {}'.format(index))

def launchergen(*args, logfilebase='{ind}.out'):
    '''
    not that useful; just use 
    echo -e "commands ${index} >> logfile_${index}" >> launchfile_name 
    in the batch script

    Parameters:
    -----------
    args: indexable of integers
        the indices to call fromcommandline with, one for each launched
        process
    logfilebase: string, formattable with argument 'ind'
        where to write the logfiles. {ind} is replaced by the index in each 
        line
    Returns:
    --------
    prints the launcher file lines. Direct output to a file to generate one.
    '''
    
    fillline = 'python ./fire_maps.py {ind} > ' + logfilebase + ' 2>&1'
    for arg in args:
        print(fillline.format(ind=arg))

if __name__ == '__main__':
    #print('fire_maps.py script started')
    if len(sys.argv) > 1:
        # generate launcher file for frontera
        # arguments: 
        #   --launchergen : generate a launcher file instead of 
        #                   default 'run with this index'
        #   --logfilebase=<string> : write log files for each launcher 
        #                   process to a file like this. Must contain a
        #                   '{ind}' part, since this is where the script
        #                   will fill in the index each process is called 
        #                   with
        #   integers :      the indices to call this script (fire_maps.py)
        #                   with in the launcher run
        if '--launchergen' in sys.argv:
            inds = [int(arg) if '-' not in arg else None \
                    for arg in sys.argv[1:]]
            while None in inds:
                inds.remove(None)
            kw = {}
            for arg in sys.argv[1:]:
                if '--logfilebase=' in arg:
                    kw['logfilebase'] = arg.split('=')[-1]
                    break
            launchergen(*inds, **kw)
        # just run code
        else:
            print('fire_maps.py script started')
            try:
                ind = int(sys.argv[1])
            except ValueError as msg1:
                msg2 = 'Could not interpret first command-line' + \
                       ' argument {} as int'
                msg2 = msg2.format(sys.argv[1])
                raise ValueError('/n'.join([msg1, msg2]))
            fromcommandline(ind)
    else:
        raise ValueError('Please specify an integer index > 1')
    
    








