


import numpy as np

from fire_an.ionrad.ion_utils import Linetable_PS20
import fire_an.mainfunc.makehist as mh
import fire_an.simlists as sl

def run_hist_carbonions_z051(opt):
    if opt >= 0 and opt < 189:
        # z = 0.0, 0.5, 1.5
        # m11 haloes pt 1
        # 189 indices
        ind = opt - 0
        simnames = sl.m11_hr_set1 # len 9
        snaps = sl.snaps_hr_051 # len 3
    elif opt >= 189 and opt < 672:
        # z = 0.0, 0.5, 1.5
        # m11 haloes pt 2
        # 483 indices
        ind = opt - 189
        simnames = sl.m11_sr_set1 # len 23
        snaps = sl.snaps_sr_051 # len 3
    elif opt >= 672 and opt < 791:
        # m12-hr z=0
        # 119 indices
        ind = opt - 672
        simnames = sl.m12_hr_all2_z0 # len 17
        snaps = [500] # len 1
    elif opt >= 791 and opt < 1043:
        # m12-hr z=1, 0.5
        # 252 indices
        ind = opt - 791
        simnames = sl.m12_hr_all2 # len 18
        snaps = [186, 258] # len 2
    elif opt >= 1043 and opt < 1071:
        # m12-sr z=0
        # 28 indices
        ind = opt - 1043
        simnames = sl.m12_sr_all2_z0 # len 4
        snaps = [60] # len 1
    elif opt >= 1071 and opt < 1127:
        # m12-sr z=1, 0.5
        # 56 indices
        ind = opt - 1071
        simnames = sl.m12_sr_all2 # len 4
        snaps = [45, 50] # len 1
    # no m13-hr to z=0
    elif opt >= 1127 and opt < 1155:
        # m13-hr z=1, 0.5
        # 28 indices
        ind = opt - 1127
        simnames = sl.m13_hr_all2 # len 2
        snaps = [186, 258] # len 2
    elif opt >= 1155 and opt < 1225:
        # m13-sr z=0
        # 70 indices
        ind = opt - 1155
        simnames = sl.m13_sr_all2_z0 # len 10
        snaps = [60] # len 1
    elif opt >= 1225 and opt < 1435:
        # m13-sr z=1, 0.5
        # 210 indices
        ind = opt - 1225
        simnames = sl.m13_sr_all2 # len 15
        snaps = [45, 50] # len 2

    wts = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'Carbon'] # len 7
    _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
    outdir = '/scratch1/08466/tg877653/output/hists/ionseries_C/'
    outname = 'hist_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
    particle_type = 0
    simi = ind // (len(snaps) * len(wts))
    snpi = (ind % (len(snaps) * len(wts))) // (len(wts))
    wti = (ind % len(wts))
    simname = simnames[simi]
    snapnum = snaps[snpi]
    wt = wts[wti]
    axtypes = []
    axtypes_args = []
    axbins = []
    
    runit = 'Rvir'
    rbins = np.arange(0.15, 4., 0.01)
    rbins = np.append(np.arange(0., 0.15, 0.005), rbins)

    # directory is halo name + resolution 
    dp2 = '_'.join(simname.split('_')[:2])
    if dp2.startswith('m13h02_'):
        dp2 = dp2.replace('m13h02', 'm13h002')
    dirpath = '/'.join([_dirpath, dp2, simname])

    if wt == 'Carbon':
        weighttype = 'Metal'
        weighttype_args = {'element': wt, 'density': False}
    else:
        weighttype = 'ion'
        weighttype_args = {'ps20depletion': False, 'ion': wt,
                            'density': False}
        if wt == 'H1':
            weighttype_args.update({'ionfrac-method': 'sim'})
    outfilen = outdir + outname.format(wt=wt, simname=simname, 
                                        snap=snapnum)
    mh.histogram_radprof(dirpath, snapnum,
                         weighttype, weighttype_args, axtypes, axtypes_args,
                         particle_type=particle_type, 
                         center='shrinksph', rbins=rbins, runit=runit,
                         logweights=True, logaxes=True, axbins=axbins,
                         outfilen=outfilen, overwrite=False)
    
def run_hist(opt):
    if opt >= 0 and opt < 60:
        ind = opt
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [50] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                   ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 60 and opt < 96:
        ind = opt - 60
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [50] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                   ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)

        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 96 and opt < 126:
        ind = opt - 96
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [258] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
                    'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
                    ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 126 and opt < 144:
        ind = opt - 126
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [258] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
                    'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
                    ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)

        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 144 and opt < 174:
        ind = opt - 144
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [258] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31_fa0.5',
                   ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 174 and opt < 192:
        ind = opt - 174
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [258] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31_fa0.5',
                   ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    
    elif opt >= 192 and opt < 207:
        ind = opt - 192
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [50] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 207 and opt < 216:
        ind = opt - 206
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [50] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    
    # add more redshifts (set 2) to previous sample
    #----------------------------------------------
    elif opt >= 216 and opt < 516:
        ind = opt - 216
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [45, 46, 47, 48, 49] # z=1.0 - 0.6
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                   ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 516 and opt < 696:
        ind = opt - 516
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [45, 46, 47, 48, 49] # z=1.0 - 0.6
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                    'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    'm13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                   ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)

        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 696 and opt < 846:
        ind = opt - 696
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [186, 197, 210, 224, 240] # z=1.0 - 0.6
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
                    'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
                    ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 846 and opt < 936:
        ind = opt - 846
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [186, 197, 210, 224, 240] # z=1.0 - 0.6
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
                    'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
                    ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)

        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 936 and opt < 1086:
        ind = opt - 936
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [186, 197, 210, 224, 240] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31_fa0.5',
                   ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 1086 and opt < 1176:
        ind = opt - 1086
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [186, 197, 210, 224, 240] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
                    'm13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31_fa0.5',
                   ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    
    elif opt >= 1176 and opt < 1251:
        ind = opt - 1176
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['Mass', 'Volume', 'H1']
        snaps = [45, 46, 47, 48, 49] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    ]
        axtypes_opts = [['sim-direct']] * 5
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/Oxygen'}],
                             [{'field': 'ElementAbundance/Neon'}],
                             [{'field': 'ElementAbundance/Magnesium'}],
                            ]
        axqts = ['Temperature', 'Density', 'Oxygen', 'Neon', 'Magnesium']
        axbins = [0.05, 0.05] + [0.1] * 3

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}\
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
    elif opt >= 1251 and opt < 1296:
        ind = opt - 1251
        outdir = '/scratch1/08466/tg877653/output/hists/clean_set1_set2/'
        outname = 'hist_{axqt}_r3D_by_{wt}_{simname}_snap{snap}_bins1_v1.hdf5'
        particle_type = 0
        wts = ['O6', 'Ne8', 'Mg10']
        snaps = [45, 46, 47, 48, 49] # z=0.5
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    ]
        axtypes_opts = [['sim-direct']] * 3
        axtypes_args_opts = [[{'field': 'Temperature'}],
                             [{'field': 'Density'}],
                             [{'field': 'ElementAbundance/{elt}'}],
                            ]
        axqts = ['Temperature', 'Density', '{elt}']
        axbins = [0.05, 0.05, 0.1]

        simi = ind // (len(snaps) * len(wts) * len(axqts))
        snpi = (ind % (len(snaps) * len(wts) * len(axqts))) \
               // (len(wts) * len(axqts))
        wti = (ind % (len(wts) * len(axqts))) // len(axqts)
        axi = ind % len(axqts)
        simname = simnames[simi]
        snapnum = snaps[snpi]
        wt = wts[wti]
        axtypes = axtypes_opts[axi]
        axtypes_args = axtypes_args_opts[axi]
        axqt = axqts[axi]
        
        runit = 'pkpc'
        rbins = np.arange(40., 501., 20.) if simname.startswith('m12') else\
                np.arange(40., 1001., 20.)
        rbins = np.append(np.arange(0., 40., 5.), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        if wt in ['Mass', 'Volume']:
            weighttype = wt
            weighttype_args = dict()
        else:
            weighttype = 'ion'
            weighttype_args = {'ps20depletion': False, 'ion': wt,
                               'density': False}
            if wt == 'H1':
                weighttype_args.update({'ionfrac-method': 'sim'})
            else:
                dummytab = Linetable_PS20(wt, 0.0, emission=False,
                                          vol=True, lintable=True)
                parentelt = dummytab.element
                axtypes_args = \
                    [{key: (dct[key]).format(elt=parentelt) for key in dct}
                     for dct in axtypes_args]
                axqt = axqt.format(elt=parentelt)
        outfilen = outdir + outname.format(axqt=axqt, wt=wt, simname=simname, 
                                           snap=snapnum)
        

        
    elif opt >= 123456789 and opt < 123467890:
        ind = opt - 1296
        outdir = '/scratch1/08466/tg877653/output/hists/massprof/'
        outname = 'hist_r3D_by_mass-pt{pt}_{simname}_snap{snap}_bins1_v1.hdf5'
        pts = [0, 1, 2, 4, 5]
        wt = 'Mass'
        snaps = []
        _dirpath = '/scratch3/01799/phopkins/fire3_suite_done/'
        simnames = ['m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000',
                    ]

        simi = ind // (len(snaps) * len(pts))
        snpi = (ind % (len(snaps) * len(pts))) // (len(pts))
        pti = (ind % len(pts))
        simname = simnames[simi]
        snapnum = snaps[snpi]
        particle_type = pts[pti]
        axtypes = []
        axtypes_args = []
        axqt = []
        axbins = []
        
        # no black holes in some simulations
        if particle_type == 5 and 'sdp1e10' in simname:
            msg = (f'Skipping particle type {particle_type} for noBH'
                   f' simulation {simname}')
            print(msg)
            return None
        
        runit = 'Rvir'
        rbins = np.arange(0.15, 4., 0.01)
        rbins = np.append(np.arange(0., 0.15, 0.005), rbins)

        # directory is halo name + resolution 
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([_dirpath, dp2, simname])

        weighttype = wt
        weighttype_args = {}
        outfilen = outdir + outname.format(pt=particle_type, simname=simname,
                                           snap=snapnum)

    mh.histogram_radprof(dirpath, snapnum,
                         weighttype, weighttype_args, axtypes, axtypes_args,
                         particle_type=particle_type, 
                         center='shrinksph', rbins=rbins, runit=runit,
                         logweights=True, logaxes=True, axbins=axbins,
                         outfilen=outfilen, overwrite=False)
