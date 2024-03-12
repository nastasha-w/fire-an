import numpy as np
import os

import fire_an.explore.clumpiness as cpn
import fire_an.simlists as sl
import fire_an.utils.opts_locs as ol

def run_clumpiness(opt):
    outdir = ol.pre + 'output/clumps/'
    if opt >= 0 and opt < 360:
        # 360 indices
        ind = opt
        simnames = sl.m13_sr_all2 # len 15
        snapnums = sl.snaps_sr
    elif opt >= 360 and opt < 408:
        # 48 indices
        ind = opt - 360
        simnames = sl.m13_hr_all2 # len 2
        snapnums = sl.snaps_hr
    if opt >= 408 and opt < 504:
        # 96 indices
        ind = opt - 408
        simnames = sl.m12_sr_all2 # len 4
        snapnums = sl.snaps_sr
    elif opt >= 504 and opt < 936:
        # 432 indices
        ind = opt - 504
        simnames = sl.m12_hr_all2 # len 18
        snapnums = sl.snaps_hr
    elif opt >= 936 and opt < 1176:
        # 240 indices
        ind = opt - 936
        simnames = sl.m12_f2md # len 10
        snapnums = sl.snaps_f2md
    elif opt >= 1176 and opt < 1248:
        outdir = ol.pre + 'clumps/'
        # 72 indices
        ind = opt - 1176
        simnames = sl.m12plus_f3nobh # len 3 (no m12g)
        snapnums = sl.snaps_hr 
    elif opt >= 1248 and opt < 1464:
        outdir = ol.pre + 'clumps/'
        # 216 indices
        ind = opt - 1248
        simnames = sl.m12plus_f3nobh_lores # len 9
        snapnums = sl.snaps_hr # len 6
    elif opt >= 1464 and opt < 1560:
        ind = opt - 1464
        # 24 halos, 96 indices
        simnames = sl.m12_fire3x_tests
        snapnums = sl.snaps_sr
    mts = ['Ne8num_Ne8dens', 'Vol_Ne8dens', 'Mass_dens', '']
    #mts = ['Ne8num_Vol', 'Ne8num_dens', 'Vol_dens', 'Ne8dens_Ne8dens']

    simi = ind // (len(snapnums) * len(mts))
    snapi = (ind % (len(snapnums) * len(mts))) // len(mts)
    mti = ind % len(mts)
    simname = simnames[simi]
    snapnum = snapnums[snapi]
    mt = mts[mti]
    
    dirpath = sl.dirpath_from_simname(simname)
    parttype = 0
    rbins = np.arange(0., 1.32, 0.05)
    outname = f'clumpines_measure_v1_{mt}_{simname}_{snapnum}'
    outname = outdir + outname + '.hdf5'
    
    if mt == 'Ne8num_Vol':
        maptype1 = 'Volume'
        maptype_args1 = {}
        maptype2 = 'ion'
        maptype_args2 = {'ion': 'Ne8', 'ps20depletion': False, 
                        'density': False, 'lintable': True}
    elif mt == 'Ne8num_dens':
        maptype1 = 'sim-direct'
        maptype_args1 = {'field': 'Density'}
        maptype2 = 'ion'
        maptype_args2 = {'ion': 'Ne8', 'ps20depletion': False, 
                        'density': False, 'lintable': True}
    elif mt == 'Vol_dens':
        maptype1 = 'sim-direct'
        maptype_args1 = {'field': 'Density'}
        maptype2 = 'Volume'
        maptype_args2 = {}
    elif mt == 'Ne8dens_Ne8dens':
        maptype1 = 'ion'
        maptype_args1 = {'ion': 'Ne8', 'ps20depletion': False, 
                        'density': True, 'lintable': True}
        maptype2 = 'ion'
        maptype_args2 = {'ion': 'Ne8', 'ps20depletion': False, 
                        'density': True, 'lintable': True}
    elif mt == 'Ne8num_Ne8dens':
        maptype1 = 'ion'
        maptype_args1 = {'ion': 'Ne8', 'ps20depletion': False, 
                        'density': False, 'lintable': True}
        maptype2 = 'ion'
        maptype_args2 = {'ion': 'Ne8', 'ps20depletion': False, 
                        'density': True, 'lintable': True}
    elif mt == 'Vol_Ne8dens':
        maptype1 = 'Volume'
        maptype_args1 = {}
        maptype2 = 'ion'
        maptype_args2 = {'ion': 'Ne8', 'ps20depletion': False, 
                        'density': True, 'lintable': True}
    elif mt == 'Mass_dens':
        maptype1 = 'sim-direct'
        maptype_args1 = {'field': 'Density'}
        maptype2 = 'Mass'
        maptype_args2 = {}
    elif mt == '':
        return None

    if os.path.isfile(outname):
        print(f'File already exists: {outname}; skipping')
        return None
    cpn.wtdavnorm(dirpath, snapnum, parttype, 
              maptype1, maptype_args1, maptype2, maptype_args2,
              rbins=rbins, rbins_unit='Rvir',
              savefile=outname)