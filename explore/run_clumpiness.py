import numpy as np
import os

import fire_an.explore.clumpiness as cpn
import fire_an.simlists as sl
import fire_an.utils.opts_locs as ol

def run_clumpiness(opt):
    if opt >= 0 and opt < 360:
        # 360 indices
        simnames = sl.m13_sr_all2 # len 15
        snapnums = sl.snaps_sr
    elif opt >= 360 and opt < 408:
        # 48 indices
        simnames = sl.m13_hr_all2 # len 2
        snapnums = sl.snaps_hr
    if opt >= 408 and opt < 504:
        # 96 indices
        simnames = sl.m12_sr_all2 # len 4
        snapnums = sl.snaps_sr
    elif opt >= 504 and opt < 936:
        # 432 indices
        simnames = sl.m12_hr_all2 # len 18
        snapnums = sl.snaps_hr
    elif opt >= 936 and opt < 1176:
        # 240 indices
        simnames = sl.m12_f2md # len 10
        snapnums = sl.snaps_f2md
    mts = ['Ne8num_Vol', 'Ne8num_dens', 'Vol_dens', 'Ne8dens_Ne8dens']

    simi = opt // (len(snapnums) * len(mts))
    snapi = (opt % (len(snapnums) * len(mts))) // len(mts)
    mti = opt % len(mts)
    simname = simnames[simi]
    snapnum = snapnums[snapi]
    mt = mts[mti]
    
    dirpath = sl.dirpath_from_simname(simname)
    parttype = 0
    rbins = np.arange(0., 1.32, 0.05)
    outdir = ol.pre + 'output/clumps/'
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

    if os.path.isfile(outname):
        print(f'File already exists: {outname}; skipping')
        return None
    cpn.wtdavnorm(dirpath, snapnum, parttype, 
              maptype1, maptype_args1, maptype2, maptype_args2,
              rbins=rbins, rbins_unit='Rvir',
              savefile=outname)