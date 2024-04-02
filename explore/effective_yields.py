'''
Calculate total stellar masses and metal masses in
zoom simulation volumes
'''

import h5py
import numpy as np

import fire_an.mainfunc.get_qty as gq
import fire_an.readfire.readin_fire_data as rfd
import fire_an.simlists as sl
import fire_an.utils.h5utils as h5u

def get_totals(simname, snapnum, outname):
    metals = ['total', 'Oxygen', 'Neon', 'Carbon', 'Nitrogen', 'Iron',
              'Magnesium', 'Sulfur']
    maptypes = ['Mass'] + ['Metal'] * len(metals)
    maptype_argss = [{}] + [{'element': val, 'density': False}
                            for val in metals]

    dirpath = sl.dirpath_from_simname(simname)
    snap = rfd.get_Firesnap(dirpath, snapnum)
    with h5py.File(outname, 'a') as f:
        hed = f.create_group('Header')
        hed.attrs.create('simname', np.string_(simname))
        hed.attrs.create('dirpath', np.string_(dirpath))
        hed.attrs.create('snapnum', snapnum)
        cosmopars = snap.cosmopars.getdct()
        csm = hed.create_group('cosmopars')
        h5u.savedict_hdf5(cosmopars, csm)
        f.create_group('gas')
        f.create_group('stars')

        for parttype, ptlabel in [(0, 'gas'), (4, 'stars')]:
            for maptype, maptype_args in zip(maptypes, maptype_argss):
                qty, toCGS, todoc = gq.get_qty(snap, parttype, 
                                               maptype, maptype_args,
                                               filterdct=None)
                tot = np.sum(qty)

                grpname = ('Mass' if maptype == 'Mass' 
                           else 'MetalMass_' + maptype_args['element'])
                cur = f[ptlabel].create_group(grpname)
                cur.attrs.create('total mass', tot)
                cur.attrs.create('toCGS', toCGS)
                cdoc = cur.create_group('doc')
                h5u.savedict_hdf5(cdoc, todoc)

def run_totals(index):
    outdir = '/projects/b1026/nastasha/hists/gas_stars_metals/'
    # leaving out the fire3_m12plus halos for npw
    if index >= 0 and index < 108:
        ind = index - 0
        simnames = sl.m12_hr_all2 # 18
        snapnums = sl.snaps_hr # 6
    elif index >= 108 and index < 132:
        ind = index - 108
        simnames = sl.m12_sr_all2 # 4
        snapnums = sl.snaps_sr # 6 
    elif index >= 132 and index < 144:
        ind = index - 132
        simnames = sl.m13_hr_all2 # 2
        snapnums = sl.snaps_hr # 6
    elif index >= 144 and index < 234:
        ind = index - 144
        simnames = sl.m13_sr_all2 # 15
        snapnums = sl.snaps_sr # 6
    elif index >= 234 and index < 294:
        ind = index - 234
        simnames = sl.m12_f2md # 10
        snapnums = sl.snaps_f2md #6
    elif index >= 294 and index < 318:
        ind = index - 294
        simnames = sl.m12_fire3x_tests # 4
        snapnums = sl.snaps_sr #6
    
    nmi = ind // len(snapnums)
    sni = ind % len(snapnums)
    simname = simnames[nmi]
    snapnum = snapnums[sni]
    
    outname = outdir + f'total_MassZ_stars_gas_{simname}_{snapnum}.hdf5'

    get_totals(simname, snapnum, outname)
            

