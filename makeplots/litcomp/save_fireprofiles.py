'''
Save Ne8 data for the observation comparison plots
'''

import h5py
import numpy as np

import fire_an.makeplots.get_2dprof as gpr
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu
import fire_an.utils.opts_locs as ol

savedir = '/projects/b1026/nastasha/plotdata/'
filen_save = savedir + 'coldens_radprof_Ne8_opt3.hdf5'
# opt1: buglist1, no m12plus haloes
# opt2: buglist2, hi, lo-res m12plus combined
# opt3: buglist2, includes FIRE-3.x test halos


def saveprof_single(h5grp, simnames, snapnums):
    profdir = '/projects/b1026/nastasha/maps/vdopmaps_all2/'
    filen_temp = ('vdoplos_by_coldens_Ne8_{simname}_snap{snapnum}'
                  '_shrink-sph-cen_BN98_depth_2.0rvir_{pax}-proj_v3.hdf5')
    otherfills = [{'pax': 'x'}, {'pax': 'y'}, {'pax': 'z'}]
    profiles = ['perc-0.1', 'perc-0.5', 'perc-0.9',
                'perc-0.2', 'perc-0.05', 'perc-0.02', 'perc-0.01',
                'perc-0.8', 'perc-0.95', 'perc-0.98', 'perc-0.99',
                'av-lin']
    rbins_pkpc = np.linspace(0., 450., 50)

    filens = [profdir + filen_temp.format(simname=simn, snapnum=snap, **ofill)
              for simn, snap in zip(simnames, snapnums)
              for ofill in otherfills]
    
    outprof = gpr.get_profile_massmap(filens, rbins_pkpc,
                                      rbin_units='pkpc',
                                      profiles=profiles,
                                      weightmap=True,
                                      absvals=False,
                                      weightrange=None)
    for proftype, profile in zip(profiles, outprof):
        h5grp.create_dataset(proftype, data=profile)
    h5grp.create_dataset('rbins', data=rbins_pkpc)
    h5grp.attrs.create('rbin_units', np.string_('pkpc'))
    h5grp.attrs.create('weightmap', True)
    h5grp.attrs.create('absvals', False)
    h5grp.create_dataset('mapfiles', 
                         data=np.array([np.string_(fn) for fn in filens]))
    h5grp.create_dataset('simnames', 
                         data=np.array([np.string_(sn) for sn in simnames]))
    h5grp.create_dataset('snapnums', 
                         data=np.array(snapnums))

def calcsave_coldensprof():
    # save profiles for:
    # individual sim/snap, 
    # all snaps for one sim (z = 0.5--1.0, 0.5--0.7)
    # all snaps/sims for one physics model and halo mass class
    # (same snaps categories)
    sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    sims_m12plus = sl.m12plus_f3nobh + sl.m12plus_f3nobh_lores
    sims_f2md = sl.m12_f2md
    sims_f3x = sl.m12_fire3x_tests
    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    snaps_f2md = sl.snaps_f2md

    snapsels = {'z0.5-0.7': slice(3, 6, None),
                'z0.5-1.0': slice(None, None, None)}
    sims_all = sims_sr + sims_hr + sims_f2md + sims_m12plus + sims_f3x
    simcats = {}
    for simn in sims_all:
        if simn in sl.buglist2:
            continue
        ic = sl.ic_from_simname(simn)
        masscat = ic[:3] # m12 or m13
        physlabel = sl.physlabel_from_simname(simn)
        catkey = masscat + '_' + physlabel
        if catkey in simcats:
            simcats[catkey].append(simn)
        else:
            simcats[catkey] = [simn]
    
    with h5py.File(filen_save, 'a') as fs:
        for sim in sims_all:
            snaps = (snaps_hr if sim in sims_hr
                     else snaps_sr if sim in sims_sr
                     else snaps_f2md if sim in sims_f2md
                     else snaps_hr if sim in sims_m12plus
                     else snaps_sr if sim in sims_f3x
                     else None) 
            for snap in snaps:
                grpname = f'{sim}_snap{snap}'
                grp = fs.create_group(grpname)
                saveprof_single(grp, [sim], [snap])
            for key in snapsels:
                grpname = f'{sim}_{key}'
                grp = fs.create_group(grpname)
                _snaps = snaps[snapsels[key]]
                saveprof_single(grp, [sim] * len(_snaps), _snaps)
        for simkey in simcats:
            sims = simcats[simkey]
            for snapkey in snapsels:
                grpname = f'{simkey}_{snapkey}'
                grp = fs.create_group(grpname)
                snapsel = snapsels[snapkey]
                snaps = [(snaps_hr[snapsel] if sim in sims_hr
                          else snaps_sr[snapsel] if sim in sims_sr
                          else snaps_f2md[snapsel] if sim in sims_f2md
                          else snaps_hr[snapsel] if sim in sims_m12plus
                          else snaps_sr[snapsel] if sim in sims_f3x
                          else None)
                          for sim in sims]
                sims_in = [sim for _snaps, sim in zip(snaps, sims) 
                           for _ in _snaps]
                snaps_in = [snap for l1 in snaps for snap in l1]
                saveprof_single(grp, sims_in, snaps_in)



    

