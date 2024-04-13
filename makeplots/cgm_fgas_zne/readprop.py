import pandas as pd
import numpy as np

import fire_an.mainfunc.cengalprop as cgp
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c

ddir = '/projects/b1026/nastasha/hists/r_vr_all2/'
totfilen = ddir + 'gas_Neon_Ne8_masses_rTcuts_v2.dat'
## just get from total masses; mostly got those to check totals anyway
#avfilen = ddir + 'mean_ZNe_by_mass_volume_rcuts.dat'
totZfilen = ddir + 'mean_Ztot_by_mass_rcuts_tcuts_v2.dat'
mstarfilen = ddir + 'stellarmass_rcuts_v2.dat'

def getbasedata(massset, inclm12plus=False, f3xset=False,
                f3xset_large=False):
    '''
    set up a dataframe with simnames, snapnums for a mass set 
    ('m12' or 'm13'), excluding runs with bugs
    '''
    snaps_f2 = sl.snaps_f2md
    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    snaps_m12plus = sl.snaps_hr
    snaps_f3x = snaps_sr
    sims_f2 = sl.m12_f2md
    sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    sims_m12plus = sl.m12plus_f3nobh + sl.m12plus_f3nobh_lores
    sims_f3x = sl.m12_fire3x_tests
    
    if f3xset:
        simnames_all = sims_f2 + sims_f3x \
                       + sl.m12_nobh_clean2 + sl.m12_nobh_rest2
    elif f3xset_large:
        simnames_all = sims_f2 + sims_f3x + sims_sr + sims_hr 
    else:
        simnames_all = sims_f2 + sims_sr + sims_hr 
        if inclm12plus:
            simnames_all = simnames_all + sims_m12plus
    simnames_all = [sn for sn in simnames_all if sn not in sl.buglist2]
    snapnums = [snaps_f2 if sn in sims_f2 else
                snaps_sr if sn in sims_sr else
                snaps_hr if sn in sims_hr else
                snaps_m12plus if sn in sims_m12plus else
                snaps_f3x if sn in sims_f3x else
                None
                for sn in simnames_all]
    snapnums = [snap for l1 in snapnums for snap in l1]
    simnames_all = [sn for sn in simnames_all for i in range(6)]

    data = pd.DataFrame({'simname': simnames_all}) #, dtype='string'
    data['snapnum'] = snapnums

    data['ic'] = np.array([sl.ic_from_simname(simname)
                           for simname in data['simname']])
    filter2 = np.array([ic.startswith(massset) for ic in data['ic']])
    data = data[filter2]
    data['physmodel'] = np.array([sl.physlabel_from_simname(simname)
                                  for simname in data['simname']])
    simpaths = [sl.dirpath_from_simname(simname) 
                for simname in data['simname']]
    _dcts = np.array([cgp.getcengalcen(_simp, _snap)[2]
                      for _simp, _snap in zip(simpaths, data['snapnum'])])
    data['Mstarcen_g'] = np.array([_dct['mstar_gal_g'] for _dct in _dcts])
    data['Mvir_g'] = np.array([_dct['halodata']['Mvir_g'] for _dct in _dcts])
    data['Rvir_cm'] = np.array([_dct['halodata']['Rvir_cm']
                                for _dct in _dcts])
    data['redshift'] = np.array([_dct['halodata_doc']['cosmopars']['z']
                                 for _dct in _dcts])
    data['Omega_b'] = np.array([_dct['halodata_doc']['cosmopars']['omegab']
                                for _dct in _dcts])
    data['Omega_m'] = np.array([_dct['halodata_doc']['cosmopars']['omegam']
                                for _dct in _dcts])
    data = data.set_index(['simname', 'snapnum'], drop=True)
    return data

def readin_all_data(rrange_rvir=(0.1, 1.0), trange_logk=(-np.inf, np.inf),
                    massset='m12', inclm12plus=False, f3xset=False,
                    f3xset_large=False):
    '''
    reads useful info for each simname, snapnum in the massset
    note: stellar mass selection is independent of the gas temperature 
          range
    '''
    data = getbasedata(massset, inclm12plus=inclm12plus, f3xset=f3xset,
                       f3xset_large=f3xset_large)

    # 1st pass total mass: gas mass, halo gas fraction
    _data = pd.read_csv(totfilen, sep='\t')
    filter = np.isclose(_data['rmin_rvir'], rrange_rvir[0])
    filter &= np.isclose(_data['rmax_rvir'], rrange_rvir[1])
    filter &= np.isclose(_data['tmin_logk'], trange_logk[0])
    filter &= np.isclose(_data['tmax_logk'], trange_logk[1])
    filter &= _data['weight'] == 'gasmass'
    _data = _data[filter]
    _data['fgas'] = _data['total [g or num. part.]'] \
                   / (_data['Mvir_g'] * _data['Omega_b'] / _data['Omega_m'])
    _data['gasmass_g'] = _data['total [g or num. part.]']
    _data = _data.set_index(['simname', 'snapnum'], drop=True)
    data['fgas'] = _data['fgas']
    data['gasmass_g'] = _data['gasmass_g']

    # 2nd pass total mass: Ne
    _data = pd.read_csv(totfilen, sep='\t')
    filter = np.isclose(_data['rmin_rvir'], rrange_rvir[0])
    filter &= np.isclose(_data['rmax_rvir'], rrange_rvir[1])
    filter &= np.isclose(_data['tmin_logk'], trange_logk[0])
    filter &= np.isclose(_data['tmax_logk'], trange_logk[1])
    filter &= _data['weight'] == 'Neon'
    _data = _data[filter]
    _data['Neon_numpart'] = _data['total [g or num. part.]']
    _data = _data.set_index(['simname', 'snapnum'], drop=True)
    data['Neon_numpart'] = _data['Neon_numpart']

    # 3rd pass total mass: Ne8
    _data = pd.read_csv(totfilen, sep='\t')
    filter = np.isclose(_data['rmin_rvir'], rrange_rvir[0])
    filter &= np.isclose(_data['rmax_rvir'], rrange_rvir[1])
    filter &= np.isclose(_data['tmin_logk'], trange_logk[0])
    filter &= np.isclose(_data['tmax_logk'], trange_logk[1])
    filter &= _data['weight'] == 'Ne8'
    _data = _data[filter]
    _data['Ne8_numpart'] = _data['total [g or num. part.]']
    _data = _data.set_index(['simname', 'snapnum'], drop=True)
    data['Ne8_numpart'] = _data['Ne8_numpart']

    # gas metallicity (total metals)
    _data = pd.read_csv(totZfilen, sep='\t')
    filter = np.isclose(_data['rmin_rvir'], rrange_rvir[0])
    filter &= np.isclose(_data['rmax_rvir'], rrange_rvir[1])
    filter &= np.isclose(_data['tmin_logk'], trange_logk[0])
    filter &= np.isclose(_data['tmax_logk'], trange_logk[1])
    _data = _data[filter]
    _data = _data.set_index(['simname', 'snapnum'], drop=True)
    data['Ztot gasmass-wtd (mass fraction)'] = \
        _data['wtd. mean total metallicity (mass fraction)']
    
    # stellar mass (in the halo)
    _data = pd.read_csv(mstarfilen, sep='\t')
    filter = np.isclose(_data['rmin_rvir'], rrange_rvir[0])
    filter &= np.isclose(_data['rmax_rvir'], rrange_rvir[1])
    _data = _data[filter]
    _data = _data.set_index(['simname', 'snapnum'], drop=True)
    data['Mstar_current_g'] = \
        _data['total current stellar mass [g]']

    #data = data.reset_index() # simname, snapnum back to regular columns
    return data