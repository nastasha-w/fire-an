'''
Calculate and store data: 
for each sim (IC/phys) and snapshot, based on previous histograms
- gas mass (for fgas): 3 radial sections (0.15-0.25, 0.45-0.55,
                                          0.9-1.0 Rvir),
                       CGM, halo, ISM:  0.1 - 1, 0 - 1, 0.1 - 1 Rvir
                       and for all T, T > 5, T > 5.5
- gas metallicity: 3 radial sections, CGM, halo, ISM
                   and for all T, T > 5, T > 5.5
                   mass- and volume-weighted
- store halo mass from haloprop calc

format: csv files
line for each simname, snapnum (+ store phys. label and IC)
files for gas mass, stellar mass, Z
different lines (can group later) for radial ranges, T ranges, 
                                      weighting
store cosmopars for each line too (Omega_b / Omega_m normalization)

_v2 files: includes FIRE-3.x data

'''
import h5py
import numpy as np

import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.math_utils as mu

datadir_rvcen = '/projects/b1026/nastasha/hists/r_vr_all2/'
filetemp_rvcen = ('hist_rcen_vcen_{axis2}_by_{weight}_{simname}'
                  '_snap{snapnum}_bins1_v1_hvcen.hdf5')
datadir_rcen = '/projects/b1026/nastasha/hists/r_wtd/'
filetemp_rcen_mstar = ('hist_rcen__by_stellarmass_{simname}_snap{snapnum}'
                       '_bins1_v1_hvcen.hdf5')
filetemp_rcen_Ztot = ('hist_rcen_Metallicity_Temperature_by_gasmass_'
                      '{simname}_snap{snapnum}_bins1_v1_hvcen.hdf5')

simnames_sr = sl.m12_sr_all2 + sl.m13_sr_all2
simnames_hr = sl.m12_hr_all2 + sl.m13_hr_all2
simnames_f2 = sl.m12_f2md
simnames_m12plus = sl.m12plus_f3nobh + sl.m12plus_f3nobh_lores
simnames_f3x = sl.m12_fire3x_tests
snaps_sr = sl.snaps_sr
snaps_hr = sl.snaps_hr
snaps_f2 = sl.snaps_f2md
snaps_m12plus = snaps_hr
snaps_f3x = snaps_sr
simnames_all = simnames_sr + simnames_hr + simnames_f2 + simnames_m12plus \
               + simnames_f3x

tranges_logk = [(-np.inf, np.inf), (5., np.inf), (5.5, np.inf)]
rsels_cat = [(0.1, 1.0), (0.0, 1.0), (0.0, 0.1)]
rsels_rad = [(0.15, 0.25), (0.45, 0.55), (0.90, 1.00)]

def getsnaps(simname):
    if simname in simnames_sr:
        return snaps_sr
    elif simname in simnames_hr:
        return snaps_hr
    elif simname in simnames_f2:
        return snaps_f2
    elif simname in simnames_m12plus:
        return snaps_m12plus
    elif simname in simnames_f3x:
        return snaps_f3x
    else:
        raise ValueError(f'No snapshots lsited for {simname}') 


def getmass_rbins(simname, snapnum, rrange_rvir=(0.1, 1.),
                  trange_logk=None, weight='gasmass'):
    axis2 = 'temperature'
    filen = datadir_rvcen + filetemp_rvcen.format(simname=simname, 
                                                  snapnum=snapnum,
                                                  weight=weight,
                                                  axis2=axis2)
    
    with h5py.File(filen, 'r') as f:
        hist = f['histogram/histogram'][:]
        islog = bool(f['histogram'].attrs['log'])
        if islog:
            hist = 10**hist
        mvir_g = f['Header/inputpars/halodata'].attrs['Mvir_g']
        #vescvir_kmps = np.sqrt(2. * c.gravity * mvir_g / rvir_cm) \
        #               * 1e-5
        hsel = [slice(None, None, None)] * 3
        #sumaxes = (0, 1)
        cosmopars = {key: val for key, val 
                     in f['Header/cosmopars'].attrs.items()}
        
        if rrange_rvir is not None:
            rbins_rvir = f['axis_0/bins'][:]
            rimin = np.where(np.isclose(rrange_rvir[0], rbins_rvir))[0]
            rimax = np.where(np.isclose(rrange_rvir[-1], rbins_rvir))[0]
            if rrange_rvir[0] == -np.inf:
                rimin = 0
            else:
                rimin = rimin[0]
            if rrange_rvir[-1] == np.inf:
                rimax = len(rbins_rvir)
            else:
                rimax = rimax[0]
            hsel[0] = slice(rimin, rimax, None)
        if trange_logk is not None:
            tbins_logk = f['axis_2/bins'][:]
            timin = np.where(np.isclose(trange_logk[0], tbins_logk))[0]
            timax = np.where(np.isclose(trange_logk[-1], tbins_logk))[0]
            if trange_logk[0] == -np.inf:
                timin = 0
            else:
                timin = timin[0]
            if trange_logk[-1] == np.inf:
                timax = len(tbins_logk)
            else:
                timax = timax[0]
            hsel[2] = slice(timin, timax, None)
        # _hist = np.sum(hist[tuple(hsel)], axis=sumaxes)
        #qtedges = f['axis_2/bins'][:]
        #perc = mu.percentiles_from_histogram(_hist, qtedges, axis=0,
        #                                     percentiles=np.array(perc))
        total = np.sum(hist[tuple(hsel)])
    out = {'cosmopars': cosmopars,
           'mvir_g': mvir_g,
           'total': total}
    return out

def getneonmassfrac_rbins(simname, snapnum, rrange_rvir=(0.1, 1.),
                          weight='gasmass'):
    axis2 = 'NeonAbundance'
    filen = datadir_rvcen + filetemp_rvcen.format(simname=simname, 
                                                  snapnum=snapnum,
                                                  weight=weight,
                                                  axis2=axis2)
    with h5py.File(filen, 'r') as f:
        hist = f['histogram/histogram'][:]
        islog = bool(f['histogram'].attrs['log'])
        if islog:
            hist = 10**hist
        mvir_g = f['Header/inputpars/halodata'].attrs['Mvir_g']
        hsel = [slice(None, None, None)] * 3
        cosmopars = {key: val for key, val 
                     in f['Header/cosmopars'].attrs.items()}
        
        if rrange_rvir is not None:
            rbins_rvir = f['axis_0/bins'][:]
            rimin = np.where(np.isclose(rrange_rvir[0], rbins_rvir))[0]
            rimax = np.where(np.isclose(rrange_rvir[-1], rbins_rvir))[0]
            if rrange_rvir[0] == -np.inf:
                rimin = 0
            else:
                rimin = rimin[0]
            if rrange_rvir[-1] == np.inf:
                rimax = len(rbins_rvir)
            else:
                rimax = rimax[0]
            hsel[0] = slice(rimin, rimax, None)
        _hist = np.sum(hist[tuple(hsel)], axis=(0, 1))
        metbins = f['axis_2/bins'][:]
        metcens = 0.5 * (metbins[:-1] + metbins[1:])
        if bool(f['axis_2'].attrs['log']):
            metcens = 10**metcens # linear average, and avoid issues from Z=0
        av = np.average(metcens, weights=_hist)
        
    out = {'cosmopars': cosmopars,
           'mvir_g': mvir_g,
           'wtdav': av}
    return out

def gettotZmassfrac_rbins(simname, snapnum, rrange_rvir=(0.1, 1.),
                          trange_logk=None):
    filen = datadir_rcen + filetemp_rcen_Ztot.format(simname=simname, 
                                                     snapnum=snapnum)
    with h5py.File(filen, 'r') as f:
        hist = f['histogram/histogram'][:]
        islog = bool(f['histogram'].attrs['log'])
        if islog:
            hist = 10**hist
        mvir_g = f['Header/inputpars/halodata'].attrs['Mvir_g']
        rvir_cm = f['Header/inputpars/halodata'].attrs['Rvir_cm']
        hsel = [slice(None, None, None)] * 3
        cosmopars = {key: val for key, val 
                     in f['Header/cosmopars'].attrs.items()}
        
        if rrange_rvir is not None:
            rbins_rvir = f['axis_0/bins'][:]
            rimin = np.where(np.isclose(rrange_rvir[0], rbins_rvir))[0]
            rimax = np.where(np.isclose(rrange_rvir[-1], rbins_rvir))[0]
            if rrange_rvir[0] == -np.inf:
                rimin = 0
            else:
                rimin = rimin[0]
            if rrange_rvir[-1] == np.inf:
                rimax = len(rbins_rvir)
            else:
                rimax = rimax[0]
            hsel[0] = slice(rimin, rimax, None)
        if trange_logk is not None:
            tbins_logk = f['axis_2/bins'][:]
            timin = np.where(np.isclose(trange_logk[0], tbins_logk))[0]
            timax = np.where(np.isclose(trange_logk[-1], tbins_logk))[0]
            if trange_logk[0] == -np.inf:
                timin = 0
            else:
                timin = timin[0]
            if trange_logk[-1] == np.inf:
                timax = len(tbins_logk)
            else:
                timax = timax[0]
            hsel[2] = slice(timin, timax, None)
        _hist = np.sum(hist[tuple(hsel)], axis=(0, 2))
        metbins = f['axis_1/bins'][:]
        metcens = 0.5 * (metbins[:-1] + metbins[1:])
        if bool(f['axis_2'].attrs['log']):
            metcens = 10**metcens # linear average, and avoid issues from Z=0
        av = np.average(metcens, weights=_hist)
        
    out = {'cosmopars': cosmopars,
           'mvir_g': mvir_g,
           'rvir_cm': rvir_cm, 
           'wtdav': av,
           'weight': 'gasmass',
           }
    return out

def getmstellar_rbins(simname, snapnum, rrange_rvir=(0.1, 1.)):
    filen = datadir_rcen + filetemp_rcen_mstar.format(simname=simname, 
                                                      snapnum=snapnum)
    with h5py.File(filen, 'r') as f:
        hist = f['histogram/histogram'][:]
        islog = bool(f['histogram'].attrs['log'])
        if islog:
            hist = 10**hist
        mvir_g = f['Header/inputpars/halodata'].attrs['Mvir_g']
        rvir_cm = f['Header/inputpars/halodata'].attrs['Rvir_cm']
        hsel = [slice(None, None, None)] * 2
        cosmopars = {key: val for key, val 
                     in f['Header/cosmopars'].attrs.items()}
        
        if rrange_rvir is not None:
            rbins_rvir = f['axis_0/bins'][:]
            rimin = np.where(np.isclose(rrange_rvir[0], rbins_rvir))[0]
            rimax = np.where(np.isclose(rrange_rvir[-1], rbins_rvir))[0]
            if rrange_rvir[0] == -np.inf:
                rimin = 0
            else:
                rimin = rimin[0]
            if rrange_rvir[-1] == np.inf:
                rimax = len(rbins_rvir)
            else:
                rimax = rimax[0]
            hsel[0] = slice(rimin, rimax, None)
        total = np.sum(hist[tuple(hsel)])
        
    out = {'cosmopars': cosmopars,
           'mvir_g': mvir_g,
           'rvir_cm': rvir_cm, 
           'total': total,
           'weight': 'gasmass',
           }
    return out

def getline_masses(outdct, simname=None, snapnum=None, rrange_rvir=None, 
                   trange_logk=None, weight=None, first=False):
    head = '\t'.join(['simname', 'snapnum', 
                      'rmin_rvir', 'rmax_rvir',
                      'tmin_logk', 'tmax_logk', 
                      'weight',
                      'total [g or num. part.]',
                      'Mvir_g',
                      'Omega_b', 
                      'Omega_m', 
                      'redshift']) + '\n'
    if first:
        return head
    values = [simname, f'{snapnum:d}', 
              f'{rrange_rvir[0]:.3f}', f'{rrange_rvir[-1]:.3f}',
              f'{trange_logk[0]:.3f}', f'{trange_logk[-1]:.3f}',
              weight,
              f'{outdct["total"]:.8e}',
              f'{outdct["mvir_g"]:.8e}',
              f'{outdct["cosmopars"]["omegab"]:.6f}',
              f'{outdct["cosmopars"]["omegam"]:.6f}',
              f'{outdct["cosmopars"]["z"]:.6f}',
              ]
    return '\t'.join(values) + '\n'
    
def getdata_masses():
    out = getline_masses(None, simname=None, snapnum=None, 
                        rrange_rvir=None, trange_logk=None, 
                        weight=None, first=True)
    for simname in simnames_all:
        snapnums = getsnaps(simname)
        for snapnum in snapnums:
            for rrange_rvir in rsels_cat + rsels_rad:
                for trange_logk in tranges_logk:
                    for weight in ['gasmass', 'Neon', 'Ne8']:
                        outdct = getmass_rbins(simname, snapnum, 
                                               rrange_rvir=rrange_rvir,
                                               trange_logk=trange_logk,
                                               weight=weight)
                        line = getline_masses(outdct, simname=simname, 
                                              snapnum=snapnum, 
                                              rrange_rvir=rrange_rvir, 
                                              trange_logk=trange_logk, 
                                              weight=weight, first=False)
                        out = out + line
    #outfile = datadir_rvcen + 'gas_Neon_Ne8_masses_rTcuts.dat'
    # with FIRE-3.x test data:
    outfile = datadir_rvcen + 'gas_Neon_Ne8_masses_rTcuts_v2.dat'
    with open(outfile, 'w') as f:
        f.write(out)

def getline_wtdav(outdct, simname=None, snapnum=None, rrange_rvir=None, 
                  weight=None, first=False):
    head = '\t'.join(['simname', 'snapnum', 
                      'rmin_rvir', 'rmax_rvir',
                      'weight',
                      'wtd. mean Ne abundance (mass fraction)',
                      'Mvir_g',
                      'Omega_b', 
                      'Omega_m', 
                      'redshift']) + '\n'
    if first:
        return head
    values = [simname, f'{snapnum:d}', 
              f'{rrange_rvir[0]:.3f}', f'{rrange_rvir[-1]:.3f}',
              weight,
              f'{outdct["wtdav"]:.8e}',
              f'{outdct["mvir_g"]:.8e}',
              f'{outdct["cosmopars"]["omegab"]:.6f}',
              f'{outdct["cosmopars"]["omegam"]:.6f}',
              f'{outdct["cosmopars"]["z"]:.6f}',
              ]
    return '\t'.join(values) + '\n'
    
def getdata_wtdav():
    out = getline_wtdav(None, simname=None, snapnum=None, 
                        rrange_rvir=None, weight=None, first=True)
    for simname in simnames_all:
        snapnums = getsnaps(simname)
        for snapnum in snapnums:
            for rrange_rvir in rsels_cat + rsels_rad:
                for weight in ['gasmass', 'gasvol']:
                    outdct = getneonmassfrac_rbins(simname, snapnum, 
                                                   rrange_rvir=rrange_rvir,
                                                   weight=weight)
                    line = getline_wtdav(outdct, simname=simname, 
                                         snapnum=snapnum, 
                                         rrange_rvir=rrange_rvir, 
                                         weight=weight, first=False)
                    out = out + line
    #outfile = datadir_rvcen + 'mean_ZNe_by_mass_volume_rcuts.dat'
    # with FIRE-3.x test data
    outfile = datadir_rvcen + 'mean_ZNe_by_mass_volume_rcuts_v2.dat'            
    with open(outfile, 'w') as f:
        f.write(out)

def getline_avZtot(outdct, simname=None, snapnum=None, rrange_rvir=None, 
                   trange_logk=None, first=False):
    head = '\t'.join(['simname', 'snapnum', 
                      'rmin_rvir', 'rmax_rvir',
                      'tmin_logk', 'tmax_logk',
                      'weight',
                      'wtd. mean total metallicity (mass fraction)',
                      'Mvir_g',
                      'Rvir_cm',
                      'Omega_b', 
                      'Omega_m', 
                      'redshift']) + '\n'
    if first:
        return head
    values = [simname, f'{snapnum:d}', 
              f'{rrange_rvir[0]:.3f}', f'{rrange_rvir[-1]:.3f}',
              f'{trange_logk[0]:.3f}', f'{trange_logk[-1]:.3f}',
              outdct['weight'],
              f'{outdct["wtdav"]:.8e}',
              f'{outdct["mvir_g"]:.8e}',
              f'{outdct["rvir_cm"]:.8e}',
              f'{outdct["cosmopars"]["omegab"]:.6f}',
              f'{outdct["cosmopars"]["omegam"]:.6f}',
              f'{outdct["cosmopars"]["z"]:.6f}',
              ]
    return '\t'.join(values) + '\n'
    
def getdata_avZtot():
    out = getline_avZtot(None, simname=None, snapnum=None,
                         rrange_rvir=None, 
                         trange_logk=None, first=True)
    for simname in simnames_all:
        snapnums = getsnaps(simname)
        for snapnum in snapnums:
            for rrange_rvir in rsels_cat + rsels_rad:
                for trange_logk in tranges_logk:
                    outdct = gettotZmassfrac_rbins(simname, snapnum, 
                                                   rrange_rvir=rrange_rvir,
                                                   trange_logk=trange_logk)
                    line = getline_avZtot(outdct, simname=simname, 
                                          snapnum=snapnum, 
                                          rrange_rvir=rrange_rvir, 
                                          trange_logk=trange_logk, 
                                          first=False)
                    out = out + line
    #outfile = datadir_rvcen + 'mean_Ztot_by_mass_rcuts_tcuts.dat'
    # with FIRE-3.x test data
    outfile = datadir_rvcen + 'mean_Ztot_by_mass_rcuts_tcuts_v2.dat'
    with open(outfile, 'w') as f:
        f.write(out)


def getline_mstellar(outdct, simname=None, snapnum=None, rrange_rvir=None, 
                     first=False):
    head = '\t'.join(['simname', 'snapnum', 
                      'rmin_rvir', 'rmax_rvir',
                      'total current stellar mass [g]',
                      'Mvir_g',
                      'Rvir_cm',
                      'Omega_b', 
                      'Omega_m', 
                      'redshift']) + '\n'
    if first:
        return head
    values = [simname, f'{snapnum:d}', 
              f'{rrange_rvir[0]:.3f}', f'{rrange_rvir[-1]:.3f}',
              f'{outdct["total"]:.8e}',
              f'{outdct["mvir_g"]:.8e}',
              f'{outdct["rvir_cm"]:.8e}',
              f'{outdct["cosmopars"]["omegab"]:.6f}',
              f'{outdct["cosmopars"]["omegam"]:.6f}',
              f'{outdct["cosmopars"]["z"]:.6f}',
              ]
    return '\t'.join(values) + '\n'
    
def getdata_mstellar():
    out = getline_mstellar(None, simname=None, snapnum=None,
                           rrange_rvir=None, first=True)
    for simname in simnames_all:
        snapnums = getsnaps(simname)
        for snapnum in snapnums:
            for rrange_rvir in rsels_cat:
                    outdct = getmstellar_rbins(simname, snapnum, 
                                               rrange_rvir=rrange_rvir)
                    line = getline_mstellar(outdct, simname=simname, 
                                            snapnum=snapnum, 
                                            rrange_rvir=rrange_rvir,
                                            first=False)
                    out = out + line
    #outfile = datadir_rvcen + 'stellarmass_rcuts.dat'
    # with FIRE-3.x data
    outfile = datadir_rvcen + 'stellarmass_rcuts_v2.dat'
    with open(outfile, 'w') as f:
        f.write(out)