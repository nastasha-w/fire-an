import h5py
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np

import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.math_utils as mu
import fire_an.utils.opts_locs as ol

histdir = '/projects/b1026/nastasha/hists/r_vr_all2/'
savedir = '/projects/b1026/nastasha/plotdata/'

def get_histdata(filen):
    with h5py.File(filen, 'r') as f:
        hdatpath = 'Header/inputpars/halodata'
        rvir_cm = f[hdatpath].attrs['Rvir_cm']
        rv = f['axis_0/bins'][:]
        _runits = f['axis_0'].attrs['units'].decode()
        yv = f['axis_2/bins'][:]
        ykey = f['axis_2'].attrs['qty_type'].decode()
        if ykey == 'sim-direct':
            _path = 'axis_2/qty_type_args'
            _key = 'field'
            ykey = (ykey, f[_path].attrs[_key].decode())
        elif ykey == 'ion':
            _path = 'axis_2/qty_type_args'
            _key = 'ion'
            ykey = (ykey, f[_path].attrs[_key].decode())
        elif ykey == 'Metal':
            _path = 'axis_2/qty_type_args'
            _key = 'element'
            ykey = (ykey, f[_path].attrs[_key].decode())
        hist_raw = f['histogram/histogram'][:]
        islog = bool(f['histogram'].attrs['log'])
        if islog:
            hist_raw = 10**hist_raw
        wkey = f['histogram'].attrs['weight_type'].decode()
        if wkey == 'ion':
            _path = 'histogram/weight_type_args'
            _key = 'ion'
            wkey = (wkey, f[_path].attrs[_key].decode())
        elif wkey == 'Metal':
            _path = 'histogram/weight_type_args'
            _key = 'element'
            wkey = (wkey, f[_path].attrs[_key].decode())
        redshift = f['Header/cosmopars'].attrs['z']
        hist = np.sum(hist_raw, axis=1) # just keep r, yqty distribution
    out = {'Rvir_cm': rvir_cm,
           'runits': _runits,
           'z': redshift,
           'hist': hist,
           'rbins': rv,
           'ybins': yv,
           'ykey': ykey,
           'wkey': wkey}
    return out

def calc_3dprof(simname, snapnum, yqty, weight, percvals=(0.1, 0.5, 0.9)):
    if yqty == 'hdens' and ('crheatfix' in simname 
                            or simname in sl.m12plus_f3nobh \
                                + sl.m12plus_f3nobh_lores):
        # renaming issue; both files store H density
        _yqty = 'density'
    else:
        _yqty = yqty
    filen = histdir + (f'hist_rcen_vcen_{_yqty}_by_{weight}_{simname}'
                       f'_snap{snapnum}_bins1_v1_hvcen.hdf5')
    histdata = get_histdata(filen)
    yvals = mu.percentiles_from_histogram(histdata['hist'], histdata['ybins'], axis=1, 
                                          percentiles=np.array(percvals))
    return histdata, yvals, filen
    
def save_3dprofs():
    outfile = savedir + 'radprof3d_nH_T_ZNe_by_vol_Ne8_opt2.hdf5'
    weights = ['Ne8', 'gasvol']
    yqtys = ['temperature', 'hdens', 'NeonAbundance']
    percvals = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.98, 0.99]
    sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    sims_f2md = sl.m12_f2md
    sims_m12plus = sl.m12plus_f3nobh + sl.m12plus_f3nobh_lores
    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    snaps_f2md = sl.snaps_f2md
    sims_all = sims_sr + sims_hr + sims_f2md + sims_m12plus
    
    with h5py.File(outfile, 'a') as fo:
        for simn in sims_all:
            g1 = fo.create_group(simn)
            snaps = (snaps_sr if simn in sims_sr 
                     else snaps_hr if simn in sims_hr
                     else snaps_f2md if simn in sims_f2md
                     else snaps_hr if simn in sims_m12plus
                     else None)
            for snap in snaps:
                g2 = g1.create_group(f'snap_{snap}')
                for weight in weights:
                    g3 = g2.create_group(f'{weight}_weighted')
                    for yqty in yqtys:
                        grp = g3.create_group(yqty)
                        hdata, yv, filen = calc_3dprof(simn, snap, 
                                                       yqty, weight,
                                                       percvals=(percvals))
                        for perc, _yv in zip(percvals, yv):
                            grp.create_dataset(f'perc_{perc:.2f}', data=_yv)
                        grp.attrs.create('histogram_file', np.string_(filen))
                        grp.attrs.create('Rvir_cm', hdata['Rvir_cm'])
                        grp.attrs.create('z', hdata['z'])
                        grp.attrs.create('runits', np.string_(hdata['runits']))
                        grp.create_dataset('rbins', data=hdata['rbins'])


