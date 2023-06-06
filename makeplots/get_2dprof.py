
import h5py
import numpy as np

import fire_an.utils.constants_and_units as c


def get_rval_massmap(filen, units='pkpc', weightmap=False,
                     absvals=False, weightrange=None):
    '''
    get radius and map quantity matched arrays

    Input:
    ------
    filen: str
        file name including full path
    units: {'pkpc', 'Rvir', 'cm'}
        units for the radius (impact parameter)
    weightmap: bool
        get the 'weightmap' array instead of the 'map' array
    absvals: bool
        get absolute values
    Returns:
    --------
    radii: float array (1D)
        impact parameters
    map values: float array (1D) 
        map values, matching impact parameters 
    '''
    with h5py.File(filen, 'r') as f:
        if weightmap:
            _map = f['weightmap'][:]
        else:
            _map = f['map'][:]
        if weightrange is not None:
            _smap = f['weightmap'][:]
            mapsel = _smap >= weightrange[0]
            mapsel &= _smap < weightrange[1]
            del _smap
        else:
            mapsel = (slice(None, None, None),) * 2
        if absvals:
            _map = np.abs(_map)
        shape = _map.shape
        xinds, yinds = np.indices(shape).astype(np.float32)
        # centered on halo
        xcen = 0.5 * float(shape[0])
        ycen = 0.5 * float(shape[1])
        dpix2 = (xinds + 0.5 - xcen)**2 + (yinds + 0.5 - ycen)**2
        dpix = np.sqrt(dpix2)
        
        dpix = dpix[mapsel]
        _map = _map[mapsel]
        
        pixssize_pkpc = f['Header/inputpars'].attrs['pixsize_pkpc']
        if units == 'pkpc':
            dpix *= pixssize_pkpc
        elif units == 'Rvir':
            # using Imran's shrinking spheres method
            rvir_cm = f['Header/inputpars/halodata'].attrs['Rvir_cm']
            dpix *= pixssize_pkpc * c.cm_per_mpc * 1e-3 / rvir_cm
        elif units == 'cm':
            dpix *= pixssize_pkpc * c.cm_per_mpc * 1e-3
    return dpix.flatten(), _map.flatten()

def get_profile_massmap(filen, rbins, rbin_units='pkpc',
                        profiles=[], weightmap=False,
                        absvals=False, weightrange=None):
    '''
    get values with impact parameter from maps. If multiple files
    are given, averages, percentiles, etc. are taken over the full
    sample, weighted by map pixel 

    Parameters:
    -----------
    filen: str or iterable of strings
        name(s) of the file(s) containing the map (with full path)
    rbins: array-like of floats
        impact parameter bin edges
    rbin_units: {'pkpc', 'Rvir', 'cm'}
        units used for the rbins
    profiles: iterable of strings
        which profiles to extract. Returned in order of input
        options are:
            'av-lin' (linear average), 
            'av-log' (log-space average),
            'perc-<float>' (percentile, values 0 -- 1)
            'min' (minimum)
            'max' (maximum)
            'fcov-<float>' (fraction of values >= given value)
    Returns:
    --------
    profiles: each an array of floats
        the profile values in each radial bin
    '''
    
    # single filename (string) or iterable of file names
    if isinstance(filen, type('')):
        filens = [filen]
    else:
        filens = filen
    first = True
    for _filen in filens:
        with h5py.File(_filen, 'r') as f:
            if weightmap:
                _islog = bool(f['weightmap'].attrs['log'])
            else:
                _islog = bool(f['map'].attrs['log'])
        rvals, mvs = get_rval_massmap(_filen, units=rbin_units,
                                      weightmap=weightmap,
                                      absvals=absvals,
                                      weightrange=weightrange)
        
        # searchsorted index 0 means value < rbins[0], 
        # index len(rbins) means values > rbins[-1] 
        rinds = np.searchsorted(rbins, rvals) - 1 
        _mvs_by_bin = [mvs[rinds == i] for i in range(len(rbins))]
        _mvs_by_bin = _mvs_by_bin[:-1]
        if first:
            mvs_by_bin = _mvs_by_bin
            islog = _islog
            first = False
        else:
            mvs_by_bin = [np.append(mvs, _mvs) \
                          for mvs, _mvs in zip(mvs_by_bin, _mvs_by_bin)]
            # not handling this for now since all maps should be log
            if _islog != islog:
                msg = 'Different files have log vs. non-log map values'
                raise RuntimeError(msg)

    out = []
    for prof in profiles:
        if prof == 'av-lin':
            if islog:
                _pr = np.log10(np.array([np.average(10**_mvs) \
                                         for _mvs in mvs_by_bin]))
            else:
                _pr = np.array([np.average(_mvs) \
                                for _mvs in mvs_by_bin])
        elif prof == 'av-log':
            if islog:
                _pr = np.array([np.average(_mvs) \
                                for _mvs in mvs_by_bin])
            else:
                _pr = 10**np.array([np.average(np.log10(_mvs)) \
                                    for _mvs in mvs_by_bin])
        elif prof == 'max':
            _pr = np.array([np.max(_mvs) for _mvs in mvs_by_bin])
        elif prof == 'max':
            _pr = np.array([np.min(_mvs) for _mvs in mvs_by_bin])
        elif prof.startswith('perc-'):
            pv = float(prof.split('-')[-1])
            _pr = np.array([np.quantile(_mvs, pv) if len(_mvs) > 0 \
                            else np.NaN for _mvs in mvs_by_bin])
        elif prof.startswith('fcov-'):
            pv = float(prof.split('-')[-1])
            _pr = np.array([np.sum(_mvs >= pv) / float(len(_mvs)) \
                            for _mvs in mvs_by_bin])
        out.append(_pr)
    return out