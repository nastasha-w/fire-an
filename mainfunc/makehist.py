

import h5py
import numbers as num
import numpy as np
import os

import fire_an.mainfunc.get_qty as gq
import fire_an.mainfunc.haloprop as hp
import fire_an.readfire.readin_fire_data as rf
import fire_an.utils.constants_and_units as c
import fire_an.utils.h5utils as h5u


def getaxbins(minfinite, maxfinite, bin, extendmin=True, extendmax=True):
    if isinstance(bin, int):
        bins = np.linspace(minfinite, maxfinite, bin + 1)
    elif isinstance(bin, float):
        minbin = np.floor(minfinite / bin) * bin
        maxbin = np.ceil(maxfinite / bin) * bin
        bins = np.arange(minbin, maxbin + 0.5 * bin, bin)
    else:
        bins = np.array(bin)
        if minfinite < bins[0]:
            extendmin = True
        if maxfinite >= bins[1]:
            extendmax = True
    if extendmin:
        bins = np.append(-np.inf, bins)
    if extendmax:
        bins = np.append(bins, np.inf)
    return bins


def histogram_radprof(dirpath, snapnum,
                      weighttype, weighttype_args, axtypes, axtypes_args,
                      particle_type=0, 
                      center='shrinksph', rbins=(0., 1.), runit='Rvir',
                      logweights=True, logaxes=True, axbins=0.1,
                      outfilen=None, overwrite=True):
    '''
    make a weightype, weighttype_args weighted histogram of 
    axtypes, axtypes_args.

    Parameters:
    -----------
    dirpath: str
        path to the directory containing the 'output' directory with the
        snapshots
    snapnum: int
        snapshot number
    weightype: float
        what to weight the histogram by. Options are maptype options in 
        get_qty
    weighttype_args: dict
        additional arguments for what to weight the histogram by. Options 
        are maptype_args options in get_qty, with some extra options from
        process_typeargs_coords
    axtypes: list
        list of what to histogram; each entry is one histogram dimension.
        Options are maptype options in get_qty.
    axtypes_args: list of dicts
        list of additional arguments for the histogram dimensions. Options 
        are maptype_args options in get_qty, with some extra options from
        process_typeargs_coords. Mind the 'density' option for
        ions and metals. These are matched to axtypes by list index.
    particle_type: int
        particle type to project (follows FIRE format)
    center: str
        how to find the halo center.
        'AHFsmooth': use halo_00000_smooth.dat from AHF 
        'rockstar-maxmass': highest mass halo at snapshot from Rockstar
        'rockstar-mainprog': main progenitor of most massive halo at
                           final snapshot from Rockstar
        'rockstar-<int>': halo with snapshot halo catalogue index <int>
                          from Rockstar 
        'shrinksph': Imran's shrinking spheres method
    rbins: array-like of floats
        bin edges in 3D distance from the halo center. Ignored if center is 
        None.
    runit: {'Rvir', 'pkpc'}
        unit to use for the rbins. These bins are never log values. 
    logweights: bool
        save log of the weight sum in each bin instead of the linear sum.
    logaxes: bool or list of bools
        save and process the histogram axis quantities in log units instead
        of linear. If a list, this is applied to each dimension by matching
        list index.
    axbins: int, float, array, or list of those.
        int: use that many bins between whatever min and maxfinite values are
             present in the data
        float: use bins of that size, in (log) cgs units. The edges are chosen
             so zero is an edge if the value range includes zero.
             A useful option to allow stacking/comparison without using too 
             much storage.
        array: just the bin edges directly. Monotonically increasing, (log) 
             cgs units.
        Note that the range are always extended with -np.inf and np.inf if the
        values include non-finite ones or values outside the specified range.
        If a list is given, the options are specified per histogram axis, 
        matched by list index. Note that if giving an array option, it should 
        be enclosed in a list to ensure it is not interpreted as a per-axis
        option list.
        Units are always (log) cgs. 
    outfilen: str
        file to save to output histogram to. None means no file is saved.
        The file must include the full path.
    overwrite: bool
        If a file with name outfilen already exists, overwrite it (True) or
        raise a ValueError (False)
    Output:
    -------
    file with saved histogram data, if a file is specified
    otherwise, the histogram and bins

    '''
    if outfilen is not None:
        if os.path.isfile(outfilen) and not overwrite:
            raise ValueError('File {} already exists.'.format(outfilen))

    todoc_gen = {}
    basepath = 'PartType{}/'.format(particle_type)
    _axvals = []
    _axbins = []
    _axdoc = []
    _axbins_outunit = []
    _logaxes = []
    _axtypes = []
    _axtypes_args = []
    if not hasattr(axbins, '__len__'):
        axbins = [axbins] * len(axtypes)
    if not hasattr(logaxes, '__len__'):
        logaxes = [logaxes] * len(axtypes)

    snap = rf.get_Firesnap(dirpath, snapnum)
    
    if center is not None:
        todoc_cen = {}
        todoc_cen['center_method'] = center
        if center == 'AHFsmooth':
            halodat = hp.mainhalodata_AHFsmooth(dirpath, snapnum)
            cen = np.array([halodat['Xc_ckpcoverh'], 
                            halodat['Yc_ckpcoverh'], 
                            halodat['Zc_ckpcoverh']])
            cen_cm = cen * snap.cosmopars.a * 1e-3 * c.cm_per_mpc \
                    / snap.cosmopars.h
            rvir_cm = halodat['Rvir_ckpcoverh'] * snap.cosmopars.a \
                    * 1e-3 * c.cm_per_mpc / snap.cosmopars.h
        elif center.startswith('rockstar'):
            select = center.split('-')[-1]
            if select not in ['maxmass', 'mainprog']:
                try:
                    select = int(select)
                except ValueError:
                    msg = 'invalid option for center: {}'.format(center)
                    raise ValueError(msg)
            halodat, _csm_halo = hp.halodata_rockstar(dirpath, snapnum,
                                                      select=select)
            cen = np.array([halodat['Xc_ckpc'], 
                            halodat['Yc_ckpc'], 
                            halodat['Zc_ckpc']])
            cen_cm = cen * snap.cosmopars.a * 1e-3 * c.cm_per_mpc 
            rvir_cm = halodat['Rvir_cm'] 
        elif center == 'shrinksph':
            halodat, _ = hp.gethalodata_shrinkingsphere(dirpath, snapnum,
                                                        meandef='BN98')
            cen_cm = np.array([halodat['Xc_cm'], 
                               halodat['Yc_cm'], 
                               halodat['Zc_cm']])
            rvir_cm = halodat['Rvir_cm']
            todoc_cen['Rvir_def'] = 'BN98'
        else:
            raise ValueError('Invalid center option {}'.format(center))
        todoc_cen['center_cm'] = cen_cm
        todoc_cen['Rvir_cm'] = rvir_cm
        todoc_cen['units'] = runit
        todoc_cen['log'] = False

        coords = snap.readarray_emulateEAGLE(basepath + 'Coordinates')
        coords_toCGS = snap.toCGS
        coords -= cen_cm / coords_toCGS
        
        if runit == 'Rvir':
            rbins_simu = np.array(rbins) * rvir_cm / coords_toCGS
            simu_to_runit = coords_toCGS / rvir_cm 
        elif runit == 'pkpc':
            rbins_simu = np.array(rbins) * c.cm_per_mpc * 1e-3 / coords_toCGS
            simu_to_runit = coords_toCGS / (c.cm_per_mpc * 1e-3)
        else:
            raise ValueError('Invalid runit option: {}'.format(runit))
        rbins2_simu = rbins_simu**2
        r2vals = np.sum(coords**2, axis=1)
        del coords
        filter = r2vals <= rbins2_simu[-1]
        r2vals = r2vals[filter]

        _axvals.append(r2vals)
        _axbins.append(rbins2_simu)
        _axdoc.append(todoc_cen)
        _axbins_outunit.append(np.sqrt(rbins2_simu) * simu_to_runit)
        _logaxes.append(False)
        _axtypes.append('halo_3Dradius')
        _axtypes_args.append({})
        filterdct = {'filter': filter}
    else:
        todoc_gen['info_halo'] = 'no halo particle selection applied'
        filterdct = {'filter': slice(None, None, None)}
        halodat = None
    
    for axt, axarg, logax, axb in zip(axtypes, axtypes_args, logaxes, axbins):
        if axt == 'coords':
            axarg, todoc = gq.process_typeargs_coords(dirpath, snapnum, axarg)
        else:
            todoc = {}
        qty, toCGS, _todoc = gq.get_qty(snap, particle_type, axt, axarg, 
                                        filterdct=filterdct)
        todoc.update(_todoc)
        if logax:
            qty = np.log10(qty)
        qty_good = np.isfinite(qty)
        minq = np.min(qty[qty_good])
        maxq = np.max(qty[qty_good])
        needext = not np.all(qty_good)
        if hasattr(axb, '__len__') \
                or (isinstance(axb, num.Number) and not isinstance(axb, int)):
            if logax:
                _axb = axb
                # log differences independent of units. smh.
                #_axb = axb - np.log10(toCGS)
            else:
                _axb = axb / toCGS
        else:
            _axb = axb
        usebins_simu = getaxbins(minq, maxq, _axb, extendmin=needext, 
                                 extendmax=needext)

        _axvals.append(qty)
        _axbins.append(usebins_simu)
        _axdoc.append(todoc)
        _logaxes.append(logax)
        _axtypes.append(axt)
        _axtypes_args.append(axarg)
        if logax:
            _bins_doc = usebins_simu + np.log10(toCGS)
        else:
            _bins_doc = usebins_simu * toCGS
        _axbins_outunit.append(_bins_doc)
    if weighttype == 'coords':
        weighttype_args, wt_todoc = gq.process_typeargs_coords(
            dirpath, snapnum, weighttype_args)
    else:
        wt_todoc = {}
    wt, wt_toCGS, _wt_todoc = gq.get_qty(snap, particle_type, weighttype,
                                         weighttype_args, 
                                         filterdct=filterdct)
    wt_todoc.update(_wt_todoc)
    #print(_axbins)
    maxperloop = 752**3 // 8
    if len(wt) <= maxperloop:
        hist, edges = np.histogramdd(_axvals, weights=wt, bins=_axbins)
    else:
        lentot = len(wt)
        numperloop = maxperloop
        slices = [slice(i * numperloop, 
                        min((i + 1) * numperloop, lentot), 
                        None) \
                  for i in range((lentot - 1) // numperloop + 1)]
        for slind in range(len(slices)):
            axdata_temp = [data[slices[slind]] for data in _axvals]
            hist_temp, edges_temp = np.histogramdd(axdata_temp, 
                                                   weights=wt[slices[slind]], 
                                                   bins=_axbins)
            if slind == 0 :
                hist = hist_temp
                edges = edges_temp
            else:
                hist += hist_temp
                if not np.all(np.array([np.all(edges[i] == edges_temp[i]) \
                              for i in range(len(edges))])):
                    msg = 'Error: edges mismatch in histogramming'+\
                          ' loop (slind = {})'.format(slind)
                    raise RuntimeError(msg)
    if logweights:
        hist = np.log10(hist)
        hist += np.log10(wt_toCGS)
    else:
        hist *= wt_toCGS

    if outfilen is not None:
        with h5py.File(outfilen, 'w') as f:
            # cosmopars (emulate make_maps format)
            hed = f.create_group('Header')
            cgrp = hed.create_group('cosmopars')
            csm = snap.cosmopars.getdct()
            for key in csm:
                cgrp.attrs.create(key, csm[key])

            # histogram and weight
            hgrp = f.create_group('histogram')
            hgrp.create_dataset('histogram', data=hist)
            hgrp.attrs.create('log', logweights)
            h5u.savedict_hdf5(hgrp, wt_todoc)
            hgrp.attrs.create('weight_type', np.string_(weighttype))
            wagrp = hgrp.create_group('weight_type_args')
            h5u.savedict_hdf5(wagrp, weighttype_args)
            
            # histogram axes
            for i in range(0, len(_axbins)):
                agrp = f.create_group('axis_{}'.format(i))
                _bins = _axbins_outunit[i]
                agrp.create_dataset('bins', data=_bins)
                agrp.attrs.create('log', _logaxes[i])
                if center is None:
                    agrp.attrs.create('bin_input', axbins[i])
                elif i == 0:
                    agrp.attrs.create('bin_input', rbins)
                else:
                    agrp.attrs.create('bin_input', axbins[i - 1])
                _todoc = _axdoc[i]
                h5u.savedict_hdf5(agrp, _todoc)
                agrp.attrs.create('qty_type', np.string_(_axtypes[i]))
                aagrp = agrp.create_group('qty_type_args')
                h5u.savedict_hdf5(aagrp, _axtypes_args[i])
            
            # direct input parameters
            igrp = hed.create_group('inputpars')
            _snf =  np.array([np.string_(fn) for fn in snap.filens])
            igrp.attrs.create('snapfiles', _snf)
            igrp.attrs.create('dirpath', np.string_(dirpath))
            igrp.attrs.create('particle_type', particle_type)
            igrp.attrs.create('outfilen', np.string_(outfilen))
            h5u.savedict_hdf5(igrp, todoc_gen)

            if halodat is not None:
                _grp = igrp.create_group('halodata')
                h5u.savedict_hdf5(_grp, halodat)