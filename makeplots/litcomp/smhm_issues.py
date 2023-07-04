
import h5py
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fire_an.makeplots.tol_colors as tc
import fire_an.makeplots.plot_utils as pu
import fire_an.mstar_mhalo.analytical as smhman
import fire_an.mstar_mhalo.loader_smdpl_sfr as smhmld
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu
import fire_an.utils.math_utils as mu
import fire_an.utils.opts_locs as ol

#datadir_b18 = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/'
datadir_b18 = '/projects/b1026/nastasha/extdata/'
dfilen_b18 = datadir_b18 + 'data_burchett_etal_2019_table1.txt'

def getdata_b19():
    #TODO CHECK: R200c or R200m! -> seems to be M200c
    # assuming impact parameters are physical/proper kpc
    #TODO ask: table footnote f says 2 systems for one Ne VIII absorber
    #          but only one line with that footnote and N(Ne VIII) value
    data_bur = pd.read_csv(dfilen_b18, comment='#', sep='\t')
    ## calculate halo masses
    # from Burchett et al. (2019):
    cosmopars_bur = {'h': 0.677, 'omegam': 0.31, 'omegalambda': 0.69}
    def hmfunc(x):
        csm = cosmopars_bur.copy()
        csm.update({'z': x.zgal, 'a': 1. / (1. + x.zgal)})
        mv = cu.mvir_from_rvir(x.rvir_kpc * 1e-3 * c.cm_per_mpc, 
                               csm, meandef='200c')
        return np.log10(mv / c.solar_mass)
    data_bur = data_bur.assign(log_Mvir_Msun=lambda x: hmfunc(x))
    return data_bur

def readin_cengaldata(simnames):
    datafn = ol.dir_halodata + 'pvcengal.hdf5'
    sms = []
    zs = []
    with h5py.File(datafn, 'r') as f:
        for sfn in simnames:
            _lsms = []
            _lz = []
            smgrp = f[sfn]
            sngrpns = [key for key in smgrp.keys() if key.startswith('snap_')]
            for sngrpn in sngrpns:
                sngrp = smgrp[sngrpn]
                csmpath = 'pv0/doc/halodata_doc_dict/cosmopars_dict'
                _lz.append(sngrp[csmpath].attrs['z'])
                _lsms.append(sngrp['pv0/doc'].attrs['mstar_gal_g'])
            zo = np.argsort(_lz)
            _lsz = np.array(_lz)[zo]
            _lsm = np.array(_lsms)[zo]
            _lsm /= c.solar_mass
            _lsm = np.log10(_lsm)
            sms.append(_lsm)
            zs.append(_lsz)
    return sms, zs

def readin_halodata(simnames):
    datafn = ol.filen_halocenrvir
    hms_bn98 = []
    hms_200c = []
    zs = []
    with h5py.File(datafn, 'r') as f:
        for sfn in simnames:
            _lbn98 = []
            _l200c = []
            _lz = []
            smgrp = f[sfn]
            sngrpns = [key for key in smgrp.keys() if key.startswith('snap_')]
            for sngrpn in sngrpns:
                sngrp = smgrp[sngrpn]
                _z = sngrp['cosmopars'].attrs['z']
                if _z < 0.4:
                    continue
                _lz.append(_z)
                _lbn98.append(sngrp['cen0/Rvir_BN98'].attrs['Mvir_g'])
                _l200c.append(sngrp['cen0/Rvir_200c'].attrs['Mvir_g'])
            zo = np.argsort(_lz)
            _lz = np.array(_lz)[zo]
            _lbn98 = np.log10(np.array(_lbn98)[zo] / c.solar_mass)
            _l200c = np.log10(np.array(_l200c)[zo] / c.solar_mass)
            hms_bn98.append(_lbn98)
            hms_200c.append(_l200c)
            zs.append(_lz)
    return hms_bn98, hms_200c, zs

def getddata_fire():
    simnames = sl.m12_hr_all2 + sl.m12_sr_all2 \
               + sl.m13_hr_all2 + sl.m13_sr_all2
    for sn in sl.buglist1:
        if sn in simnames:
            simnames.remove(sn)
    hms_bn98, hms_200c, zs = readin_halodata(simnames)
    sms, zs = readin_cengaldata(simnames)
    return simnames, hms_bn98, hms_200c, sms, zs


def plot_sm_hm_b19_vs_fire():
    data_b19 = getdata_b19()
    z_bur = data_b19['zgal']
    ms_bur = data_b19['log_Mstar_Msun']
    ms_bur_err = data_b19['log_Mstar_Msun_err']
    mh_bur = data_b19['log_Mvir_Msun']
    isul = data_b19['log_N_Ne8_isUL']

    simnames, hms_bn98, hms_200c, sms, zs = getddata_fire()

    physcolors = sl.physcolors
    icclasses = [sn[:3] for sn in simnames]
    physclasses = ['noBH' if '_sdp1e10_' in sn
                   else 'AGN-CR' if '_MHDCRspec1_' in sn
                   else 'AGN-noCR'
                   for sn in simnames]
    _b19colors = tc.tol_cset('vibrant')
    color_allb19 = _b19colors[0]
    color_detb19 = _b19colors[1]
    alpha = 0.3

    fig = plt.figure(figsize=(11., 6.))
    grid = gsp.GridSpec(nrows=3, ncols=3, hspace=0., wspace=0.4)
    fontsize = 12

    ## stellar masses
    axfire = fig.add_subplot(grid[0, 0])
    axdist = fig.add_subplot(grid[1, 0])
    axindiv = fig.add_subplot(grid[2, 0])

    axindiv.set_xlabel('$\\log_{10}\\, \\mathrm{M}_{\\star} \\;'
                       '[\\mathrm{M}_{\\odot}]$',
                       fontsize=fontsize)
    axindiv.set_ylabel('B+19 pdfs', fontsize=fontsize)
    axdist.set_ylabel('sum B+19 pdfs', fontsize=fontsize)
    axfire.set_ylabel('FIRE counts', fontsize=fontsize)
    axindiv.tick_params(which='both', direction='in', 
                        labelsize=fontsize - 1., top=True, 
                        right=True, labelbottom=True, labelleft=True) 
    axdist.tick_params(which='both', direction='in', 
                       labelsize=fontsize - 1., top=True, 
                       right=True, labelbottom=False, labelleft=True)
    axfire.tick_params(which='both', direction='in', 
                       labelsize=fontsize - 1., top=True, 
                       right=True, labelbottom=False, labelleft=True)  
    binsize = 0.2
    _xmin = min(np.min(ms_bur - 2. * ms_bur_err), np.min(sms) - binsize)
    _xmax = max(np.max(ms_bur + 2. * ms_bur_err), np.max(sms) + binsize)
    xmin = np.floor(_xmin / binsize) * binsize
    xmax = np.ceil(_xmax / binsize) * binsize
    xbins = np.arange(xmin, xmax + 0.5 * binsize, binsize)
    xcens = 0.5 * (xbins[:-1] + xbins[1:])

    b19_uldist = np.zeros(len(xcens))
    b19_detdist = np.zeros(len(xcens))
    for ms, mserr, _isul in zip(ms_bur, ms_bur_err, isul):
        _pcbin = smhman.cumulgauss((xbins - ms) / mserr)
        _pdf = np.diff(_pcbin) / binsize
        if _isul:
            color = color_allb19
            b19_uldist += _pdf
        else:
            color = color_detb19
            b19_detdist += _pdf
        axindiv.plot(xcens, _pdf, color=color, alpha=alpha)
    axdist.bar(xcens, b19_detdist, width=binsize, align='center',
               color=color_detb19)
    axdist.bar(xcens, b19_uldist, bottom=b19_detdist, width=binsize, 
               align='center', color=color_allb19)
    
    firehists = {'m12': {'noBH': np.zeros(len(xcens)),
                         'AGN-noCR': np.zeros(len(xcens)),
                         'AGN-CR': np.zeros(len(xcens))},
                 'm13': {'noBH': np.zeros(len(xcens)),
                         'AGN-noCR': np.zeros(len(xcens)),
                         'AGN-CR': np.zeros(len(xcens))}
                 }
    for physclass, icclass, sm in zip(physclasses, icclasses, sms):
        _hist, _ = np.histogram(sm, bins=xbins)
        firehists[icclass][physclass] += _hist
    bottom = np.zeros(len(xcens))
    for icclass in ['m12', 'm13']:
        if icclass == 'm13':
            kwa = {'hatch': 'x'}
        else:
            kwa = {}
        for physclass in ['noBH', 'AGN-noCR', 'AGN-CR']:
            color = physcolors[physclass]
            _hist = firehists[icclass][physclass]
            axfire.bar(xcens, _hist, width=binsize, align='center',
                       color=color, bottom=bottom, **kwa)
            bottom += _hist
    
    axfire.set_xlim((xmin, xmax))
    axdist.set_xlim((xmin, xmax))
    axindiv.set_xlim((xmin, xmax))
    axfire.set_title('stellar mass', fontsize=fontsize)
    

    ## halo masses (200c/B19)
    axfire = fig.add_subplot(grid[0, 1])
    axdist = fig.add_subplot(grid[1, 1])
    axindiv = fig.add_subplot(grid[2, 1])

    axindiv.set_xlabel('$\\log_{10}\\, \\mathrm{M}_{\\mathrm{200c}} \\;'
                       '[\\mathrm{M}_{\\odot}]$',
                       fontsize=fontsize)
    axindiv.set_ylabel('B+19 pdfs', fontsize=fontsize)
    axdist.set_ylabel('sum B+19 pdfs', fontsize=fontsize)
    axfire.set_ylabel('FIRE counts', fontsize=fontsize)
    axindiv.tick_params(which='both', direction='in', 
                        labelsize=fontsize - 1., top=True, 
                        right=True, labelbottom=True, labelleft=True) 
    axdist.tick_params(which='both', direction='in', 
                       labelsize=fontsize - 1., top=True, 
                       right=True, labelbottom=False, labelleft=True)
    axfire.tick_params(which='both', direction='in', 
                       labelsize=fontsize - 1., top=True, 
                       right=True, labelbottom=False, labelleft=True)  
    binsize = 0.2
    _xmin = min(np.min(mh_bur) - 0.4, np.min(hms_200c) - binsize, 10.4)
    _xmax = max(np.max(mh_bur) + 0.4, np.max(hms_200c) + binsize, 14.6)
    xmin = np.floor(_xmin / binsize) * binsize
    xmax = np.ceil(_xmax / binsize) * binsize
    xbins = np.arange(xmin, xmax + 0.5 * binsize, binsize)
    xcens = 0.5 * (xbins[:-1] + xbins[1:])
    binsize_fine_fac = 4
    xbins_fine = np.arange(xmin, xmax + 0.5 * binsize / binsize_fine_fac, 
                           binsize / binsize_fine_fac)
    xcens_fine = 0.5 * (xbins_fine[:-1] + xbins_fine[1:])

    b19_uldist = np.zeros(len(xcens))
    b19_detdist = np.zeros(len(xcens))
    print('mh (B19), mh (conversion B19 method)')
    for mh, ms, mserr, _isul, _zgal in zip(mh_bur, ms_bur, ms_bur_err, 
                                           isul, z_bur):
        mhvals = xbins_fine
        msfrommh = np.log10(smhman.mstar_burchett_etal_2019(10**xbins_fine,
                                                            _zgal))
        _Pcbin_ms = smhman.cumulgauss((msfrommh - ms) / mserr)
        _Pbin_ms = np.diff(_Pcbin_ms)
        mh_mhtoms = mu.linterpsolve(msfrommh, xbins_fine, ms)
        print(mh, mh_mhtoms)
        mh_pdist = _Pbin_ms / np.diff(xbins_fine)
        mh_pdist_add = np.sum(_Pbin_ms.reshape(len(xcens), binsize_fine_fac),
                              axis=1) / np.diff(xbins)
        if _isul:
            color = color_allb19
            b19_uldist += mh_pdist_add
        else:
            color = color_detb19
            b19_detdist += mh_pdist_add
        axindiv.plot(xcens_fine, mh_pdist, color=color, alpha=alpha)
        y_cenest = mu.linterpsolve(xcens_fine, mh_pdist, mh_mhtoms)
        axindiv.scatter([mh_mhtoms], [y_cenest],
                        color=color, marker='o', s=5,
                        alpha=alpha)
    axdist.bar(xcens, b19_detdist, width=binsize, align='center',
               color=color_detb19)
    axdist.bar(xcens, b19_uldist, bottom=b19_detdist, width=binsize, 
               align='center', color=color_allb19)
    
    firehists = {'m12': {'noBH': np.zeros(len(xcens)),
                         'AGN-noCR': np.zeros(len(xcens)),
                         'AGN-CR': np.zeros(len(xcens))},
                 'm13': {'noBH': np.zeros(len(xcens)),
                         'AGN-noCR': np.zeros(len(xcens)),
                         'AGN-CR': np.zeros(len(xcens))}
                 }
    for physclass, icclass, sm in zip(physclasses, icclasses, hms_200c):
        _hist, _ = np.histogram(sm, bins=xbins)
        firehists[icclass][physclass] += _hist
    bottom = np.zeros(len(xcens))
    for icclass in ['m12', 'm13']:
        if icclass == 'm13':
            kwa = {'hatch': 'x'}
        else:
            kwa = {}
        for physclass in ['noBH', 'AGN-noCR', 'AGN-CR']:
            color = physcolors[physclass]
            _hist = firehists[icclass][physclass]
            axfire.bar(xcens, _hist, width=binsize, align='center',
                       color=color, bottom=bottom, **kwa)
            bottom += _hist
    
    axfire.set_xlim((xmin, xmax))
    axdist.set_xlim((xmin, xmax))
    axindiv.set_xlim((xmin, xmax))
    axfire.set_title('halo mass (B+19, 200c)', fontsize=fontsize)

    # halo masses (UM/BN98)
    axfire = fig.add_subplot(grid[0, 2])
    axdist = fig.add_subplot(grid[1, 2])
    axindiv = fig.add_subplot(grid[2, 2])

    axindiv.set_xlabel('$\\log_{10}\\, \\mathrm{M}_{\\mathrm{vir}} \\;'
                       '[\\mathrm{M}_{\\odot}]$',
                       fontsize=fontsize)
    axindiv.set_ylabel('B+19 pdfs', fontsize=fontsize)
    axdist.set_ylabel('sum B+19 pdfs', fontsize=fontsize)
    axfire.set_ylabel('FIRE counts', fontsize=fontsize)
    axindiv.tick_params(which='both', direction='in', 
                        labelsize=fontsize - 1., top=True, 
                        right=True, labelbottom=True, labelleft=True) 
    axdist.tick_params(which='both', direction='in', 
                       labelsize=fontsize - 1., top=True, 
                       right=True, labelbottom=False, labelleft=True)
    axfire.tick_params(which='both', direction='in', 
                       labelsize=fontsize - 1., top=True, 
                       right=True, labelbottom=False, labelleft=True)  
    binsize = 0.1
    _xmin = min(np.min(mh_bur) - 0.8, np.min(hms_200c) - 3. * binsize,
                10.4)
    _xmax = max(np.max(mh_bur) + 2.0, np.max(hms_200c) + 4. * binsize,
                14.6)
    xmin = np.floor(_xmin / binsize) * binsize
    xmax = np.ceil(_xmax / binsize) * binsize
    xbins = np.arange(xmin, xmax + 0.5 * binsize, binsize)
    xcens = 0.5 * (xbins[:-1] + xbins[1:])

    b19_uldist = None
    b19_detdist = None
    histobj = smhmld.SMHMhists(np.array(z_bur), binsize=binsize)
    for mh, ms, mserr, _isul, _zgal in zip(mh_bur, ms_bur, ms_bur_err, 
                                           isul, z_bur):
        msbins_hist = histobj.getbins(_zgal, mode='ms')
        _Pcbin_ms = smhman.cumulgauss((msbins_hist - ms) / mserr)
        _Pbin_ms = np.diff(_Pcbin_ms)
        mhp, _mhbins = histobj.matrixconv(_Pbin_ms, _zgal, 
                                          mode='mstomh')
        mhpd = mhp / np.diff(_mhbins)
        mhcens = 0.5 * (_mhbins[:-1] + _mhbins[1:])
        if _isul:
            color = color_allb19
            if b19_uldist is None:
                b19_uldist = mhpd
                b19_uldist_bins = [np.array(_mhbins)]
            else:
                b19_uldist, b19_uldist_bins = mu.combine_hists(
                    b19_uldist, mhpd, b19_uldist_bins, [_mhbins], 
                    rtol=1e-5, atol=1e-8, add=True)
        else:
            color = color_detb19
            if b19_detdist is None:
                b19_detdist = mhpd
                b19_detdist_bins = [np.array(_mhbins)]
            else:
                b19_detdist, b19_detdist_bins = mu.combine_hists(
                    b19_detdist, mhpd, b19_detdist_bins, [_mhbins], 
                    rtol=1e-5, atol=1e-8, add=True)
        axindiv.plot(mhcens, mhpd, color=color, alpha=alpha)
    b19_uldist, b19_uldist_bins = mu.combine_hists(
                    b19_uldist, np.zeros((len(b19_detdist),)), 
                    b19_uldist_bins, b19_detdist_bins, 
                    rtol=1e-5, atol=1e-8, add=True)
    b19_detdist, b19_detdist_bins = mu.combine_hists(
                    b19_detdist, np.zeros((len(b19_uldist),)), 
                    b19_detdist_bins, b19_uldist_bins, 
                    rtol=1e-5, atol=1e-8, add=True)
    b19_uldist_cens = 0.5 * (b19_uldist_bins[0][:-1] 
                             +  b19_uldist_bins[0][1:])
    b19_detdist_cens = 0.5 * (b19_detdist_bins[0][:-1] 
                              +  b19_detdist_bins[0][1:])
    axdist.bar(b19_detdist_cens, b19_detdist, width=binsize, align='center',
               color=color_detb19)
    axdist.bar(b19_uldist_cens, b19_uldist, bottom=b19_detdist, width=binsize, 
               align='center', color=color_allb19)
    
    firehists = {'m12': {'noBH': np.zeros(len(xcens)),
                         'AGN-noCR': np.zeros(len(xcens)),
                         'AGN-CR': np.zeros(len(xcens))},
                 'm13': {'noBH': np.zeros(len(xcens)),
                         'AGN-noCR': np.zeros(len(xcens)),
                         'AGN-CR': np.zeros(len(xcens))}
                 }
    for physclass, icclass, sm in zip(physclasses, icclasses, hms_bn98):
        _hist, _ = np.histogram(sm, bins=xbins)
        firehists[icclass][physclass] += _hist
    bottom = np.zeros(len(xcens))
    for icclass in ['m12', 'm13']:
        if icclass == 'm13':
            kwa = {'hatch': 'x'}
        else:
            kwa = {}
        for physclass in ['noBH', 'AGN-noCR', 'AGN-CR']:
            color = physcolors[physclass]
            _hist = firehists[icclass][physclass]
            axfire.bar(xcens, _hist, width=binsize, align='center',
                       color=color, bottom=bottom, **kwa)
            bottom += _hist
    
    axfire.set_xlim((xmin, xmax))
    axdist.set_xlim((xmin, xmax))
    axindiv.set_xlim((xmin, xmax))
    axfire.set_title('halo mass (UM, BN98)', fontsize=fontsize)

    outdir = '/projects/b1026/nastasha/imgs/datacomp/smhm/'
    outname = 'mstar_mh_distcomp_B19_methods_fire.pdf'
    plt.savefig(outdir + outname, bbox_inches='tight')

def plot_hm_uncertainty_um_b19():
    data_b19 = getdata_b19()
    z_bur = data_b19['zgal']
    ms_bur = data_b19['log_Mstar_Msun']
    ms_bur_err = data_b19['log_Mstar_Msun_err']
    isul = data_b19['log_N_Ne8_isUL']
    binsize = 0.1

    _b19colors = tc.tol_cset('vibrant')
    color_allb19 = _b19colors[0]
    color_detb19 = _b19colors[1]

    fig = plt.figure(figsize=(5.5, 4.))
    ax = fig.add_subplot(1, 1, 1)
    fontsize = 12

    ax.set_xlabel('$\\log_{10}\\, \\mathrm{M}_{\\mathrm{vir}} \\;'
                  '[\\mathrm{M}_{\\odot}]$',
                  fontsize=fontsize)
    ax.set_ylabel('probability density', fontsize=fontsize)

    ax.tick_params(which='both', direction='in', 
                   labelsize=fontsize - 1., top=True, 
                   right=True, labelbottom=True, labelleft=True)  
    
    histobj = smhmld.SMHMhists(np.array(z_bur), binsize=binsize)
    ullabeldone = False
    detlabeldone = False
    for ms, mserr, _isul, _zgal in zip(ms_bur, ms_bur_err, 
                                       isul, z_bur):
        msbins_hist = histobj.getbins(_zgal, mode='ms')
        _Pcbin_ms = smhman.cumulgauss((msbins_hist - ms) / mserr)
        _Pbin_ms = np.diff(_Pcbin_ms)
        mhp, _mhbins = histobj.matrixconv(_Pbin_ms, _zgal, 
                                          mode='mstomh')
        mhpd = mhp / np.diff(_mhbins)
        mhcens = 0.5 * (_mhbins[:-1] + _mhbins[1:])
        if _isul:
            color = color_allb19
            alpha = 0.3
            if not ullabeldone:
                _label = 'B+19 Ne VIII UL'
                ullabeldone = True
            else:
                _label = None
        else:
            color = color_detb19
            alpha = 0.6
            if not detlabeldone:
                _label = 'B+19 Ne VIII det.'
                detlabeldone = True
            else:
                _label = None
        ax.plot(mhcens, mhpd, color=color, alpha=alpha, label=_label)
    
    ax.set_xlim((10.4, 14.8))
    ax.legend(fontsize=fontsize, loc='upper right')

    outdir = '/projects/b1026/nastasha/imgs/datacomp/smhm/'
    outname = 'mh_pdists_b19_umcalc.pdf'
    plt.savefig(outdir + outname, bbox_inches='tight')
    
    




    



