
import h5py
import numpy as np

import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt

from fire_an.ionrad.ion_utils import Linetable_PS20
import fire_an.mainfunc.get_qty as gq
import fire_an.makeplots.plot_utils as pu
import fire_an.readfire.readin_fire_data as rf
import fire_an.utils.constants_and_units as c



def test_ionbal_calc(dirpath, snapnum, ion, target_Z=0.01, delta_Z=0.001,
                     ps20depletion=False, outfilen='ionbal_test.hdf5',
                     lintable=False):
    snap = rf.get_Firesnap(dirpath, snapnum)
    cosmopars = snap.cosmopars.getdct()
    
    # filter sim. particles and calculate ion balances, rho, T, Z
    metallicity = snap.readarray_emulateEAGLE('PartType0/Metallicity')
    zfilter = metallicity >= target_Z - delta_Z
    zfilter &= metallicity <= target_Z + delta_Z
    metallicity = metallicity[zfilter]
    indct = {'filter': zfilter}
    ionbals = gq.get_ionfrac(snap, ion, indct=indct, table='PS20', 
                             simtype='fire',
                             ps20depletion=ps20depletion, lintable=lintable)
    temperature = snap.readarray_emulateEAGLE('PartType0/Temperature')[zfilter]
    temperature *= snap.toCGS
    hdens = snap.readarray_emulateEAGLE('PartType0/Density')[zfilter]
    hconv = snap.toCGS
    hdens *= snap.readarray_emulateEAGLE('PartType0/ElementAbundance/Hydrogen')[zfilter]
    hconv *= snap.toCGS
    hconv /= (c.atomw_H * c.u)
    hdens *= hconv
    
    # get corresponding ion balance table
    # for table read-in only, lin/log shouldn't matter; problems weren't there
    iontab = Linetable_PS20(ion, cosmopars['z'], emission=False, vol=True,
                            lintable=False)
    iontab.findiontable()
    tab_logT = iontab.logTK
    tab_lognH = iontab.lognHcm3
    tab_logZ = iontab.logZsol + np.log10(iontab.solarZ)
    tab_ionbal_T_Z_nH = iontab.iontable_T_Z_nH.copy()
    if ps20depletion:
        tab_ionbal_T_Z_nH = 10**tab_ionbal_T_Z_nH
        iontab.finddepletiontable()
        tab_depletion_T_Z_nH = iontab.depletiontable_T_Z_nH.copy()
        tab_ionbal_T_Z_nH *= (1. - 10**tab_depletion_T_Z_nH)
        tab_ionbal_T_Z_nH = np.log10(tab_ionbal_T_Z_nH)
        tab_ionbal_T_Z_nH[tab_ionbal_T_Z_nH < -50.] = -50.

    interpvalZ = np.log10(target_Z)
    iZhi = np.where(tab_logZ >= interpvalZ)[0][0]
    iZlo = np.where(tab_logZ <= interpvalZ)[0][-1]
    if iZlo == iZhi:
        tab_ionbal_T_nH = tab_ionbal_T_Z_nH[:, iZlo, :] 
    else:
        hiZ = tab_logZ[iZhi]
        loZ = tab_logZ[iZlo]
        tab_ionbal_T_nH = (hiZ - interpvalZ) / (hiZ - loZ) * tab_ionbal_T_Z_nH[:, iZlo, :] +\
                          (interpvalZ - loZ) / (hiZ - loZ) * tab_ionbal_T_Z_nH[:, iZhi, :]
    tab_ionbal_T_nH = 10**tab_ionbal_T_nH

    # save data
    with h5py.File(outfilen, 'w') as f:
        hed = f.create_group('Header')
        cgrp = hed.create_group('cosmopars')
        cosmopars = snap.cosmopars.getdct()
        for key in cosmopars:
            cgrp.attrs.create(key, cosmopars[key])
        hed.attrs.create('snapnum', snapnum)
        hed.attrs.create('filepath_first', np.string_(snap.firstfilen))
        _info = 'FIRE calculated ion balances and the underlying ion balance table'
        hed.attrs.create('info', np.string_(_info))
        hed.attrs.create('target_Z', target_Z)
        hed.attrs.create('delta_Z', delta_Z)
        hed.attrs.create('ion', np.string_(ion))
        hed.attrs.create('ps20depletion', ps20depletion)
        hed.attrs.create('lintable', lintable)
        
        gsim = f.create_group('simulation_data')
        gsim.create_dataset('ionbal', data=ionbals)
        gsim.create_dataset('T_K', data=temperature)
        gsim.create_dataset('nH_cm**-3', data=hdens)
        gsim.create_dataset('metallicity_abs_mass_frac', data=metallicity)
        
        gtab = f.create_group('iontab_data')
        print('About to save tab_ionbal_T_nH')
        print('{} / {} NaN'.format(np.sum(np.isnan(tab_ionbal_T_nH)), 
                                   np.prod(tab_ionbal_T_nH.shape)))
        print(tab_ionbal_T_nH)
        gtab.create_dataset('ionbal_T_nH', data=tab_ionbal_T_nH)
        gtab.create_dataset('logT_K', data=tab_logT)
        gtab.create_dataset('lognH_cm**-3', data=tab_lognH)

def run_ionbal_test(opt=0):
    dirpath1 = '/projects/b1026/snapshots/fire3/m13h206_m3e5/' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000/'
    simname1 = 'm13h206_m3e5__' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'
    snaps1 = [27, 45]
    ions1 = ['O{}'.format(i) for i in range(1, 10)]

    outdir =  '/projects/b1026/nastasha/tests/start_fire/ionbal_tests/'
    outtemplate1 = outdir + 'ionbal_test_PS20_{ion}_depletion-{dp}_Z-{Z}' + \
                            '_snap{snap:03d}_{sim}.hdf5'
    outtemplate2 = outdir + 'ionbal_test_PS20_{ion}_depletion-{dp}_Z-{Z}' + \
                            '_snap{snap:03d}_lintable-{lintable}_{sim}.hdf5'
    
    if opt >= 0 and opt < 6:
        dirpath = dirpath1
        simname = simname1
        ions = ions1
        ps20depletion = bool(opt % 2)
        snapnum = snaps1[opt // 4]
        target_Z = [0.01, 0.0001][(opt  // 2) % 2]
        delta_Z = 0.1 * target_Z
        outtemplate = outtemplate1
        dolintable = False
    if opt >= 6 and opt < 18:
        _opt = opt - 6
        __opt = _opt % 6
        lintable = bool(_opt // 6)
        dirpath = dirpath1
        simname = simname1
        ions = ions1
        ps20depletion = bool(__opt % 2)
        snapnum = snaps1[__opt // 4]
        target_Z = [0.01, 0.0001][(__opt  // 2) % 2]
        delta_Z = 0.1 * target_Z
        outtemplate = outtemplate2
        dolintable = True
    else:
        raise ValueError('Invalid opt {}'.format(opt))
    for ion in ions:
        if dolintable:
            outfilen = outtemplate.format(ion=ion, dp=ps20depletion, 
                                          Z=target_Z, sim=simname, 
                                          snap=snapnum, lintable=lintable)
            test_ionbal_calc(dirpath, snapnum, ion, target_Z=target_Z, 
                             delta_Z=delta_Z, ps20depletion=ps20depletion, 
                             outfilen=outfilen, lintable=lintable)
        else:
            outfilen = outtemplate.format(ion=ion, dp=ps20depletion, 
                                          Z=target_Z, sim=simname, 
                                          snap=snapnum)
            test_ionbal_calc(dirpath, snapnum, ion, target_Z=target_Z, 
                             delta_Z=delta_Z, ps20depletion=ps20depletion, 
                             outfilen=outfilen)
            
def ionbal_test(filens, simlabel=None):
    '''
    for a series of files containing the tables and data for a full
    ionisation series: plot the table and interpolated ion balances
    and get a sense of whether the differences are ok
    Also get the total element fraction, to check if it makes sense
    with depletion.
    '''
    if simlabel is None:
        simlabel = 'testhalo1-m13h206_m3e5'

    tot_iontab = None
    logTtab = None
    lognHtab = None

    logTsim = None
    lognHsim = None
    Zsim = None
    tot_ionsim = None

    target_Z = None
    delta_Z = None
    
    redshift = None
    ions = []

    imgdir = '/'.join(filens[0].split('/')[:-1]) + '/'
    _imgname = 'ionbal-test_{ion}_depletion-{dep}_Z-{met}_z-{z:.1f}{lt}_{simname}.pdf'

    for filen in filens:
        with h5py.File(filen, 'r') as f:
            # title info
            cosmopars = {key: val for key, val in \
                         f['Header/cosmopars'].attrs.items()}
            if redshift is None:
                redshift = cosmopars['z']
            elif not redshift == cosmopars['z']:
                msg = 'Input files have different redshifts; found in {}'
                raise ValueError(msg.format(filen))
            _delta_Z = f['Header'].attrs['delta_Z']
            if delta_Z is None:
                delta_Z = _delta_Z
            elif not delta_Z == _delta_Z:
                msg = 'Input files have different delta_Z; found in {}'
                raise ValueError(msg.format(filen))
            _target_Z = f['Header'].attrs['target_Z']
            if target_Z is None:
                target_Z = _target_Z
            elif not target_Z == _target_Z:
                msg = 'Input files have different target_Z; found in {}'
                raise ValueError(msg.format(filen))
            ion = f['Header'].attrs['ion'].decode()
            ps20depletion = bool(f['Header'].attrs['ps20depletion'])
            ions.append(ion)
            if 'lintable' in f['Header'].attrs:
                dolintable = True
                lintable = bool(f['Header'].attrs['lintable']) 
                ltlabel = '_lintable-{}'.format(lintable)
            else:
                dolintable = False
                ltlabel = ''
  
            # table data
            if logTtab is None:
                logTtab = f['iontab_data/logT_K'][:]
            elif not np.all(logTtab == f['iontab_data/logT_K'][:]):
                msg = 'Input files have different table log T bins; found in {}'
                raise ValueError(msg.format(filen))
            if lognHtab is None:
                lognHtab = f['iontab_data/lognH_cm**-3'][:]
            elif not np.all(lognHtab == f['iontab_data/lognH_cm**-3'][:]):
                msg = 'Input files have different table log nH bins; found in {}'
                raise ValueError(msg.format(filen))
            iontab = f['iontab_data/ionbal_T_nH'][:].T
            if tot_iontab is None:
                tot_iontab = iontab
            else:
                tot_iontab += iontab 
            
            # sim data
            if logTsim is None:
                logTsim = np.log10(f['simulation_data/T_K'])
            elif not np.all(logTsim == np.log10(f['simulation_data/T_K'])):
                msg = 'Input files have different simulation log T values; found in {}'
                raise ValueError(msg.format(filen))
            if lognHsim is None:
                lognHsim = np.log10(f['simulation_data/nH_cm**-3'])
            elif not np.all(lognHsim == np.log10(f['simulation_data/nH_cm**-3'])):
                msg = 'Input files have different simulation log nH values; found in {}'
                raise ValueError(msg.format(filen))
            if Zsim is None:
                Zsim = f['simulation_data/metallicity_abs_mass_frac'][:]
            elif not np.all(Zsim == f['simulation_data/metallicity_abs_mass_frac'][:]):
                msg = 'Input files have different simulation Z values; found in {}'
                raise ValueError(msg.format(filen))
            ionsim = f['simulation_data/ionbal'][:]
            if tot_ionsim is None:
                tot_ionsim = ionsim
            else:
                tot_ionsim += ionsim
        
        outname = imgdir + _imgname.format(simname=simlabel, ion=ion, 
                                           dep=ps20depletion,
                                           met=target_Z, z=redshift,
                                           lt=ltlabel)
            
        title = '{ion} PS20 table at z={z:.2f}, Z={met:.1e} vs. interp. FIRE data, dust depl. {dep}' 
        title = title.format(ion=ion, dep=ps20depletion, z=redshift, met=target_Z)
        if dolintable:
            title = title + ', lin. table {}'.format(lintable)
        
        fig = plt.figure(figsize=(13., 4.))
        grid = gsp.GridSpec(nrows=1, ncols=5, hspace=0.0, wspace=0.5, 
                            width_ratios=[1., 0.1, 1., 0.1, 1.],
                            left=0.05, right=0.95)
        axes = [fig.add_subplot(grid[0, i]) for i in range(5)]
        fontsize = 12
        cmap = 'viridis'
        size = 10.
        vmin = -10
        vmax = 0.

        fig.suptitle(title, fontsize=fontsize)

        ax = axes[0]
        cax = axes[1]

        xedges = lognHtab[:-1] - 0.5 * np.diff(lognHtab)
        xend = [lognHtab[-1] - 0.5 * (lognHtab[-1] - lognHtab[-2]),
                lognHtab[-1] + 0.5 * (lognHtab[-1] - lognHtab[-2])
                ]
        xedges = np.append(xedges, xend)
        yedges = logTtab[:-1] - 0.5 * np.diff(logTtab)
        yend = [logTtab[-1] - 0.5 * (logTtab[-1] - logTtab[-2]),
                logTtab[-1] + 0.5 * (logTtab[-1] - logTtab[-2])
                ]
        yedges = np.append(yedges, yend)
        img = ax.pcolormesh(xedges, yedges, np.log10(iontab.T), cmap=cmap, 
                            rasterized=True, vmin=vmin, vmax=vmax)
        plt.colorbar(img, cax=cax)
        cax.set_ylabel('log ion fraction', fontsize=fontsize)

        ax.scatter(lognHsim, logTsim, s=size, c=np.log10(ionsim),
                   edgecolor='black', cmap=cmap, vmin=vmin, vmax=vmax,
                   rasterized=True)
        ax.set_ylabel('$\\log \\, \\mathrm{T} \\; [\\mathrm{K]}$', 
                      fontsize=fontsize)
        xlabel = '$\\log \\, \\mathrm{n}_{\\mathrm{H}} \\;' + \
                 ' [\\mathrm{cm}^{-3}]$'
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.tick_params(which='both', axis='both', labelsize=fontsize - 1)
        cax.tick_params(labelsize=fontsize - 1, labelright=False, labelleft=True,
                        right=False, left=True)
        cax.yaxis.set_label_position('left')
        
        ax = axes[2]
        cax = axes[3]
        
        Tinds = np.argmin(np.abs(logTsim[:, np.newaxis] 
                                 - logTtab[np.newaxis, :]), axis=1)
        nHinds = np.argmin(np.abs(lognHsim[:, np.newaxis] 
                                 - lognHtab[np.newaxis, :]), axis=1)
        closest_gridtosim = iontab[(nHinds, Tinds)]
        dz = Zsim - target_Z
        vmin = -1. * delta_Z
        vmax = delta_Z
        
        delta_sim = closest_gridtosim - ionsim
        xvals_sim = np.random.uniform(low=0.0, high=1.0, size=len(logTsim))

        delta_tab1 = (iontab[1:, :] - iontab[:-1, :]).flatten()
        delta_tab2 = (iontab[:, 1:] - iontab[:, :-1]).flatten()
        xvals_tab1 = np.random.uniform(low=0.0, high=1.0, size=len(delta_tab1))
        xvals_tab2 = np.random.uniform(low=0.0, high=1.0, size=len(delta_tab2))

        ax.scatter(xvals_tab1, delta_tab1, s=0.5*size, color='black', 
                   label='$\\Delta$ table grid', rasterized=True)
        ax.scatter(xvals_tab1, -1. * delta_tab1, s=0.5*size, color='black',
                   rasterized=True)
        ax.scatter(xvals_tab2, delta_tab2, s=0.5*size, color='black',
                   rasterized=True)
        ax.scatter(xvals_tab2, -1. * delta_tab2, s=0.5*size, color='black',
                   rasterized=True)

        img = ax.scatter(xvals_sim, delta_sim, s=size, c=dz,
                         edgecolor='black', cmap='RdBu', vmin=vmin, vmax=vmax,
                         label='sim - table', rasterized=True)
        plt.colorbar(img, cax=cax, extend='neither')
        ax.set_ylabel('difference with nearest table value', 
                      fontsize=fontsize)
        #cax.set_ylabel('simulation Z - table Z', fontsize=fontsize,
        #               horizontalalignment='center', verticalalignment='center',
        #               x=0.5, y=0.5)
        cax.text(0.5, 0.5, 'simulation Z - table Z', fontsize=fontsize,
                 horizontalalignment='center', verticalalignment='center',
                 rotation='vertical', transform=cax.transAxes)
        #cax.yaxis.set_label_coords(0.5, 0.5)
        cax.tick_params(labelsize=fontsize - 3, labelright=False, labelleft=True,
                        right=False, left=True)
        #cax.yaxis.set_label_position('left')
        ax.legend(fontsize=fontsize - 1)
        ax.tick_params(labelbottom=False, bottom=False)

        ax = axes[4]

        nbins = 100
        maxv = np.max(np.abs(delta_tab1))
        maxv = max(maxv, np.max(np.abs(delta_tab2)))
        maxv = max(maxv, np.max(np.abs(delta_sim)))
        bins = np.linspace(-1. * maxv, maxv, nbins)
        
        tabvals = np.append(delta_tab1, -1. * delta_tab1)
        tabvals = np.append(tabvals, delta_tab2)
        tabvals = np.append(tabvals, -1. * delta_tab2)

        ax.set_yscale('log')
        ax.hist(tabvals, bins=bins, histtype='step', color='black',
                label='$\\Delta$ table grid', align='mid', density=True)
        ax.hist(delta_sim, bins=bins, histtype='step', color='blue',
                label='sim - table', linestyle='dashed', align='mid', 
                density=True)

        ax.set_xlabel('difference with nearest table value', 
                      fontsize=fontsize)
        ax.set_ylabel('probability density', fontsize=fontsize)
        ax.legend(fontsize=fontsize - 1)
        
        plt.savefig(outname, bbox_inches='tight')
    
    outname = imgdir + _imgname.format(simname=simlabel, ion='-'.join(ions), 
                                       dep=ps20depletion,
                                       met=target_Z, z=redshift,
                                       lt=ltlabel)
        
    title = 'sum of {ion} PS20 tables at z={z:.2f}, Z={met:.1e} vs. ' + \
            'interp. FIRE data, dust depl. {dep}' 
    title = title.format(ion=', '.join(ions), dep=ps20depletion, 
                         z=redshift, met=target_Z)
    if dolintable:
        title = title + ', lin. table {}'.format(lintable)
    
    fig = plt.figure(figsize=(13., 4.))
    grid = gsp.GridSpec(nrows=1, ncols=5, hspace=0.0, wspace=0.5, 
                        width_ratios=[1., 0.1, 1., 0.1, 1.],
                        left=0.05, right=0.95)
    axes = [fig.add_subplot(grid[0, i]) for i in range(5)]
    fontsize = 12
    cmap = 'RdBu'
    size = 10.
        
    vmax = np.max(np.abs(tot_ionsim - 1.))
    vmax = max(vmax, np.max(np.abs(tot_iontab - 1.)))
    vmin = 1. - vmax
    vmax = 1. + vmax

    fig.suptitle(title, fontsize=fontsize)

    ax = axes[0]
    cax = axes[1]

    xedges = lognHtab[:-1] - 0.5 * np.diff(lognHtab)
    xend = [lognHtab[-1] - 0.5 * (lognHtab[-1] - lognHtab[-2]),
            lognHtab[-1] + 0.5 * (lognHtab[-1] - lognHtab[-2])
            ]
    xedges = np.append(xedges, xend)
    yedges = logTtab[:-1] - 0.5 * np.diff(logTtab)
    yend = [logTtab[-1] - 0.5 * (logTtab[-1] - logTtab[-2]),
            logTtab[-1] + 0.5 * (logTtab[-1] - logTtab[-2])
            ]
    yedges = np.append(yedges, yend)
    img = ax.pcolormesh(xedges, yedges, tot_iontab.T, cmap=cmap, 
                        rasterized=True, vmin=vmin, vmax=vmax)
    plt.colorbar(img, cax=cax)
    cax.set_ylabel('ion fraction sum', fontsize=fontsize)

    ax.scatter(lognHsim, logTsim, s=size, c=tot_ionsim,
                edgecolor='black', cmap=cmap, vmin=vmin, vmax=vmax,
                rasterized=True)
    ax.set_ylabel('$\\log \\, \\mathrm{T} \\; [\\mathrm{K]}$', 
                    fontsize=fontsize)
    xlabel = '$\\log \\, \\mathrm{n}_{\\mathrm{H}} \\;' + \
                ' [\\mathrm{cm}^{-3}]$'
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.tick_params(which='both', axis='both', labelsize=fontsize - 1)
    cax.tick_params(labelsize=fontsize - 1, labelright=False, labelleft=True,
                    right=False, left=True)
    cax.yaxis.set_label_position('left')
    
    ax = axes[2]
    cax = axes[3]
    
    Tinds = np.argmin(np.abs(logTsim[:, np.newaxis] 
                                - logTtab[np.newaxis, :]), axis=1)
    nHinds = np.argmin(np.abs(lognHsim[:, np.newaxis] 
                                - lognHtab[np.newaxis, :]), axis=1)
    closest_gridtosim = tot_iontab[(nHinds, Tinds)]
    dz = Zsim - target_Z
    vmin = -1. * delta_Z
    vmax = delta_Z
    
    delta_sim = closest_gridtosim - tot_ionsim
    xvals_sim = np.random.uniform(low=0.0, high=1.0, size=len(logTsim))

    delta_tab1 = (tot_iontab[1:, :] - tot_iontab[:-1, :]).flatten()
    delta_tab2 = (tot_iontab[:, 1:] - tot_iontab[:, :-1]).flatten()
    xvals_tab1 = np.random.uniform(low=0.0, high=1.0, size=len(delta_tab1))
    xvals_tab2 = np.random.uniform(low=0.0, high=1.0, size=len(delta_tab2))

    ax.scatter(xvals_tab1, delta_tab1, s=0.5*size, color='black', 
                label='$\\Delta$ table grid', rasterized=True)
    ax.scatter(xvals_tab1, -1. * delta_tab1, s=0.5*size, color='black',
                rasterized=True)
    ax.scatter(xvals_tab2, delta_tab2, s=0.5*size, color='black',
                rasterized=True)
    ax.scatter(xvals_tab2, -1. * delta_tab2, s=0.5*size, color='black',
                rasterized=True)

    img = ax.scatter(xvals_sim, delta_sim, s=size, c=dz,
                     edgecolor='black', cmap='RdBu', vmin=vmin, vmax=vmax,
                     label='sim - table', rasterized=True)
    plt.colorbar(img, cax=cax, extend='neither')
    ax.set_ylabel('difference with nearest table value', 
                fontsize=fontsize)
    #cax.set_ylabel('simulation Z - table Z', fontsize=fontsize,
    #               horizontalalignment='center', verticalalignment='center',
    #               x=0.5, y=0.5)
    cax.text(0.5, 0.5, 'simulation Z - table Z', fontsize=fontsize,
                horizontalalignment='center', verticalalignment='center',
                rotation='vertical', transform=cax.transAxes)
    #cax.yaxis.set_label_coords(0.5, 0.5)
    cax.tick_params(labelsize=fontsize - 3, labelright=False, labelleft=True,
                    right=False, left=True)
    #cax.yaxis.set_label_position('left')
    ax.legend(fontsize=fontsize - 1)
    ax.tick_params(labelbottom=False, bottom=False)

    ax = axes[4]

    nbins = 100
    maxv = np.max(np.abs(delta_tab1))
    maxv = max(maxv, np.max(np.abs(delta_tab2)))
    maxv = max(maxv, np.max(np.abs(delta_sim)))
    bins = np.linspace(-1. * maxv, maxv, nbins)
    
    tabvals = np.append(delta_tab1, -1. * delta_tab1)
    tabvals = np.append(tabvals, delta_tab2)
    tabvals = np.append(tabvals, -1. * delta_tab2)

    ax.set_yscale('log')
    ax.hist(tabvals, bins=bins, histtype='step', color='black',
            label='$\\Delta$ table grid', align='mid', density=True)
    ax.hist(delta_sim, bins=bins, histtype='step', color='blue',
            label='sim - table', linestyle='dashed', align='mid', 
            density=True)

    ax.set_xlabel('difference with nearest table value', 
                    fontsize=fontsize)
    ax.set_ylabel('probability density', fontsize=fontsize)
    ax.legend(fontsize=fontsize - 1)
    
    plt.savefig(outname, bbox_inches='tight')

    # sum test tables
    _imgname = 'ionbal-test_tablesum_{ion}_depletion-{dep}_Z-{met}_z-{z:.1f}{lt}.pdf'
    outname = imgdir + _imgname.format(simname=simlabel, ion=ion, 
                                       dep=ps20depletion,
                                       met=target_Z, z=redshift,
                                       lt=ltlabel)

    title = 'sum of {ion} PS20 tables at z={z:.2f}, Z={met:.1e} vs. 1.0 dust depl. {dep}' 
    title = title.format(ion=', '.join(ions), dep=ps20depletion, 
                         z=redshift, met=target_Z)
    
    fig = plt.figure(figsize=(13., 4.))
    grid = gsp.GridSpec(nrows=1, ncols=5, hspace=0.0, wspace=0.5, 
                        width_ratios=[1., 0.1, 1., 0.1, 1.],
                        left=0.05, right=0.95)
    axes = [fig.add_subplot(grid[0, i]) for i in range(5)]
    fontsize = 12
    cmap = 'RdBu'
    size = 10.
        
    vmax = np.max(np.abs(tot_ionsim - 1.))
    vmax = max(vmax, np.max(np.abs(tot_iontab - 1.)))
    vmin = 1. - vmax
    vmax = 1. + vmax

    fig.suptitle(title, fontsize=fontsize)

    ax = axes[0]
    cax = axes[1]

    xedges = lognHtab[:-1] - 0.5 * np.diff(lognHtab)
    xend = [lognHtab[-1] - 0.5 * (lognHtab[-1] - lognHtab[-2]),
            lognHtab[-1] + 0.5 * (lognHtab[-1] - lognHtab[-2])
            ]
    xedges = np.append(xedges, xend)
    yedges = logTtab[:-1] - 0.5 * np.diff(logTtab)
    yend = [logTtab[-1] - 0.5 * (logTtab[-1] - logTtab[-2]),
            logTtab[-1] + 0.5 * (logTtab[-1] - logTtab[-2])
            ]
    yedges = np.append(yedges, yend)
    img = ax.pcolormesh(xedges, yedges, tot_iontab.T, cmap=cmap, 
                        rasterized=True, vmin=vmin, vmax=vmax)
    plt.colorbar(img, cax=cax)
    cax.set_ylabel('ion fraction sum', fontsize=fontsize)

    ax.set_ylabel('$\\log \\, \\mathrm{T} \\; [\\mathrm{K]}$', 
                    fontsize=fontsize)
    xlabel = '$\\log \\, \\mathrm{n}_{\\mathrm{H}} \\;' + \
                ' [\\mathrm{cm}^{-3}]$'
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.tick_params(which='both', axis='both', labelsize=fontsize - 1)
    cax.tick_params(labelsize=fontsize - 1, labelright=False, labelleft=True,
                    right=False, left=True)
    cax.yaxis.set_label_position('left')
    
    ax = axes[2]
    cax = axes[3]
    
    delta_tab1 = (tot_iontab - 1.).flatten()
    xvals_tab1 = np.random.uniform(low=0.0, high=1.0, size=len(delta_tab1))
    c_tab1 = np.repeat(lognHtab[:, np.newaxis], tot_iontab.shape[1], axis=1)
    vmid = 1.5
    vmin = np.min(lognHtab)
    vmax = np.max(lognHtab)
    tv = (vmid - vmin) / (vmax - vmin) 
    cmap = pu.paste_cmaps(['plasma', 'viridis'], [vmin, vmid, vmax], 
                          [(0., tv), (tv, 1.)])

    img = ax.scatter(xvals_tab1, delta_tab1, s=0.5*size, c=c_tab1, 
                     rasterized=True, cmap=cmap)

    plt.colorbar(img, cax=cax, extend='neither')
    ax.set_ylabel('table sum - 1.', fontsize=fontsize)
    #cax.set_ylabel('simulation Z - table Z', fontsize=fontsize,
    #               horizontalalignment='center', verticalalignment='center',
    #               x=0.5, y=0.5)
    clabel = '$\\log_{10}\\,\\mathrm{n}_{\\mathrm{H}}\\;[\\mathrm{cm}^{-3}]$'
    cax.text(0.5, 0.5, clabel, fontsize=fontsize,
                horizontalalignment='center', verticalalignment='center',
                rotation='vertical', transform=cax.transAxes)
    #cax.yaxis.set_label_coords(0.5, 0.5)
    cax.tick_params(labelsize=fontsize - 3, labelright=False, labelleft=True,
                    right=False, left=True)
    #cax.yaxis.set_label_position('left')
    ax.legend(fontsize=fontsize - 1)
    ax.tick_params(labelbottom=False, bottom=False)

    ax = axes[4]

    nbins = 100
    maxv = np.max(np.abs(delta_tab1))
    bins = np.linspace(-1. * maxv, maxv, nbins)
    
    tabvals = delta_tab1

    ax.set_yscale('log')
    ax.hist(tabvals, bins=bins, histtype='step', color='black',
            label=None, align='mid', density=True)

    ax.set_xlabel('table sum - 1.', 
                    fontsize=fontsize)
    ax.set_ylabel('probability density', fontsize=fontsize)
    ax.legend(fontsize=fontsize - 1)
    
    plt.savefig(outname, bbox_inches='tight')

    outname = imgdir + _imgname.format(simname=simlabel, ion='-'.join(ions), 
                                       dep=ps20depletion,
                                       met=target_Z, z=redshift,
                                       lt=ltlabel)

    # sum test interp
    _imgname = 'ionbal-test_interp-sum_{ion}_depletion-{dep}_Z-{met}'+\
               '_z-{z:.1f}{lt}_{simname}.pdf'
    outname = imgdir + _imgname.format(simname=simlabel, ion='-'.join(ions), 
                                       dep=ps20depletion,
                                       met=target_Z, z=redshift,
                                       lt=ltlabel)
        
    title = 'sum of interpolated {ion} PS20 tables at z={z:.2f}, '+ \
            'Z={met:.1e} from FIRE data, dust depl. {dep}' 
    title = title.format(ion=', '.join(ions), dep=ps20depletion, 
                         z=redshift, met=target_Z)
    if dolintable:
        title = title + ', lin. table {}'.format(lintable)
    
    fig = plt.figure(figsize=(13., 4.))
    grid = gsp.GridSpec(nrows=1, ncols=5, hspace=0.0, wspace=0.5, 
                        width_ratios=[1., 0.1, 1., 0.1, 1.],
                        left=0.05, right=0.95)
    axes = [fig.add_subplot(grid[0, i]) for i in range(5)]
    fontsize = 12
    cmap = 'RdBu'
    size = 10.
        
    vmax = np.max(np.abs(tot_ionsim - 1.))
    vmax = max(vmax, np.max(np.abs(tot_iontab - 1.)))
    vmin = 1. - vmax
    vmax = 1. + vmax

    fig.suptitle(title, fontsize=fontsize)

    ax = axes[0]
    cax = axes[1]

    img = ax.scatter(lognHsim, logTsim, s=size, c=tot_ionsim,
                     edgecolor='black', cmap=cmap, vmin=vmin, vmax=vmax,
                     rasterized=True)
    plt.colorbar(img, cax=cax)
    cax.set_ylabel('ion fraction sum', fontsize=fontsize)

    ax.set_ylabel('$\\log \\, \\mathrm{T} \\; [\\mathrm{K]}$', 
                    fontsize=fontsize)
    xlabel = '$\\log \\, \\mathrm{n}_{\\mathrm{H}} \\;' + \
                ' [\\mathrm{cm}^{-3}]$'
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.tick_params(which='both', axis='both', labelsize=fontsize - 1)
    cax.tick_params(labelsize=fontsize - 1, labelright=False, labelleft=True,
                    right=False, left=True)
    cax.yaxis.set_label_position('left')
    
    ax = axes[2]
    cax = axes[3]
    
    delta_sim = tot_ionsim - 1.
    xvals_sim = np.random.uniform(low=0.0, high=1.0, size=len(logTsim))
    vmid = 1.5
    vmin = np.min(lognHsim)
    vmax = np.max(lognHsim)
    tv = (vmid - vmin) / (vmax - vmin) 
    if vmid >= vmax:
        cmap = 'plasma'
    else:
        cmap = pu.paste_cmaps(['plasma', 'viridis'], [vmin, vmid, vmax], 
                              [(0., tv), (tv, 1.)])

    img = ax.scatter(xvals_sim, delta_sim, s=size, c=lognHsim,
                     edgecolor='black', cmap=cmap,
                     label='interp. sum - 1.', rasterized=True)
    plt.colorbar(img, cax=cax, extend='neither')
    ax.set_ylabel('interp. sum - 1.', fontsize=fontsize)
    #cax.set_ylabel('simulation Z - table Z', fontsize=fontsize,
    #               horizontalalignment='center', verticalalignment='center',
    #               x=0.5, y=0.5)
    cax.text(0.5, 0.5, clabel, fontsize=fontsize,
                horizontalalignment='center', verticalalignment='center',
                rotation='vertical', transform=cax.transAxes)
    #cax.yaxis.set_label_coords(0.5, 0.5)
    cax.tick_params(labelsize=fontsize - 3, labelright=False, labelleft=True,
                    right=False, left=True)
    #cax.yaxis.set_label_position('left')
    ax.legend(fontsize=fontsize - 1)
    ax.tick_params(labelbottom=False, bottom=False)

    ax = axes[4]

    nbins = 100
    maxv = np.max(np.abs(delta_sim))
    bins = np.linspace(-1. * maxv, maxv, nbins)

    ax.set_yscale('log')
    ax.hist(delta_sim, bins=bins, histtype='step', color='blue',
            label=None, linestyle='dashed', align='mid', 
            density=True)

    ax.set_xlabel('interp. sim - 1.', fontsize=fontsize)
    ax.set_ylabel('probability density', fontsize=fontsize)
    ax.legend(fontsize=fontsize - 1)
    
    plt.savefig(outname, bbox_inches='tight')


def run_ionbal_tests(index):
    # laptop
    ddir = '/Users/nastasha/ciera/tests/fire_start/ionbal_tests/'
    simname = 'm13h206_m3e5__m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1'+\
              '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'
    pre = 'ionbal_test_PS20'
    ftmp = ['{pre}_{ion}_depletion-False_Z-0.01_snap045{lt}_{simname}.hdf5',
            '{pre}_{ion}_depletion-False_Z-0.0001_snap027{lt}_{simname}.hdf5',
            '{pre}_{ion}_depletion-False_Z-0.01_snap027{lt}_{simname}.hdf5',
            '{pre}_{ion}_depletion-True_Z-0.0001_snap027{lt}_{simname}.hdf5',
            '{pre}_{ion}_depletion-True_Z-0.01_snap027{lt}_{simname}.hdf5',
            '{pre}_{ion}_depletion-True_Z-0.01_snap045{lt}_{simname}.hdf5',
            ]
    _index = index % 6
    lti = index // 6
    if lti == 0:
        lt = ''
    elif lti == 1:
        lt = '_lintable-False'
    elif lti == 2:
        lt = '_lintable-True'
    ions = ['O{}'.format(i) for i in range(1, 10)]
    filens = [ddir + ftmp[_index].format(pre=pre, simname=simname, 
                                         ion=ion, lt=lt) \
              for ion in ions]
    ionbal_test(filens)

def test_tablesum_direct(ztargets = [0., 1., 3.], element='oxygen'):
    '''
    forget interpolations, do the maps look right on their own.
    '''
    fn = '/Users/nastasha/phd/tables/ionbal/lines_sp20/'+\
         'UVB_dust1_CR1_G1_shield1.hdf5'
    outdir = '/Users/nastasha/ciera/tests/fire_start/ionbal_tests/'
    savename = outdir + 'sumtest-nodepl_z-{z:.3f}_{elt}_tables-direct.pdf'

    with h5py.File(fn, 'r') as f:
        mgrp = f['Tdep/IonFractions']
        eltkeys = list(mgrp.keys())
        eltgrp = [key if element.lower() in key else None\
                  for key in eltkeys]
        eltgrp = list(set(eltgrp))
        eltgrp.remove(None)
        eltgrp = mgrp[eltgrp[0]]
        # redshift, T, Z, nH, ion
        ztab = f['TableBins/RedshiftBins'][:]
        logT = f['TableBins/TemperatureBins'][:]
        met = f['TableBins/MetallicityBins'][:]
        lognH = f['TableBins/DensityBins'][:]
        nHcut = 1.0
        nHcuti = np.min(np.where(lognH >= nHcut)[0])
        fontsize = 12
        nHedges = lognH[:-1] - 0.5 * np.diff(lognH)
        nHend = [lognH[-1] - 0.5 * (lognH[-1] - lognH[-2]),
                 lognH[-1] + 0.5 * (lognH[-1] - lognH[-2])
                 ]
        nHedges = np.append(nHedges, nHend)
        Tedges = logT[:-1] - 0.5 * np.diff(logT)
        Tend = [logT[-1] - 0.5 * (logT[-1] - logT[-2]),
                logT[-1] + 0.5 * (logT[-1] - logT[-2])
                ]
        Tedges = np.append(Tedges, Tend)
        
        for ztar in ztargets:
            Zstart = 0 if element in ['hydrogen', 'helium'] else 1

            ncols = 4
            nz = len(met) - Zstart
            nrows = ((nz - 1) // ncols + 1) * 2
            width = 11.
            width_ratios = [1.] * ncols + [0.1]
            wspace = 0.3
            panelwidth = width / (sum(width_ratios) + ncols * wspace)
            hspace = 0.4
            height = (hspace * (nrows - 1.) + nrows) * panelwidth 

            fig = plt.figure(figsize=(width, height))
            grid = gsp.GridSpec(nrows=nrows, ncols=ncols + 1, hspace=hspace,
                                wspace=wspace, 
                                width_ratios=width_ratios,
                                top=0.95, bottom=0.05)
            cax = fig.add_subplot(grid[:2, ncols])
            pdaxes = [fig.add_subplot(grid[2 * (i // 4), i % 4]) \
                      for i in range(nz)]
            haxes =  [fig.add_subplot(grid[2 * (i // 4) + 1, i % 4]) \
                       for i in range(nz)]

            zi = np.argmin(np.abs(ztab - ztar))
            zval = ztab[zi]
            _savename = savename.format(z=zval, elt=element)
            title = 'Sum of {elt} ions at z={z:.3f}'
            fig.suptitle(title.format(elt=element, z=zval), fontsize=fontsize)
            # T, Z, nH
            tabsum = np.sum(10**eltgrp[zi, :, Zstart:, :, :], axis=3)
            vmax = np.max(np.abs(tabsum - 1.))
            vmin = 1. - vmax
            vmax = 1. + vmax
            cmap = 'RdBu'
            #hbins = np.linspace(vmin, vmax, 100)

            nHlabel = '$\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}} \\;' + \
                      '[\\mathrm{cm}^{-3}]$'
            Tlabel = '$\\log_{10} \\, \\mathrm{T} \\; [\\mathrm{K}]$'
            flabel = 'sum of ion fractions'

            # first value is Z = 0.0 -> no metal ions 
            for _Zi, _Z in enumerate(met[Zstart:]):
                label = '$ \\log_{{10}}Z/Z_{{\\odot}}$\n={:.1f}'.format(_Z)
                pdax = pdaxes[_Zi]
                hax = haxes[_Zi]
                pdax.text(0.05, 0.95, label, fontsize=fontsize - 2, 
                          horizontalalignment='left', verticalalignment='top',
                          transform=pdax.transAxes)
                hax.text(0.05, 0.95, label, fontsize=fontsize - 2, 
                          horizontalalignment='left', verticalalignment='top',
                          transform=hax.transAxes)
                pdax.set_xlabel(nHlabel, fontsize=fontsize)
                hax.set_xlabel(flabel, fontsize=fontsize)
                if _Zi % ncols == 0: 
                    pdax.set_ylabel(Tlabel, fontsize=fontsize)
                    hax.set_ylabel('probability density', fontsize=fontsize)

                img = pdax.pcolormesh(nHedges, Tedges, 
                                      tabsum[:, _Zi, :], 
                                      vmin=vmin, vmax=vmax, cmap=cmap)
                pdax.axvline(nHcut, color='black', linestyle='dashed')

                nHhi = (tabsum[:, _Zi, nHcuti:]).flatten()
                nHlo = (tabsum[:, _Zi, :nHcuti]).flatten()
                
                hax.set_yscale('log')
                dmax = np.max(np.abs(nHhi - 1.))
                dmax = max(dmax, np.max(np.abs(nHlo - 1.)))
                hmin = 1. - dmax
                hmax = 1 + dmax
                hbins = np.linspace(hmin, hmax, 100)
                
                histlo, _ = np.histogram(nHlo, bins=hbins, density=True)
                histhi, _ = np.histogram(nHhi, bins=hbins, density=True)
                hax.set_yscale('log')
                hax.step(hbins[:-1], histlo, where='post', 
                         label='log nH >= {:.1f}'.format(nHcut),
                         color='blue')
                hax.step(hbins[:-1], histhi, where='post', 
                         label='log nH >= {:.1f}'.format(nHcut),
                         color='black', linestyle='dashed')

                #hv1, _, _ = hax.hist(nHlo, bins=hbins, density=True, 
                #                     label='log nH >= {:.1f}'.format(nHcut),
                #                     color='blue', align='mid', 
                #                     histtype='step')
                #hv2, _, _ = hax.hist(nHhi, bins=hbins, density=True, 
                #                     label='log nH < {:.1f}'.format(nHcut),
                #                     color='black', linestyle='dashed', 
                #                     align='mid', histtype='step')
                ymin = min(np.min(histlo[histlo > 0.]), 
                           np.min(histhi[histhi > 0.]))
                _ymin, _ymax = hax.get_ylim()
                hax.set_ylim((0.7 * ymin, _ymax))

                if _Zi == 0:
                    hax.legend(fontsize=fontsize - 2)
            
            plt.colorbar(img, cax=cax)
            cax.set_ylabel(flabel, fontsize=fontsize)
            plt.savefig(_savename, bbox_inches='tight')

def test_lin_option_ionbal(filen_old, filen_new_log, filen_new_lin):

    imgdir = '/'.join(filen_old.split('/')[:-1]) + '/'
    imgname_oldnew = 'ionbal_comp_test_ps20-logtable_ps20-old_{ion}' + \
                     '_{sim}_z-{z:.2f}_logZ-{met:.2f}_depletion-{dep}.pdf'
    imgname_linlog = 'ionbal_comp_test_ps20-log-lintable' + \
                     '_{sim}_z-{z:.2f}_logZ-{met:.2f}_depletion-{dep}.pdf'
    simlabel = 'testhalo1-m13h206_m3e5'

    logTsim = None
    lognHsim = None
    Zsim = None

    target_Z = None
    delta_Z = None
    ion = None
    ps20depletion = None
    
    redshift = None

    dct_old = {}
    dct_new_log = {}
    dct_new_lin = {}

    for filen, dct in zip([filen_old, filen_new_log, filen_new_lin],
                          [dct_old, dct_new_log, dct_new_lin]):
        with h5py.File(filen, 'r') as f:
            # title info
            cosmopars = {key: val for key, val in \
                        f['Header/cosmopars'].attrs.items()}
            if redshift is None:
                redshift = cosmopars['z']
            elif not redshift == cosmopars['z']:
                msg = 'Input files have different redshifts; found in {}'
                raise ValueError(msg.format(filen))
            _delta_Z = f['Header'].attrs['delta_Z']
            if delta_Z is None:
                delta_Z = _delta_Z
            elif not delta_Z == _delta_Z:
                msg = 'Input files have different delta_Z; found in {}'
                raise ValueError(msg.format(filen))
            _target_Z = f['Header'].attrs['target_Z']
            if target_Z is None:
                target_Z = _target_Z
            elif not target_Z == _target_Z:
                msg = 'Input files have different target_Z; found in {}'
                raise ValueError(msg.format(filen))
            _ion = f['Header'].attrs['ion'].decode()
            if ion is None:
                ion = _ion
            elif not ion == _ion:
                msg = 'Input files have different ion; found in {}'
                raise ValueError(msg.format(filen))
            _ps20depletion = bool(f['Header'].attrs['ps20depletion'])
            if ps20depletion is None:
                ps20depletion = _ps20depletion
            elif not ps20depletion == _ps20depletion:
                msg = 'Input files have different ps20depletion; found in {}'
                raise ValueError(msg.format(filen))
            if 'lintable' in f['Header'].attrs:
                dct['lintable'] = bool(f['Header'].attrs['lintable']) 
            else:
                dct['lintable'] = None
            
            # sim data
            if logTsim is None:
                logTsim = np.log10(f['simulation_data/T_K'])
            elif not np.all(logTsim == np.log10(f['simulation_data/T_K'])):
                msg = 'Input files have different simulation log T values; found in {}'
                raise ValueError(msg.format(filen))
            if lognHsim is None:
                lognHsim = np.log10(f['simulation_data/nH_cm**-3'])
            elif not np.all(lognHsim == np.log10(f['simulation_data/nH_cm**-3'])):
                msg = 'Input files have different simulation log nH values; found in {}'
                raise ValueError(msg.format(filen))
            if Zsim is None:
                Zsim = f['simulation_data/metallicity_abs_mass_frac'][:]
            elif not np.all(Zsim == f['simulation_data/metallicity_abs_mass_frac'][:]):
                msg = 'Input files have different simulation Z values; found in {}'
                raise ValueError(msg.format(filen))
            dct['ionsim'] = f['simulation_data/ionbal'][:]
    
    # new log vs. old version
    title = 'interpolated {ion} PS20 tables at z={z:.2f}, '+ \
            'Z={met:.1e} from FIRE data, dust depl. {dep}' 
    title = title.format(ion=ion, dep=ps20depletion, 
                         z=redshift, met=target_Z)
    
    fig = plt.figure(figsize=(13., 4.))
    grid = gsp.GridSpec(nrows=1, ncols=5, hspace=0.0, wspace=0.5, 
                        width_ratios=[1., 0.1, 1., 0.1, 1.],
                        left=0.05, right=0.95)
    axes = [fig.add_subplot(grid[0, i]) for i in range(5)]
    fontsize = 12
    cmap = 'viridis'
    size = 10.
    
    ionbal_old = dct_old['ionsim']
    ionbal_new_log = dct_new_log['ionsim']
    ionbal_new_lin = dct_new_lin['ionsim']
    vmax = 0.
    vmin = -10.

    fig.suptitle(title, fontsize=fontsize)

    ax = axes[0]
    cax = axes[1]
    img = ax.scatter(lognHsim, logTsim, s=size, c=np.log10(ionbal_old),
                     edgecolor='black', cmap=cmap, vmin=vmin, vmax=vmax,
                     rasterized=True)
    plt.colorbar(img, cax=cax)
    cax.set_ylabel('$\\log_{10}$ ion fraction', fontsize=fontsize)

    ax.set_ylabel('$\\log \\, \\mathrm{T} \\; [\\mathrm{K]}$', 
                    fontsize=fontsize)
    xlabel = '$\\log \\, \\mathrm{n}_{\\mathrm{H}} \\;' + \
                ' [\\mathrm{cm}^{-3}]$'
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.tick_params(which='both', axis='both', labelsize=fontsize - 1)
    cax.tick_params(labelsize=fontsize - 1, labelright=False, labelleft=True,
                    right=False, left=True)
    cax.yaxis.set_label_position('left')
    ax.text(0.05, 0.05, 'Old table interp.,\nlog space',
            horizontalalignment='left', verticalalignment='bottom',
            transform=ax.transAxes)

    ax = axes[2]
    cax = axes[3]
    delta_oldnew = np.log10(ionbal_new_log) - np.log10(ionbal_old)
    if np.all(ionbal_new_log == ionbal_old):
        print('Old, new log maps are same to fp precision')
        vmax = 1e-49
        vmin = -1. * vmax
        dohist = False
    else:
        vmax = np.max(np.abs(delta_oldnew))
        vmin = -1. * vmax
        dohist = True

    img = ax.scatter(lognHsim, logTsim, s=size, 
                     c=delta_oldnew,
                     edgecolor='black', cmap='RdBu', vmin=vmin, vmax=vmax,
                     rasterized=True)
    plt.colorbar(img, cax=cax)
    clabel = '$\\Delta_{\\mathrm{new, old}} \\, \\log_{10}$ ion fraction'
    cax.set_ylabel(clabel, fontsize=fontsize)

    ax.set_ylabel('$\\log \\, \\mathrm{T} \\; [\\mathrm{K]}$', 
                    fontsize=fontsize)
    xlabel = '$\\log \\, \\mathrm{n}_{\\mathrm{H}} \\;' + \
                ' [\\mathrm{cm}^{-3}]$'
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.tick_params(which='both', axis='both', labelsize=fontsize - 1)
    cax.tick_params(labelsize=fontsize - 1, labelright=False, labelleft=True,
                    right=False, left=True)
    cax.yaxis.set_label_position('left')
 
    ax = axes[4]
    
    if dohist: 
        nbins = 100
        maxv = np.max(np.abs(delta_oldnew))
        bins = np.linspace(-1. * maxv, maxv, nbins)

        ax.set_yscale('log')
        ax.hist(delta_oldnew, bins=bins, histtype='step', color='blue',
                label=None, linestyle='dashed', align='mid', 
                density=True)

        ax.set_xlabel(clabel, fontsize=fontsize)
        ax.set_ylabel('probability density', fontsize=fontsize)
        ax.legend(fontsize=fontsize - 1)
    else:
        ax.text(0.5, 0.5, 'Old, new log maps are\nthe same to fp precision',
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes)
        ax.tick_params(left=False, bottom=False, labelleft=False, 
                       labelbottom=False, which='both')
        
    
    outname = imgdir + imgname_oldnew.format(sim=simlabel, dep=ps20depletion,
                                             z=redshift, met=np.log10(target_Z),
                                             ion=ion)
    plt.savefig(outname, bbox_inches='tight')

    # new log vs. new lin
    title = 'interpolated {ion} PS20 tables at z={z:.2f}, '+ \
            'Z={met:.1e} from FIRE data, dust depl. {dep}' 
    title = title.format(ion=ion, dep=ps20depletion, 
                         z=redshift, met=target_Z)
    
    fig = plt.figure(figsize=(13., 4.))
    grid = gsp.GridSpec(nrows=1, ncols=5, hspace=0.0, wspace=0.5, 
                        width_ratios=[1., 0.1, 1., 0.1, 1.],
                        left=0.05, right=0.95)
    axes = [fig.add_subplot(grid[0, i]) for i in range(5)]
    fontsize = 12
    cmap = 'viridis'
    size = 10.
    
    vmax = 0.
    vmin = -10.

    fig.suptitle(title, fontsize=fontsize)

    ax = axes[0]
    cax = axes[1]
    img = ax.scatter(lognHsim, logTsim, s=size, c=np.log10(ionbal_new_lin),
                     edgecolor='black', cmap=cmap, vmin=vmin, vmax=vmax,
                     rasterized=True)
    plt.colorbar(img, cax=cax)
    cax.set_ylabel('$\\log_{10}$ ion fraction', fontsize=fontsize)

    ax.set_ylabel('$\\log \\, \\mathrm{T} \\; [\\mathrm{K]}$', 
                    fontsize=fontsize)
    xlabel = '$\\log \\, \\mathrm{n}_{\\mathrm{H}} \\;' + \
                ' [\\mathrm{cm}^{-3}]$'
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.tick_params(which='both', axis='both', labelsize=fontsize - 1)
    cax.tick_params(labelsize=fontsize - 1, labelright=False, labelleft=True,
                    right=False, left=True)
    cax.yaxis.set_label_position('left')
    ax.text(0.05, 0.05, 'New table interp.,\nlin space',
            horizontalalignment='left', verticalalignment='bottom',
            transform=ax.transAxes)

    ax = axes[2]
    cax = axes[3]
    delta_linlog = np.log10(ionbal_new_lin) - np.log10(ionbal_new_log)
    if np.all(ionbal_new_log == ionbal_new_lin):
        print('New log, lin maps are the same to fp precision')
        vmax = 1e-49
        dohist = False
    else:
        vmax = np.max(np.abs(delta_linlog[np.isfinite(delta_linlog)]))
        extend = 'neither'
        if vmax > 1.:
            vmax = 1.
            extend = 'both'
        vmin = -1. * vmax
        dohist = True
    img = ax.scatter(lognHsim, logTsim, s=size, 
                     c=delta_linlog,
                     edgecolor='black', cmap='RdBu', vmin=vmin, vmax=vmax,
                     rasterized=True)
    plt.colorbar(img, cax=cax, extend=extend)
    clabel = '$\\Delta_{\\mathrm{lin, log}} \\, \\log_{10}$ ion fraction'
    cax.set_ylabel(clabel, fontsize=fontsize)

    ax.set_ylabel('$\\log \\, \\mathrm{T} \\; [\\mathrm{K]}$', 
                    fontsize=fontsize)
    xlabel = '$\\log \\, \\mathrm{n}_{\\mathrm{H}} \\;' + \
                ' [\\mathrm{cm}^{-3}]$'
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.tick_params(which='both', axis='both', labelsize=fontsize - 1)
    cax.tick_params(labelsize=fontsize - 1, labelright=False, labelleft=True,
                    right=False, left=True)
    cax.yaxis.set_label_position('left')
 
    ax = axes[4]
    
    if dohist:
        nbins = 100
        maxv = np.max(np.abs(delta_linlog[np.isfinite(delta_linlog)]))
        bins = np.linspace(-1. * maxv, maxv, nbins)

        ax.set_yscale('log')
        ax.hist(delta_linlog, bins=bins, histtype='step', color='blue',
                label=None, linestyle='dashed', align='mid', 
                density=True)

        ax.set_xlabel(clabel, fontsize=fontsize)
        ax.set_ylabel('probability density', fontsize=fontsize)
        ax.legend(fontsize=fontsize - 1)
    else:
        ax.text(0.5, 0.5, 'New lin, log maps are\nthe same to fp precision',
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes)
        ax.tick_params(left=False, bottom=False, labelleft=False, 
                       labelbottom=False, which='both')

    outname = imgdir + imgname_linlog.format(sim=simlabel, dep=ps20depletion,
                                             z=redshift, met=np.log10(target_Z),
                                             ion=ion)
    plt.savefig(outname, bbox_inches='tight')
    
def run_test_lin_option_ionbal(index):
     # laptop
    ddir = '/Users/nastasha/ciera/tests/fire_start/ionbal_tests/'
    simname = 'm13h206_m3e5__m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1'+\
              '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'
    pre = 'ionbal_test_PS20'
    ftmp = ['{pre}_{ion}_depletion-False_Z-0.01_snap045{lt}_{simname}.hdf5',
            '{pre}_{ion}_depletion-False_Z-0.0001_snap027{lt}_{simname}.hdf5',
            '{pre}_{ion}_depletion-False_Z-0.01_snap027{lt}_{simname}.hdf5',
            '{pre}_{ion}_depletion-True_Z-0.0001_snap027{lt}_{simname}.hdf5',
            '{pre}_{ion}_depletion-True_Z-0.01_snap027{lt}_{simname}.hdf5',
            '{pre}_{ion}_depletion-True_Z-0.01_snap045{lt}_{simname}.hdf5',
            ]
    ions = ['O{}'.format(i) for i in range(1, 10)]
    lts = ['', '_lintable-False', '_lintable-True']

    ioni = index % len(ions)
    vari = index // len(ions)
    
    filen_old, filen_new_log, filen_new_lin = \
        (ddir + ftmp[vari].format(pre=pre, simname=simname, 
                                  ion=ions[ioni], lt=lt) \
         for lt in lts)
    test_lin_option_ionbal(filen_old, filen_new_log, filen_new_lin)


def test_tablesum_interpolate_to_tabulated(ztargets = [0., 1., 3.], 
                                           element='oxygen'):
    '''
    Ok, can it get the right values when interpolating to the
    listed values then.
    '''
    fn = '/Users/nastasha/phd/tables/ionbal/lines_sp20/'+\
         'UVB_dust1_CR1_G1_shield1.hdf5'
    outdir = '/Users/nastasha/ciera/tests/fire_start/ionbal_tests/'
    savename = outdir + 'sumtest-nodepl_z-{z:.3f}_{elt}_tables-direct.pdf'

    with h5py.File(fn, 'r') as f:
        mgrp = f['Tdep/IonFractions']
        eltkeys = list(mgrp.keys())
        eltgrp = [key if element.lower() in key else None\
                  for key in eltkeys]
        eltgrp = list(set(eltgrp))
        eltgrp.remove(None)
        eltgrp = mgrp[eltgrp[0]]
        # redshift, T, Z, nH, ion
        ztab = f['TableBins/RedshiftBins'][:]
        logT = f['TableBins/TemperatureBins'][:]
        met = f['TableBins/MetallicityBins'][:]
        lognH = f['TableBins/DensityBins'][:]
        nHcut = 1.0
        nHcuti = np.min(np.where(lognH >= nHcut)[0])
        fontsize = 12
        nHedges = lognH[:-1] - 0.5 * np.diff(lognH)
        nHend = [lognH[-1] - 0.5 * (lognH[-1] - lognH[-2]),
                 lognH[-1] + 0.5 * (lognH[-1] - lognH[-2])
                 ]
        nHedges = np.append(nHedges, nHend)
        Tedges = logT[:-1] - 0.5 * np.diff(logT)
        Tend = [logT[-1] - 0.5 * (logT[-1] - logT[-2]),
                logT[-1] + 0.5 * (logT[-1] - logT[-2])
                ]
        Tedges = np.append(Tedges, Tend)
        
        for ztar in ztargets:
            Zstart = 0 if element in ['hydrogen', 'helium'] else 1

            ncols = 4
            nz = len(met) - Zstart
            nrows = ((nz - 1) // ncols + 1) * 2
            width = 11.
            width_ratios = [1.] * ncols + [0.1]
            wspace = 0.3
            panelwidth = width / (sum(width_ratios) + ncols * wspace)
            hspace = 0.4
            height = (hspace * (nrows - 1.) + nrows) * panelwidth 

            fig = plt.figure(figsize=(width, height))
            grid = gsp.GridSpec(nrows=nrows, ncols=ncols + 1, hspace=hspace,
                                wspace=wspace, 
                                width_ratios=width_ratios,
                                top=0.95, bottom=0.05)
            cax = fig.add_subplot(grid[:2, ncols])
            pdaxes = [fig.add_subplot(grid[2 * (i // 4), i % 4]) \
                      for i in range(nz)]
            haxes =  [fig.add_subplot(grid[2 * (i // 4) + 1, i % 4]) \
                       for i in range(nz)]

            zi = np.argmin(np.abs(ztab - ztar))
            zval = ztab[zi]
            _savename = savename.format(z=zval, elt=element)
            title = 'Sum of {elt} ions at z={z:.3f}'
            fig.suptitle(title.format(elt=element, z=zval), fontsize=fontsize)
            # T, Z, nH
            tabsum = np.sum(10**eltgrp[zi, :, Zstart:, :, :], axis=3)
            vmax = np.max(np.abs(tabsum - 1.))
            vmin = 1. - vmax
            vmax = 1. + vmax
            cmap = 'RdBu'
            #hbins = np.linspace(vmin, vmax, 100)

            nHlabel = '$\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}} \\;' + \
                      '[\\mathrm{cm}^{-3}]$'
            Tlabel = '$\\log_{10} \\, \\mathrm{T} \\; [\\mathrm{K}]$'
            flabel = 'sum of ion fractions'

            # first value is Z = 0.0 -> no metal ions 
            for _Zi, _Z in enumerate(met[Zstart:]):
                label = '$ \\log_{{10}}Z/Z_{{\\odot}}$\n={:.1f}'.format(_Z)
                pdax = pdaxes[_Zi]
                hax = haxes[_Zi]
                pdax.text(0.05, 0.95, label, fontsize=fontsize - 2, 
                          horizontalalignment='left', verticalalignment='top',
                          transform=pdax.transAxes)
                hax.text(0.05, 0.95, label, fontsize=fontsize - 2, 
                          horizontalalignment='left', verticalalignment='top',
                          transform=hax.transAxes)
                pdax.set_xlabel(nHlabel, fontsize=fontsize)
                hax.set_xlabel(flabel, fontsize=fontsize)
                if _Zi % ncols == 0: 
                    pdax.set_ylabel(Tlabel, fontsize=fontsize)
                    hax.set_ylabel('probability density', fontsize=fontsize)

                img = pdax.pcolormesh(nHedges, Tedges, 
                                      tabsum[:, _Zi, :], 
                                      vmin=vmin, vmax=vmax, cmap=cmap)
                pdax.axvline(nHcut, color='black', linestyle='dashed')

                nHhi = (tabsum[:, _Zi, nHcuti:]).flatten()
                nHlo = (tabsum[:, _Zi, :nHcuti]).flatten()
                
                hax.set_yscale('log')
                dmax = np.max(np.abs(nHhi - 1.))
                dmax = max(dmax, np.max(np.abs(nHlo - 1.)))
                hmin = 1. - dmax
                hmax = 1 + dmax
                hbins = np.linspace(hmin, hmax, 100)
                
                histlo, _ = np.histogram(nHlo, bins=hbins, density=True)
                histhi, _ = np.histogram(nHhi, bins=hbins, density=True)
                hax.set_yscale('log')
                hax.step(hbins[:-1], histlo, where='post', 
                         label='log nH >= {:.1f}'.format(nHcut),
                         color='blue')
                hax.step(hbins[:-1], histhi, where='post', 
                         label='log nH >= {:.1f}'.format(nHcut),
                         color='black', linestyle='dashed')

                #hv1, _, _ = hax.hist(nHlo, bins=hbins, density=True, 
                #                     label='log nH >= {:.1f}'.format(nHcut),
                #                     color='blue', align='mid', 
                #                     histtype='step')
                #hv2, _, _ = hax.hist(nHhi, bins=hbins, density=True, 
                #                     label='log nH < {:.1f}'.format(nHcut),
                #                     color='black', linestyle='dashed', 
                #                     align='mid', histtype='step')
                ymin = min(np.min(histlo[histlo > 0.]), 
                           np.min(histhi[histhi > 0.]))
                _ymin, _ymax = hax.get_ylim()
                hax.set_ylim((0.7 * ymin, _ymax))

                if _Zi == 0:
                    hax.legend(fontsize=fontsize - 2)
            
            plt.colorbar(img, cax=cax)
            cax.set_ylabel(flabel, fontsize=fontsize)
            plt.savefig(_savename, bbox_inches='tight')
