import h5py
import numpy as np

import matplotlib.gridspec as gsp
import matplotlib.patheffects as mppe 
import matplotlib.pyplot as plt

import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c


def readin_data(filen):
    with h5py.File(filen) as f:
        hist = 10**f['histogram/histogram'][:]
        rbins_rvir = f['axis_0/bins'][:]
        cosmopars = {key: val for key, val
                     in f['Header/cosmopars'].attrs.items()}
        halomass = f['Header/inputpars/halodata'].attrs['Mvir_g']
        halomass_msun = halomass / c.solar_mass
    return rbins_rvir, hist, cosmopars, halomass_msun

def plotfracs_haloes(filen_template, fills_sim, title=None, outname=None,
                     rmin_rvir=0.1, rmax_rvir=1.0):
    ions = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'Carbon']
    ionnames = ['C I', 'C II', 'C III', 'C IV', 'C V', 'C VI']
    colors = tc.tol_cset('bright')[:len(ions) - 1]
    m11list = []
    m12list = []
    m13list = []
    
    for sfill in fills_sim:
        masses = {}
        for ion in ions:
            filen = filen_template.format(ion=ion, **sfill)
            rbins_rvir, hist, cosmopars, halomass_msun = readin_data(filen)
            print(cosmopars['z'])
            stag = sfill['simname']
            stag = stag.split('_')[0]
            imin = np.where(np.isclose(rmin_rvir, rbins_rvir))[0]
            if len(imin) == 0:
                print(f'No bin match found for radius {rmin_rvir} Rvir')
            else:
                imin = imin[0]
            imax = np.where(np.isclose(rmax_rvir, rbins_rvir))[0]
            if len(imax) == 0:
                print(f'No bin match found for radius {rmax_rvir} Rvir')
            else:
                imax = imax[0]
            mass = np.sum(hist[imin : imax])
            masses[ion] = mass
        if stag.startswith('m11'):
            _ls = m11list
        elif stag.startswith('m12'):
            _ls = m12list
        elif stag.startswith('m13'):
            _ls = m13list
        #print(stag)
        _ls.append({'masses': masses, 'halomass_msun': halomass_msun, 
                    'stag': stag})
    m11list.sort(key=lambda x: x['halomass_msun'])
    m12list.sort(key=lambda x: x['halomass_msun'])
    m13list.sort(key=lambda x: x['halomass_msun'])
    # None: skip leave a gap between m11/m12/m13 haloes
    halolist = m11list + [None] + m12list + [None] + m13list

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(11., 3.))
    fontsize = 12
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    
    tickposs = []
    ticklabels = []
    fracss = []
    for xpos, halodat in enumerate(halolist):
        if halodat is None:
            continue
        xpos = xpos + 0.5
        tickposs.append(xpos)
        tl = halodat['stag'] + f' ({np.log10(halodat["halomass_msun"]):.1f})'
        ticklabels.append(tl)

        fracs = [halodat['masses'][ion] / halodat['masses'][ions[-1]]
                 for ion in ions[:-1]]
        #print(fracs)
        fracss.append(fracs)
    fracss = np.array(fracss)
    for ii in range(len(ions) - 1):
        ax.bar(tickposs, fracss[:, ii], color=colors[ii], edgecolor='black',
               bottom=np.sum(fracss[:, :ii], axis=1))
    ax.set_xticks(tickposs, labels=ticklabels, ha='right')
    nnames = len(ions) - 1
    textoutline = [mppe.Stroke(linewidth=1.5, foreground='black'),
                   mppe.Normal()]
    for ii in range(nnames):
        ax.text(tickposs[-1] + 0.5, (ii + 0.5) / (nnames + 1), ionnames[ii],
                fontsize=fontsize, color=colors[ii],
                horizontalalignment='left', verticalalignment='center',
                path_effects=textoutline)
    ax.tick_params(which='both', axis='x', labelsize=fontsize, direction='out',
                   labelrotation=45.)
    ax.tick_params(which='both', axis='y', labelsize=fontsize - 1, 
                   direction='in', right=True)
    ax.set_ylim((0., 1.))
    ylab = (f'fraction of {ions[-1]}, '
            f'${rmin_rvir:.1f} \\endash {rmax_rvir:.1f}'
            f' \\, \\mathrm{{R}}_{{\\mathrm{{vir}}}}$')
    ax.set_ylabel(ylab, fontsize=fontsize)

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def phys_from_simn(simname):
    if '_sdp1e10_' in simname:
        return 'noBH'
    elif '_MHDCRspec1_' in simname:
        return 'AGN-CR'
    else:
        return 'AGN-noCR'

def plotsetfracs_haloes():
    ddir = '/Users/nastasha/ciera/profiles/fire/ionseries_C/'
    ftemp = 'hist_r3D_by_{ion}_{simname}_snap{snap}_bins1_v1.hdf5'
    filen_template = ddir + ftemp
    outdir = '/Users/nastasha/ciera/projects_co/hsiao-wen_carbon_ions/'
    rmin_rvir=0.1
    rmax_rvir=1.0
    for zi, redshift in enumerate([1.0, 0.5, 0.0]):
        snap_hr = sl.snaps_hr_051[zi]
        snap_sr = sl.snaps_sr_051[zi]
        
        simns_m11_hr = sl.m11_hr_agncr_set1 +\
                       sl.m11_hr_agnnocr_set1 +\
                       sl.m11_hr_nobh_set1 
        simns_m11_sr = sl.m11_sr_agncr_set1 \
                       + sl.m11_sr_agnnocr_set1 \
                       + sl.m11_sr_nobh_set1
        if redshift == 0.:
            simns_hi_hr = sl.m12_hr_all2_z0 + sl.m13_hr_all2_z0
            simns_hi_sr = sl.m12_sr_all2_z0 + sl.m13_sr_all2_z0
        else:
            simns_hi_hr = sl.m12_hr_all2 + sl.m13_hr_all2
            simns_hi_sr = sl.m12_sr_all2 + sl.m13_sr_all2
        simns_hr = simns_hi_hr + simns_m11_hr
        simns_sr = simns_hi_sr + simns_m11_sr
        
        for phys in ['noBH', 'AGN-noCR', 'AGN-CR']:
            msg = f'plotting z={redshift} (snap {snap_hr}/{snap_sr}), {phys}'
            print(msg)
            _simns_hr = [sim for sim in simns_hr 
                         if phys_from_simn(sim) == phys]
            _simns_sr = [sim for sim in simns_sr 
                         if phys_from_simn(sim) == phys]
            fills_sim = [{'snap': snap_sr, 'simname': simn} 
                         for simn in _simns_sr]
            fills_sim += [{'snap': snap_hr, 'simname': simn} 
                          for simn in _simns_hr]
            title = f'CGM carbon ion fractions: {phys}, $z={redshift:.1f}$'
            outname = (f'ionfrac_C_CGM_{rmin_rvir:.2f}_to_{rmax_rvir:.2f}'
                       f'_RBN98_{phys}_z{redshift:.1f}.pdf')
            outname = outdir + outname
            plotfracs_haloes(filen_template, fills_sim, title=title, 
                             outname=outname, rmin_rvir=rmin_rvir, 
                             rmax_rvir=rmax_rvir)

def plotradfracs_haloes(filen_template, fills_sim, ics_av,
                        title=None, outname=None, rmax_rvir=2.0):
    ions = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6'][::-1] + ['Carbon']
    ionnames = ['C I', 'C II', 'C III', 'C IV', 'C V', 'C VI'][::-1]
    colors = tc.tol_cset('bright')[len(ions) - 2 : : -1]
    print(colors)
    
    data = []
    rcens_all = None
    for sfill in fills_sim:
        masses = {}
        for ion in ions:
            filen = filen_template.format(ion=ion, **sfill)
            rbins_rvir, hist, cosmopars, halomass_msun = readin_data(filen)
            stag = sfill['simname']
            stag = stag.split('_')[0]
            imax = np.where(np.isclose(rmax_rvir, rbins_rvir))[0]
            if len(imax) == 0:
                print(f'No bin match found for radius {rmax_rvir} Rvir')
            else:
                imax = imax[0]
            mass = hist[: imax]
            masses[ion] = mass
            rbins_rvir = rbins_rvir[:imax + 1]
            rcens_rvir = 0.5 * (rbins_rvir[:-1] + rbins_rvir[1:])
            if rcens_all is None:
                rcens_all = rcens_rvir
            elif not np.allclose(rcens_all, rcens_rvir):
                print(f'Issue in file {filen}')
                raise RuntimeError('Radial bins do not match between files')
        data.append({'masses': masses, 'halomass_msun': halomass_msun, 
                     'stag': stag})

    fig = plt.figure(figsize=(11., 5.5))
    grid1 = gsp.GridSpec(nrows=1, ncols=2,
                         hspace=0., wspace=0.2,
                         width_ratios=[1., 1.])
    avax = fig.add_subplot(grid1[0, 0])
    frameax = fig.add_subplot(grid1[0, 1])
    grid = gsp.GridSpecFromSubplotSpec(nrows=3, ncols=2, hspace=0., wspace=0.,
                                       subplot_spec=grid1[0, 1])
    axes = [fig.add_subplot(grid[i // 2, i % 2]) for i in range(6)]

    fontsize = 12
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    xlab = '$r_{\\mathrm{3D}} \\,/\\, \\mathrm{R}_{\\mathrm{vir}}$'
    ylab = f'ion mass / {ions[-1]} mass'
    frameax.tick_params(bottom=False, left=False, labelbottom=False, 
                        labelleft=False)
    frameax.spines['top'].set_visible(False)
    frameax.spines['bottom'].set_visible(False)
    frameax.spines['left'].set_visible(False)
    frameax.spines['right'].set_visible(False)
    frameax.set_xlabel(xlab, fontsize=fontsize, labelpad=20)
    frameax.set_ylabel(ylab, fontsize=fontsize, labelpad=28)
    #avax.set_xscale('log')
    avax.set_xlabel(xlab, fontsize=fontsize)
    avax.set_ylabel(ylab, fontsize=fontsize)
    avax.tick_params(which='both', labelsize=fontsize-1, 
                     top=True, right=True, direction='in',
                     labelleft=True, labelbottom=True)

    bottom = np.zeros(len(rcens_all))
    for ii, ion in enumerate(ions[:-1]):
        ax = axes[ii]
        #ax.set_xscale('log')
        ax.set_yscale('log')
        ax.tick_params(which='both', labelsize=fontsize-1, 
                       top=True, right=True, direction='in',
                       labelleft=(ii % 2 == 0), 
                       labelbottom=(ii // 2 >= 2))
        color = colors[ii]
        print(ii, ion)
        
        _avparts = []
        for hd in data:
            ionf = hd['masses'][ion] / hd['masses'][ions[-1]]
            if hd['stag'] in ics_av:
                _avparts.append(ionf)
                lw = 1.0
                _c = 'black'
            else:
                lw = 1.0
                _c = tc.tol_cset('bright').grey
            ax.plot(rcens_all, ionf, linewidth=lw, color=_c,
                    linestyle='solid')
        _av = np.average(_avparts, axis=0)
        _min = np.min(_avparts, axis=0)
        _max = np.max(_avparts, axis=0)
        
        pe = pu.getoutline(1.5)
        ax.errorbar(rcens_all, _av, yerr=(_av - _min, _max - _av),
                    color=color, linewidth=1.5, linestyle='solid',
                    errorevery=13, capsize=2., path_effects=pe)
        avax.plot(rcens_all, bottom + _av,
                  color=color, linewidth=1.5, linestyle='solid')
        #avax.errorbar(rcens_all, bottom + _av, yerr=(_av - _min, _max - _av),
        #              color=color, linewidth=2., linestyle='solid')
        avax.fill_between(rcens_all, bottom, bottom + _av,
                          color=color, alpha=0.3)
        bottom += _av
            
        if _av[-1] <= 0.5:
            ytpos = 0.95
            va = 'top'
        else:
            ytpos = 0.05
            va = 'bottom'
        textoutline = [mppe.Stroke(linewidth=1.5, foreground='black'),
                       mppe.Normal()]
        ax.text(0.95, ytpos, ionnames[ii], fontsize=fontsize,
                horizontalalignment='right', verticalalignment=va,
                transform=ax.transAxes, path_effects=textoutline, color=color)
    ylims = [ax.get_ylim() for ax in axes]
    ymax = np.max([yl[1] for yl in ylims])
    ymax = min(ymax, 1.3)
    ymin = max(1e-5, ymax * 3e-5)
    [ax.set_ylim((ymin, ymax)) for ax in axes]
    avlim = avax.get_ylim()
    avax.set_ylim(0., avlim[1])

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotsetradfracs_haloes():
    ddir = '/Users/nastasha/ciera/profiles/fire/ionseries_C/'
    ftemp = 'hist_r3D_by_{ion}_{simname}_snap{snap}_bins1_v1.hdf5'
    filen_template = ddir + ftemp
    outdir = '/Users/nastasha/ciera/projects_co/hsiao-wen_carbon_ions/'
    rmax_rvir=1.0
    for zi, redshift in enumerate([1.0, 0.5, 0.0]):
        snap_hr = sl.snaps_hr_051[zi]
        snap_sr = sl.snaps_sr_051[zi]
        
        simns_m11_hr = sl.m11_hr_agncr_set1 +\
                       sl.m11_hr_agnnocr_set1 +\
                       sl.m11_hr_nobh_set1 
        simns_m11_sr = sl.m11_sr_agncr_set1 \
                       + sl.m11_sr_agnnocr_set1 \
                       + sl.m11_sr_nobh_set1
        if redshift == 0.:
            simns_hi_hr = sl.m12_hr_all2_z0 + sl.m13_hr_all2_z0
            simns_hi_sr = sl.m12_sr_all2_z0 + sl.m13_sr_all2_z0
        else:
            simns_hi_hr = sl.m12_hr_all2 + sl.m13_hr_all2
            simns_hi_sr = sl.m12_sr_all2 + sl.m13_sr_all2
        simns_hr = simns_hi_hr + simns_m11_hr
        simns_sr = simns_hi_sr + simns_m11_sr
        
        for phys in ['noBH', 'AGN-noCR', 'AGN-CR']:
            for icset in ['m11', 'm12', 'm13']:
                msg = (f'plotting z={redshift} (snap {snap_hr}/{snap_sr})'
                       f', {phys}, {icset}')
                print(msg)
                _simns_hr = [sim for sim in simns_hr 
                             if phys_from_simn(sim) == phys 
                             and (sim.split('/')[0]).startswith(icset)]
                _simns_sr = [sim for sim in simns_sr 
                             if phys_from_simn(sim) == phys
                             and (sim.split('/')[0]).startswith(icset)]
                fills_sim = [{'snap': snap_sr, 'simname': simn} 
                              for simn in _simns_sr]
                fills_sim += [{'snap': snap_hr, 'simname': simn} 
                              for simn in _simns_hr]
                title = (f'carbon ion fraction profiles: {icset}'
                         f', {phys}, $z={redshift:.1f}$')
                outname = (f'ionfracprof_C_to_{rmax_rvir:.2f}'
                           f'_RBN98_{icset}_{phys}_z{redshift:.1f}.pdf')
                if icset == 'm11':
                    ics_av = ['m11a', 'm11b', 'm11d', 'm11e', 'm11f', 'm11g',
                              'm11h', 'm11i', 'm11q', 'm11v']
                elif icset == 'm12':
                    if redshift == 0.:
                        ics_av = ['m12f', 'm12i', 'm12m']
                    else:
                        ics_av = ['m12f', 'm12i', 'm12m', 'm12q']
                elif icset == 'm13':
                    if redshift == 0.:
                        # no m13 AGN-noCR halos reached z=0, so only
                        # require AGN-CR and noBH runs
                        ics_av = ['m13h007', 'm13h206']
                    else:
                        ics_av = ['m13h113', 'm13h206']
                if len(fills_sim) == 0:
                    print('skipping; no haloes')
                    continue
                outname = outdir + outname
                plotradfracs_haloes(filen_template, fills_sim, ics_av,
                                    title=title, outname=outname, rmax_rvir=2.0)
