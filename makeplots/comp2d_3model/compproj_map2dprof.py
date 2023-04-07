
import h5py
import matplotlib.collections as mcol
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import numpy as np

import fire_an.makeplots.get_2dprof as g2d
import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.math_utils as mu


def plotcomp_projax(filen_template, qtyfills, paxfills,
                    qtylabels=None, qtyclabels=None, title=None,
                    outname=None):
    '''
    map/profile comp: for each snapshot+phys+ic
    maps of gas, Ne, 3 ions (columns) in 3 projections (rows)
    and at the bottom, profiles for 3 projections and total

    filen_template: str
        a string useable with .format to get the file names.
        should include the full path.
    qtyfills: list of dicts
        what to fill in filen_template to get the different ion,
        mass, and metal profiles. keywords should match filen_template.
    paxfills: list of dicts
        like qtyfills, but here it's what to fill in to get the 
        different axis projections.
    qtylabels: list of strings
        labels to use for the columns (projected quantities)
    qtyclabels: list of strings
        color bar labels for the columns
    title: str or None
        figure title (none if not given)
    '''
    
    rvirs_all = []
    mvirs_all = []
    maps_all = []
    extents_all = []
    vmins_all = []
    vmaxs_all = []
    paxlabels_all = []
    filens_all = []
    xlabels_all = []
    ylabels_all = []
    for qfill in qtyfills:
        _kwa = qfill.copy()
        _rv = []
        _mv = []
        _mp = []
        _xt = []
        _vn = []
        _vx = []
        _pl = []
        _fn = []
        _xl = []
        _yl = []
        for afill in paxfills:
            kwa = _kwa.copy()
            kwa.update(afill)
            filen = filen_template.format(**kwa)
            _fn.append(filen)
            with h5py.File(filen, 'r') as f:
                map = f['map'][:]
                vmin = f['map'].attrs['minfinite']
                vmax = f['map'].attrs['max']

                box_cm = f['Header/inputpars'].attrs['diameter_used_cm']
                cosmopars = {key: val for key, val in \
                            f['Header/inputpars/cosmopars'].attrs.items()}
                #print(cosmopars)
                if 'Rvir_ckpcoverh' in f['Header/inputpars/halodata'].attrs:
                    h5path = 'Header/inputpars/halodata'
                    rvir_ckpcoverh = f[h5path].attrs['Rvir_ckpcoverh']
                    rvir_pkpc = rvir_ckpcoverh * cosmopars['a'] \
                                / cosmopars['h']
                elif 'Rvir_cm' in f['Header/inputpars/halodata'].attrs:
                    rvir_cm = f['Header/inputpars/halodata'].attrs['Rvir_cm']
                    rvir_pkpc = rvir_cm / (c.cm_per_mpc * 1e-3)
                xax = f['Header/inputpars'].attrs['Axis1']
                yax = f['Header/inputpars'].attrs['Axis2']
                zax = f['Header/inputpars'].attrs['Axis3']
                mvir_msun = f['Header/inputpars/halodata'].attrs['Mvir_g']
            mvir_msun /= c.solar_mass
            box_pkpc = box_cm / (1e-3 * c.cm_per_mpc)
            extent = (-0.5 * box_pkpc[xax], 0.5 * box_pkpc[xax],
                      -0.5 * box_pkpc[yax], 0.5 * box_pkpc[yax])
            paxlabel = 'XYZ'[zax] + '-proj.'
            xlabel = 'XYZ'[xax] + ' [pkpc]'
            ylabel = 'XYZ'[yax] + ' [pkpc]'
            
            _mp.append(map)
            _vn.append(vmin)
            _vx.append(vmax)
            _xt.append(extent)
            _rv.append(rvir_pkpc)
            _mv.append(mvir_msun)
            _pl.append(paxlabel)
            _xl.append(xlabel)
            _yl.append(ylabel)

            #maptype = f['Header/inputpars'].attrs['maptype'].decode()
            #if maptype == 'Mass':
            #    ion_used = 'Mass'
            #elif maptype == 'Metal':
            #    pathn = 'Header/inputpars/maptype_args_dict'
            #    ion_used = f[pathn].attrs['element'].decode()
            #elif maptype == 'ion':
            #    pathn = 'Header/inputpars/maptype_args_dict'
            #    ion_used = f[pathn].attrs['ion'].decode()
        rvirs_all.append(_rv)
        mvirs_all.append(_mv)
        maps_all.append(_mp)
        extents_all.append(_xt)
        vmins_all.append(_vn)
        vmaxs_all.append(_vx)
        paxlabels_all.append(_pl)
        filens_all.append(_fn)
        xlabels_all.append(_xl)
        ylabels_all.append(_yl)
    #print(extents_all)
    #print(rvirs_all)
    minrvir = np.min([np.min(rvl) for rvl in rvirs_all])
    rbins = np.linspace(0., 2. * minrvir, 50.)
    rc = 0.5 * (rbins[:-1] + rbins[1:])

    vmaxs = [np.max(_vx) for _vx in vmaxs_all]
    vmins = [np.max(_vn) for _vn in vmins_all]
    vmins = [max(vmins[i], vmaxs[i] - 5.) for i in range(len(vmins))]
    cmap = 'plasma'
    fontsize = 12

    ncols = len(qtyfills)
    nrows = len(paxfills)
    panelsize = 2.5
    hspace = 0.4
    wspace = 0.4
    cheight = 0.4

    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize] * nrows + [cheight, panelsize]
    width = sum(width_ratios) \
            * (1. + (len(width_ratios) - 1.) / len(width_ratios) * wspace)
    height = sum(height_ratios) \
             * (1. + (len(height_ratios) - 1.) / len(height_ratios) * hspace)
    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(nrows=nrows + 2, ncols=ncols, hspace=hspace, 
                        wspace=wspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    maxes = [[fig.add_subplot(grid[i, j]) \
              for i in range(nrows)] for j in range(ncols)]
    caxes = [fig.add_subplot(grid[nrows, j]) \
              for j in range(ncols)]
    paxes = [fig.add_subplot(grid[nrows + 1, j]) \
              for j in range(ncols)]
    
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    cset = tc.tol_cset('vibrant')
    colors = [cset.blue, cset.teal, cset.orange, cset.red, cset.magenta]
    for col in range(ncols):
        _rvs = rvirs_all[col]
        _mvs = mvirs_all[col]
        _xts = extents_all[col]
        _pls = paxlabels_all[col]
        _xls = xlabels_all[col]
        _yls = ylabels_all[col]
        _fns = filens_all[col]
        _mps = maps_all[col]
        _axs = maxes[col]
        _pax = paxes[col]
        _cax = caxes[col]

        vmin = vmins[col]
        vmax = vmaxs[col]
        qtylab = qtylabels[col]
        qtyclab = qtyclabels[col]

        for i, (_ax, _mp, _xt) in enumerate(zip(_axs, _mps, _xts)):
            img = _ax.imshow(_mp.T, origin='lower', interpolation='nearest',
                             cmap=cmap, vmin=vmin, vmax=vmax,
                             rasterized=True, extent=_xt)
            
            patches = [mpatch.Circle((0., 0.), _rvs[i])]
            collection = mcol.PatchCollection(patches)
            collection.set(edgecolor=['green'], facecolor='none', 
                           linewidth=1.5)
            _ax.add_collection(collection)

            _ax.set_xlabel(_xls[i], fontsize=fontsize)
            _ax.set_ylabel(_yls[i], fontsize=fontsize)
            _ax.tick_params(labelsize=fontsize - 1)
            if i == 0 and col == 0:
                _ax.text(1.05 * 2**-0.5 * _rvs[i], 
                         1.05 * 2**-0.5 * _rvs[i], 
                         '$R_{\\mathrm{vir}}$',
                         color='green', fontsize=fontsize)
                mvir = np.log10(_mvs[i])
                mvlab = ('$\\log_{10} \\, \\mathrm{M}_{\\mathrm{vir}}'
                         '\\, / \\, \\mathrm{M}_{\\mathrm{\\odot}} = '
                         f'{mvir:.1f}$')
                _ax.text(0.02, 0.98, mvlab, fontsize=fontsize - 1, 
                         color='green',
                         horizontalalignment='left', verticalalignment='top',
                         transform=_ax.transAxes)
                _ax.text(0.02, 0.02, f'$z={cosmopars["z"]:.1f}$',
                         fontsize=fontsize - 1, 
                         color='green',
                         horizontalalignment='left', 
                         verticalalignment='bottom',
                         transform=_ax.transAxes)
            if i == 0:
                _ax.set_title(qtylab, fontsize=fontsize)
            if col == ncols - 1:
                _ax.text(1.05, 0.5, _pls[i], fontsize=fontsize,
                         transform=_ax.transAxes, horizontalalignment='left',
                         verticalalignment='center', rotation=90)
        plt.colorbar(img, cax=_cax, orientation='horizontal',
                     aspect=0.1)
        _cax.set_xlabel(qtyclab, fontsize=fontsize)
        
        plo, pmed, phi = g2d.get_profile_massmap(_fns, rbins, 
            rbin_units='pkpc', profiles=['perc-0.1', 'perc-0.5', 'perc-0.9'])
        _pax.plot(rc, pmed, color='black', linewidth=2., label='all')
        _pax.plot(rc, plo, color='black', linewidth=2., linestyle='dashed')
        _pax.plot(rc, phi, color='black', linewidth=2., linestyle='dashed')
        for fi, filen in enumerate(_fns):
            plo, pmed, phi = g2d.get_profile_massmap(filen, rbins, 
                rbin_units='pkpc', 
                profiles=['perc-0.1', 'perc-0.5', 'perc-0.9'])
            _pax.plot(rc, pmed, color=colors[fi], linewidth=1.5,
                      label=_pls[fi])
            _pax.fill_between(rc, plo, phi, color=colors[fi],
                              alpha=0.3)
        _pax.set_xlabel('$\\mathrm{r}_{\\perp} \\, [\\mathrm{pkpc}]$',
                        fontsize=fontsize)
        _pax.set_ylabel(qtyclab, fontsize=fontsize)
        if col == 0:
            _pax.legend(fontsize=fontsize - 1, handlelength=1.)
    
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')


def plotcompset_projax(mapset='set3_model3'):
    
    if mapset == 'set3_model3':
        # quest
        outdir = '/projects/b1026/nastasha/imgs/2dcomp_set3_model3/'
        fdir = '/projects/b1026/nastasha/maps/set3_model3/'
        ftemp = ('coldens_{{qty}}_{simname}_snap{snapnum}_'
                 'shrink-sph-cen_BN98_2rvir_{{pax}}-proj_v3.hdf5')
        simnames = [sim for sim in sl.m12_hr_all1 
                    for i in range(len(sl.snaplists['m12_hr']))] \
                   + [sim for sim in sl.m12_sr_all1 
                      for i in range(len(sl.snaplists['m12_sr']))] \
                   + [sim for sim in sl.m13_hr_all1 
                      for i in range(len(sl.snaplists['m13_hr']))] \
                   + [sim for sim in sl.m13_sr_all1 
                      for i in range(len(sl.snaplists['m13_sr']))]

        snapnums = sl.snaplists['m12_hr'] * len(sl.m12_hr_all1) \
                   + sl.snaplists['m12_sr'] * len(sl.m12_sr_all1) \
                   + sl.snaplists['m13_hr'] * len(sl.m13_hr_all1) \
                   + sl.snaplists['m13_sr'] * len(sl.m13_sr_all1)
        _qfills = ['gas-mass', 'Neon', 'Ne8', 'O6', 'Mg10']
        qtyfillss = [[{'qty': val} for val in _qfills]] * len(simnames)
        _afills = ['x', 'y', 'z']
        paxfillss = [[{'pax': val} for val in _afills]] * len(simnames)
        _qtylab = ['Gas', 'Neon', 'Ne VIII', 'O VI', 'Mg X']
        qtylabelss = [_qtylab] * len(simnames)
        _qtyclab = [('$\\log_{10} \\, \\Sigma(\\mathrm{gas})'
                     ' \\; [\\mathrm{g}\\,\\mathrm{cm}^{-2}]$'),
                    ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Neon})'
                     ' \\; [\\mathrm{cm}^{-2}]$'),
                    ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\, VIII})'
                     ' [\\mathrm{cm}^{-2}]$'),
                    ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{O\\, VI})'
                     ' [\\mathrm{cm}^{-2}]$'),
                    ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Mg\\, X})'
                     ' [\\mathrm{cm}^{-2}]$')
                   ]
        qtyclabelss = [_qtyclab] * len(simnames)
        outname = ('comp_map_2dprof_projax_gas_Ne_Ne8_O6_Mg10_{ic}'
                   '_{phys}_{snap}.pdf')
    
    for simname, snapnum, qtyfills, paxfills, qtyclabels, qtylabels \
            in zip(simnames, snapnums, qtyfillss, paxfillss, 
                   qtyclabelss, qtylabelss):
        filen_template = fdir + ftemp.format(simname=simname, snapnum=snapnum)
        ic = simname.split('_')[0]
        physmodel = 'AGN-CR' if 'MHDCRspec1' in simname \
                    else 'noBH' if 'sdp1e10' in simname \
                    else 'AGN-noCR'
        title = f'{ic} {physmodel}, snapshot {snapnum}'
        if simname in sl.buglist1:
            title = title + ', possible bug'
        title = title + '\n' + simname
        _outname = outdir + outname.format(ic=ic, phys=physmodel, 
                                           snap=snapnum)
        
        plotcomp_projax(filen_template, qtyfills, paxfills,
                        qtylabels=qtylabels, qtyclabels=qtyclabels, 
                        title=title, outname=_outname)
        plt.close() # avoid using too much memory

# to avoid the previous 'annoying long nested loop' issue
def gethalodct_z_pax_phys(temp, afill, zfill, pfill):
    fills = afill.copy()
    fills.update(zfill)
    fills.update(pfill)
    filen = temp.format(**fills)
    with h5py.File(filen, 'r') as f:
        cosmopars = {key: val for key, val in \
            f['Header/inputpars/cosmopars'].attrs.items()}
        hdpath = 'Header/inputpars/halodata'
        mvir_msun = f[hdpath].attrs['Mvir_g']
        if 'Rvir_ckpcoverh' in f[hdpath].attrs:
            rvir_ckpcoverh = f[hdpath].attrs['Rvir_ckpcoverh']
            rvir_pkpc = rvir_ckpcoverh * cosmopars['a'] \
                        / cosmopars['h']
        elif 'Rvir_cm' in f[hdpath].attrs:
            rvir_cm = f[hdpath].attrs['Rvir_cm']
            rvir_pkpc = rvir_cm / (c.cm_per_mpc * 1e-3)
        pax = f['Header/inputpars'].attrs['Axis3']
    mvir_msun /= c.solar_mass
    outdct = {'mv': mvir_msun, 'rv': rvir_pkpc,
              'a3': pax, 'z': cosmopars['z']}
    return outdct

def plotcomp_zev_projax_phys(filen_template, paxfills, zfills,
                             physfills, physlabels, title=None, 
                             outname=None, ylabel=None):
    '''
    in each plot, for one IC and type of column, show the profile
    with impact parameter for different redshifts (line colors),
    projections axes (columns; last is all), and physics models
    (different rows). The bottom panels show redshift ensemble
    comparisons between models.

    Parameters:
    -----------
    filen_templates_phys: str
        file name, used to make individual file names with .format.
        Should include the full path.
    paxfills: list of dicts
        the keyword/value sets usable with .format to make the file
        names for different projection axes.
    physfills: list of dicts
        like paxfills, but the keywords and values produce file names
        for different physics models.
    zfills: list of list of dicts
        like paxfills, but these keywords and values give files for 
        different redshifts. The outer list layer is index-matched
        to the physfills, the inner lists should give redshifts for 
        the different physics models which are the same or very close. 
    physlabels: list of str
        the labels for the different physics models. Index-matched to
        physfills.
    title: str or None
        figure title, if any.
    outname: str
        file to save the image to. Should include the full path.
    
    '''
    
    # axes outer to inner: projection axis, physics, redshift
    mapprops = [[[gethalodct_z_pax_phys(filen_template, afill, zfill, pfill)
                  for zfill in zflist]
                 for zflist, pfill in zip(zfills, physfills)]
                for afill in paxfills]
    rvirs = [[[l3['rv'] for l3 in l2] for l2 in l1] for l1 in mapprops]
    mvirs = [[[l3['mv'] for l3 in l2] for l2 in l1] for l1 in mapprops]
    axis3s = [[[l3['a3'] for l3 in l2] for l2 in l1] for l1 in mapprops]
    redshifts = [[[l3['z'] for l3 in l2] for l2 in l1] for l1 in mapprops]
    ## axes outer to inner: projection axis, physics, redshift
    # should not depend on projection axis
    if not np.all([[[rvirs[0][j][k] == rvirs[i][j][k]
                     for k in range(len(zfills))]
                    for j in range(len(physfills))]
                   for i in range(len(paxfills))]):
        msg = ('Different projection axes (outermost index) list different'
               f'virial radii [pkpc]:\n{rvirs}')
        raise ValueError(msg)
    if not np.all([[[mvirs[0][j][k] == mvirs[i][j][k]
                     for k in range(len(zfills))]
                    for j in range(len(physfills))]
                   for i in range(len(paxfills))]):
        msg = ('Different projection axes (outermost index) list different'
               f'virial masses [Msun]:\n{mvirs}')
        raise ValueError(msg)
    if not np.all([[[axis3s[i][0][0] == axis3s[i][j][k]
                     for k in range(len(zfills))]
                    for j in range(len(physfills))]
                   for i in range(len(paxfills))]):
        msg = ('Different phys/redshifts list different'
               f'projection axes (inner 2 indices):\n{axis3s}')
        raise ValueError(msg)
    if not np.all([[[np.isclose(redshifts[0][0][k], redshifts[i][j][k],
                                rtol=1e-2, atol=1e-3)
                     for k in range(len(zfills))]
                    for j in range(len(physfills))]
                   for i in range(len(paxfills))]):
        msg = ('Different phys/proj. ax list different'
               f'redshifts (outer 2 indices):\n{redshifts}')
        raise ValueError(msg)
    redshifts = redshifts[0][0]
    inds_zslotohi = np.argsort(redshifts)
    
    fontsize = 12
    ncols = len(paxfills)
    nrows = len(physfills)
    panelsize = 2.5
    hspace = 0.
    wspace = 0.
    laxheight = 1.
    width_ratios = [panelsize] * (ncols + 1)
    height_ratios = [panelsize] * (nrows + 1) + [laxheight]
    width = sum(width_ratios) \
            * (1. + (len(width_ratios) - 1.) / len(width_ratios) * wspace)
    height = sum(height_ratios) \
             * (1. + (len(height_ratios) - 1.) / len(height_ratios) * hspace)
    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(nrows=nrows + 2, ncols=ncols + 1, hspace=hspace, 
                        wspace=wspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    maxes = [[fig.add_subplot(grid[i, j]) \
              for j in range(ncols)] for i in range(nrows)]
    allphaxes = [fig.add_subplot(grid[nrows, j]) \
                 for j in range(ncols)]
    alla3axes = [fig.add_subplot(grid[i, ncols]) \
                 for i in range(nrows)]
    allcombax = fig.add_subplot(grid[nrows, ncols])
    alax = fig.add_subplot(grid[nrows + 1, -1])
    plax = fig.add_subplot(grid[nrows + 1, 0])
    zlax = fig.add_subplot(grid[nrows + 1, 1:-1])
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    nz = len(zfills[0])
    _colors = tc.tol_cmap('rainbow_discrete', nz)
    zcolors = _colors(np.linspace(1. / nz, 1. - 1. / nz, nz))[::-1]
    pcolors = sl.physcolors
    alpha = 0.3
    _cset = tc.tol_cset('vibrant')
    acolors = [_cset.blue, _cset.teal, _cset.orange, _cset.red, _cset.magenta]

    rvmin = np.min(rvirs)
    rbins = np.linspace(0., 2. * rvmin, 50.)
    rc = 0.5 * (rbins[:-1] + rbins[1:])
    pv = ['perc-0.1', 'perc-0.5', 'perc-0.9']

    for pi, pfill in enumerate(physfills):
        __fills = pfill.copy()
        fns_all = []
        for a3i, afill in enumerate(paxfills):
            _fills = __fills.copy()
            _fills.update(afill)
            ax = maxes[pi][a3i]
            ax.tick_params(labelsize=fontsize - 1., direction='in',
                           top=True, left=True, which='both',
                           labelbottom=False, labelleft=(a3i == 0))
            ax.grid(True)
            if a3i == 0 and ylabel is not None:
                ax.set_ylabel(ylabel, fontsize=fontsize)
            if pi == 0:
                collabel = 'XYZ'[axis3s[a3i][pi][0]] + '-proj.'
                ax.set_title(collabel, fontsize=fontsize)
            if a3i == 0:
                mvmin = np.min(mvirs[0][pi])
                mvmax = np.max(mvirs[0][pi])
                mlabel = ('$\\log \\, \\mathrm{M}_{\\mathrm{vir}}'
                          '/ \\mathrm{M}_{\\odot} ='
                          f'{np.log10(mvmin):.1f} \\endash '
                          f'{np.log10(mvmax):.1f}$')
                ax.text(0.98, 0.98, mlabel, fontsize=fontsize - 1, 
                        transform=ax.transAxes, horizontalalignment='right',
                        verticalalignment='top')
            fns = []
            for zi in inds_zslotohi:
                zfill = zfills[pi][zi]
                color = zcolors[zi]
                fills = _fills.copy()
                fills.update(zfill)
                fn = filen_template.format(**fills)
                fns.append(fn)
                plo, pmed, phi = g2d.get_profile_massmap(fn, rbins, 
                                                         rbin_units='pkpc',
                                                         profiles=pv)
                ax.plot(rc, pmed, color=color, linestyle='solid',
                        linewidth=1.5)
                ax.plot(rc, plo, color=color, linestyle='dashed',
                        linewidth=1.5)
                ax.plot(rc, phi, color=color, linestyle='dashed',
                        linewidth=1.5)
                rv = rvirs[a3i][pi][zi]
                yp_med = mu.linterpsolve(rc, pmed, rv)
                ax.scatter([rv], [yp_med], marker='o', c=[color], s=60)
            
            plo, pmed, phi = g2d.get_profile_massmap(fns, rbins, 
                                                     rbin_units='pkpc',
                                                     profiles=pv)
            ax.plot(rc, pmed, color='black', linestyle='solid',
                    linewidth=1.5)
            ax.fill_between(rc, plo, phi, color='black', alpha=alpha)
            
            _ax = alla3axes[pi]
            color = acolors[a3i]
            _ax.plot(rc, pmed, color=color, linestyle='solid',
                     linewidth=1.5)
            _ax.plot(rc, plo, color=color, linestyle='dashed',
                     linewidth=1.5)
            _ax.plot(rc, phi, color=color, linestyle='dashed',
                     linewidth=1.5)
            fns_all = fns_all + fns
            
        # compare all z proj. ax differences
        _ax.tick_params(labelsize=fontsize - 1., direction='in',
                        top=True, left=True, which='both',
                        labelbottom=False, labelleft=False)
        _ax.grid(True)
        _ax.text(0.98, 0.98, 'all z', fontsize=fontsize,
                 transform=_ax.transAxes, horizontalalignment='right',
                 verticalalignment='top')
        if pi == 0:
            collabel = 'all proj.'
            _ax.set_title(collabel, fontsize=fontsize)
        _ax.text(1.05, 0.5, physlabels[pi], fontsize=fontsize,
                 transform=_ax.transAxes, horizontalalignment='left',
                 verticalalignment='center', rotation=90.)
        plo, pmed, phi = g2d.get_profile_massmap(fns_all, rbins, 
                                                 rbin_units='pkpc',
                                                 profiles=pv)
        _ax.plot(rc, pmed, color='black', linestyle='solid',
                 linewidth=1.5)
        _ax.fill_between(rc, plo, phi, color='black', alpha=alpha)
    
    fns_phys = {key: [] for key in physlabels}
    for a3i, afill in enumerate(paxfills):
        ax = allphaxes[a3i]
        ax.set_xlabel('$\\mathrm{r}_{\\perp}\\;[\\mathrm{pkpc}]$',
                      fontsize=fontsize)
        __fills = afill.copy()
        ax.tick_params(labelsize=fontsize - 1., direction='in',
                       top=True, left=True, which='both',
                       labelbottom=True, labelleft=a3i == 0)
        ax.grid(True)
        ax.text(0.98, 0.98, 'all z', fontsize=fontsize,
                transform=ax.transAxes, horizontalalignment='right',
                verticalalignment='top')
        if a3i == 0 and ylabel is not None:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        for pi, pfill in enumerate(physfills):
            _fills = __fills.copy()
            _fills.update(pfill)
            color = pcolors[physlabels[pi]]
            fns = []
            for zi in inds_zslotohi:
                zfill = zfills[pi][zi]
                fills = _fills.copy()
                fills.update(zfill)
                fn = filen_template.format(**fills)
                fns.append(fn)
            plo, pmed, phi = g2d.get_profile_massmap(fns, rbins, 
                                                     rbin_units='pkpc',
                                                     profiles=pv)
            ax.plot(rc, pmed, color=color, linestyle='solid',
                    linewidth=1.5)
            ax.fill_between(rc, plo, phi, color=color, alpha=alpha)

            fns_phys[physlabels[pi]] = fns_phys[physlabels[pi]] + fns
    
    ax = allcombax
    ax.set_xlabel('$\\mathrm{r}_{\\perp}\\;[\\mathrm{pkpc}]$',
                  fontsize=fontsize)
    ax.tick_params(labelsize=fontsize - 1., direction='in',
                    top=True, left=True, which='both',
                    labelbottom=True, labelleft=False)
    ax.grid(True)
    ax.text(1.05, 0.5, 'all phys.', fontsize=fontsize,
            transform=ax.transAxes, horizontalalignment='left',
            verticalalignment='center', rotation=90.)
    ax.text(0.98, 0.98, 'all z, proj.', fontsize=fontsize,
            transform=ax.transAxes, horizontalalignment='right',
            verticalalignment='top')
    for plabel in physlabels:
        color = pcolors[plabel]
        fns = fns_phys[plabel]
        plo, pmed, phi = g2d.get_profile_massmap(fns, rbins, 
                                                 rbin_units='pkpc',
                                                 profiles=pv)
        ax.plot(rc, pmed, color=color, linestyle='solid',
                linewidth=1.5)
        ax.fill_between(rc, plo, phi, color=color, alpha=alpha)
    
    allcombax.set_xscale('log')
    [ax.set_xscale('log') for ax in alla3axes]
    [ax.set_xscale('log') for ax in allphaxes]
    [ax.set_xscale('log') for l1 in maxes for ax in l1]
    #sync ax ranges
    ymin, ymax = allcombax.get_ylim()
    xmin, xmax = allcombax.get_xlim()
    ylims = [ax.get_ylim() for ax in alla3axes]
    xlims = [ax.get_xlim() for ax in alla3axes]
    ylims += [ax.get_ylim() for ax in allphaxes]
    xlims += [ax.get_xlim() for ax in allphaxes]
    ylims += [ax.get_ylim() for l1 in maxes for ax in l1]
    xlims += [ax.get_xlim() for l1 in maxes for ax in l1]
    ymin = min(ymin, np.min([yl[0] for yl in ylims]))
    ymax = max(ymax, np.max([yl[1] for yl in ylims]))
    ymin = max(ymin, ymax - 5.)
    xmin = min(xmin, np.min([xl[0] for xl in xlims]))
    xmax = max(xmax, np.max([xl[1] for xl in xlims]))
    allcombax.set_ylim((ymin, ymax))
    allcombax.set_xlim((xmin, xmax))
    [ax.set_ylim(ymin, ymax) for ax in alla3axes]
    [ax.set_xlim(xmin, xmax) for ax in alla3axes]
    [ax.set_ylim(ymin, ymax) for ax in allphaxes]
    [ax.set_xlim(xmin, xmax) for ax in allphaxes]
    [ax.set_ylim(ymin, ymax) for l1 in maxes for ax in l1]
    [ax.set_xlim(xmin, xmax) for l1 in maxes for ax in l1]

    zlax.axis('off')
    line = [[(0, 0)]]
    lcs = mcol.LineCollection(line * len(zcolors), linestyle='solid', 
                              linewidth=1.5, colors=zcolors)
    lcd = mcol.LineCollection(line * len(zcolors), linestyle='dashed', 
                              linewidth=1.5, colors=zcolors)
    zhandles = [lcs, lcd,
                mlines.Line2D((), (), linestyle='solid', linewidth=1.5,
                              label='all med.', color='black'),
                mpatch.Patch(label='all perc. 10-90', linewidth=0.5, 
                             color='black', alpha=alpha),
                mlines.Line2D((), (), linestyle=None, marker='o',
                              label='$\\mathrm{R}_{\\mathrm{vir}}$', 
                              color='gray', markersize=5)
                ]
    zlabels = ['median', '10, 90%', 'all z med.', 'all z 10-90%', 
               '$\\mathrm{R}_{\\mathrm{vir}}$']
    zhandles += [mlines.Line2D((), (), linewidth=1.5, linestyle='solid',
                              label=f'$z={redshifts[zi]:.1f}$', 
                              color=zcolors[zi])
                 for zi in inds_zslotohi]
    zlabels += [f'$z={redshifts[zi]:.1f}$' for zi in inds_zslotohi]
    zlax.legend(zhandles, zlabels, 
                fontsize=fontsize, ncol=(ncols - 1),
                handler_map={type(lcs): pu.HandlerDashedLines()},
                bbox_to_anchor=(0.5, 0.2), loc='upper center',
                title='redshift comp.', title_fontsize=fontsize)
    
    alax.axis('off')
    line = [[(0, 0)]]
    lcs = mcol.LineCollection(line * len(paxfills), linestyle='solid', 
                              linewidth=1.5, colors=acolors[:len(paxfills)])
    lcd = mcol.LineCollection(line * len(paxfills), linestyle='dashed', 
                              linewidth=1.5, colors=acolors[:len(paxfills)])
    ahandles = [lcs, lcd,
                mlines.Line2D((), (), linestyle='solid', linewidth=1.5,
                              label='all med.', color='black'),
                mpatch.Patch(label='all perc. 10-90', linewidth=0.5, 
                             color='black', alpha=alpha)
                ]
    alabels = ['median', '10, 90%', 'all pr. med.', 'all pr. 10-90%']
    ahandles += [mlines.Line2D((), (), linewidth=1.5, linestyle='solid',
                              label=f'{"XYZ"[a3i]}-proj.', 
                              color=acolors[a3i])
                 for a3i in np.arange(3)]
    alabels += [f'{"XYZ"[a3i]}-proj.' for a3i in np.arange(3)]
    alax.legend(ahandles, alabels, 
                fontsize=fontsize, ncol=1,
                handler_map={type(lcs): pu.HandlerDashedLines()},
                bbox_to_anchor=(1.0, 0.2), loc='upper right',
                title='proj. ax comp.', title_fontsize=fontsize)
    
    plax.axis('off')
    line = [[(0, 0)]]
    cvals = [pcolors[key] for key in pcolors]
    lcs = mcol.LineCollection(line * len(pcolors), linestyle='solid', 
                              linewidth=1.5, colors=cvals)
    phandles = [lcs,
                mpatch.Patch(label='all perc. 10-90', linewidth=0.5, 
                             color='gray', alpha=alpha)
                ]
    plabels = ['median', '10-90%']
    phandles += [mlines.Line2D((), (), linewidth=1.5, linestyle='solid',
                              label=physlab, 
                              color=pcolors[physlab])
                 for physlab in physlabels]
    plabels += physlabels
    plax.legend(phandles, plabels, 
                fontsize=fontsize, ncol=1,
                handler_map={type(lcs): pu.HandlerDashedLines()},
                bbox_to_anchor=(0.0, 0.2), loc='upper left',
                title='physics comp.', title_fontsize=fontsize)
    
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')


def plotsetcomp_zev_projax_phys(fileset='set3_model3'):

    if fileset == 'set3_model3':
         # quest
        outdir = '/projects/b1026/nastasha/imgs/2dcomp_set3_model3/'
        fdir = '/projects/b1026/nastasha/maps/set3_model3/'
        ftemp = ('coldens_{qty}_{{simname}}_snap{{snapnum}}_'
                 'shrink-sph-cen_BN98_2rvir_{{pax}}-proj_v3.hdf5')
        simnames_all = sl.m12_hr_all1 + sl.m12_sr_all1 \
                       + sl.m13_hr_all1 + sl.m13_sr_all1
        simlists = []
        iclab = []
        for simname in simnames_all:
            ic = simname.split('_')[0]
            if ic in iclab:
                si = np.where([ic == _ic for _ic in iclab])[0][0]
                simlists[si].append(simname)
            else:
                iclab.append(ic)
                simlists.append([simname])
        all_hr = sl.m12_hr_all1 + sl.m13_hr_all1
        all_sr = sl.m12_sr_all1 + sl.m13_sr_all1
        zfillsets = [[[{'snapnum': snap} for snap in sl.snaps_hr]
                      if simname in all_hr else 
                      [{'snapnum': snap} for snap in sl.snaps_sr]
                      if simname in all_sr else 
                      None
                      for simname in simlist] 
                     for simlist in simlists]
        physlabels = [['AGN-CR' if 'MHDCRspec1' in simname
                        else 'noBH' if 'sdp1e10' in simname
                        else 'AGN-noCR'
                        for simname in simlist]
                       for simlist in simlists]
        _afills = ['x', 'y', 'z']
        paxfills = [[{'pax': val} for val in _afills]]
        qtys = ['gas-mass', 'Neon', 'Ne8', 'O6', 'Mg10']
        qtylabs = ['Gas', 'Neon', 'Ne VIII', 'O VI', 'Mg X']
        qtyylabs = [('$\\log_{10} \\, \\Sigma(\\mathrm{gas})'
                     ' \\; [\\mathrm{g}\\,\\mathrm{cm}^{-2}]$'),
                    ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Neon})'
                     ' \\; [\\mathrm{cm}^{-2}]$'),
                    ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\, VIII})'
                     ' [\\mathrm{cm}^{-2}]$'),
                    ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{O\\, VI})'
                     ' [\\mathrm{cm}^{-2}]$'),
                    ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Mg\\, X})'
                     ' [\\mathrm{cm}^{-2}]$')
                   ]
        outname = ('comp_2dprof_z_phys_projax_{qty}_{ic}.pdf')
        simlists_ext = [[{'simname': simn} for simn in siml] 
                        for siml in simlists] * len(qtys)
        zfills_ext = zfillsets * len(qtys)
        physlabels_ext = physlabels * len(qtys)
        ics_ext = iclab * len(qtys)
        paxfills_ext = paxfills * len(zfills_ext)
        qtylabels_ext = [lab for lab in qtylabs for i in range(len(simlists))]
        ylabels_ext = [lab for lab in qtyylabs for i in range(len(simlists))]
        qtys_ext = [lab for lab in qtys for i in range(len(simlists))]

    for (simfills, qty, ylabel, qtylabel, paxfills, zfills, physlabels, ic) \
            in zip(simlists_ext, qtys_ext, ylabels_ext, qtylabels_ext,
                   paxfills_ext, zfills_ext, physlabels_ext, ics_ext):
        filen_template = fdir + ftemp.format(qty=qty)
        title = f'{ic}, z=0.5-1.0, {qtylabel}'
        print('running ', title)
        if np.any([sd['simname'] in sl.buglist1 for sd in simfills]):
            title = title + ', possible bug'
        title = title + '\n' 
        title += '\n'.join([f'{plab}: {sn["simname"]}' 
                            for plab, sn in zip(physlabels, simfills)])
        _outname = outdir + outname.format(ic=ic, qty=qty)

        plotcomp_zev_projax_phys(filen_template, paxfills, zfills,
                                simfills, physlabels, title=title, 
                                outname=_outname, ylabel=ylabel)
        plt.close() # don't overload memory with too many plots

# to avoid the previous 'annoying long nested loop' issue
def gethalodct_phys_haloset(temp, afill, zfill, pfill):
    fills = afill.copy()
    fills.update(zfill)
    fills.update(pfill)
    filen = temp.format(**fills)
    with h5py.File(filen, 'r') as f:
        cosmopars = {key: val for key, val in \
            f['Header/inputpars/cosmopars'].attrs.items()}
        hdpath = 'Header/inputpars/halodata'
        mvir_msun = f[hdpath].attrs['Mvir_g']
        if 'Rvir_ckpcoverh' in f[hdpath].attrs:
            rvir_ckpcoverh = f[hdpath].attrs['Rvir_ckpcoverh']
            rvir_pkpc = rvir_ckpcoverh * cosmopars['a'] \
                        / cosmopars['h']
        elif 'Rvir_cm' in f[hdpath].attrs:
            rvir_cm = f[hdpath].attrs['Rvir_cm']
            rvir_pkpc = rvir_cm / (c.cm_per_mpc * 1e-3)
        pax = f['Header/inputpars'].attrs['Axis3']
    mvir_msun /= c.solar_mass
    outdct = {'mv': mvir_msun, 'rv': rvir_pkpc,
              'a3': pax, 'z': cosmopars['z'],
              'fn': filen}
    return outdct

def plotcomp_phys_haloset(filen_template, paxfills, zfills,
                          physfills, iclabels, physlabels, title=None,
                          outname=None, ylabel=None):
    '''
    in each panel, redshift- and projection-axis-combined profiles
    are shown. Different panels show different ICs, different plots
    are for different halo sets (m12/m13) and projected quantities.

    Parameters:
    -----------
    filen_template_phys: str
        file name, used to make individual file names with .format.
        Should include the full path.
    paxfills: list of dicts
        the keyword/value sets usable with .format to make the file
        names for different projection axes.
    physfills: list of lists of dicts
        like paxfills, but the keywords and values produce file names
        for different physics models. The outer list corresponds to
        different panels (ICs), inner lists are different physics for
        the same ICs.
    zfills: list of lists of lists of dicts (I know...)
        like paxfills, but these keywords and values give files for 
        different redshifts. The outer list layer is index-matched
        to the panels the middle lists are index-matched to the 
        different physics in each panel, and the innermost lists of
        dicts should give redshifts for each different physics model
        which are the same or very close. 
    iclabels: list of str
        the labels for the different ICs (panels). Index-matched to
        the outer list of the physfills.
    physlabels: list of lists of str
        the labels for the different physics models. Index-matched to
        physfills.
    title: str or None
        figure title, if any.
    ylabel: str or None
        y axis label, if any.
    outname: str
        file to save the image to. Should include the full path.
    
    '''
    
    # axes outer to inner: projection axis, physics, redshift
    mapprops = [[[[gethalodct_phys_haloset(filen_template, afill, zfill, 
                                           pfill)
                   for afill in paxfills] # proj. axis
                  for zfill in zflist2] # redshift
                 for zflist2, pfill in zip(zflist1, pflist)] # phys.
                for zflist1, pflist in zip(zfills, physfills)] # ICS
                
    rvirs = [[[[l4['rv'] for l4 in l3]for l3 in l2] for l2 in l1] 
             for l1 in mapprops]
    mvirs = [[[[l4['mv'] for l4 in l3] for l3 in l2] for l2 in l1] 
             for l1 in mapprops]
    axis3s = [[[[l4['a3'] for l4 in l3] for l3 in l2] for l2 in l1] 
              for l1 in mapprops]
    redshifts = [[[[l4['z'] for l4 in l3] for l3 in l2] for l2 in l1] 
                 for l1 in mapprops]
    filens = [[[[l4['fn'] for l4 in l3] for l3 in l2] for l2 in l1] 
                 for l1 in mapprops]
    ## axes outer to inner: projection axis, physics, redshift
    # should not depend on projection axis
    if not np.all([np.all([np.all([[rvirs[i][j][k][0] == rvirs[i][j][k][l]
                      for l in range(len(paxfills))]
                     for k in range(len(zfills[i][j]))])
                    for j in range(len(physfills[i]))])
                   for i in range(len(physfills))]):
        msg = ('Different projection axes (innermost index) list different'
               f'virial radii [pkpc]:\n{rvirs}')
        raise ValueError(msg)
    if not np.all([np.all([np.all([[mvirs[i][j][k][0] == mvirs[i][j][k][l]
                      for l in range(len(paxfills))]
                     for k in range(len(zfills[i][j]))])
                    for j in range(len(physfills[i]))])
                   for i in range(len(physfills))]):
        msg = ('Different projection axes (innermost index) list different'
               f'virial masses [Msun]:\n{mvirs}')
        raise ValueError(msg)
    if not np.all([np.all([np.all([[axis3s[0][0][0][l] == axis3s[i][j][k][l]
                      for l in range(len(paxfills))]
                     for k in range(len(zfills[i][j]))])
                    for j in range(len(physfills[i]))])
                   for i in range(len(physfills))]):
        msg = ('Different phys/redshifts list different'
               f'projection axes (all but innermost axis):\n{axis3s}')
        raise ValueError(msg)
    if not np.all([np.all([np.all([[np.isclose(redshifts[0][0][k][0],
                                               redshifts[i][j][k][l],
                                               rtol=1e-2, atol=1e-3)
                      for l in range(len(paxfills))]
                     for k in range(len(zfills[i][j]))])
                    for j in range(len(physfills[i]))])
                   for i in range(len(physfills))]):
        msg = ('Different ics/phys/proj. ax list different'
               f'redshifts (all but 3rd index):\n{redshifts}')
        raise ValueError(msg)
    
    npanels = len(physfills)
    ncols_max = 4
    if npanels <= ncols_max:
        ncols = npanels
        nrows = 1
        lax_below = True
    elif npanels % ncols_max == 0:
        ncols = ncols_max
        nrows = npanels // ncols
        lax_below = True
    else:
        ncols = ncols_max
        nrows = (npanels - 1) // ncols + 1
        lax_below = False
    panelsize = 3.
    laxheight = panelsize 
    hspace = 0.
    wspace = 0.
    ncols_legend = (ncols - npanels % ncols)
    if lax_below:
        width_ratios = [panelsize] * (ncols)
        height_ratios = [panelsize] * (nrows) + [laxheight]
        width = sum(width_ratios) \
                * (1. + (len(width_ratios) - 1.) / len(width_ratios) * wspace)
        height = sum(height_ratios) \
                * (1. + (len(height_ratios) - 1.) / len(height_ratios) * hspace)
        fig = plt.figure(figsize=(width, height))
        grid = gsp.GridSpec(nrows=nrows + 1, ncols=ncols, hspace=hspace, 
                            wspace=wspace, width_ratios=width_ratios,
                            height_ratios=height_ratios)
        axes = [fig.add_subplot(grid[i // ncols, i % ncols]) \
                for i in range(npanels)]
        lax = fig.add_subplot(grid[nrows, :])
    else:
        width_ratios = [panelsize] * (ncols)
        height_ratios = [panelsize] * (nrows)
        width = sum(width_ratios) \
                * (1. + (len(width_ratios) - 1.) / len(width_ratios) * wspace)
        height = sum(height_ratios) \
                * (1. + (len(height_ratios) - 1.) / len(height_ratios) * hspace)
        fig = plt.figure(figsize=(width, height))
        grid = gsp.GridSpec(nrows=nrows, ncols=ncols, hspace=hspace, 
                            wspace=wspace, width_ratios=width_ratios,
                            height_ratios=height_ratios)
        axes = [fig.add_subplot(grid[i // ncols, i % ncols]) \
                for i in range(npanels)]
        lax = fig.add_subplot(grid[nrows - 1, npanels % ncols:])
    
    phys_used = set()
    pcolors = sl.physcolors.copy()
    pv = ['perc-0.1', 'perc-0.5', 'perc-0.9']
    alpha = 0.3
    fontsize = 12
    for ici in range(npanels):
        lby = ici % ncols == 0
        lbx = ici >= npanels - ncols
        ax = axes[ici]
        ax.set_xscale('log')
        ax.tick_params(labelsize=fontsize - 1, direction='in',
                       top=True, right=True, labelleft=lby,
                       labelbottom=lbx, which='both')
        ax.grid(True)
        ax.text(0.98, 0.98, iclabels[ici], fontsize=fontsize,
                transform=ax.transAxes, horizontalalignment='right',
                verticalalignment='top')
        if lby and ylabel is not None:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        if lbx:
            ax.set_xlabel('$\\mathrm{r}_{\\perp} \\; [\\mathrm{pkpc}]$',
                          fontsize=fontsize)
        _rvs = [rv for l1 in rvirs[ici] for l2 in l1 for rv in l2]
        rvmin = np.min(_rvs)
        rvmax = np.max(_rvs)
        ax.axvspan(rvmin, rvmax, alpha=0.2, color='gray')
        for pi, plab in enumerate(physlabels[ici]):
            color = pcolors[plab]
            phys_used.add(plab)
            
            _fns = [fn for l1 in filens[ici][pi] for fn in l1]
            rbins = np.linspace(0., 2. * rvmin, 50)
            rc = 0.5 * (rbins[:-1] + rbins[1:])
            plo, pmed, phi = g2d.get_profile_massmap(_fns, rbins, 
                                                     rbin_units='pkpc',
                                                     profiles=pv)
            ax.plot(rc, pmed, color=color, linestyle='solid',
                    linewidth=1.5)
            ax.fill_between(rc, plo, phi, color=color, alpha=alpha)

    # sync axis ranges
    xlims = [ax.get_xlim() for ax in axes]
    xmin = np.min(xlims)
    xmax = np.max(xlims)
    [ax.set_xlim((xmin, xmax)) for ax in axes]
    ylims = [ax.get_ylim() for ax in axes]
    ymin = np.min(ylims)
    ymax = np.max(ylims)
    ymin = max(ymin, ymax - 5.)
    [ax.set_ylim((ymin, ymax)) for ax in axes]

    lax.axis('off')
    line = [[(0, 0)]]
    phys_used = sorted(list(phys_used))
    cvals = [pcolors[key] for key in phys_used]
    lcs = mcol.LineCollection(line * len(pcolors), linestyle='solid', 
                              linewidth=1.5, colors=cvals)
    phandles = [lcs,
                mpatch.Patch(label='perc. 10-90', linewidth=0.5, 
                             color='brown', alpha=alpha)
                ]
    plabels = ['median', '10-90%']
    phandles += [mlines.Line2D((), (), linewidth=2., linestyle='solid',
                              label=physlab, 
                              color=pcolors[physlab])
                 for physlab in phys_used]
    plabels += phys_used
    lax.legend(phandles, plabels, 
               fontsize=fontsize, ncol=ncols_legend,
               handler_map={type(lcs): pu.HandlerDashedLines()},
               bbox_to_anchor=(0.5, 0.7), loc='upper center')
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotsetcomp_phys_haloset(fileset='set3_model3'):

    if fileset == 'set3_model3':
         # quest
        outdir = '/projects/b1026/nastasha/imgs/2dcomp_set3_model3/'
        fdir = '/projects/b1026/nastasha/maps/set3_model3/'
        ftemp = ('coldens_{qty}_{{simname}}_snap{{snapnum}}_'
                 'shrink-sph-cen_BN98_2rvir_{{pax}}-proj_v3.hdf5')
        simnames_all = [sl.m12_hr_all1 + sl.m12_sr_all1] \
                       + [sl.m13_hr_all1 + sl.m13_sr_all1]
        simgrouplabels = ['m12_haloes', 'm13_haloes']
        simgroups = []
        iclabelss = []
        for simgroup in simnames_all:
            simlists = []
            iclab = []
            for simname in simgroup:
                ic = simname.split('_')[0]
                if ic in iclab:
                    si = np.where([ic == _ic for _ic in iclab])[0][0]
                    simlists[si].append(simname)
                else:
                    iclab.append(ic)
                    simlists.append([simname])
            simgroups.append(simlists)
            iclabelss.append(iclab)

        all_hr = sl.m12_hr_all1 + sl.m13_hr_all1
        all_sr = sl.m12_sr_all1 + sl.m13_sr_all1
        zfillsets = [[[[{'snapnum': snap} for snap in sl.snaps_hr]
                       if simname in all_hr else 
                       [{'snapnum': snap} for snap in sl.snaps_sr]
                       if simname in all_sr else 
                       None
                       for simname in simlist] 
                      for simlist in simlists]
                     for simlists in simgroups]
        physlabels = [[['AGN-CR' if 'MHDCRspec1' in simname
                         else 'noBH' if 'sdp1e10' in simname
                         else 'AGN-noCR'
                         for simname in simlist]
                        for simlist in simlists]
                       for simlist in simgroups]
        _afills = ['x', 'y', 'z']
        paxfills = [[{'pax': val} for val in _afills]]
        qtys = ['gas-mass', 'Neon', 'Ne8', 'O6', 'Mg10']
        qtylabs = ['Gas', 'Neon', 'Ne VIII', 'O VI', 'Mg X']
        qtyylabs = [('$\\log_{10} \\, \\Sigma(\\mathrm{gas})'
                     ' \\; [\\mathrm{g}\\,\\mathrm{cm}^{-2}]$'),
                    ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Neon})'
                     ' \\; [\\mathrm{cm}^{-2}]$'),
                    ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Ne\\, VIII})'
                     ' [\\mathrm{cm}^{-2}]$'),
                    ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{O\\, VI})'
                     ' [\\mathrm{cm}^{-2}]$'),
                    ('$\\log_{10} \\, \\mathrm{N}(\\mathrm{Mg\\, X})'
                     ' [\\mathrm{cm}^{-2}]$')
                   ]
        outname = ('comp_2dprof_phys_{qty}_{sglab}.pdf')
        simlists_ext = [[[{'simname': simn} for simn in siml] 
                         for siml in simlists] for simlists in simgroups] \
                       * len(qtys)
        iclab_ext = [[(l2[0]['simname'].split('/')[-1]).split('_')[0]
                       for l2 in l1] for l1 in simlists_ext]
        zfills_ext = zfillsets * len(qtys)
        physlabels_ext = physlabels * len(qtys)
        icset_ext = simgrouplabels * len(qtys)
        paxfills_ext = paxfills * len(zfills_ext)
        qtylabels_ext = [lab for lab in qtylabs 
                         for i in range(len(simgroups))]
        ylabels_ext = [lab for lab in qtyylabs for i in range(len(simgroups))]
        qtys_ext = [lab for lab in qtys for i in range(len(simgroups))]

    for (simfills, qty, ylabel, qtylabel, paxfills, zfills, physlabels, 
            icset, iclabs) \
            in zip(simlists_ext, qtys_ext, ylabels_ext, qtylabels_ext,
                   paxfills_ext, zfills_ext, physlabels_ext, icset_ext,
                   iclab_ext):
        filen_template = fdir + ftemp.format(qty=qty)
        title = f'{icset}, z=0.5-1.0, {qtylabel}'
        _outname = outdir + outname.format(sglab=icset, qty=qty)

        plotcomp_phys_haloset(filen_template, paxfills, zfills,
                              simfills, iclabs, physlabels, title=title,
                              outname=_outname, ylabel=ylabel)
        plt.close() # don't overload memory with too many plots