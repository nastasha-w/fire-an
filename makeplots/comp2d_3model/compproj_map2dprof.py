
import h5py
import matplotlib.collections as mcol
import matplotlib.gridspec as gsp
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import numpy as np

import fire_an.makeplots.get_2dprof as g2d
import fire_an.makeplots.tol_colors as tc
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c


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
    if not np.all([[[np.isclose(redshifts[0][0][k], redshifts[i][j][k])
                     for k in range(len(zfills))]
                    for j in range(len(physfills))]
                   for i in range(len(paxfills))]):
        msg = ('Different phys/proj. ax list different'
               f'redshifts (outer 2 indices):\n{redshifts}')
        raise ValueError(msg)
    redshifts = redshifts[0][0]
    inds_zslotihi = np.argsort(redshifts)
    
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
    maxes = [[fig.add_subplot(grid[j, i]) \
              for i in range(nrows)] for j in range(ncols)]
    allphaxes = [fig.add_subplot(grid[nrows, j]) \
                 for j in range(ncols)]
    alla3axes = [fig.add_subplot(grid[i, ncols]) \
                 for i in range(nrows)]
    allcombax = fig.add_subplot(grid[nrows, ncols])
    alax = fig.add_subplot(grid[nrows + 1, -1])
    plax = fig.add_subplot(grid[nrows + 1, 0])
    zlax = fig.add_subplot(grid[nrows + 1, 1:-1])

    nlines = len(zfills)
    _colors = tc.tol_cmap('rainbow_discrete', nlines)
    zcolors = _colors(np.linspace(1. / nlines, 1. - 1. / nlines, nlines))
    pcolors = sl.physcolors
    alpha = 0.3
    _cset = tc.tol_cset('vibrant')
    acolors = [cset.blue, cset.teal, cset.orange, cset.red, cset.magenta]

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
                           top=True, left=True, which='both')
            ax.grid(True)
            if a3i == 0 and ylabel is not None:
                ax.set_ylabel(ylabel, fontsize=fontsize)
            if phi == 0:
                collabel = 'XYZ'[axis3s[a3i][phi][0]] '-proj.'
                ax.set_title(collabel, fontsize=fontsize)
            if a3i == 0:
                mvmin = np.min(mvir[0][pi])
                mvmax = np.max(mvir[0][pi])
                mlabel = ('$\\log_{10} \\, \\mathrm{M}_{\\mathmr{vir}}'
                          '\\, / \\, \\mathrm{M}_{\\odot} ='
                          f'{np.log10(mvmin):.1f} \\endash '
                          f'{np.log10(mvmax):.1f}$')
                ax.text(0.98, 0.98, mlabel, fontsize=fontsize - 1, 
                        transform=ax.transAxes, horizontalalignment='top',
                        verticalalignment='right')
            fns = []
            for zi in inds_zslotihi:
                zfill = zfills[zi]
                color = zcolors[zi]
                fill = _fill.copy()
                fill.update(zfill)
                fn = filen_template(**fills)
                fns.append(fn)
                plo, pmed, phi = get_profile_massmap(fn, rbins, 
                                                     rbin_units='pkpc',
                                                     profiles=pv)
                ax.plot(rc, pmed, color=color, linestyle='solid',
                        linewidth=1.5)
                ax.plot(rc, plo, color=color, linestyle='dashed',
                        linewidth=1.5)
                ax.plot(rc, phi, color=color, linestyle='dashed',
                        linewidth=1.5)
                rv = rvirs[a3i][pi][zi]
                yp_med = pu.linterpsolve(rc, pmed, rv)
                ax.scatter([rv], [yp_med], marker='o', c=color, s=60)
            
            plo, pmed, phi = get_profile_massmap(fns, rbins, 
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
                 transform=ax.transAxes, horizontalalignment='right',
                 verticalalignment='top')
        if phi == 0:
            collabel = 'all proj.'
            _ax.set_title(collabel, fontsize=fontsize)
        _ax.text(1.05, 0.05, physlabels[pi], fontsize=fontsize,
                 transform=_ax.transAxes, horizontalalignment='left',
                 verticalalignment='center', rotation=90.)
        plo, pmed, phi = get_profile_massmap(fns_all, rbins, 
                                             rbin_units='pkpc',
                                             profiles=pv)
        _ax.plot(rc, pmed, color='black', linestyle='solid',
                 linewidth=1.5)
        _ax.fill_between(rc, plo, phi, color='black', alpha=alpha)
    
    fns_phys = {key: [] for key in physlabels}
    for a3i, afill in enumerate(paxfills):
        ax = alla3xes[a3i]
        ax.set_xlabel('$\\mathrm{r}_{\\perp}\\;[\\mathrm{pkpc}]$',
                      fontsize=fontsize)
        __fills = afills.copy()
        ax.tick_params(labelsize=fontsize - 1., direction='in',
                       top=True, left=True, which='both',
                       labelbottom=True, labelleft=a3i == 0)
        ax.grid(True)
        if a3i == 0 and ylabel is not None:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        for pi, pfill in enumerate(physfills):
            _fills = __fills.copy()
            _fills.update(pfill)
            color = pcolors[physlabels[pi]]
            fns = []
            for zi in inds_zslotihi:
                zfill = zfills[zi]
                fill = _fill.copy()
                fill.update(zfill)
                fn = filen_template(**fills)
                fns.append(fn)
            plo, pmed, phi = get_profile_massmap(fns, rbins, 
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
    ax.text(1.05, 0.05, physlabels[phi], fontsize=fontsize,
            transform=ax.transAxes, horizontalalignment='left',
            verticalalignment='center', rotation=90.)
    ax.text(0.98, 0.98, 'all z, proj.', fontsize=fontsize,
            transform=ax.transAxes, horizontalalignment='right',
            verticalalignment='top')
    for plabel in physlabels:
        color = pcolors[plabel]
        fns = fns_phys[plabel]
        plo, pmed, phi = get_profile_massmap(fns, rbins, 
                                             rbin_units='pkpc',
                                             profiles=pv)
        ax.plot(rc, pmed, color=color, linestyle='solid',
                linewidth=1.5)
        ax.fill_between(rc, plo, phi, color=color, alpha=alpha)
    
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
                 for zi in inds_zslotihi]
    zlabels += [f'$z={redshifts[zi]:.1f}$' for zi in inds_zslotih]
    zlax.legend(zhandles, zlabels, 
                fontsize=fontsize, ncols=(ncols - 1),
                handler_map={type(lc): pu.HandlerDashedLines()},
                bbox_to_anchor=(0.5, 1.0), loc='upper center',
                title='redshift comp.', title_fontsize=fontsize)
    
    alax.axis('off')
    line = [[(0, 0)]]
    lcs = mcol.LineCollection(line * len(acolors), linestyle='solid', 
                              linewidth=1.5, colors=acolors)
    lcd = mcol.LineCollection(line * len(acolors), linestyle='dashed', 
                              linewidth=1.5, colors=acolors)
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
                fontsize=fontsize, ncols=1,
                handler_map={type(lc): pu.HandlerDashedLines()},
                bbox_to_anchor=(1.0, 0.65), loc='upper right',
                title='proj. ax comp.', title_fontsize=fontsize)
    
    plax.axis('off')
    line = [[(0, 0)]]
    lcs = mcol.LineCollection(line * len(pcolors), linestyle='solid', 
                              linewidth=1.5, colors=pcolors)
    phandles = [lcs,
                mpatch.Patch(label='all perc. 10-90', linewidth=0.5, 
                             color='gray', alpha=alpha)
                ]
    plabels = ['median', '10-90%']
    phandles += [mlines.Line2D((), (), linewidth=1.5, linestyle='solid',
                              label=physlab, 
                              color=pcolor[physlab])
                 for physlab in physlabels]
    plabels += physlabels
    plax.legend(phandles, plabels, 
                fontsize=fontsize, ncols=1,
                handler_map={type(lc): pu.HandlerDashedLines()},
                bbox_to_anchor=(0.0, 0.65), loc='upper left',
                title='physics comp.', title_fontsize=fontsize)
    
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotsetcomp_zev_projax_phys(fileset='set3_model3'):

    if fileset == 'set3_model3':
         # quest
        outdir = '/projects/b1026/nastasha/imgs/2dcomp_set3_model3/'
        fdir = '/projects/b1026/nastasha/maps/set3_model3/'
        ftemp = ('coldens_{{qty}}_{simname}_snap{snapnum}_'
                 'shrink-sph-cen_BN98_2rvir_{pax}-proj_v3.hdf5')
        simnames_all = sl.m12_hr_all1 + sl.m12_sr_all1 \
                       + sl.m13_hr_all1 + sl.m13_sr_all1

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



    plotcomp_zev_projax_phys(filen_template, paxfills, zfills,
                             physfills, physlabels, title=None, 
                             outname=None, ylabel=None)