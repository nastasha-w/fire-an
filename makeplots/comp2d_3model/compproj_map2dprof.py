
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

