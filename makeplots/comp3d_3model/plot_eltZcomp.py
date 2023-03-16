'''
In each panel, plot an individual snapshot's abundance median and 
scatter, for different elements, as a function of radius. A second
set of panels shows the element ratios compared to solar. Plot columns
vary the redshift, and rows show different weights. There are different
plots for different ICs and physics models.
'''

import h5py
import numpy as np

import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt

import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
import fire_an.utils.constants_and_units as c
import fire_an.utils.opts_locs as ol


def plot3dprof_eltZcomp_weightzgrid(filen_template, fillkw_weights, 
                                    fillkw_vals, fillkw_z,
                                    figtitle=None, outname=None):
    zolab = '$\\mathrm{Z}_{\\mathrm{O}}$'
    znelab = '$\\mathrm{Z}_{\\mathrm{Ne}}$'
    zmglab = '$\\mathrm{Z}_{\\mathrm{Mg}}$'
    vlabopts = {('sim-direct', 'ElementAbundance/Oxygen'): zolab,
                ('sim-direct', 'ElementAbundance/Neon'): znelab,
                ('sim-direct', 'ElementAbundance/Magnesium'): zmglab,
               }
    ylabel = ('$\\log_{10}\\,\\mathrm{Z}_{\\mathrm{elt}} \\;'
              ' [\\mathrm{Z}_{\\mathrm{elt}, \\odot}]$')
    zonorm = ol.solar_abunds_ea['oxygen']
    znenorm = ol.solar_abunds_ea['neon']
    zmgnorm = ol.solar_abunds_ea['magnesium']
    ynormopts = {('sim-direct', 'ElementAbundance/Oxygen'): zonorm,
                 ('sim-direct', 'ElementAbundance/Neon'): znenorm,
                 ('sim-direct', 'ElementAbundance/Magnesium'): zmgnorm,
                }
    wlabopts = {'Mass': 'Mass-wtd',
                'Volume': 'Volume-wtd',
                ('ion', 'O6'): 'O VI-wtd',
                ('ion', 'Ne8'): 'Ne VIII-wtd',
                ('ion', 'Mg10'): 'Mg X-wtd',
                ('ion', 'H1'): 'H I-wtd',
                }
    cset = tc.tol_cset('vibrant')
    vcolors = {('sim-direct', 'ElementAbundance/Oxygen'): cset.orange,
               ('sim-direct', 'ElementAbundance/Neon'): cset.teal,
               ('sim-direct', 'ElementAbundance/Magnesium'): cset.cyan,
               }
    hists_raw = []
    rvals = []
    yvals = []
    ykeys = []
    wkeys = []
    pkeys = []
    rvirs = []
    zvals = []
    runits = None
    for wdict in fillkw_weights:
        tdict = wdict.copy()
        _rvals = []
        _yvals = []
        _ykeys = []
        _wkeys = []
        _pkeys = []
        _rvirs = []
        _zvals = []
        _hists_raw = []
        for zdict in fillkw_z:
            _tdict = tdict.copy()
            _tdict.update(zdict)
            __rvals = []
            __yvals = []
            __ykeys = []
            __wkeys = []
            __pkeys = []
            __rvirs = []
            __zvals = []
            __hists_raw = []
            for vdict in fillkw_vals:
                __tdict = _tdict.copy()
                __tdict.update(vdict)
                filen = filen_template.format(**__tdict)
                with h5py.File(filen, 'r') as f:
                    hdatpath = 'Header/inputpars/halodata'
                    rvir_cm = f[hdatpath].attrs['Rvir_cm']
                    rv = f['axis_0/bins'][:]
                    _runits = f['axis_0'].attrs['units'].decode()
                    yv = f['axis_1/bins'][:]
                    ykey = f['axis_1'].attrs['qty_type'].decode()
                    if ykey == 'sim-direct':
                        _path = 'axis_1/qty_type_args'
                        _key = 'field'
                        ykey = (ykey, f[_path].attrs[_key].decode())
                    elif ykey == 'ion':
                        _path = 'axis_1/qty_type_args'
                        _key = 'ion'
                        ykey = (ykey, f[_path].attrs[_key].decode())
                    elif ykey == 'Metal':
                        _path = 'axis_1/qty_type_args'
                        _key = 'element'
                        ykey = (ykey, f[_path].attrs[_key].decode())
                    hist_raw = f['histogram/histogram'][:]
                    wkey = f['histogram'].attrs['weight_type'].decode()
                    if wkey == 'ion':
                        _path = 'histogram/weight_type_args'
                        _key = 'ion'
                        wkey = (wkey, f[_path].attrs[_key].decode())
                    elif wkey == 'Metal':
                        _path = 'histogram/weight_type_args'
                        _key = 'element'
                        wkey = (wkey, f[_path].attrs[_key].decode())
                    pkey = 'noBH' if '_sdp1e10_' in filen else \
                        'AGN-CR' if '_MHDCRspec1_' in filen else \
                        'AGN-noCR' if '_MHD_' in filen else \
                        None
                    redshift = f['Header/cosmopars'].attrs['z']
                    
                    # adjust for plotting purposes
                    if yv[0] == -np.inf:
                        yv[0] = 2. * yv[1] - yv[2]
                    if yv[-1] == np.inf:
                        yv[-1] = 2. * yv[-2] - yv[-3]
                    
                    if runits is not None:
                        if _runits != runits:
                            raise ValueError('mismatch in r3D units')
                    if _runits == 'pkpc':
                        rvir = rvir_cm / (c.cm_per_mpc * 1e-3)
                    elif _runits == 'Rvir':
                        rvir = 1.
                __hists_raw.append(hist_raw)
                __rvals.append(rv)
                __yvals.append(yv)
                __ykeys.append(ykey)
                __wkeys.append(wkey)
                __pkeys.append(pkey)
                __rvirs.append(rvir)
                __zvals.append(redshift)
            _hists_raw.append(__hists_raw)
            _rvals.append(__rvals)
            _yvals.append(__yvals)
            _ykeys.append(__ykeys)
            _wkeys.append(__wkeys)
            _pkeys.append(__pkeys)
            _rvirs.append(__rvirs)
            _zvals.append(__zvals)
        hists_raw.append(_hists_raw)
        rvals.append(_rvals)
        yvals.append(_yvals)
        ykeys.append(_ykeys)
        wkeys.append(_wkeys)
        pkeys.append(_pkeys)
        rvirs.append(_rvirs)
        zvals.append(_zvals)
    #print(ykeys)
    #print(wkeys)

    if runits == 'Rvir':
        xlabel = '$\\mathrm{r} \\; [\\mathrm{R}_{\\mathrm{vir}}]$'
        markrvir = False
    else:
        xlabel = '$\\mathrm{r} \\; [\\mathrm{pkpc}]$'
        markrvir = True

    panelsize = 2.
    diffpanelheight = 1.5
    legheight = 2.
    ncols = len(fillkw_z)
    nrows = len(fillkw_weights)
    wspace = 0.
    hspace = 0. 
    nZ = len(fillkw_vals)
    width_ratios = [panelsize] * ncols
    height_ratios = [panelsize, diffpanelheight] * nrows + [legheight]
    width = sum(width_ratios) \
            * (1. + (len(width_ratios) - 1.) / len(width_ratios) * wspace)
    height = sum(height_ratios) \
             * (1. + (len(height_ratios) - 1.) / len(height_ratios) * hspace)

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(nrows=nrows * 2 + 1, ncols=ncols, hspace=hspace, 
                        wspace=wspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [[fig.add_subplot(grid[2 * i, j]) \
             for j in range(ncols)] for i in range(nrows)]
    daxes = [[fig.add_subplot(grid[2 * i + 1, j]) \
              for j in range(ncols)] for i in range(nrows)]
    lax = fig.add_subplot(grid[2 * nrows, :])
    fontsize = 12

    if not np.all([np.all([np.all([ykeys[i][j][k] == ykeys[0][0][k]\
                                   for k in range(nZ)]) \
                           for j in range(ncols)])\
                   for i in range(nrows)]):
        print('ykeys:\n', ykeys)
        raise RuntimeError('profile quantities differ between panels')
    if not np.all([np.all([np.all([wkeys[i][j][k] == wkeys[i][0][0]\
                                   for k in range(nZ)]) \
                           for j in range(ncols)])\
                   for i in range(nrows)]):
        print('wkeys:\n', wkeys)
        raise RuntimeError('weights differ in same row')
    if not np.all([np.all([np.all([np.isclose(zvals[0][j][0], zvals[i][j][zi],
                                              atol=3e-3, rtol=3e-3)\
                                   for zi in range(nZ)])\
                           for j in range(ncols)])\
                   for i in range(nrows)]):
        print('zvals:\n', zvals)
        raise RuntimeError('redshift values differ between rows')
    
    styleargs_lines = {'med': {'linewidth': 2., 'linestyle': 'solid'},
                       'sct': {'linewidth': 1., 'linestyle': 'dashed'},
                       'scthi': {'linewidth': 1., 'linestyle': 'dashed'},
                       'sctlo': {'linewidth': 1., 'linestyle': 'dashdot'},
                       'sct_med': {'elinewidth': 1., 'linestyle': 'none'},
                       }
    styleargs_shade = {'med': {'linewidth': 2., 'linestyle': 'solid'},
                       'sct': {'alpha': 0.2},
                       'sct_med': {'elinewidth': 1., 'linestyle': 'none'},
                      }

    for ri in range(nrows):
        for ci in range(ncols):
            ax = axes[ri][ci]
            dax = daxes[ri][ci]
            _ykeys = ykeys[ri][ci] 
            zval = zvals[ri][ci][0] # already checked equal in columns
            wkey = wkeys[ri][ci][0] # checked same in rows
            _rvirs = rvirs[ri][ci]
            _rvals = rvals[ri][ci]
            _yvals = yvals[ri][ci]
            _hists_raw = hists_raw[ri][ci]
            
            _yvals = [yv - np.log10(ynormopts[ykey]) \
                      for yv, ykey in zip(_yvals, _ykeys)]
            if ri == nrows - 1:
                dax.set_xlabel(xlabel, fontsize=fontsize)
                labelx = True
            else:
                labelx = False
            if ri == 0:
                ax.set_title(f'$z={zval:.1f}$')
            if ci == 0:
                ax.set_ylabel(ylabel, fontsize=fontsize)
                dax.set_ylabel('$\\Delta \\, \\log_{10} \\mathrm{Z}$',
                               fontsize=fontsize)
                labely = True
                ax.text(0.05, 0.05, wlabopts[wkey], fontsize=fontsize,
                        transform=ax.transAxes,
                        horizontalalignment='left', 
                        verticalalignment='bottom')
            else:
                labely = False
            ax.tick_params(which='both', direction='in', top=True,
                           left=True, bottom=True, right=True,
                           labelbottom=False, labelleft=labely,
                           labelsize=fontsize - 1.)
            dax.tick_params(which='both', direction='in', top=True,
                            left=True, bottom=True, right=True,
                            labelbottom=labelx, labelleft=labely,
                            labelsize=fontsize - 1.)
            ax.set_xscale('log')
            ax.grid(visible=True)
            dax.set_xscale('log')
            dax.grid(visible=True)

            pv = np.array([0.1, 0.5, 0.9])
            
            plos = []
            phis = []
            pmeds = []
            for rvir, rv, yv, hist_raw, ykey in \
                    zip(_rvirs, _rvals, _yvals, _hists_raw, _ykeys):
                    
                plo, pmed, phi = pu.percentiles_from_histogram(
                    10**hist_raw, yv, axis=1, percentiles=pv)
                rc = 0.5 * (rv[1:] + rv[:-1])
                plos.append(plo)
                phis.append(phi)
                pmeds.append(pmed)

                pstyle = 'lines'
                color = vcolors[ykey]
                vlab = vlabopts[ykey]
                medlab = 'med. of ' + vlab
                sctlab_lines = '10, 90% of ' + vlab
                sctlab_shade = '10-90% of ' + vlab
                if pstyle == 'lines':
                    ax.plot(rc, pmed, color=color, label=medlab,
                            **styleargs_lines['med'])
                    ax.plot(rc, phi, color=color, label=sctlab_lines,
                            **styleargs_lines['sct'])
                    ax.plot(rc, plo, color=color, 
                            **styleargs_lines['sct'])
                elif pstyle == 'shade':
                    ax.plot(rc, pmed, color=color, label=medlab,
                            **styleargs_shade['med'])
                    ax.fill_between(rc, plo, phi, color=color, 
                                    **styleargs_shade['sct'],
                                    label=sctlab_shade)
            if markrvir:
                yp_med = pu.linterpsolve(rc, pmed, rvir)
                ax.scatter([rvir], [yp_med], marker='o', c=color, s=15)
            
            compelt = 'Oxygen'
            ykey_base = ('sim-direct', f'ElementAbundance/{compelt}')
            base_ind = np.where([ykey == ykey_base for ykey in _ykeys])[0][0]
            for phi, plo, pmed, ykey in zip(phis, plos, pmeds, _ykeys):
                if ykey == ykey_base:
                    continue
                medlabel = (f'med. $\\log_{{10}}${vlabopts[ykey]}'
                            f' - med. $\\log_{{10}}${vlabopts[ykey_base]}')
                sctlolabel = (f'10% $\\log_{{10}}${vlabopts[ykey]} - '
                             f'10% $\\log_{{10}}${vlabopts[ykey_base]}')
                scthilabel = (f'90% $\\log_{{10}}${vlabopts[ykey]} - '
                              f'90% $\\log_{{10}}${vlabopts[ykey_base]}')
                dax.plot(rc, pmed - pmeds[base_ind], color=vcolors[ykey],
                         label=medlabel, **styleargs_lines['med'])
                dax.plot(rc, phi - phis[base_ind], color=vcolors[ykey],
                         label=scthilabel, **styleargs_lines['scthi'])
                dax.plot(rc, plo - plos[base_ind], color=vcolors[ykey],
                         label=sctlolabel, **styleargs_lines['sctlo']) 

            if ykey[0] == 'sim-direct' and \
                    ykey[1].startswith('ElementAbundance'):
                ymin, ymax = ax.get_ylim()
                ymin = max(ymin, -10.)
                ax.set_ylim(ymin, ymax)

        ylims = [axes[ri][i].get_ylim() for i in range(ncols)]
        ymin = min([ylim[0] for ylim in ylims])
        ymax = max([ylim[1] for ylim in ylims])
        [axes[ri][i].set_ylim((ymin, ymax)) for i in range(ncols)]

        ylims = [daxes[ri][i].get_ylim() for i in range(ncols)]
        ymin = min([ylim[0] for ylim in ylims])
        ymax = max([ylim[1] for ylim in ylims])
        [daxes[ri][i].set_ylim((ymin, ymax)) for i in range(ncols)]

    handles, _ = axes[0][0].get_legend_handles_labels()
    _handles, _ = daxes[0][0].get_legend_handles_labels()
    lax.set_axis_off()
    lax.legend(handles=handles + _handles, fontsize=fontsize - 1., 
               ncol=int(np.floor(0.5 * ncols)),
               loc='upper center', bbox_to_anchor=(0.5, 0.65))

    if figtitle is not None:
        fig.suptitle(figtitle, fontsize=fontsize)

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotset3dprof_eltZcomp_weightzgrid(fileset):
    if fileset == 'clean_set1-2':
        filedir = '/Users/nastasha/ciera/profiles/fire/clean_set1_set2/'
        filen_template = ('hist_{{yq}}_r3D_by_{{wq}}_{simname}'
                          '_snap{{snap}}_bins1_v1.hdf5')
        simnames = [('m13h113_m3e5_MHDCRspec1_fire3_fireBH'
                     '_fireCR1_Oct252021_crdiffc1_sdp1e-4'
                     '_gacc31_fa0.5_fcr1e-3_vw3000'),
                     ('m13h113_m3e4_MHD_fire3_fireBH'
                      '_Sep182021_hr_crdiffc690_sdp1e-4'
                      '_gacc31_fa0.5'),
                     ('m13h113_m3e5_MHD_fire3_fireBH'
                      '_Sep182021_crdiffc690_sdp1e10'
                      '_gacc31_fa0.5'),
                     ('m13h206_m3e5_MHDCRspec1_fire3_fireBH'
                      '_fireCR1_Oct252021_crdiffc1_sdp1e-4'
                      '_gacc31_fa0.5_fcr1e-3_vw3000'),
                     ('m13h206_m3e4_MHD_fire3_fireBH'
                      '_Sep182021_hr_crdiffc690_sdp3e-4'
                      '_gacc31_fa0.5'),
                     ('m13h206_m3e5_MHD_fire3_fireBH'
                      '_Sep182021_crdiffc690_sdp1e10'
                      '_gacc31_fa0.5'),
                     ('m12f_m6e4_MHDCRspec1_fire3_fireBH'
                      '_fireCR1_Oct252021_crdiffc1_sdp1e-4'
                      '_gacc31_fa0.5_fcr1e-3_vw3000'),
                     ('m12f_m7e3_MHD_fire3_fireBH'
                      '_Sep182021_hr_crdiffc690_sdp2e-4'
                      '_gacc31_fa0.5'),
                     ('m12f_m7e3_MHD_fire3_fireBH'
                      '_Sep182021_hr_crdiffc690_sdp1e10'
                      '_gacc31_fa0.5'),
                     ]
        simlabs = ['m13h113 AGN-CR', 
                   'm13h113 AGN-noCR',
                   'm13h113 noBH',
                   'm13h206 AGN-CR',
                   'm13h206 AGN-noCR',
                   'm13h206 noBH',
                   'm12f AGN-CR',
                   'm12f AGN-noCR',
                   'm12f noBH',
                   ]
        srst = ['m13h113_m3e5', 'm13h206_m3e5', 'm12f_m6e4']
        hrst = ['m13h113_m3e4', 'm13h206_m3e4', 'm12f_m7e3']
        snaps_sr = [50, 49, 48, 47, 46, 45]
        snaps_hr = [258, 240, 224, 210, 197, 186]
        #redshifts = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        simlen = len(simnames)
        fillkw_valss = [[{'yq': 'Neon'},
                         {'yq': 'Oxygen'},
                         {'yq': 'Magnesium'},
                        ]] * simlen
        fillkw_weightss = [[{'wq': 'Mass'},
                            {'wq': 'Volume'},
                            {'wq': 'H1'},
                           ]] * simlen
        fillkw_zs = [[{'snap': snhr} \
                           if np.any([sim.startswith(st) for st in hrst])\
                       else {'snap': snsr} \
                           if np.any([sim.startswith(st) for st in srst])\
                       else None
                      for snhr, snsr in zip(snaps_hr, snaps_sr)] \
                     for sim in simnames]
        filen_templates = \
            [filedir + filen_template.format(simname=sim)\
             for sim in simnames]
        outdir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/3dprof/'
        outtemplate = 'prof3D_Zcomp_M_V_HI_{simlab}.pdf'
        outnames = [outdir + outtemplate.format(
                        simlab=simlab.replace(' ', '_'))\
                    for simlab in simlabs]
        figtitles = simlabs
    
    for filen_template, fillkw_weights, fillkw_vals, fillkw_z, figtitle, \
            outname in zip(filen_templates, fillkw_weightss, fillkw_valss,
                           fillkw_zs, figtitles, outnames):
        plot3dprof_eltZcomp_weightzgrid(filen_template, fillkw_weights, 
                                        fillkw_vals, fillkw_z,
                                        figtitle=figtitle, outname=outname)