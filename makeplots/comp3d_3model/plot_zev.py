'''
In each panel, plot an individual IC+physics model at different
snapshots. Colored lines show different redshifts, and the black line
shows an ensemble profile. Plot columns vary the weight, and rows show
the properties density, temperature, and metallicity.
There are different plots for different weight/property sets, ICs, and
physics models.
'''

import h5py
import numpy as np

import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt

import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
import fire_an.utils.constants_and_units as c
import fire_an.utils.opts_locs as ol

def plot3dprof_zev_weightvalgrid(filen_template, fillkw_weights, fillkw_vals,
                                 fillkw_z, figtitle=None, outname=None):
    
    rholab = ('$\\log_{10}\\,\\rho \\; '
              '[\\mathrm{H}(X=0.752) \\, \\mathrm{cm}^{-3}]$')
    tlab = '$\\log_{10}\\,\\mathrm{T} \\; [\\mathrm{K}]$'
    zolab = ('$\\log_{10}\\,\\mathrm{Z}_\\mathrm{O} \\;'
             ' [\\mathrm{Z}_{\\mathrm{O}, \\odot}]$')
    znelab = ('$\\log_{10}\\,\\mathrm{Z}_\\mathrm{Ne} \\;'
              ' [\\mathrm{Z}_{\\mathrm{Ne}, \\odot}]$')
    zmglab = ('$\\log_{10}\\,\\mathrm{Z}_\\mathrm{Mg} \\;'
              ' [\\mathrm{Z}_{\\mathrm{Mg}, \\odot}]$')
    ylabopts = {('sim-direct', 'Density'): rholab,
                ('sim-direct', 'Temperature'): tlab,
                ('sim-direct', 'ElementAbundance/Oxygen'): zolab,
                ('sim-direct', 'ElementAbundance/Neon'): znelab,
                ('sim-direct', 'ElementAbundance/Magnesium'): zmglab,
               }
    rhonorm = c.atomw_H * c.u / 0.752
    tnorm = 1.
    zonorm = ol.solar_abunds_ea['oxygen']
    znenorm = ol.solar_abunds_ea['neon']
    zmgnorm = ol.solar_abunds_ea['magnesium']
    ynormopts = {('sim-direct', 'Density'): rhonorm,
                 ('sim-direct', 'Temperature'): tnorm,
                 ('sim-direct', 'ElementAbundance/Oxygen'): zonorm,
                 ('sim-direct', 'ElementAbundance/Neon'): znenorm,
                 ('sim-direct', 'ElementAbundance/Magnesium'): zmgnorm,
                }
    wlabopts = {'Mass': 'Mass-weighted',
                'Volume': 'Volume-weighted',
                ('ion', 'O6'): 'O VI-weighted',
                ('ion', 'Ne8'): 'Ne VIII-weighted',
                ('ion', 'Mg10'): 'Mg X-weighted',
                ('ion', 'H1'): 'H I-weighted',
                }
    nz = len(fillkw_z)
    hists_raw = []
    rvals = []
    yvals = []
    ykeys = []
    wkeys = []
    rvirs = []
    zvals = []
    runits = None
    for vdict in fillkw_vals:
        _rvals = []
        _yvals = []
        _ykeys = []
        _wkeys = []
        _rvirs = []
        _zvals = []
        _hists_raw = []
        for wdict in fillkw_weights:
            _tdict = wdict.copy()
            _tdict.update(vdict)
            __rvals = []
            __yvals = []
            __ykeys = []
            __wkeys = []
            __rvirs = []
            __zvals = []
            __hists_raw = []
            for zdict in fillkw_z:
                __tdict = _tdict.copy()
                __tdict.update(zdict)
                filen = filen_template.format(**__tdict)
                with h5py.File(filen, 'r') as f:
                    rvir_cm = f['Header/inputpars/halodata'].attrs['Rvir_cm']
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
                __rvirs.append(rvir)
                __zvals.append(redshift)
            _hists_raw.append(__hists_raw)
            _rvals.append(__rvals)
            _yvals.append(__yvals)
            _ykeys.append(__ykeys)
            _wkeys.append(__wkeys)
            _rvirs.append(__rvirs)
            _zvals.append(__zvals)
        hists_raw.append(_hists_raw)
        rvals.append(_rvals)
        yvals.append(_yvals)
        ykeys.append(_ykeys)
        wkeys.append(_wkeys)
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
    caxwidth = 0.5
    ncols = len(fillkw_weights)
    nrows = len(fillkw_vals)
    wspace = 0.
    hspace = 0. 
    width_ratios = [panelsize] * ncols + [caxwidth]
    height_ratios = [panelsize] * nrows
    width = sum(width_ratios) \
            * (1. + (len(width_ratios) - 1.) / len(width_ratios) * wspace)
    height = sum(height_ratios) \
             * (1. + (len(height_ratios) - 1.) / len(height_ratios) * hspace)

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(nrows=nrows, ncols=ncols + 1, hspace=hspace, 
                        wspace=wspace, width_ratios=width_ratios,
                        height_ratios=height_ratios)
    axes = [[fig.add_subplot(grid[i, j]) \
             for j in range(ncols)] for i in range(nrows)]
    cax = fig.add_subplot(grid[0, ncols])
    fontsize = 12

    if not np.all([np.all([np.all([key == rowkeys[0][0] \
                                   for key in key_list])\
                           for key_list in rowkeys])\
                   for rowkeys in ykeys]):
        print('ykeys:\n', ykeys)
        raise RuntimeError('y-axis quantities differ in same row')
    if not np.all([np.all([np.all([wkeys[0][i][0] == colkeys[i][j] \
                                   for j in range(nz)])\
                           for i in range(ncols)])\
                   for colkeys in wkeys]):
        print('wkeys:\n', wkeys)
        raise RuntimeError('weights differ in same column')
    if not np.all([np.all([np.allclose(zvals[0][0], zvals[j][i])\
                           for i in range(ncols)])\
                   for j in range(nrows)]):
        print('zvals:\n', zvals)
        raise RuntimeError('redshift values differ between panels')
    
    _colors = tc.tol_cmap('rainbow_discrete', nz)
    zticks = np.linspace(0.5 / nz, 1. - 0.5 / nz, nz)
    colors = _colors(zticks)
    # zip/list comprehension issues or something when using 
    # tol_cmap outputs directly
    colors = [mcolors.to_rgb(col) for col in colors]
    
    norm = mcolors.Normalize(vmin=0., vmax=1.)
    plt.colorbar(mcm.ScalarMappable(cmap=_colors, norm=norm), cax=cax, 
                 extend='neither', orientation='vertical')
    cax.set_ylabel('redshift', fontsize=fontsize)
    cax.set_yticks(zticks)
    cax.set_yticklabels([f'{z:.1f}' for z in zvals[0][0]])
    cax.tick_params(labelsize=fontsize - 1.)

    for ri in range(nrows):
        for ci in range(ncols):
            ax = axes[ri][ci]
            ykey = ykeys[ri][ci][0] # already checked equal in panel
            wkey = wkeys[ri][ci][0] # already checked equal in panel
            _rvirs = rvirs[ri][ci]
            _rvals = rvals[ri][ci]
            _yvals = yvals[ri][ci]
            _hists_raw = hists_raw[ri][ci]
            
            _yvals = [yv - np.log10(ynormopts[ykey]) for yv in _yvals]
            if ri == nrows - 1:
                ax.set_xlabel(xlabel, fontsize=fontsize)
                labelx = True
            else:
                labelx = False
            if ri == 0:
                ax.set_title(wlabopts[wkey], fontsize=fontsize)
            if ci == 0:
                ax.set_ylabel(ylabopts[ykey], fontsize=fontsize)
                labely = True
            else:
                labely = False
            ax.tick_params(which='both', direction='in', top=True,
                           left=True, bottom=True, right=True,
                           labelbottom=labelx, labelleft=labely,
                           labelsize=fontsize - 1.)
            
            pv = np.array([0.1, 0.5, 0.9])
            rc_all = None
            plos = []
            pmeds = []
            phis = []
            for rvir, rv, yv, hist_raw, color in\
                    zip(_rvirs, _rvals, _yvals, _hists_raw, colors):
                plo, pmed, phi = pu.percentiles_from_histogram(10**hist_raw, 
                                     yv, axis=1, percentiles=pv)
                rc = 0.5 * (rv[1:] + rv[:-1])
                ax.plot(rc, pmed, color=color, linestyle='solid', 
                        linewidth=2.)
                ax.plot(rc, plo, color=color, linestyle='dashed', 
                        linewidth=0.7)
                ax.plot(rc, phi, color=color, linestyle='dashed', 
                        linewidth=0.7)

                if markrvir:
                    yp_med = pu.linterpsolve(rc, pmed, rvir)
                    ax.scatter([rvir], [yp_med], marker='o', c=color, s=15)
                if rc_all is None:
                    rc_all = rc
                elif not np.allclose(rc_all, rc):
                    msg = 'radial bins mismatched between redshifts'
                    raise RuntimeError(msg)
                plos.append(plo)
                pmeds.append(pmed)
                phis.append(phi)
            plos = np.array(plos)
            pmeds = np.array(pmeds)
            phis = np.array(phis)
            plo_med = np.quantile(plos, 0.5, axis=0)
            phi_med = np.quantile(phis, 0.5, axis=0)
            pmed_med = np.quantile(pmeds, 0.5, axis=0)

            ax.plot(rc_all, pmed_med, color='black', linestyle='solid', 
                    linewidth=2.)
            ax.fill_between(rc_all, plo_med, phi_med, color='black', 
                            alpha=0.3)
            ax.set_xscale('log')

        ylims = [axes[ri][i].get_ylim() for i in range(ncols)]
        ymin = min([ylim[0] for ylim in ylims])
        ymax = max([ylim[1] for ylim in ylims])
        [axes[ri][i].set_ylim((ymin, ymax)) for i in range(ncols)]

    if figtitle is not None:
        fig.suptitle(figtitle, fontsize=fontsize)

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotset3dprof_zev_weightvalgrid(fileset):
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
        fillkw_valss = [[{'yq': 'Density'},
                         {'yq': 'Temperature'},
                         {'yq': 'Neon'},
                         {'yq': 'Oxygen'},
                         {'yq': 'Magnesium'},
                        ],
                        [{'yq': 'Density'},
                         {'yq': 'Temperature'},
                         {'yq': 'Neon'},
                        ],
                        [{'yq': 'Density'},
                         {'yq': 'Temperature'},
                         {'yq': 'Oxygen'},
                        ],
                        [{'yq': 'Density'},
                         {'yq': 'Temperature'},
                         {'yq': 'Magnesium'},
                        ]
                       ] * simlen
        fillkw_weightss = [[{'wq': 'Mass'},
                            {'wq': 'Volume'},
                            {'wq': 'H1'},
                           ],
                           [{'wq': 'Mass'},
                            {'wq': 'Volume'},
                            {'wq': 'Ne8'},
                           ],
                           [{'wq': 'Mass'},
                            {'wq': 'Volume'},
                            {'wq': 'O6'},
                           ],
                           [{'wq': 'Mass'},
                            {'wq': 'Volume'},
                            {'wq': 'Mg10'},
                           ],
                          ] * simlen
        weightsvals = ['M_V_H1', 'M_V_Ne8', 'M_V_O6', 'M_V_Mg10']
        fillkw_zs = [[{'snap': snhr} \
                           if np.any([sim.startswith(st) for st in hrst])\
                       else {'snap': snsr} \
                           if np.any([sim.startswith(st) for st in srst])\
                       else None
                      for snhr, snsr in zip(snaps_hr, snaps_sr)] \
                     for sim in simnames for i in range(len(weightsvals))]
        filen_templates = \
            [filedir + filen_template.format(simname=sim)\
             for sim in simnames for i in range(len(weightsvals))]
        outdir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/3dprof/'
        outtemplate = 'prof3D_zev_{wv}_{simlab}.pdf'
        outnames = [outdir + outtemplate.format(
                        simlab=simlab.replace(' ', '_'),
                        wv=wv)\
                    for simlab in simlabs for wv in weightsvals]
        figtitles = [simlab for simlab in simlabs \
                     for i in range(len(weightsvals))]
    
    for filen_template, fillkw_weights, fillkw_vals, fillkw_z, figtitle, \
            outname in zip(filen_templates, fillkw_weightss, fillkw_valss,
                           fillkw_zs, figtitles, outnames):
        plot3dprof_zev_weightvalgrid(filen_template, fillkw_weights, 
                                     fillkw_vals, fillkw_z, figtitle=figtitle, 
                                     outname=outname)