'''
In each panel, plot an individual snapshot's property/radius 
histograms, and the weighted property median and scatter as a 
function or radius. Plot columns vary the weight, and rows show the
properties density, temperature, and metallicity.
There are different plots for different weight/property sets, 
redshifts, ICs, and physics models.
'''

import h5py
import numpy as np

import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt

import fire_an.makeplots.plot_utils as pu
import fire_an.utils.constants_and_units as c
import fire_an.utils.opts_locs as ol

def plot3dprof_weightvalgrid(filen_template, fillkw_weights, fillkw_vals,
                             figtitle=None, outname=None):
    
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
    hists = []
    hists_raw = []
    rvals = []
    yvals = []
    ykeys = []
    wkeys = []
    rvirs = []
    runits = None
    for vdict in fillkw_vals:
        _hists = []
        _rvals = []
        _yvals = []
        _ykeys = []
        _wkeys = []
        _rvirs = []
        _hists_raw = []
        for wdict in fillkw_weights:
            _tdict = wdict.copy()
            _tdict.update(vdict)
            filen = filen_template.format(**_tdict)
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
               
                deltax = np.diff(rv)
                if yv[0] == -np.inf:
                    yv[0] = 2. * yv[1] - yv[2]
                if yv[-1] == np.inf:
                    yv[-1] = 2. * yv[-2] - yv[-3]
                deltay = np.diff(yv)
                hist_norm = hist_raw - np.log10(np.sum(10**hist_raw))
                hist_norm = hist_norm - np.log10(deltax)[:, np.newaxis]
                hist_norm = hist_norm - np.log10(deltay)[np.newaxis, :]
                
            if runits is not None:
                if _runits != runits:
                    raise ValueError('mismatch in r3D units')
            if _runits == 'pkpc':
                rvir = rvir_cm / (c.cm_per_mpc * 1e-3)
            elif _runits == 'Rvir':
                rvir = 1.
            
            _hists_raw.append(hist_raw)
            _hists.append(hist_norm)
            _rvals.append(rv)
            _yvals.append(yv)
            _ykeys.append(ykey)
            _wkeys.append(wkey)
            _rvirs.append(rvir)
        hists_raw.append(_hists_raw)
        hists.append(_hists)
        rvals.append(_rvals)
        yvals.append(_yvals)
        ykeys.append(_ykeys)
        wkeys.append(_wkeys)
        rvirs.append(_rvirs)
    print(ykeys)
    print(wkeys)

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
    cax = fig.add_subplot(grid[:2, ncols])
    fontsize = 12

    cmap = 'gist_yarg'
    vmin = min([min([np.min(hist[np.isfinite(hist)]) \
                     for hist in _])\
                for _ in hists])
    vmax = max([max([np.max(hist) for hist in _]) for _ in hists])

    if not np.all([np.all([key == colkeys[0] for key in colkeys])\
                   for colkeys in ykeys]):
        print('ykeys:\n', ykeys)
        raise RuntimeError('y-axis quantities differ in same row')
    if not np.all([np.all([wkeys[0][i] == colkeys[i] for i in range(ncols)])\
                   for colkeys in wkeys]):
        print('wkeys:\n', wkeys)
        raise RuntimeError('weights differ in same column')
    
    for ri in range(nrows):
        for ci in range(ncols):
            hist_raw = hists_raw[ri][ci]
            hist = hists[ri][ci]
            ykey = ykeys[ri][ci]
            wkey = wkeys[ri][ci]
            rvir = rvirs[ri][ci]
            rv = rvals[ri][ci]
            yv = yvals[ri][ci]
            ax = axes[ri][ci]
            
            yv = yv - np.log10(ynormopts[ykey])
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
            img = ax.pcolormesh(rv, yv, hist.T, shading='flat', 
                                rasterized=True, cmap=cmap, vmin=vmin,
                                vmax=vmax)
            
            pv = np.array([0.1, 0.5, 0.9])
            plo, pmed, phi = pu.percentiles_from_histogram(10**hist_raw, yv, 
                                                           axis=1,
                                                           percentiles=pv)
            rc = 0.5 * (rv[1:] + rv[:-1])
            ax.plot(rc, pmed, color='red', linestyle='dotted', linewidth=1.)
            ax.plot(rc, plo, color='red', linestyle='dashed', linewidth=1.)
            ax.plot(rc, phi, color='red', linestyle='dashed', linewidth=1.)

            if markrvir:
                yp_med = pu.linterpsolve(rc, pmed, rvir)
                ax.scatter([rvir], [yp_med], marker='o', c='red', s=10)
        
        ylims = [axes[ri][i].get_ylim() for i in range(ncols)]
        ymin = min([ylim[0] for ylim in ylims])
        ymax = max([ylim[1] for ylim in ylims])
        [axes[ri][i].set_ylim((ymin, ymax)) for i in range(ncols)]
    
    plt.colorbar(img, cax=cax)
    clabel = ('$\\log_{10} \\; \\Delta(\\mathrm{weight}) \\, / \\,'
              '\\mathrm{total}\\,\\mathrm{weight} \\, / \\,'
              '\\Delta(\\mathrm{r}) \\, / \\, \\Delta(\\log_{10}\\, y)$')
    cax.set_ylabel(clabel, fontsize=fontsize)

    if figtitle is not None:
        fig.suptitle(figtitle, fontsize=fontsize)

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')
            
def plotset3dprof_weightvalgrid(fileset):
    if fileset == 'test1':
        filedir = '/Users/nastasha/ciera/profiles/fire/clean_set1_set2/'
        filen_template = ('hist_{yq}_r3D_by_{wq}_m13h206_m3e5_MHD_fire3'
                          '_fireBH_Sep182021_crdiffc690_sdp1e10'
                          '_gacc31_fa0.5_snap50_bins1_v1.hdf5')
        filen_templates = [filedir + filen_template]
        fillkw_valss = [[{'yq': 'Density'},
                         {'yq': 'Temperature'},
                         {'yq': 'Neon'},
                         {'yq': 'Oxygen'},
                         ]]
        fillkw_weightss = [[{'wq': 'Mass'}]]
        outdir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/3dprof/'
        outnames = [outdir + 'test1_plot_3dprof.pdf']
        figtitles = ['m13h206 noBH, z=0.5']
    if fileset == 'clean_set1':
        filedir = '/Users/nastasha/ciera/profiles/fire/clean_set1_set2/'
        filen_template = ('hist_{{yq}}_r3D_by_{{wq}}_{simname}'
                          '_snap{snap}_bins1_v1.hdf5')
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
        snaps_sr = [45, 46, 47, 48, 49, 50]
        snaps_hr = [186, 197, 210, 224, 240, 258]
        redshifts = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
        simsnaplen = len(snaps_hr + snaps_hr) * len(simnames)
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
                       ] * simsnaplen
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
                          ] * simsnaplen
        weightsvals = ['M_V_H1', 'M_V_Ne8', 'M_V_O6', 'M_V_Mg10']
        filen_templates = \
            [filedir + filen_template.format(
                simname=sim,
                snap=snhr if np.any([sim.startswith(st) for st in hrst])\
                     else snsr if np.any([sim.startswith(st) for st in srst])\
                     else None)\
             for sim in simnames for snhr, snsr in zip(snaps_hr, snaps_sr)\
             for i in range(len(weightsvals))]
        outdir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/3dprof/'
        outtemplate = 'prof3D_indivdist_{wv}_{simlab}_z{z}.pdf'
        outnames = [outdir + outtemplate.format(
                        simlab=simlab.replace(' ', '_'),
                        z=f'{z:.1f}'.replace('.', 'p'),
                        wv=wv)\
                    for simlab in simlabs for z in redshifts \
                    for wv in weightsvals]
        figtitles = [f'{simlab}, z={z:.1f}'\
                     for simlab in simlabs for z in redshifts\
                     for i in range(len(weightsvals))]
    
    for filen_template, fillkw_weights, fillkw_vals, figtitle, outname \
            in zip(filen_templates, fillkw_weightss, fillkw_valss,
                   figtitles, outnames):
        print(fillkw_weights)
        print(fillkw_vals)
        plot3dprof_weightvalgrid(filen_template, fillkw_weights, fillkw_vals,
                                 figtitle=figtitle, outname=outname)
 

