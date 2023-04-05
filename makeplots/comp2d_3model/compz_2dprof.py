
import h5py
import matplotlib.collections as mcol
import matplotlib.colors as mcolors
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

def plotcomp_zev_ic(filen_template, paxfills, zfills,
                    physfills, iclabels, title=None,
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
        the physfills.
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
                  for zfill in zlist2] # snapshot
                 for zlist2 in zlist1] # redshift bin
                for pfill, zlist1 in zip(physfills, zfills)] # ICS
                
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
    ## axes outer to inner: 
    # projection axis, snapshot, redshift set, proj. ax.
    # should not depend on projection axis
    if not np.all([np.all([np.all([[rvirs[i][j][k][0] == rvirs[i][j][k][l]
                      for l in range(len(paxfills))]
                     for k in range(len(zfills[i][j]))])
                    for j in range(len(zfills[i]))])
                   for i in range(len(physfills))]):
        msg = ('Different projection axes (innermost index) list different'
               f'virial radii [pkpc]:\n{rvirs}')
        raise ValueError(msg)
    if not np.all([np.all([np.all([[mvirs[i][j][k][0] == mvirs[i][j][k][l]
                      for l in range(len(paxfills))]
                     for k in range(len(zfills[i][j]))])
                    for j in range(len(zfills[i]))])
                   for i in range(len(physfills))]):
        msg = ('Different projection axes (innermost index) list different'
               f'virial masses [Msun]:\n{mvirs}')
        raise ValueError(msg)
    if not np.all([np.all([np.all([[axis3s[0][0][0][l] == axis3s[i][j][k][l]
                      for l in range(len(paxfills))]
                     for k in range(len(zfills[i][j]))])
                    for j in range(len(zfills[i]))])
                   for i in range(len(physfills))]):
        msg = ('Different ics/redshifts list different'
               f'projection axes (all but innermost axis):\n{axis3s}')
        raise ValueError(msg)
    if not np.all([np.all([np.all([[np.isclose(redshifts[0][j][k][0],
                                               redshifts[i][j][k][l],
                                               rtol=1e-2, atol=1e-3)
                      for l in range(len(paxfills))]
                     for k in range(len(zfills[i][j]))])
                    for j in range(len(zfills[i]))])
                   for i in range(len(physfills))]):
        msg = ('Different ics/proj. ax list different'
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
                * (1. + (len(height_ratios) - 1.) / len(height_ratios) \
                * hspace)
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
    
    pv = ['perc-0.1', 'perc-0.5', 'perc-0.9']
    maxscattershow = 3
    alpha = 0.3
    fontsize = 12
    nlines = len(redshifts[0])
    _colors = tc.tol_cmap('rainbow_discrete', nlines)
    zcolors = _colors(np.linspace(1. - 0.5 / nlines, 0.5 / nlines, nlines))
    # zip/list comprehension issues or something when using 
    # tol_cmap outputs directly
    zcolors = [mcolors.to_rgb(col) for col in zcolors]
    zvals = [set() for _ in range(nlines)] 
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
        for zsi in range(len(zfills[ici])):
            color = zcolors[zsi]
            zvals[zsi] |= set([val for _ in redshifts[ici][zsi] for val in _])
            _fns = [fn for l1 in filens[ici][zsi] for fn in l1]
            _rvs = [rv for l1 in rvirs[ici][zsi] for rv in l1]
            rvmin = np.min(_rvs)
            rvmax = np.max(_rvs)
            if rvmin < rvmax * 0.99:
                ax.axvspan(rvmin, rvmax, alpha=0.2, color=color)
            else:
                rvav = 0.5 * (rvmin + rvmax)
                ax.axvline(rvav, color=color, linestyle='dotted')
            rbins = np.linspace(0., 2. * rvmin, 50)
            rc = 0.5 * (rbins[:-1] + rbins[1:])
            plo, pmed, phi = g2d.get_profile_massmap(_fns, rbins, 
                                                     rbin_units='pkpc',
                                                     profiles=pv)
            ax.plot(rc, pmed, color=color, linestyle='solid',
                    linewidth=1.5)
            if (len(zfills[ici]) <= maxscattershow 
                or zsi % (len(zfills) // maxscattershow) == 0):
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
    cvals = [zcolors[i] for i in range(len(zfills[0]))]
    lcs = mcol.LineCollection(line * len(zcolors), linestyle='solid', 
                              linewidth=1.5, colors=cvals)
    zhandles = [lcs,
                mpatch.Patch(label='perc. 10-90', linewidth=0.5, 
                             color='brown', alpha=alpha)
                ]
    plabels = ['median', '10-90%']
    zlabels = [f'$z={min(_zv):.1f} \endash {max(_zv):.1f}$' for _zv in zvals]
    zhandles += [mlines.Line2D((), (), linewidth=2., linestyle='solid',
                               label=zlab, 
                               color=zcolors[zsi])
                 for zsi, zlab in enumerate(zlabels)]
    plabels += zlabels
    lax.legend(zhandles, plabels, 
               fontsize=fontsize, ncol=ncols_legend,
               handler_map={type(lcs): pu.HandlerDashedLines()},
               bbox_to_anchor=(0.5, 0.7), loc='upper center')
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotsetcomp_zev_ic(fileset='set3_model3', numzsets=3):

    if fileset == 'set3_model3':
        # quest
        outdir = '/projects/b1026/nastasha/imgs/2dcomp_set3_model3/'
        fdir = '/projects/b1026/nastasha/maps/set3_model3/'
        ftemp = ('coldens_{qty}_{{simname}}_snap{{snapnum}}_'
                 'shrink-sph-cen_BN98_2rvir_{{pax}}-proj_v3.hdf5')
        simnames_all = sl.m12_hr_all1 + sl.m12_sr_all1 + \
                       sl.m13_hr_all1 + sl.m13_sr_all1
        simgroups = []
        iclabelss = []
        sgtitles = []
        # sort into lists by halo mass group, physics model
        for simname in simnames_all:
            ic = simname.split('_')[0]
            physlab = ('AGN-CR' if 'MHDCRspec1' in simname
                       else 'noBH' if 'sdp1e10' in simname
                       else 'AGN-noCR')
            hmlab = ic[:3] # e.g. m12f -> m12, m13h007 -> m13
            title = hmlab + ', ' + physlab
            if title in sgtitles:
                si = np.where([title == _tl for _tl in sgtitles])[0][0]
                simgroups[si].append(simname)
                iclabelss[si].append(ic)
            else:
                simgroups.append([simname])
                iclabelss.append([ic])
                sgtitles.append(title)
        sgfills = [title.replace(', ', '_') for title in sgtitles]
        sgfills = [title.replace('-', '_') for title in sgfills]

        all_hr = sl.m12_hr_all1 + sl.m13_hr_all1
        all_sr = sl.m12_sr_all1 + sl.m13_sr_all1
        nsn_hr = len(sl.snaps_hr)
        numz_hr = [nsn_hr // numzsets + 1 if i < nsn_hr % numzsets 
                   else nsn_hr // numzsets
                   for i in range(numzsets)]
        slices_hr = [slice(sum(numz_hr[:i]), sum(numz_hr[:i + 1]), None)
                     for i in range(numzsets)]
        nsn_sr = len(sl.snaps_sr)
        numz_sr = [nsn_sr // numzsets + 1 if i < nsn_sr % numzsets 
                   else nsn_sr // numzsets
                   for i in range(numzsets)]
        slices_sr = [slice(sum(numz_sr[:i]), sum(numz_sr[:i + 1]), None)
                     for i in range(numzsets)]
        zfillsets = [[[[{'snapnum': snap} 
                        for snap in sl.snaps_hr[slices_hr[zsi]]]
                       if simname in all_hr else 
                       [{'snapnum': snap} 
                        for snap in sl.snaps_sr[slices_sr[zsi]]]
                       if simname in all_sr else 
                       None
                       for zsi in range(numzsets)]
                      for simname in simlist] 
                     for simlist in simgroups]

        _afills = ['x', 'y', 'z']
        paxfills = [{'pax': val} for val in _afills]
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
        outname = (f'comp_2dprof_zev_{numzsets}bins_{{qty}}_{{sglab}}.pdf')
        simlists_ext = [[{'simname': simn} for simn in siml] 
                         for siml in simgroups] \
                       * len(qtys)
        iclabels_ext = iclabelss * len(qtys)
        zfills_ext = zfillsets * len(qtys)
        sgtitles_ext = sgtitles * len(qtys)
        sglabels_ext = sgfills * len(qtys)
        paxfills_ext = paxfills
        qtylabels_ext = [lab for lab in qtylabs 
                         for i in range(len(simgroups))]
        ylabels_ext = [lab for lab in qtyylabs for i in range(len(simgroups))]
        qtys_ext = [lab for lab in qtys for i in range(len(simgroups))]

    for (simfills, qty, ylabel, qtylabel, zfills, iclabs, sgtitle, sglab) \
            in zip(simlists_ext, qtys_ext, ylabels_ext, qtylabels_ext,
                   zfills_ext, iclabels_ext, sgtitles_ext, sglabels_ext):
        filen_template = fdir + ftemp.format(qty=qty)
        title = f'{qtylabel}, {sgtitle}'
        _outname = outdir + outname.format(sglab=sglab, 
                                           qty=qty.replace('-', '_'))

        plotcomp_zev_ic(filen_template, paxfills_ext, zfills,
                        simfills, iclabs, title=title,
                        outname=_outname, ylabel=ylabel)
        #plt.show()
        plt.close() # don't overload memory with too many plots