import h5py
import matplotlib as mpl
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np


import fire_an.makeplots.ionmech.plotpds as ppd
import fire_an.makeplots.plot_utils as pu
import fire_an.makeplots.tol_colors as tc
from fire_an.ionrad.ion_utils import Linetable_PS20
import fire_an.simlists as sl 

def readpd(filen, rrange_rvir=(0.1, 1.)):
    with h5py.File(filen, 'r') as f:
        linpd = f['histogram/histogram'][:]
        if bool(f['histogram'].attrs['log']):
            linpd = 10**linpd
        rbins_rvir = f['axis_0/bins'][:]
        ri0 = np.where(np.isclose(rbins_rvir, rrange_rvir[0]))[0][0]
        ri1 = np.where(np.isclose(rbins_rvir, rrange_rvir[1]))[0][0]
        linpd = np.sum(linpd[ri0 : ri1, :, :], axis=0)
        norm = np.sum(linpd)
        linpd *= 1. / norm # fraction of ions in each bin
        tbins = rbins_rvir = f['axis_1/bins'][:]
        nHbins = rbins_rvir = f['axis_2/bins'][:]
        linpddens = linpd / np.diff(tbins)[:, np.newaxis] \
                    /  np.diff(nHbins)[np.newaxis, :]
    return {'linpd': linpd, 'linpddens': linpddens,
            'Tbins': tbins, 'nHbins': nHbins}

def checkNe8pds_physmodel(physlabels, ics, redshift=1.0,
                          rrange_rvir_main=(0.1, 1.),
                          rranges_rvir_sub=((0.15, 0.25), 
                                            (0.45, 0.55),
                                            (0.9, 1.0))):
    '''
    For the m12s, in the halo overall and at 0.9 -- 1.0 Rvir,
    Ne8 often seems to be close to the CIE/PIE/P+CIE split point.
    Here, getting the 'corner' right will be better, so putting the
    CIE/mixed cutoff at Clayton Strawn's transition density, if it 
    exists, at given T, is probably a good idea. Interpolating 
    between the points to classify general rho/T will be a pain though.

    For the m13s, it probably doesn't matter as much, since higher T Ne8
    is more typical there, including at large radii.
    '''
    ddir = '/projects/b1026/nastasha/hists/phasediagrams_all2/'
    fbase = ('hist_rcen_temperature_density_by_{wtstr}_'
             '{simname}_snap{snapnum}'
             '_bins1_v1.hdf5')
    filetemp = ddir + fbase
    
    mdir = '/projects/b1026/nastasha/imgs/pcie/cpie_summaries/'
    simslabel = '_'.join([f'{ic}_{mod}' for ic, mod in zip(ics, physlabels)])
    outname = (f'phasediag_Ne8_Neon_mass_{simslabel}'
               f'_z{redshift:.1f}')
    outname = mdir + outname.replace('.', 'p') + '.pdf'

    zopts = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
    snapi = np.where(np.isclose(zopts, redshift))[0][0]
    snap_sr = sl.snaps_sr[snapi]
    snap_hr = sl.snaps_hr[snapi]
    snap_f2 = sl.snaps_f2md[snapi]
    sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2
    sims_f2 = sl.m12_f2md

    allsims = sl.m12_sr_all2 + sl.m12_hr_all2 + sl.m12_f2md + \
              sl.m13_sr_all2 + sl.m13_hr_all2
    simnames_base = [sn for sn in allsims if (sn not in sl.buglist1)]
    simnames = [[sn for sn in simnames_base 
                 if sl.ic_from_simname(sn) == ic 
                    and sl.physlabel_from_simname(sn) == pl][0]
                for ic, pl in zip(ics, physlabels)]
    snaps = [snap_f2 if sn in sims_f2
             else snap_hr if sn in sims_hr
             else snap_sr if sn in sims_sr
             else None
             for sn in simnames]
    
    ion = 'Ne8'
    wts = ['Ne8', 'gasmass'] 
    # 'Neon': very similar to total gas in warm-hot CGM
    wtlabels = {'gasmass': 'gas',
                'Neon': 'Ne',
                'Ne8': 'Ne8',
                'gasvol': 'Vol.'}
    colors = tc.tol_cset('vibrant')
    #wtcolors = {wt: col for wt, col in zip(wts, colors)}
    wtcolors = {'gasmass': colors[0], 'Ne8': colors[1]}
    wtlinestyles = {wt: 'solid' for wt in wts}
    ibcolor = colors[3]
    iblinestyle = 'dashed'
    encllevels = [0.98, 0.9, 0.5]
    iblevels = [0.001, 0.01, 0.1]
    linewidths = [1.5, 1., 0.7]
    refZ = 0. # log solar mass fraction units

    iontab = Linetable_PS20(ion, redshift, emission=False, 
                            vol=True, lintable=True)    
    iontab.findiontable()
    logTK = iontab.logTK
    lognH = iontab.lognHcm3
    logZsol = iontab.logZsol
    tab_T_Z_nH = iontab.iontable_T_Z_nH 
    iZ = np.where(np.isclose(refZ, logZsol))[0][0]
    tab_T_nH = tab_T_Z_nH[:, iZ, :] 

    ionpds = {(simname, rrange): readpd(filetemp.format(simname=simname,
                                                        snapnum=snap,
                                                        wtstr=ion), 
                              rrange_rvir=rrange)
              for simname, snap in zip(simnames, snaps)
              for rrange in [rrange_rvir_main] + list(rranges_rvir_sub)}
    vmax = max([np.max(ionpds[(simname, rrange)]['linpddens']) 
                for simname in simnames
                for rrange in [rrange_rvir_main] + list(rranges_rvir_sub)])
    vmax = np.log10(vmax)
    vmin = vmax - 6.

    npanels = len(simnames)
    ncols = min(npanels, 4)
    nrows = (npanels - 1) // ncols + 1
    panelwidth = 2.5
    panelheight = 3.5
    cbarwidth = 0.2
    laxheight = 1.
    width_ratios = [panelwidth] * ncols + [cbarwidth]
    height_ratios = [panelheight] * nrows + [laxheight]
    figsize = (sum(width_ratios), sum(height_ratios))

    fig = plt.figure(figsize=figsize)
    grid = gsp.GridSpec(ncols=ncols + 1, nrows=nrows + 1,
                        width_ratios=width_ratios, 
                        height_ratios=height_ratios,
                        hspace=0.1, wspace=0.1)
    cax = fig.add_subplot(grid[:-1, -1])
    lax = fig.add_subplot(grid[-1, :])
    lax.axis('off')
    nsubcols = len(rranges_rvir_sub)
    
    clabel = ('$\\log_{10} \\, \\left( \\mathrm{f}(\\mathrm{ions})'
              ' \\,/\\, \\Delta \\, \\log_{10} \\mathrm{n}_{\\mathrm{H}}'
              ' \\,/\\, \\Delta \\, \\log_{10} \\mathrm{T} \\right)$')
    xlabel = ('$\\log_{10} \\, \\mathrm{n}_{\\mathrm{H}}'
              '\\; [\\mathrm{cm}^{-3}]$')
    ylabel = '$\\log_{10} \\, \\mathrm{T} \\; [\\mathrm{K}]$'
    _cmap = mpl.cm.get_cmap('gist_yarg')
    cmap = pu.truncate_colormap(_cmap, minval=0., maxval=0.6)
    fontsize = 12

    mainaxes = []
    smallaxes = []
    for supi, (simname, snap) in enumerate(zip(simnames, snaps)):
        xi = supi % ncols
        yi = supi // ncols
        doleft = xi == 0
        dobottom = supi >= npanels - ncols
        
        subgrid = gsp.GridSpecFromSubplotSpec(ncols=nsubcols, nrows=2,
            hspace=0.2, wspace=0., subplot_spec=grid[yi, xi],
            height_ratios=(panelheight - panelwidth, panelwidth))
        _mainax = fig.add_subplot(subgrid[1, :])
        _smallaxes = [fig.add_subplot(subgrid[0, i]) for i in range(nsubcols)]

        mainaxes.append(_mainax)
        smallaxes.append(_smallaxes)
        _mainax.tick_params(which='both', direction='in', 
                            labelsize=fontsize - 1,
                            top=True, right=True, labelbottom=dobottom,
                            labelleft=doleft)
        [_sax.tick_params(which='both', direction='in', 
                          labelsize=fontsize - 2,
                          top=True, right=True, labelbottom=si % 2 == 0,
                          labelleft=(si == 0 and doleft))
         for si, _sax in enumerate(_smallaxes)]
        if dobottom:
            _mainax.set_xlabel(xlabel, fontsize=fontsize)
        if doleft:
            _mainax.set_ylabel(ylabel, fontsize=fontsize)
        
        iondat_main = ionpds[(simname, rrange_rvir_main)]
        toplot = np.log10(iondat_main['linpddens'])
        extent = (iondat_main['nHbins'][0], iondat_main['nHbins'][-1],
                  iondat_main['Tbins'][0], iondat_main['Tbins'][-1])
        img = _mainax.imshow(toplot, extent=extent,
                             cmap=cmap, vmin=vmin, vmax=vmax, 
                             interpolation='nearest', origin='lower',
                             aspect='auto')
        for wt in wts[1:]:
            #for _sax, _rrange_rvir in zip(_smallaxes, rranges_rvir_sub):
            #    dat = readpd(filetemp.format(simname=simname, snapnum=snap,
            #                                 wtstr=wt), 
            #                 rrange_rvir=_rrange_rvir)
            #    pu.add_2dhist_contours(_sax, dat['linpd'].T, 
            #                           [dat['nHbins'], dat['Tbins']],
            #                           (0, 1),
            #                           histlegend=False, fraclevels=True, 
            #                           levels=encllevels, colors=wtcolors[wt],
            #                           linestyles=wtlinestyles[wt], 
            #                           linewidths=linewidths)
            dat = readpd(filetemp.format(simname=simname, snapnum=snap,
                                         wtstr=wt), 
                         rrange_rvir=rrange_rvir_main)
            pu.add_2dhist_contours(_mainax, dat['linpd'].T, 
                                   [dat['nHbins'], dat['Tbins']],
                                   (0, 1),
                                   histlegend=False, fraclevels=True, 
                                   levels=encllevels, colors=wtcolors[wt],
                                   linestyles=wtlinestyles[wt], 
                                   linewidths=linewidths)
        pu.add_2dhist_contours(_mainax, iondat_main['linpd'].T, 
                               [iondat_main['nHbins'], 
                                iondat_main['Tbins']],
                                (0, 1),
                                histlegend=False, fraclevels=True, 
                                levels=encllevels, colors=wtcolors[ion],
                                linestyles=wtlinestyles[ion], 
                                linewidths=linewidths)
        for _sax, _rrange_rvir in zip(_smallaxes, rranges_rvir_sub):
            iondat = ionpds[(simname, _rrange_rvir)]
            toplot = np.log10(iondat['linpddens'])
            extent = (iondat['nHbins'][0], iondat['nHbins'][-1],
                      iondat['Tbins'][0], iondat['Tbins'][-1])
            img = _sax.imshow(toplot, extent=extent,
                              cmap=cmap, vmin=vmin, vmax=vmax, 
                              interpolation='nearest', origin='lower',
                              aspect='auto')
            pu.add_2dhist_contours(_sax, iondat['linpd'].T, 
                              [iondat['nHbins'], iondat['Tbins']],
                              (0, 1),
                              histlegend=False, fraclevels=True, 
                              levels=[encllevels[0]], colors=wtcolors[ion],
                              linestyles=wtlinestyles[ion], 
                              linewidths=[linewidths[0]])
            rl = (f'${_rrange_rvir[0]:.2f} \\, \\endash$ \n'
                  f'${_rrange_rvir[1]:.2f}'
                  '\\, \\mathrm{R}_{\\mathrm{vir}}$')
            _sax.text(0.5, 1.05, rl, color='black',
                      fontsize=fontsize - 2, horizontalalignment='center',
                      verticalalignment='bottom', transform=_sax.transAxes)
        if supi == -1:
            posy = 0.05
            va = 'bottom'
        else:
            posy = 0.95
            va = 'top'
        supl = sl.ic_from_simname(simname) + ' ' + \
               sl.plotlabel_from_physlabel_short[
                   sl.physlabel_from_simname(simname)]
        _mainax.text(0.05, posy, supl, color='black',
                     fontsize=fontsize - 2, horizontalalignment='left',
                     verticalalignment=va, transform=_mainax.transAxes)
        rl = (f'${rrange_rvir_main[0]:.1f} \\endash {rrange_rvir_main[1]:.1f}'
              '\\mathrm{R}_{\\mathrm{vir}}$')
        _mainax.text(0.05, 0.03, rl, color='black',
                     fontsize=fontsize - 2, horizontalalignment='left',
                     verticalalignment='bottom', 
                     transform=_mainax.transAxes)
        
    plt.colorbar(img, cax=cax, orientation='vertical', aspect=15.,
                 extend='min')
    cax.set_ylabel(clabel, fontsize=fontsize)

    axes = mainaxes + [ax for sub in smallaxes for ax in sub]
    xlims = [ax.get_xlim() for ax in axes]
    xlim = [min([l[0] for l in xlims]), max([l[1] for l in xlims])]
    xlim[0] = max([xlim[0], -5.4])
    xlim[1] = min([xlim[1], -1.7])
    ylims = [ax.get_ylim() for ax in axes]
    ylim = [min([l[0] for l in ylims]), max([l[1] for l in ylims])]
    ylim[0] = max([ylim[0], 4.05])
    ylim[1] = min([ylim[1], 7.4])
    for axi, ax in enumerate(axes):
        #xlim = ax.get_xlim()
        #ylim = ax.get_ylim()
        ax.contour(lognH, logTK, tab_T_nH, origin='lower',
                   levels=iblevels, linewidths=linewidths, colors=ibcolor,
                   linestyles=iblinestyle)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    
    legtitle = (f'encl. frac.')
    handles0 = [mlines.Line2D((), (), linewidth=2., 
                              label=wtlabels[wt],
                              color=wtcolors[wt], 
                              linestyle='solid',)
                for wt in wts]
    handles1 =  [mlines.Line2D((), (), linewidth=lw, 
                              label=f'{100. * cl:.0f}%',
                              color='gray', 
                              linestyle='solid')
                for cl, lw in zip(encllevels, linewidths)]     
    leg0 = lax.legend(handles=handles0 + handles1, loc='upper left',
                      bbox_to_anchor=(0.0, 0.5),
                      fontsize=fontsize - 1, 
                      title=legtitle,
                      title_fontsize=fontsize - 1,
                      ncol=int(np.ceil(1.1 * npanels)),
                      columnspacing=0.8, handlelength=1.2,
                      handletextpad=0.6)

    legtitle = (f'ion frac.')
    handles = [mlines.Line2D((), (), linewidth=lw, 
                             label=f'{ibl:.0e}'.replace('e-0', 'e-'),
                             color=ibcolor, 
                             linestyle=iblinestyle)
               for ibl, lw in zip(iblevels, linewidths)]
    lax.legend(handles=handles, loc='upper right',
               bbox_to_anchor=(1.0, 0.5),
               fontsize=fontsize - 1, 
               title=legtitle,
               title_fontsize=fontsize - 1,
               ncol=int(np.floor(1.1 * npanels)),
               columnspacing=0.8, handletextpad=0.6,
               handlelength=1.6)
               #handlelength=1.5, handletextpad=0.4
    lax.add_artist(leg0)

    plt.savefig(outname, bbox_inches='tight')
