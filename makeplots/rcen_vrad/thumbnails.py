import h5py
import matplotlib as mpl
import matplotlib.collections as mcol
import matplotlib.colors as mcolors
import matplotlib.gridspec as gsp
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import numpy as np

import fire_an.makeplots.plot_utils as pu
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c


def getmap(filen, weightmap=False):
    with h5py.File(filen, 'r') as f:
        halodata = {key: val for key, val 
                    in f['Header/inputpars/halodata'].attrs.items()}
        cosmopars = {key: val for key, val 
                     in f['Header/inputpars/cosmopars'].attrs.items()}
        rvir_pkpc = halodata['Rvir_cm'] / (c.cm_per_mpc * 1e-3)
        axis1 = f['Header/inputpars'].attrs['Axis1']
        axis2 = f['Header/inputpars'].attrs['Axis2']
        axis3 = f['Header/inputpars'].attrs['Axis3']
        boxdiam_cm = f['Header/inputpars'].attrs['diameter_used_cm']
        extent = (-0.5 * boxdiam_cm[axis1] / (c.cm_per_mpc * 1e-3),
                  0.5 * boxdiam_cm[axis1] / (c.cm_per_mpc * 1e-3),
                  -0.5 * boxdiam_cm[axis2] / (c.cm_per_mpc * 1e-3),
                  0.5 * boxdiam_cm[axis2] / (c.cm_per_mpc * 1e-3),
                  )
        depth_cm = boxdiam_cm[axis3]
        pax = ['x', 'y', 'z'][axis3]
        if weightmap:
            _map = f['weightmap'][:]
            # TODO: divide by depth in cm for weightmap
            _minf = f['weightmap'].attrs['minfinite']
            _maxf = f['weightmap'].attrs['max']
            if bool(f['weightmap'].attrs['log']):
                fct = np.log10(depth_cm)
                _map -= fct
                _minf -= fct
                _maxf -= fct
            else:
                fct = 1. / depth_cm
                _map *= fct
                _minf *= fct
                _maxf *= fct
        else:
            _map = f['map'][:]
            _minf = f['map'].attrs['min']
            _maxf = f['map'].attrs['max']

    # depth in cm (for col. dens. -> dens.)
    out = {'map': _map, 'minv': _minf, 'maxv': _maxf,
           'extent': extent, 'rvir_pkpc': rvir_pkpc,
           'depth_cm': depth_cm, 'pax': pax,
           'cosmopars': cosmopars}
    return out

def getpv(filen, weight='Volume'):
    # !! only works for the one run I'm writing this for !!
    weightis = ['Mass', 'Volume', 'Neon', 'Oxygen',
                'O6', 'O7', 'O8', 'Ne8']
    wi = np.where([wt == weight for wt in weightis])[0][0]
    selpath = f'selection_{wi}'
    with h5py.File(filen, 'r') as f:
        pos_pkpc = f[selpath]['pos_norot_pkpc'][:]
        vel_kmps = f[selpath]['vel_norot_kmps'][:]
        csmpath = 'Header/halodata/doc/cosmopars_dict'
        cosmopars = {key: val for key, val in f[csmpath].attrs.items()}
        halodat = {key: val for key, val 
                   in f['Header/halodata'].attrs.items()}
        # just to check we got the right weights
        print('selqty: ', f[selpath].attrs['selqty'].decode())
        #print(f[selpath].attrs['selqty_args'].decode())
        if f[selpath].attrs['selqty_args'].decode() == 'dict':
            print('selqty_args: ', 
                  list(f[selpath]['selqty_args_dict'].attrs.items()))
    out = {'pos_pkpc': pos_pkpc, 'vel_kmps': vel_kmps, 
           'cosmopars': cosmopars, 'halodat': halodat}
    return out

def plotgrid_velweight(filen_temp_map, filen_temp_pv, 
                       simnames, colfills, weight='Ne8',
                       vmark=100., vmax=100., wmax=None, wmin=None,
                       clabelw=None, collabels=None, outname=None):
    sims_hr = sl.m13_hr_all2 + sl.m12_hr_all2
    sims_sr = sl.m13_sr_all2 + sl.m12_sr_all2

    snapshots = [sl.snaps_hr if simname in sims_hr
                 else sl.snaps_sr if simname in sims_sr
                 else None
                 for simname in simnames]
    vlosmaps = [[getmap(filen_temp_map.format(snapnum=snap, 
                                              simname=simname,
                                              **colfill),
                        weightmap=False)
                 for snap in _snapshots]
                for colfill, simname, _snapshots in zip(colfills, simnames, 
                                                        snapshots)]
    totmaps = [[getmap(filen_temp_map.format(snapnum=snap, 
                                             simname=simname,
                                             **colfill),
                        weightmap=True)
                  for snap in _snapshots]
                for colfill, simname, _snapshots in zip(colfills, simnames, 
                                                        snapshots)]
    pvs = [[getpv(filen_temp_pv.format(snapnum=snap, simname=simname,
                                       **colfill),
                  weight=weight)
            for snap in _snapshots]
           for colfill, simname, _snapshots in zip(colfills, simnames, 
                                                   snapshots)]
    ncfill = len(colfills)
    ncols = 2 * len(colfills)
    nrows = len(snapshots[0])
    panelsize = 1.7
    caxspace = 0.5
    height_ratios = [panelsize] * nrows + [caxspace]
    width_ratios = [panelsize] * ncols
    hspace = 0.
    wspace = 0.
    height = sum(height_ratios)
    width = sum(width_ratios)

    fontsize = 12

    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(ncols=ncols, nrows=nrows + 1,
                        hspace=hspace, wspace=wspace, 
                        height_ratios=height_ratios,
                        width_ratios=width_ratios)
    #print('nrows: ', nrows)
    #print('ncols: ', ncols)
    #print(grid)

    waxes = [[fig.add_subplot(grid[i, 2 * j]) for j in range(ncfill)]
             for i in range(nrows)]
    vaxes = [[fig.add_subplot(grid[i, 2 * j + 1]) 
              for j in range(ncfill)]
             for i in range(nrows)]
    caxw = fig.add_subplot(grid[nrows, :ncols // 3])
    #caxv = fig.add_subplot(grid[nrows, ncols // 3 : 2 * (ncols // 3)])
    caxv = fig.add_subplot(grid[nrows, 2 * (ncols // 3):])
    
    rvmin = np.min([np.min([elt['rvir_pkpc'] for elt in row]) 
                    for row in totmaps])
    if wmax is None or wmin is None:
        wmin = np.min([np.min([elt['minv'] for elt in col]) 
                       for col in totmaps])
        wmax = np.max([np.max([elt['maxv'] for elt in col])
                       for col in totmaps])
        wmin = max(wmin, wmax - 5.)
    vmin = -vmax
    _cmapv = 'RdBu'
    _cmapw = 'viridis'
    cmapv = mpl.cm.get_cmap(_cmapv)
    cmapw = mpl.cm.get_cmap(_cmapw)
    normv = mcolors.Normalize(vmin=vmin, vmax=vmax)
    normw = mcolors.Normalize(vmin=wmin, vmax=wmax)
    cmapv.set_under(cmapv(0.))
    cmapv.set_over(cmapv(1.))
    cmapw.set_under(cmapw(0.))
    cmapw.set_over(cmapw(1.))
    clabelv = '$v_{\\mathrm{los}} \\; [\\mathrm{km}\\,\\mathrm{s}^{-1}]$'
    if clabelw is None:
        clabelw = '$\\log_{10}\\,$' + weight + ' density'
    quiverscale = 5. # uniform across plot instead of per panel
    
    rvmar = 1.7
    xlim = (-rvmar * rvmin, rvmar * rvmin)
    ylim = (-rvmar * rvmin, rvmar * rvmin)
    for ri in range(nrows):
        for ai in range(ncfill):
            wax = waxes[ri][ai]
            vax = vaxes[ri][ai]
            wax.axis('off')
            vax.axis('off')
            vdata = vlosmaps[ai][ri]
            wdata = totmaps[ai][ri]
            pvdata = pvs[ai][ri]
            pax = wdata['pax']
            if pax == 'z':
                axis1 = 0
                axis2 = 1
                axis3 = 2
            elif pax == 'x':
                axis1 = 1
                axis2 = 2
                axis3 = 0
            elif pax == 'y':
                axis1 = 2
                axis2 = 0
                axis3 = 1
            
            wimg = wax.imshow(wdata['map'].T, extent=wdata['extent'],
                              cmap=cmapw, norm=normw, origin='lower',
                              interpolation='nearest', rasterized=True)
            vimg = vax.imshow(vdata['map'].T * 1e-5,
                              extent=vdata['extent'],
                              cmap=cmapv, norm=normv, origin='lower',
                              interpolation='nearest', rasterized=True)
            patches = [mpatch.Circle((0., 0.), wdata['rvir_pkpc'])]
            wcollection = mcol.PatchCollection(patches)
            #pe = pu.getoutline(1.5)
            wcollection.set(edgecolor=['red'], facecolor='none', 
                            linewidth=1.5, linestyle='dashed')
            vcollection = mcol.PatchCollection(patches)
            #pe = pu.getoutline(1.5)
            vcollection.set(edgecolor=['black'], facecolor='none', 
                            linewidth=1.5, linestyle='dashed')
            wax.add_collection(wcollection)
            vax.add_collection(vcollection)
            
            _pos = (pvdata['pos_pkpc'][:, axis1], 
                    pvdata['pos_pkpc'][:, axis2])
            _vel = (pvdata['vel_kmps'][:, axis1], 
                    pvdata['vel_kmps'][:, axis2])
            _c = pvdata['vel_kmps'][:, axis3]
            quiver = wax.quiver(*((_pos + _vel) + (_c,)), cmap=cmapv,
                                norm=normv, edgecolor='black',
                                angles='xy', pivot='tail', 
                                scale=quiverscale,
                                scale_units='xy')
            if ai == 0 and ri == 0:
                # Rvir label on circle
                vax.text(1.05 * 2**-0.5 * wdata['rvir_pkpc'], 
                         1.05 * 2**-0.5 * wdata['rvir_pkpc'], 
                         '$R_{\\mathrm{vir}}$',
                         color='black', fontsize=fontsize)
                # quiver scale indicators
                vunits = '$\\, \\mathrm{km}/\\mathrm{s}$'
                wax.quiverkey(quiver, 0.7, 0.95, vmark, 
                              color='red',
                              label=f'{vmark:.0f}' + vunits,
                              angle=0, coordinates='axes', labelpos='S',
                              fontproperties={'size': fontsize - 2.},
                              labelcolor='red', labelsep=0.07)
            if ri == 0 and collabels is not None:
                # hacky version of title across w/v colums
                wax.text(1.0, 1.05, collabels[ai], 
                         fontsize=fontsize - 1, color='black',
                         horizontalalignment='center', 
                         verticalalignment='bottom',
                         transform=wax.transAxes)
            if ai == 0:
                wax.text(-0.05, 0.5, f'$z={wdata["cosmopars"]["z"]:.1f}$',
                         color='black', fontsize=fontsize,
                         horizontalalignment='right',
                         verticalalignment='center', 
                         transform=wax.transAxes,
                         rotation=90.)
            if ri == 0 and ai == 1:
                # scale marker
                scalelen_pkpc = 100.
                x0 = xlim[0]
                xr = xlim[1] - xlim[0]
                y0 = ylim[0]
                yr = ylim[1] - ylim[0]
                axcoords = (0.07, 0.07)
                xplot = [x0 + axcoords[0] * xr,
                         x0 + axcoords[0] * xr + scalelen_pkpc]
                yplot = [y0 + axcoords[1] * yr] * 2
                vax.plot(xplot, yplot, color='black', linestyle='solid',
                         linewidth=2.)
                vax.text(xplot[0], yplot[0] + 0.02 * yr,
                         f'{scalelen_pkpc:.0f} pkpc', 
                         fontsize=fontsize - 2.,
                         horizontalalignment='left', 
                         verticalalignment='bottom')
            
            wax.set_xlim(xlim)
            wax.set_ylim(ylim)
            vax.set_xlim(xlim)
            vax.set_ylim(ylim)
    
    plt.colorbar(wimg, cax=caxw, extend='min', 
                 orientation='horizontal', aspect=15.)
    caxw.set_xlabel(clabelw, fontsize=fontsize)
    plt.colorbar(vimg, cax=caxv, extend='both', orientation='horizontal',
                 aspect=15.)
    caxv.set_xlabel(clabelv, fontsize=fontsize)
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotset_clean_samephys_xyz(hset='m12-pax', weight='Ne8'):
    ddir_maps = '/projects/b1026/nastasha/maps/clean2_vlos/'
    ftemp_maps = ('vlos_by_coldens_{weight}_{simname}_snap{snapnum}'
                  '_shrink-sph-cen_BN98_depth_2.0rvir_{pax}-proj_v3.hdf5')
    ddir_pvsel = '/projects/b1026/nastasha/pvsel/clean2/'
    ftemp_pvsel = ('weightedsel_150_{simname}_snap{snapnum}'
                   '_shrink-sph-cen_BN98_depth_0.1rvir_{pax}-slice_v1.hdf5')
    filen_temp_map = ddir_maps + ftemp_maps
    filen_temp_pv = ddir_pvsel + ftemp_pvsel
    outdir = '/projects/b1026/nastasha/imgs/thumbnails/'
    outfile = 'maps_{weight}dens_vlos_{comp}_{simset}.pdf'

    if weight == 'Mass':
        clabelw = ('$\\log_{10} \\, \\rho \\; '
                   '[\\mathrm{g}\\, \\mathrm{cm}^{-3}]$')
        wmin = -30.
        wmax = -26.
    elif weight == 'Oxygen':
        clabelw = ('$\\log_{10} \\, \\mathrm{n}(\\mathrm{O}) \\; '
                   '[\\mathrm{cm}^{-3}]$')
        wmin = -12.
        wmax = -7.
    elif weight == 'Neon':
        clabelw = ('$\\log_{10} \\, \\mathrm{n}(\\mathrm{Ne}) \\; '
                   '[\\mathrm{cm}^{-3}]$')
        wmin = -12.
        wmax = -7.
    elif weight == 'Ne8':
        clabelw = ('$\\log_{10} \\, \\mathrm{n}(\\mathrm{Ne\\, VIII}) \\; '
                   '[\\mathrm{cm}^{-3}]$')
        wmin = -13.
        wmax = -9.5
    elif weight == 'O6':
        clabelw = ('$\\log_{10} \\, \\mathrm{n}(\\mathrm{O\\,VI}) \\; '
                   '[\\mathrm{cm}^{-3}]$')
        wmin = -13.
        wmax = -9.5
    elif weight == 'O7':
        clabelw = ('$\\log_{10} \\, \\mathrm{n}(\\mathrm{O\\,VII}) \\; '
                   '[\\mathrm{cm}^{-3}]$')
        wmin = -11.5
        wmax = -8.
    elif weight == 'O8':
        clabelw = ('$\\log_{10} \\, \\mathrm{n}(\\mathrm{O\\,VIII}) \\; '
                   '[\\mathrm{cm}^{-3}]$')
        wmin = -11.5
        wmax = -8.
 
    if hset == 'm12-pax':
        _simnames = [
            ('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
             '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
            ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
             '_sdp2e-4_gacc31_fa0.5'),
            ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
             '_sdp1e10_gacc31_fa0.5'),
            ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
             '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
            ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
             '_sdp2e-4_gacc31_fa0.5'),
            ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
             '_sdp1e10_gacc31_fa0.5')]
        simlabels = ['m12q_agncr', 'm12q_agnnocr', 'm12q_nobh',
                     'm12f_agncr', 'm12f_agnnocr', 'm12f_nobh']
        simsets = [[simname] * 3 for simname in _simnames]
        vmax = 100.
        vmark = 100.
        colfills = [{'weight': ('gas-mass' if weight == 'Mass' 
                                else 'gas-vol' if weight == 'Volume'
                                else weight), 
                     'pax': pax}
                    for pax in ['x', 'y', 'z']]
        collabels = [pax + '-projection' for pax in ['x', 'y', 'z']]
        comp = 'pax_comp'
    elif hset == 'm13-pax':
        _simnames = [
            ('m13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
             '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
            ('m13h113_m3e4_MHD_fire3_fireBH_Sep182021'
             '_hr_crdiffc690_sdp1e-4_gacc31_fa0.5'),
            ('m13h113_m3e5_MHD_fire3_fireBH_Sep182021'
             '_crdiffc690_sdp1e10_gacc31_fa0.5'),
            ('m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
             '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
            ('m13h206_m3e4_MHD_fire3_fireBH_Sep182021'
             '_hr_crdiffc690_sdp3e-4_gacc31_fa0.5'),
            ('m13h206_m3e5_MHD_fire3_fireBH_Sep182021'
             '_crdiffc690_sdp1e10_gacc31_fa0.5'),
            ]
        simlabels = ['m13h113_agncr', 'm13h113_agnnocr', 'm13h113_nobh',
                     'm13h206_agncr', 'm13h206_agnnocr', 'm13h206_nobh']
        simsets = [[simname] * 3 for simname in _simnames]
        vmax = 150.
        vmark = 150.
        colfills = [{'weight': ('gas-mass' if weight == 'Mass' 
                                else 'gas-vol' if weight == 'Volume'
                                else weight), 
                     'pax': pax}
                    for pax in ['x', 'y', 'z']]
        collabels = [pax + '-projection' for pax in ['x', 'y', 'z']]
        comp = 'pax_comp'
    elif hset == 'm12-phys':
        simsets = [
            [('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
              '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
             ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
              '_sdp2e-4_gacc31_fa0.5'),
             ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
              '_sdp1e10_gacc31_fa0.5')],
            [('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
              '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
             ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
              '_sdp2e-4_gacc31_fa0.5'),
             ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
              '_sdp1e10_gacc31_fa0.5')]
              ]
        simlabels = ['m12q', 'm12f']
        vmax = 100.
        vmark = 100.
        colfills = [{'weight': ('gas-mass' if weight == 'Mass' 
                                else 'gas-vol' if weight == 'Volume'
                                else weight), 
                     'pax': 'z'}] * 3
        collabels = ['AGN-CR', 'AGN-noCR', 'noBH']
        comp = 'phys_comp_zproj'
    elif hset == 'm13-phys':
        simsets = [
            [('m13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
              '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
             ('m13h113_m3e4_MHD_fire3_fireBH_Sep182021'
              '_hr_crdiffc690_sdp1e-4_gacc31_fa0.5'),
             ('m13h113_m3e5_MHD_fire3_fireBH_Sep182021'
              '_crdiffc690_sdp1e10_gacc31_fa0.5')],
            [('m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021'
              '_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
             ('m13h206_m3e4_MHD_fire3_fireBH_Sep182021'
              '_hr_crdiffc690_sdp3e-4_gacc31_fa0.5'),
             ('m13h206_m3e5_MHD_fire3_fireBH_Sep182021'
              '_crdiffc690_sdp1e10_gacc31_fa0.5')],
            ]
        simlabels = ['m13h113', 'm13h206']
        vmax = 150.
        vmark = 150.
        colfills = [{'weight': ('gas-mass' if weight == 'Mass' 
                                else 'gas-vol' if weight == 'Volume'
                                else weight), 
                     'pax': 'z'}] * 3
        collabels = ['AGN-CR', 'AGN-noCR', 'noBH']
        comp = 'phys_comp_zproj'

    for (simnames, simlabel) in zip(simsets, simlabels):
        outname = outdir + outfile.format(weight=weight, comp=comp,
                                          simset=simlabel)
        plotgrid_velweight(filen_temp_map, filen_temp_pv, 
                        simnames, colfills, weight=weight,
                        vmark=vmark,
                        vmax=vmax, wmax=wmin, wmin=wmax,
                        clabelw=clabelw, collabels=collabels, 
                        outname=outname)







    
    
    

