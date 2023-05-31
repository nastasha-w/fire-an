import h5py
from matplotlib import patheffects
import matplotlib.collections as mcol
import matplotlib.colors as mcolors
import matplotlib.gridspec as gsp
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import numpy as np

import fire_an.makeplots.plot_utils as pu
import fire_an.simlists as sl


def getmap(filen, weightmap=False):
    with h5py.File(filen, 'r') as f:
        _map = f['map'][:]
        # TODO: divide by depth in cm for weightmap
        _minf = f['map'].attrs['minfinite']
        _maxf = f['map'].attrs['max']
    # depth in cm (for col. dens. -> dens.)
    out = {'map': _map, 'minv': _minf, 'maxv': _maxf,
           'extent_pkpc': None, 'rvir_pkpc': None,
           'depth_cm': None, 'pax': None}
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
        print('selqty_args: ', 
              list(f[selpath]['selqty_args'].attrs.items()))
    out = {'pos_pkpc': pos_pkpc, 'vel_kmps': vel_kmps, 
           'cosmopars': cosmopars, 'halodat': halodat}
    return out

def plotgrid_velweight(filen_temp_map, filen_temp_pv, 
                       simnames, colfills, weight='Ne8',
                       vmarks=(-100., -50., 50., 100.),
                       vmax=100.,
                       clabelw=None, collabels=None, outname=None):
    sims_hr = sl.m13_hr_all2 + sl.m12_hr_all2
    sims_sr = sl.m13_sr_all2 + sl.m12_sr_all2

    snapshots = [sl.snaps_hr if simname in sims_hr
                 else sl.snaps_sr if simname in sims_sr
                 else None
                 for simname in simnames]
    vlosmaps = [[getmap(filen_temp_map.format(snapnum=snap,
                                              **colfill),
                        weightmap=False)
                 for snap in _snapshots]
                for colfill, _snapshots in zip(colfills, snapshots)]
    totmaps = [[getmap(filen_temp_map.format(snapnum=snap, **colfill),
                        weightmap=True)
                  for snap in _snapshots]
                for colfill, _snapshots in zip(colfills, snapshots)]
    pvs = [[getpv(filen_temp_pv.format(snapnum=snap, **colfill),
                  weight=weight)
            for snap in _snapshots]
           for colfill, _snapshots in zip(colfills, snapshots)]

    ncols = 2 * len(colfills)
    nrows = len(snapshots)
    panelsize = 1.7
    caxspace = 1.
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
    waxes = [[fig.add_subplot(grid[i, 2 * j]) for j in range(ncols)]
             for i in range(nrows)]
    vaxes = [[fig.add_subplot(grid[i, 2 * j + 1]) 
              for j in range(ncols)]
             for i in range(nrows)]
    caxw = fig.add_subplot(grid[nrows, :ncols // 3])
    caxv = fig.add_subplot(grid[nrows, ncols // 3 : 2 * (ncols // 3)])
    qvax = fig.add_subplot(grid[nrows, 2 * (ncols // 3):])
    
    rvmin = np.min([np.min([elt['rvir_pkpc'] for elt in row]) 
                    for row in totmaps])
    wmin = np.min([np.min([elt['minv'] for elt in col]) for col in totmaps])
    wmax = np.max([np.max([elt['maxv'] for elt in col]) for col in totmaps])
    wmin = max(wmin, wmax - 5.)
    vmin = -vmax
    _cmapv = 'RdBu'
    _cmapw = 'viridis'
    cmapv = mcolors.Colormap(_cmapv)
    cmapw = mcolors.Colormap(_cmapw)
    normv = mcolors.Normalize(vmin=vmin, vmax=vmax)
    normw = mcolors.Normalize(vmin=wmin, vmax=wmax)
    cmapv.set_under(cmapv(0.))
    cmapv.set_over(cmapv(1.))
    cmapw.set_under(cmapw(0.))
    cmapw.set_over('black')
    clabelv = '$v_{\\mathrm{los}} \\; [\\mathrm{km}\\,\\mathrm{s}^{-1}]$'
    if clabelw is None:
        clabelw = '$\\log_{10}\\,$' + weight + ' density'
    quiverscale = 10. # uniform across plot instead of per panel
    
    rvmar = 1.5
    xlim = (-rvmar * rvmin, rvmar * rvmin)
    ylim = (-rvmar * rvmin, rvmar * rvmin)
    for ri in range(nrows):
        for ai in range(ncols):
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
            
            wimg = wax.imgshow(wdata['map'].T, extent=wdata['extent'],
                               cmap=cmapw, norm=normw, origin='lower',
                               interpolation='nearest', rasterized=True)
            vimg = vax.imgshow(vdata['map'].T, extent=vdata['extent'],
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
            quiver = wax.quiver(*(_pos + _vel), c=_c, cmap=cmapv,
                                norm=normv, edgecolor='black',
                                angles='xy', pivot='tail', 
                                scale=quiverscale,
                                scale_units='xy')
            if ai == 0 and ri == 0:
               vax.text(1.05 * 2**-0.5 * wdata['rvir_pkpc'], 
                         1.05 * 2**-0.5 * wdata['rvir_pkpc'], 
                         '$R_{\\mathrm{vir}}$',
                         color='red', fontsize=fontsize)
            if ri == 0 and collabels is not None:
                # hacky version of title across w/v colums
                wax.text(1.0, 1.05, collabels[ai], 
                         fontsize=fontsize - 1, color='black',
                         horizontalalignment='center', 
                         verticalalignment='bottom',
                         transform=wax.transAxes)
            if ai == 0:
                # scale marker per row
                scalelen_pkpc = 50.
                x0 = xlim[0]
                xr = xlim[1] - xlim[0]
                y0 = ylim[0]
                yr = ylim[1] - ylim[0]
                axcoords = (0.07, 0.07)
                xplot = [x0 + axcoords[0] * xr,
                         x0 + axcoords[0] * xr + scalelen_pkpc]
                yplot = [y0 + axcoords[1] * yr] * 2
                vax.plot(xplot, yplot, color='red', linestyle='solid',
                         linewidth=2.)
                vax.text(0.5 * (xplot[0] + xplot[1]), yplot[0] + 0.02 * yr,
                         f'{scalelen_pkpc:.0f} pkpc', fontsize=fontsize,
                         horizontalalignment='center', 
                         verticalalignment='bottom')
            
            wax.set_xlim(xlim)
            wax.set_ylim(ylim)
            vax.set_xlim(xlim)
            vax.set_ylim(ylim)
    
    plt.colorbar(wimg, cax=caxw, extend='min', orientation='horizontal')
    caxw.set_xlabel(clabelw, fontsize=fontsize)
    plt.colorbar(vimg, cax=caxv, extend='both', orientation='horizontal')
    caxw.set_xlabel(clabelv, fontsize=fontsize)
    nar = len(vmarks)
    xp = np.linspace(0.25, 0.75, nar)
    yp = [0.8] * nar
    vunits = '$\\, \\mathrm{km}/\\mathrm${s}'
    for i in range(nar):
        qvax.quiverkey(quiver, xp[i], yp[i], vmarks[i], 
                       label=f'{vmarks[i]:.0f}' + vunits,
                       angle=0, coordinates='axes', labelpos='S',
                       fontproperties={'fontsize': fontsize})
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

            





    
    
    

