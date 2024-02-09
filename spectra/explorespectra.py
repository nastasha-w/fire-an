import glob

import h5py
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# used for read-in, not actual fitting
import fire_an.makeplots.plot_utils as pu
import fire_an.spectra.findcomponents as fc
import fire_an.spectra.genspectra as gs
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu



logN_defaults = (-np.inf, 12.5, 13.0, 13.5, 14.0, np.inf)

def plotoverview_spectra(filepattern: str,
                         infofile: str | None = None,
                         logNbins: np.ndarray | None 
                                   = np.array(logN_defaults),
                         title: str | None = None,
                         outname: str | None = None):
    nrows = len(logN_defaults) - 1
    panelheight = 1.5
    panelwidth = 7.5
    height = panelheight * nrows
    fontsize = 12

    grid = gsp.GridSpec(ncols=2, nrows=nrows, hspace=0.,
                        wspace=0.,
                        width_ratios=[6.5, 0.5])
    fig = plt.figure(figsize=(panelwidth, height))
    axes = [fig.add_subplot(grid[ri, 0]) for ri in range(nrows)]
    cax = fig.add_subplot(grid[:, 1])
    cmapn = 'rainbow'

    if title is not None:
        fig.suptitle(title, fontsize=fontsize)

    if infofile is not None:
        cgpath = 'Header/cengal/halodata_doc_dict'
        with h5py.File(infofile, 'r') as f:
            cosmopars = {key: val for key, val in 
                         f[cgpath]['cosmopars_dict'].attrs.items()}
            pgal_cm = f['Header/cengal'].attrs['pcen_cm']
            vcom_cmps = f['Header/cengal'].attrs['vcom_cmps']
            hpar = cu.Hubble(cosmopars)
            axis = (f['Header/sample'].attrs['axis']).decode()
            xaxi, yaxi, losaxi = gs.getinds_ax(axis)
            starts = f['startpos_cm'][:]
            ends = f['endpos_cm'][:]
            if not (np.allclose(starts[:, losaxi], starts[0, losaxi]) 
                    and np.allclose(ends[:, losaxi], ends[0, losaxi])):
                print('Sightlines start and end at different positions,'
                      ' so they will not be on a common velocity grid.')
                vgal_kmps = 0.
            p0 = starts[0, losaxi]
            poff = pgal_cm[losaxi] - p0
            zgal = (-1. * vcom_cmps[losaxi] / c.c + 1.) \
                   * (1. - poff * hpar / c.c) - 1.
            # for vbins, snapshot redshift is factored out
                   #* (1. + cosmopars['z']) - 1
            vgal_kmps = zgal * c.c * 1e-5
            print(zgal)
            print(vgal_kmps)

            ipars = (starts[:, xaxi] - pgal_cm[xaxi])**2 \
                    + (starts[:, yaxi] - pgal_cm[yaxi])**2
            ipars = np.sqrt(ipars) / (c.cm_per_mpc * 1e-3)
            
    else:
        vgal_kmps = 0. # just use the output velocities
    maxipar = np.max(ipars)
    #print(maxipar)
    cmap = pu.paste_cmaps([cmapn], edges=[0., maxipar])
    clabel = 'impact parameter [kpc]'

    filens = glob.glob(filepattern)
    for filen in filens:
        #print(filen)
        sli = filen.split('.')[-2]
        sli = sli.split('_')[-1]
        sli = int(sli)
        ipar = ipars[sli]
        spec = fc.SpectrumFitFreq(fc.ne8_770, filen=filen)
        #print(spec.tau_raw)
        #print(spec.vel_kmps)
        logN = np.log10(spec.line.tau_to_coldens(spec.tau_raw, spec.vel_kmps))
        #print(logN)
        ri = np.searchsorted(logNbins, logN) - 1
        if ri < 0 or ri >= nrows:
            print(f'Skipping spectrum {filen}; out of logNbins range')
            continue
        ax = axes[ri]
        color = cmap(ipar / maxipar)
        ax.plot(spec.vel_kmps - vgal_kmps, spec.spec_raw, linestyle='solid',
                color=color, alpha=0.5, linewidth=1.)
    
    xlims = [ax.get_xlim() for ax in axes]
    xmin = min([xl[0] for xl in xlims])
    xmax = max([xl[1] for xl in xlims])
    pu.add_colorbar(cax, vmin=0., vmax=maxipar, cmap=cmap,
                    clabel=clabel, fontsize=fontsize,
                    extend='neither', orientation='vertical')

    for ri, ax in enumerate(axes):
        #ax.set_ylim(0.0, 1.05)
        ax.set_xlim(xmin, xmax)
        ax.tick_params(which='both', labelleft=True, 
                       labelbottom=(ri == nrows - 1),
                       direction='in', labelsize=fontsize - 1)
        ax.axhline(1., color='black', linestyle='dotted', linewidth=1.,
                   alpha=0.7)
        label = ('$\\log_{10} \\mathrm{N} \\,[\\mathrm{cm}^{-2}] = '
                f'{logNbins[ri]:.1f} \\endash {logNbins[ri + 1]:.1f}$')
        ax.text(0.98, 0.02, label, fontsize=fontsize,
                transform=ax.transAxes, horizontalalignment='right',
                verticalalignment='bottom')
        if ri == nrows - 1:
            ax.set_xlabel('$v \\; [\\mathrm{km}\\,\\mathrm{s}^{-1}]$',
                          fontsize=fontsize)
        if ri == nrows // 2:
            ax.set_ylabel('transmission', fontsize=fontsize)

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')


def plotoverview(testset=4):
    if testset == 4:
        simnames = sl.m12_f2md
        snapshots = [sl.snaps_f2md[0], sl.snaps_f2md[1]]
        filedir = '/projects/b1026/nastasha/spectra/test4/'
        outdir =  '/projects/b1026/nastasha/imgs/spectra/test4/'
        # the one that worked the first time
        simnames = ['crheatfix_m12i_r7100', 'crheatfix_m12f_r7100']
        snapshots = [294, 277]
    
    
    for simname in simnames:
        for snapnum in snapshots:
            filebase = f'/tridentray_{simname}_{snapnum}'
            filepattern = filedir + filebase + '_*.txt'
            infofile = filedir + filebase + '_info.hdf5'
            title = f'FIRE-2 core, {simname}, snapshot {snapnum}'
            outname = outdir + f'spectra_overview_{simname}_{snapnum}.pdf'
            plotoverview_spectra(filepattern, infofile=infofile,
                                 title=title, outname=outname)





