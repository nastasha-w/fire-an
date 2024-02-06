import glob

import h5py
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# used for read-in, not actual fitting
import fire_an.spectra.findcomponents as fc
import fire_an.spectra.genspectra as gs
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

    grid = gsp.GridSpec(ncols=1, nrows=nrows, hspace=0.,
                        wspace=0.)
    fig = plt.figure(figsize=(panelwidth, height))
    axes = [fig.add_subplot(grid[ri, 0]) for ri in range(nrows)]
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)

    if infofile is not None:
        cgpath = 'Header/cengal/halodata_doc_dict'
        with h5py.File(infofile, 'r') as f:
            cosmopars = {key: val for key, val in 
                         f[cgpath]['cosmopars_dict'].attrs.items()}
            pgal_cm = f[cgpath].attrs['pcen_cm']
            vcom_cmps = f[cgpath].attrs['vcom_cmps']
            hpar = cu.Hubble(cosmopars)
            axis = np.string(f['Header/sample'].attrs['axis'])
            _, _, losaxi = gs.getinds_ax(axis)
            starts = f['startpos_cm'][:, losaxi]
            ends = f['endpos_cm'][:, losaxi]
            if not (np.allclose(starts, starts[0]) 
                    and np.allclose(ends, ends[0])):
                print('Sightlines start and end at different positions,'
                      ' so they will not be on a common velocity grid.')
                vgal_kmps = 0.
            p0 = starts[0]
            poff = pgal_cm[losaxi] - p0
            zgal = (vcom_cmps[losaxi] / c.c + 1.) \
                   * (1. - poff * hpar / c.c) \
                   * (1. + cosmopars['z']) - 1
            vgal_kmps = zgal * c.c * 1e-5
            
    else:
        vgal_kmps = 0. # just use the output velocities

    filens = glob.glob(filepattern)
    for filen in filens:
        #print(filen)
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
        ax.plot(spec.vel_kmps - vgal_kmps, spec.spec_raw, linestyle='solid',
                color='gray', alpha=0.5, linewidth=1.)
    
    xlims = [ax.get_xlim() for ax in axes]
    xmin = min([xl[0] for xl in xlims])
    xmax = max([xl[1] for xl in xlims])
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

        







