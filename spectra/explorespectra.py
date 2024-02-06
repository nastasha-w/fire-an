import glob

import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# used for read-in, not actual fitting
import findcomponents as fc


logN_defaults = (-np.inf, 13.0, 13.5, 14.0, 14.5, np.inf)

def plotoverview_spectra(filepattern: str,
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
    fig = plt.figure(size=(panelwidth, height))
    axes = [fig.add_subplot(grid[ri, 0]) for ri in range(nrows)]
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)

    filens = glob.glob(filepattern)
    for filen in filens:
        spec = fc.SpectrumFitFreq(fc.ne8_770, filen=filen)
        logN = np.log10(spec.line.tau_to_coldens(spec.tau_raw, spec.vel_kmps))
        ri = np.searchsorted(logNbins, logN) - 1
        if ri < 0 or ri >= nrows:
            print(f'Skipping spectrum {filen}; out of logNbins range')
            continue
        ax = axes[ri]
        ax.plot(spec.vel_kmps, spec.spec_raw, linestyle='solid',
                color='gray', alpha='0.5', linewidth=1.)
    
    xlims = [ax.get_xlim() for ax in axes]
    xmin = min([xl[0] for xl in xlims])
    xmax = max([xl[1] for xl in xlims])
    for ri, ax in enumerate(axes):
        ax.set_ylim(0.0, 1.05)
        ax.set_xlim(xmin, xmax)
        ax.tick_params(which='both', labelleft=True, 
                       labelbottom=(ri == nrows - 1),
                       direction='in', labelsize=fontsize - 1)
        ax.hline(1., color='black', linestyle='dotted', linewidth=1.,
                 alpha=0.7)
        label = ('$\\log \\mathrm{N} [\\mathrm{cm}^{-2}] = '
                f'{logNbins[ri]:.1f} \\endash {logNbins[ri + 1]:.1f}$')
        ax.text(0.98, 0.02, label, fontsize=fontsize,
                transform=ax.transAxes, horizontalalignment='right',
                verticalalignment='right')
        if ri == nrows - 1:
            ax.set_xlabel('$v \\; [\\mathrm{km}\\,\\mathrm{s}^{-1}]$',
                          fontsize=fontsize)
        if ri == nrows // 2:
            ax.set_ylabel('transmission', fontsize=fontsize)

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

        







