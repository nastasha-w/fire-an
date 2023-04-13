import h5py
import numpy as np

import matplotlib.gridspec as gsp
import matplotlib.patheffects as mppe 
import matplotlib.pyplot as plt

import fire_an.makeplots.tol_colors as tc
import fire_an.utils.constants_and_units as c



def readin_data(filen):
    with h5py.File(filen) as f:
        hist = f['histpath']
        rbins_rvir = f['rbinspath']
        cosmopars = {key: val for key, val in f['cppath'].attrs.items()}
        halomass = f['halopars'].attrs['halomass_g']
        halomass_msun = halomass / c.solar_mass
    return rbins_rvir, hist, cosmopars, halomass_msun

def plotfracs_haloes(filen_template, fills_sim, title=None, outname=None,
                     rmin_rvir=0.1, rmax_rvir=1.0):
    ions = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'Carbon']
    ionnames = ['C I', 'C II', 'C III', 'C IV', 'C V', 'C VI']
    colors = tc.tol_cmap('bright')[:len(ions) - 1]
    m11list = []
    m12list = []
    m13list = []
    
    for sfill in fills_sim:
        masses = {}
        for ion in ions:
            filen = filen_template.format(ion=ion, **sfill)
            rbins_rvir, hist, cosmopars, halomass_msun = readin_data(filen)
            stag = filen.split('/')[-1]
            stag = stag.split('_')[0]
            imin = np.where(np.isclose(rmin_rvir, rbins_rvir))[0]
            if len(imin) == 0:
                print(f'No bin match found for radius {rmin_rvir} Rvir')
            else:
                imin = imin[0]
            imax = np.where(np.isclose(rmax_rvir, rbins_rvir))[0]
            if len(imax) == 0:
                print(f'No bin match found for radius {rmax_rvir} Rvir')
            else:
                imax = imax[0]
            mass = np.sum(hist[imin : imax])
            masses[ion] = mass
        if stag.startswith('m11'):
            _ls = m11list
        elif stag.startswith('m12'):
            _ls = m12list
        elif stag.startswith('m13'):
            _ls = m13list
        _ls.append({'masses': masses, 'halomass_msun': halomass_msun, 
                    'stag': stag})
    m11list.sort(key=lambda x: x['halomass_msun'])
    m12list.sort(key=lambda x: x['halomass_msun'])
    m13list.sort(key=lambda x: x['halomass_msun'])
    # None: skip leave a gap between m11/m12/m13 haloes
    halolist = m11list + [None] + m12list + [None] + m13list

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3., 11.))
    fontsize = 12
    if title is not None:
        fig.suptitle(title, fontsize=fontsize)
    
    tickposs = []
    ticklabels = []
    fracss = []
    for xpos, halodat in enumerate(halolist):
        if halodat is None:
            continue
        xpos = xpos + 0.5
        tickposs.append(xpos)
        tl = halodat['stag'] + f' ({np.log10(halodat["halomass_msun"]):.1f})'
        ticklabels.append(tl)

        fracs = [halodat['masses'][ion] / halodat['masses'][ions[-1]]
                 for ion in ions[:-1]]
        fracss.append(fracs)
    fracss = np.array(fracss)
    for ii in range(len(ions) - 1):
        ax.bar(tickposs, fracss[:, ii], color=colors[ii], edgecolor='black',
               bottom=np.sum(fracss[:, :ii], axis=1), tick_label=ticklabels)
    nnames = len(ions) - 1
    textoutline = [mppe.Stroke(linewidth=1.5, foreground='black'),
                   mppe.Normal()]
    for ii in range(nnames):
        ax.text(tickposs[-1] + 0.5, (ii + 0.5) / (nnames + 1), ionnames[ii],
                fontsize=fontsize, color=colors[ii],
                horizontalalignment='left', verticalalignment='center',
                patheffects=textoutline)
    ax.tick_params(which='both', axis='x', labelsize=fontsize, direction='out',
                   labelrotation=45.)
    ax.tick_params(which='both', axis='y', labelsize=fontsize - 1, 
                   direction='in', right=True)
    ax.set_ylim((0., 1.))
    ylab = (f'fraction of {ions[-1]}, '
            f'${rmin_rvir:.1f} \\endash ${rmax_rvir:.1f}'
            f' \\, \\mathrm{{R}}_{{\\mathrm{{vir}}}}$')
    ax.set_ylabel(ylab)

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')
    
        
        

    

