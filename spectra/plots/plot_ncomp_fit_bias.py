
import h5py
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np

import fire_an.spectra.findcomponents as fc
import fire_an.spectra.test_ncomp_bias as tnb
import fire_an.utils.constants_and_units as c

def fmtlinedesc(cd_logcm2, b_kmps, v_kmps):
    cdstr = ', '.join([f"{val:.1f}" for val in cd_logcm2])
    bstr = ', '.join([f"{val:.0f}" for val in b_kmps])
    vstr = ', '.join([f"{val:.0f}" for val in v_kmps])
    out = ("$\\mathrm{N}: " + cdstr + 
           ","
           "b: " + bstr + ", " +
           "v: " + vstr + "$")
    return out


def plotpanel_spec(ax, h5group, fontsize=12):
    ax.axhline(1., color='gray', linestyle='dotted')
    line = fc.getline(h5group['line'])
    lsf_sigma_kmps = h5group.attrs['lsf_sigma_kmps']

    nu_Hz = h5group['nu_Hz'][:]
    noisyspec = h5group['noisyspec'][:]
    wl_cm = c.c / nu_Hz
    vbins_kmps = (wl_cm / (line.wavelength_A * 1e-8) - 1.) * c.c * 1e-5
    
    b_in = h5group.attrs['b_kmps_in']
    cd_in = h5group.attrs['logN_cm2_in']
    v_in = h5group.attrs['vcen_kmps_in']
    spec_in = line.getspectrum(nu_Hz, cd_in, b_in, v_in)
    spec_nonoise = tnb.convgauss(vbins_kmps, spec_in, lsf_sigma_kmps)
    
    ax.plot(vbins_kmps, noisyspec, color='black', linestyle='solid')
    ax.plot(vbins_kmps, spec_nonoise, color='gray', linestyle='solid')
    note_in = 'In: ' + fmtlinedesc(cd_in, b_in, v_in)
    ax.text(0.02, 0.98, note_in,
            color='gray', fontsize=fontsize,
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes)

    grp_c1 = h5group['fit_1comp']
    fit_failed = bool(grp_c1.attrs['fit_failed'])
    if not fit_failed:
        b_fit = grp_c1.attrs['b_kmps']
        cd_fit = grp_c1.attrs['logN_cm2']
        v_fit = grp_c1.attrs['vcen_kmps']
        _spec_fit = line.getspectrum(nu_Hz, cd_fit, b_fit, v_fit)
        spec_fit = tnb.convgauss(vbins_kmps, _spec_fit, lsf_sigma_kmps)
        ax.plot(vbins_kmps, spec_fit, color='blue', linestyle='solid')
        note_fit = '1-comp. fit: ' + fmtlinedesc(cd_fit, b_fit, v_fit)
    else:
        note_fit = '1-comp. fit: Failed'
    ax.text(0.02, 0.20, note_fit,
            color='blue', fontsize=fontsize,
            horizontalalignment='left', verticalalignment='bottom',
            transform=ax.transAxes)
    
    grp_c2 = h5group['fit_2comp']
    fit_failed = bool(grp_c2.attrs['fit_failed'])
    if not fit_failed:
        b_fit = grp_c2.attrs['b_kmps']
        cd_fit = grp_c2.attrs['logN_cm2']
        v_fit = grp_c2.attrs['vcen_kmps']
        _spec_fit = line.getspectrum(nu_Hz, cd_fit, b_fit, v_fit)
        spec_fit = tnb.convgauss(vbins_kmps, _spec_fit, lsf_sigma_kmps)
        ax.plot(vbins_kmps, spec_fit, color='red', linestyle='solid')
        note_fit = '2-comp. fit: ' + fmtlinedesc(cd_fit, b_fit, v_fit)
    else:
        note_fit = '2-comp. fit: Failed'
    ax.text(0.02, 0.02, note_fit,
            color='red', fontsize=fontsize,
            horizontalalignment='left', verticalalignment='bottom',
            transform=ax.transAxes)


def plotfits(filen: str, subset: str = '1compin'):
    """
    
    Parameters:
    -----------
    filen:
        name of the hdf5 file with the data
    subset:
        "1compin" or "2compin_set<number>"
    """
    fontsize = 12
    panelwidth = 4.5
    panelheight = 1.5
    wspace = 0.2
    hspace = 0.2
    ncols_max = 3
    outname = filen[:-5] + f'_specplots_{subset}.pdf'
    with h5py.File(filen, 'r')  as f:
        _grpns = list(f.keys())
        if subset == '1compin':
            grpns = [grpn for grpn in _grpns if grpn.startswith('comp1')]
        elif subset.startswith('2compin'): # 2compin_set<number>
            nperset = 64
            setn = int((subset.split("_")[-1])[3:])
            _grpns = [grpn for grpn in _grpns if grpn.startswith('comp2')]
            _grpns.sort()
            sel = slice(nperset * setn, nperset * (setn + 1))
            grpns = _grpns[sel]
        npanels = len(grpns)
        ncols = min(npanels, ncols_max)
        nrows = (npanels - 1) // ncols + 1
        figsize = (ncols * panelwidth, nrows * panelheight)

        fig = plt.figure(figsize=figsize)
        title = ("$\\mathrm{N}: \\log_{10}(\\mathrm{cm}^{-2}), "
                 " b, v: \\mathrm{km}\\,\\mathrm{s}^{-1}$")
        fig.suptitle(title, fontsize=fontsize)
        grid = gsp.GridSpec(ncols=ncols, nrows=nrows, wspace=wspace,
                            hspace=hspace)
        axes = [fig.add_subplot(grid[i // ncols, i % ncols])
                for i in range(npanels)]
        for i, grpn in enumerate(grpns):
            ax = axes[i]
            dobottom = i >= npanels - ncols
            doleft = i % ncols == 0
            if dobottom:
                ax.set_xlabel('$v \\; [\\mathrm{km}\\,\\mathrm{s}^{-1}]$',
                              fontsize=fontsize)
            if doleft:
                ax.set_ylabel('transmission', fontsize=fontsize)
            ax.tick_params(which='both', labelsize=fontsize - 1, 
                           direction='in', right=True, top=True)
            plotpanel_spec(ax, f[grpn], fontsize=fontsize - 3)
    fig.savefig(outname, bbox_inches='tight')


def get_N_sigma(h5group, ncomp_fit):
    out = {}
    line = fc.getline(h5group['line'])
    nu_Hz = h5group['nu_Hz'][:]
    wl_cm = c.c / nu_Hz
    vbins_kmps = (wl_cm / (line.wavelength_A * 1e-8) - 1.) * c.c * 1e-5
    
    b_in = h5group.attrs['b_kmps_in']
    cd_in = h5group.attrs['logN_cm2_in']
    v_in = h5group.attrs['vcen_kmps_in']
    spec_in = line.getspectrum(nu_Hz, cd_in, b_in, v_in)
    out['logN_in'] = np.log10(np.sum(10**cd_in))
    weight = -1. * np.log(spec_in) # tau from transmission
    vcen = np.sum(weight * vbins_kmps) / np.sum(weight)
    dv = vbins_kmps - vcen
    sigma = np.sqrt(np.sum(weight * dv**2) / np.sum(weight))
    out['sigma_in'] = sigma
    
    grp = h5group[f'fit_{ncomp_fit}comp']
    fit_failed = bool(grp.attrs['fit_failed'])
    if fit_failed:
        out['logN_fit'] = None
        out['sigma_fit'] = None
    else:
        b_fit = grp.attrs['b_kmps']
        cd_fit = grp.attrs['logN_cm2']
        v_fit = grp.attrs['vcen_kmps']
        spec_fit = line.getspectrum(nu_Hz, cd_fit, b_fit, v_fit)
        out['logN_fit'] = np.log10(np.sum(10**cd_fit))
        weight = -1. * np.log(spec_fit) # tau from transmission
        vcen = np.sum(weight * vbins_kmps) / np.sum(weight)
        dv = vbins_kmps - vcen
        sigma = np.sqrt(np.sum(weight * dv**2) / np.sum(weight))
        out['sigma_fit'] = sigma
    return out


def plot_N_sigma_bias(filen, ncomp_in=2, ncomp_fit=1):
    outname = filen[:-5] \
              + f'_N_sigma_{ncomp_in}compin_{ncomp_fit}compfit.pdf'

    fontsize = 12
    fig = plt.figure(figsize=(5.5, 5.))
    ax = fig.add_subplot()
    ax.tick_params(which='both', direction='in', labelsize=fontsize - 1,
                   top=True, right=True)
    ax.set_xlabel("$\\log_{10}\\, \\mathrm{N} \\; [\\mathrm{cm}^{-2}]$",
                  fontsize=fontsize)
    ax.set_ylabel("$\\sigma_{\\tau}(v) \\; "
                  "[\\mathrm{km}\\,\\mathrm{s}^{-1}]$",
                  fontsize=fontsize)
    fig.suptitle(f"{ncomp_in} input components,"
                 f" fitted with {ncomp_fit} components",
                 fontsize=fontsize)
    with h5py.File(filen, 'r')  as f:
        _grpns = list(f.keys())
        gstart = f"comp{ncomp_in}"
        grpns = [grpn for grpn in _grpns if grpn.startswith(gstart)]
        first = True
        for grpn in grpns:
            dp = get_N_sigma(f[grpn], ncomp_fit)
            p_in = (dp['logN_in'], dp['sigma_in'])
            p_fit = (dp['logN_fit'], dp['sigma_fit'])
            if dp['logN_fit'] is None:
                ax.scatter([p_in[0]], [p_in[1]], s=5,
                           edgecolors='black', facecolors='None',
                           marker="o", alpha=0.5)
            else:
                if first:
                    label1 = 'input'
                    label2 = 'fit'
                else:
                    label1 = None
                    label2 = None
                ax.plot([p_in[0], p_fit[0]], [p_in[1], p_fit[1]],
                        color='gray', linestyle='solid', alpha=0.5)
                ax.scatter([p_in[0]], [p_in[1]], s=5,
                           edgecolors='black', facecolors='None',
                           marker="o",
                           label=label1)
                ax.scatter([p_fit[0]], [p_fit[1]], s=5,
                           edgecolors='magenta', facecolors='None',
                           marker="o",
                           label=label2)
                first = False
        ax.legend(fontsize=fontsize)

    fig.savefig(outname, bbox_inches='tight')
    
