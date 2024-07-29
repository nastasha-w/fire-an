import matplotlib.pyplot as plt
import numpy as np

import fire_an.spectra.findcomponents as fc
import fire_an.utils.constants_and_units as c
import fire_an.utils.math_utils as mu

def get_roughdetmask(spectrum: fc.SpectrumFitBayes,
                     lsf_sigma_kmps: float = 30.,
                     coldetlim_logN_cm2: float = 14.,
                     maxbpar_kmps: float = 100.,
                     mincutfrac: float = 1e-2):
    '''
    untested draft
    '''
    vvals = spectrum.vel_kmps
    rawspec = spectrum.spec_raw
    convflux = spectrum._convolve_gauss(rawspec, width_kmps=lsf_sigma_kmps)
    
    # calc. min flux level below which an absorber must dip
    line = spectrum.line
    nuline_Hz = c.c / (line.wavelength_A * 1e-8)
    equiv_bval_kmps = np.sqrt(maxbpar_kmps**2 + 2. * lsf_sigma_kmps**2)
    vrange_cmps = equiv_bval_kmps * 1e5 * 10.
    dv_cmps = 0.02 * equiv_bval_kmps * 1e5
    v_minspec_cmps = np.arange( -1. * vrange_cmps, 
                               vrange_cmps + 0.5 * dv_cmps,
                               dv_cmps)
    nu_minspec_Hz = nuline_Hz / (1. + v_minspec_cmps / c.c)
    minspec = line.getspectrum(nu_minspec_Hz,
                               np.array([coldetlim_logN_cm2]),
                               np.array([equiv_bval_kmps]),
                               np.array([0.]))
    minlevel = np.min(minspec)

    ## find min level crossings:
    ## candidate components lie between those
    vcross = mu.find_intercepts(convflux, vvals, minlevel)
    # crosses minlevel between index and index + 1
    crossinds = np.searchsorted(vvals, vcross) - 1
    enclinds = []
    if np.all(convflux < minlevel):
        enclinds = [[0, len(vvals)]]
    else:
        if convflux[0] < minlevel: 
            enclinds.append([0, crossinds[0]])
            crossinds = crossinds[1:]
        numi = len(crossinds) // 2 # round down to an even number
        enclinds += [[crossinds[2 * i], crossinds[2 * i + 1] + 1]
                    for i in range(numi)]
        if convflux[-1] < minlevel:
            enclinds.append([crossinds[-1], len(vvals)])
    ## extend regions to capture more of the absorbers
    cutlevel = 1. - (1. - minlevel) * mincutfrac
    cutinds = []
    vcross = mu.find_intercepts(convflux, vvals, cutlevel)
    # crosses cutnlevel between index and index + 1
    crossinds = np.searchsorted(vvals, vcross) - 1
    for ai, encli in enumerate(enclinds):
        cuti = [None, None]
        # first component: region starts at zero or first crossing
        if ai == 0: 
            if crossinds[0] < encli[0]:
                cuti[0] = crossinds[0]
            else:
                cuti[0] = 0
        # last component: region ends at array end or last crossing
        if ai == len(enclinds) - 1: 
            if crossinds[-1] > encli[1] - 1:
                cuti[1] = crossinds[-1] + 1
            else:
                cuti[1] = len(vvals)
        # if a start (end) is not set, there is a component before
        # (after) this one
        if cuti[0] is not None:
            lastend = enclinds[ai - 1][1]
            closestcross = np.max(crossinds[crossinds < encli[ai][0]])
            # cutlevel crossing closest to component start
            if closestcross > lastend: 
                cuti[0] = closestcross
            else: # minimum between components
                cuti[0] = np.argmax(convflux[lastend : enclinds[ai][0]])
                cuti[0] = cuti[0] + lastend
        if cuti[1] is not None:
            nextstart = enclinds[ai + 1][0]
            closestcross = np.mix(crossinds[crossinds >= encli[ai][1]])
            if closestcross < nextstart: 
                cuti[1] = closestcross + 1
            else:
                cuti[1] = np.argmax(convflux[enclinds[ai][1] : nextstart])
                cuti[1] = cuti[1] + enclinds[ai][1] 
        cutinds.append(cuti)
    # check each 'component' for detectability (raw spectrum)
    enclinds_det = []
    cds = []
    mask = np.zeros(len(vvals), dtype=bool)
    for imin, imax in cutinds:
        sel = slice(imin, imax, None)
        cd = line.tau_to_coldens(spectrum.tau_raw[sel], vvals[sel])
        if cd >= 10**coldetlim_logN_cm2:
            enclinds_det.append([imin, imax])
            cds.append(cd)
            mask[sel] = True
    return mask, enclinds_det, cds
    
def test_rough_detmasks(spectrum: fc.SpectrumFitBayes,
                        outname: str | None = None,
                        lsf_sigma_kmps: float = 30.,
                        coldetlim_logN_cm2: float = 13.5,
                        maxbpar_kmps: float = 100.):
    vvals = spectrum.vel_kmps
    rawspec = spectrum.spec_raw
    convflux = spectrum._convolve_gauss(rawspec, width_kmps=lsf_sigma_kmps)
    mask, enclinds, cds = get_roughdetmask(
        spectrum, lsf_sigma_kmps=lsf_sigma_kmps,
        coldetlim_logN_cm2=coldetlim_logN_cm2, maxbpar_kmps=maxbpar_kmps)
    
    fig = plt.figure(figsize=(6., 2.5))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(vvals, rawspec, linestyle='solid', color='blue', linewidth=1.5,
            label='raw spectrum')
    ax.plot(vvals, convflux, linestyle='solid', color='black', linewidth=1.5,
            label='spectrum + LSF')
    for ci, (imin, imax) in enumerate(enclinds):
        sl = slice(imin, imax + 1, None)
        label = 'detected' if ci == 0 else None
        ax.plot(vvals[sl], rawspec[sl], linestyle='dashed', color='red',
                linewidth=1.5, label=label)
        vcen = 0.5 * (vvals[imin] + vvals[imax])
        ax.text(vcen, 1.02, '$ \\log_{10} \\, \\mathrm{N} = '
                f'{np.log10(cds[ci]):.2f}$',
                horizontalalignment='center',
                verticalalignment='bottom',
                fontsize=10)    
    ax.axhline(1., linestyle='dotted', linewidth=1.5, color='gray')
    pmask = np.asarray(np.logical_not(mask), dtype=np.float32)
    pmask[mask] = min(np.min(rawspec), np.min(convflux)) - 0.02 
    ax.plot(vvals, pmask, linestyle='dotted', linewidth=1.5, color='red')

    ax.set_xlabel('$v \\;[\\mathrm{km} / \\mathrm{s}]$', fontsize=12)
    ax.set_ylabel('transmission', fontsize=12)
    ax.tick_params(which='both', direction='in', labelsize=11.,
                   top=True, right=True)
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')


def runtest_rough_detmasks():
    # some pretty much random sightlines,
    # but close to the galaxy for hopefully higher column densities
    outdir = '/projects/b1026/nastasha/tests/absspec_rough_compdet/'
    outname_temp = 'test_compdet_ne8_770_{specname}_lsf30_detlim13p0_maxb50.pdf'
    specdir = '/projects/b1026/nastasha/spectra/test4/'
    specfiles = ['tridentray_m12z_r4200_294_210.txt',
                 #'tridentray_m12z_r4200_294_208.txt',
                 #'tridentray_m12z_r4200_294_212.txt',
                 #'tridentray_m12w_r7100_294_210.txt',
                 #'tridentray_m12w_r7100_294_208.txt',
                 #'tridentray_m12w_r7100_294_212.txt',
                 #'tridentray_crheatfix_m12f_r7100_294_210.txt',
                 #'tridentray_crheatfix_m12f_r7100_294_208.txt',
                 #'tridentray_crheatfix_m12f_r7100_294_212.txt',
                 #'tridentray_m12m_r7100_294_210.txt',
                 #'tridentray_m12m_r7100_294_208.txt',
                 #'tridentray_m12m_r7100_294_212.txt',
                 ]
                 

    for specfile in specfiles:
        outname = outdir + outname_temp.format(specname=specfile)
        spectrum = fc.SpectrumFitBayes(fc.ne8_770,
                                       filen=specdir + specfile)
        test_rough_detmasks(spectrum,
                            outname=outname,
                            lsf_sigma_kmps=30.,
                            coldetlim_logN_cm2=13.0,
                            maxbpar_kmps=50.)