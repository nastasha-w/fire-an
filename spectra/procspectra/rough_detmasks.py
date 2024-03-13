import numpy as np

import fire_an.spectra.findcomponents as fc
import fire_an.utils.constants_and_units as c
import fire_an.utils.math_utils as mu

def get_roughdetmask(spectrum: fc.SpectrumFitBayes,
                     lsf_sigma_kmps: float = 30.,
                     coldetlim_logN_cm2: float = 14.,
                     maxbpar_kmps: float = 100.):
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
                               coldetlim_logN_cm2,
                               equiv_bval_kmps,
                               nuline_Hz)
    minlevel = np.min(minspec)

    # find min level crossings: candidate components lie between those
    vcross = mu.find_intercepts(convflux, vvals, minlevel)
    # crosses minlevel between index and index + 1
    crossinds = np.searchsorted(vvals, vcross) - 1
    enclinds = []
    if convflux[0] > minlevel: 
        enclinds.append([0, crossinds[0]])
        crossinds = crossinds[1:]
    numi = len(crossinds) // 2 # round down to an even number
    enclinds += [[crossinds[2 * i], crossinds[2 * i + 1] + 1]
                 for i in range(numi)]
    if convflux[-1] > minlevel:
        enclinds.append([crossinds[-1], len(vvals)])
    
    # check each 'component' for detectability (raw spectrum)
    enclinds_det = []
    cds = []
    mask = np.zeros(len(vvals), dtype=bool)
    for imin, imax in enclinds:
        sel = slice(imin, imax + 1, None)
        cd = line.tau_to_coldens(spectrum.tau_raw[sel], vvals[sel])
        if cd >= 10**coldetlim_logN_cm2:
            enclinds_det.append([imin, imax])
            cds.append(cd)
            mask[sel] = True
    
    return mask, enclinds_det, cds
    

