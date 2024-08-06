import h5py
import numpy as np
import scipy.optimize as so
from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel

import fire_an.spectra.findcomponents as fc
import fire_an.utils.constants_and_units as c


_bnorm_fit = 10.
_vnorm_fit = 10.
_bounds_logN_cm2 = (10., 16.)
_bounds_vcen_kmps = (-500., 500.)
_bounds_b_kmps = (1., 200.)

def convgauss(vbins_kmps: np.ndarray[float],
              spec: np.ndarray[float],
              sigma_kmps: float) -> np.ndarray[float]:
    dv = np.average(np.diff(vbins_kmps))
    width_bins = sigma_kmps / dv
    kernel = Gaussian1DKernel(width_bins)
    spec_smoothed = convolve(spec, kernel, 
                             boundary='fill', fill_value=1.)
    return spec_smoothed

rndgen = np.random.default_rng(seed=None)

# CUBS: sigma_v for LSF ~10 km/s,
#       typical b pars are ~30 km/s
#       snr 10-30 per 'resolution element'
def noisyspectrum(logN_cm2: np.ndarray[float],
                  b_kmps: np.ndarray[float],
                  vcen_kmps: np.ndarray[float],
                  lsf_sigma_kmps: float = 10.,
                  snr: float = 20.,
                  line: fc.Line = fc.ne8_770):
    oversamplef = 5.
    dv = lsf_sigma_kmps / oversamplef
    vmin = np.floor((np.min(vcen_kmps) - 10. * np.max(b_kmps)) / dv) * dv
    vmax = np.ceil((np.max(vcen_kmps) + 10. * np.max(b_kmps)) / dv) * dv
    vvals_kmps = np.arange(vmin, vmax + 0.5 * dv, dv)
    wls_A = line.wavelength_A * (1. + vvals_kmps * 1e5 / c.c)
    nuspec_Hz = c.c / (wls_A * 1e-8)
    
    snr_use = snr / np.sqrt(oversamplef)
    specvals = line.getspectrum(nuspec_Hz, logN_cm2, b_kmps, vcen_kmps)
    outspec = convgauss(vvals_kmps, specvals, lsf_sigma_kmps)
    noise = rndgen.normal(loc=0.0, scale=1.0 / snr_use, size=outspec.shape)
    outspec += noise
    return nuspec_Hz, outspec

def getchisqfunc(nuspec_Hz: np.ndarray[float],
                 specvals: np.ndarray[float],
                 line: fc.Line = fc.ne8_770,
                 lsf_sigma_kmps: float = 10.,
                 ):
        wls_A = (c.c / nuspec_Hz) * 1e8
        vbins_kmps = (wls_A / line.wavelength_A - 1.) * c.c * 1e-5

        def getchisq(fitpars: np.ndarray[float]):
            try_logcds_cm2 = fitpars[0::3]
            try_bvals_kmps = fitpars[1::3] * _bnorm_fit
            try_centers_kmps = fitpars[2::3] * _vnorm_fit
            # no negative widths
            if np.any(try_bvals_kmps) <= 0.:
                 return 1e30 
            else:
                try_normflux = line.getspectrum(
                    nuspec_Hz, try_logcds_cm2, try_bvals_kmps,
                    try_centers_kmps)
                try_convflux = convgauss(vbins_kmps, try_normflux,
                                         lsf_sigma_kmps)
            sumsq = np.sum((try_convflux - specvals)**2)
            #print('try_notmflux: ', self.try_normflux)
            #print('try_convflux: ', self.try_convflux)
            #print('target: ', self.target)
            #print('fitfunc chi^2: ', sumsq)
            return sumsq
        return getchisq

def fitnoisyspec(nuspec_Hz: np.ndarray[float],
                 specvals: np.ndarray[float],
                 line: fc.Line = fc.ne8_770,
                 ncomp: int = 2,
                 lsf_sigma_kmps: float = 10.,
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    fitfunc = getchisqfunc(nuspec_Hz, specvals, line=line,
                           lsf_sigma_kmps=lsf_sigma_kmps)
    rhobeg = 1.
    bounds = [_bounds_logN_cm2,
              _bounds_b_kmps / _bnorm_fit, 
              _bounds_vcen_kmps / _vnorm_fit] * ncomp
    components_guess = [[14., 30. /  _bnorm_fit, 
                         -20. + 40. * float(i) / ncomp] for i in ncomp]
    components_guess = np.array([v for l in components_guess for v in l])

    minimizer_kwa = {'args': (lsf_sigma_kmps,), 
                     'bounds': bounds,
                     'method': 'COBYLA',
                     'options': {'rhobeg': rhobeg}}
    optimizeresult = so.basinhopping(
        fitfunc, components_guess, T=1e-2, stepsize=3.,
        minimizer_kwargs=minimizer_kwa
    )

    if not optimizeresult.success:
        msg = (f'Fitting spectrum with {ncomp} '
                'components failed. scipy.optimize message:\n')
        msg = msg + str(optimizeresult.message)
        raise RuntimeError(msg)

    components = np.copy(optimizeresult.x)
    logcds_cm2 = components[0::3]
    bvals_kmps = components[1::3] * _bnorm_fit
    vcens_kmps = components[2::3] * _vnorm_fit
    return logcds_cm2, bvals_kmps, vcens_kmps

_ncomp_default = [1, 2]
def runfits(logN_cm2: np.ndarray[float],
            b_kmps: np.ndarray[float],
            vcen_kmps: np.ndarray[float],
            ncomp_fit: list[int] = _ncomp_default,
            h5group=None):

    speckw = {"lsf_sigma_kmps": 10.,
              "snr": 20.,
              "line": fc.ne8_770,
              } 
    nuspec_Hz, fitspec = noisyspectrum(logN_cm2, b_kmps, vcen_kmps, **speckw)
    comp_found = {}
    for ncomp in ncomp_fit:
        res = fitnoisyspec(nuspec_Hz, fitspec, line=speckw["line"],
                           ncomp=ncomp, 
                           lsf_sigma_kmps=speckw["lsf_sigma_kmps"],
                           )
        comp_found[ncomp] = res
    if h5group is not None:
        h5group.attrs.create("lsf_sigma_kmps", speckw["lsf_sigma_kmps"])
        h5group.attrs.create("snr", speckw["snr"])
        _grp = h5group.create_group("line")
        (speckw["line"]).save(_grp)
        h5group.create_dataset("nu_Hz", data=nuspec_Hz)
        h5group.create_dataset("noisyspec", data=fitspec)
        h5group.attrs.create("logN_cm2_in", logN_cm2)
        h5group.attrs.create("b_kmps_in", b_kmps)
        h5group.attrs.create("vcen_kmps_in", vcen_kmps)
        for ncomp in ncomp_fit:
            sgrp = h5group.create_group(f"fit_{ncomp}comp")
            sgrp.attrs.create("logN_cm2", comp_found[ncomp][0])
            sgrp.attrs.create("b_kmps", comp_found[ncomp][1])
            sgrp.attrs.create("vcen_kmps", comp_found[ncomp][2])
    return comp_found

        
def runfitgrid():
    h5filen = ""
    c1_logN_cm2 = [14.5, 14.25, 14.0, 13.75, 13.5]
    c1_b_kmps = [120., 60., 30., 15.]
    c1_v_kmps = [0.]
    
    c2_logN_cm2_diff = [-0.1, -0.2, -0.3, -0.5, -1.0]
    c2_b_kmps = [120., 60., 30., 15.]
    c2_v_kmps = [-100., -50., 20., 10.]
    
    with h5py.File(h5filen, "w") as fo:
        for lN1 in c1_logN_cm2:
            for b1 in c1_b_kmps:
                for v1 in c1_v_kmps:
                    logN_cm2 = np.array([lN1])
                    b_kmps = np.array([b1])
                    vcen_kmps = np.array([v1])
                    grpn = f"comp1_N{lN1:.2f}_b{b1:.0f}_v{v1:.0f}"
                    grp = fo.create_group(grpn)
                    runfits(logN_cm2, b_kmps, vcen_kmps,
                            ncomp_fit=[1, 2], h5group=grp)

                    for lNd in c2_logN_cm2_diff:
                        for b2 in c2_b_kmps:
                            for v2 in c2_v_kmps:
                                lN2 = lN1 + lNd
                                logN_cm2 = np.array([lN1, lN2])
                                b_kmps = np.array([b1, b2])
                                vcen_kmps = np.array([v1, v2])

                                grpn = (f"comp2_N{lN1:.2f}_{lN2:.2f}"
                                        f"_b{b1:.0f}_{b2:.0f}"
                                        f"_v{v1:.0f}_{v2:.0f}")
                                grp = fo.create_group(grpn)
                                runfits(logN_cm2, b_kmps, vcen_kmps,
                                        ncomp_fit=[1, 2], h5group=grp)

    
