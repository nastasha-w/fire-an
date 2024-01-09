from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel
from dataclasses import dataclass
import numpy as np
import pandas as pd
import scipy.optimize as so
import scipy.special as sp
import scipy.stats as st

import fire_an.utils.constants_and_units as c

@dataclass
class Line:
    name: str
    wavelength_A: float
    fosc: float
    atrans_Hz: float

    def getspectrum(self,
                    nuspec_Hz: np.array[float],
                    logcds_cm2: np.array[float],
                    bvals_kmps: np.array[float],
                    centers_kmps: np.array[float]) -> np.array[float]:
        # from wikipedia: Voigt profile
        # for a Voigt function (normalized to 1) V(x; sigma, gamma)
        # where sigma is the gaussian width sigma
        # and gamma is the lorentzian gamma: 
        #   L(x; gamma) = gamma / (pi * (x**2 + gamma**2))
        # then V(x; sigma, gamma) = Re[wofz(z)] / (sigma * sqrt(2 * pi))
        # where z = (x + i * gamma) / (sigma * sqrt(2))
        #
        # normalization:
        # https://casper.astro.berkeley.edu/astrobaki/index.php/Line_Profile_Functions
        # that was tricky to find though, and seems to depend on
        # wavelength vs. frequency space for V(x; sigma, gamma)
        nuzero = c.c / (self.wavelength_A * 1e-8)
        nucenters_Hz = nuzero * (1. - centers_kmps * 1e5 / c.c)
        cds_cm2 = 10**logcds_cm2
        sigmas_Hz = bvals_kmps * 1e5 \
                    / (self.wavelength_A * 1e-8 * 2.**0.5)
        hwhm_cauchy_Hz = self.atrans_Hz / (4. * np.pi)
        norms = (np.pi * c.electroncharge**2 / (c.electronmass * c.c) \
                 * self.fosc * cds_cm2)
        xsample = (nuspec_Hz[:, np.newaxis] - nucenters_Hz[np.newaxis, :])
        z_in = (xsample + hwhm_cauchy_Hz * 1j) \
                / (sigmas_Hz[np.newaxis, :] * 2.**0.5)   
        vps = np.real(sp.wofz(z_in)) / (sigmas_Hz * (2. * np.pi)**0.5)
        tau = np.sum(vps * norms[:, np.newaxis], axis=1)
        normflux = 1. - np.exp(-1. * tau)
        return normflux
    
    def tau_to_coldens(self, spectau, specv_kmps):
        specnu = c.c \
                 / (self.wavelength_A * 1e-8 * (1. + specv_kmps * 1e5/ c.c))
        # expecting small v range compared to c
        dnu = np.average(np.diff(specnu)) 
        norm = (c.electronmass * c.c) \
               / (np.pi * c.electroncharge**2 * self.fosc)
        coldens = norm * np.sum(spectau) * dnu
        return coldens
class Spectrum:
    
    def __init__(self, line: Line, filen: str):
        self.line = line
        self.filen = filen
        self._readin_txtspectrum()

    def _readin_txtspectrum(self):
        spec = pd.read_csv(self.filen, comment='#', 
                           columns=['velocity_kmps', 'tau', 'flux', 'flux_error'])
        self.spec_raw = np.array(spec['flux'])
        self.vel_kmps = np.array(spec['velocity_kmps'])

    def _convolve_gauss(self, spec: np.array[float], width_kmps: float=30.):
        self.dv = np.average(np.diff(self.spec_raw))
        self.width_bins = width_kmps / self.dv
        self.kernel = Gaussian1DKernel(self.width_bins)
        spec_smoothed = convolve(spec, self.kernel, 
                                 boundary='fill', fill_value=1.)
        return spec_smoothed

    def _fitfunc(self,
                 fitpars: np.array[float],
                 specres_kmps: float) -> float:
        '''
        fitpars:
            array of (log column density [cm**-2], 
                      b value [km/s], 
                      line center [kmps])
            of each component to fit, in order.
            (array.reshape(ncomp, 3) 
            would give a table of values for each component)
        specres_kmps: 
            spectral resolution in km/s; should match that of the
            target array
        '''
        logcds_cm2 = 10**fitpars[0::3]
        bvals_kmps = fitpars[1::3]
        centers_kmps = fitpars[2::3]
        normflux = self.line.getspectrum(logcds_cm2, bvals_kmps, centers_kmps)
        convflux = self._convolve_gauss(normflux, specres_kmps)
        sumsq = np.sum((convflux - self.target)**2)
        return sumsq

    def _ftest(self, 
               prev: np.array[float], 
               current: np.array[float], 
               target: np.array[float], 
               ncomp_prev: int,
               ncomp_cur: int,
               snr: float) -> bool:
        '''
        Based on:
        https://online.stat.psu.edu/stat501/lesson/6/6.2 
        https://en.wikipedia.org/wiki/F-test
        '''
        # snr = sigma: for signal = 1 (continuum level), 
        # noise = sigma = 1
        chisq_prev = np.sum((prev - target)**2) / snr**2
        chisq_curr = np.sum((current - target)**2) / snr**2
        # larger fitting spectra are gaurantueed free of extra info
        nobs = len(self.spec_raw)
        pars_prev = 3 * ncomp_prev # each component has N, b, v
        pars_curr = 3 * ncomp_cur
        fstat = ((chisq_prev - chisq_curr) / (pars_curr - pars_prev)) \
                / (chisq_curr / (nobs - pars_curr))
        pval = st.fdist.cdf(fstat, pars_curr - pars_prev, nobs - pars_curr)
        return pval > self.minpval_newcomp
        
    def findcomponents(self, 
                       resolution_kmps: float=30., 
                       snr: float=30.,
                       sigma_det: float=3.,
                       save: bool = False):
        self.minpval_newcomp = st.NormalDist().cdf(sigma_det)
        self.target = self._convolve_gauss(self.spec_raw, resolution_kmps)

        self.fit_prev = np.ones(len(self.target)) # 0 components
        self.components_prev = np.array((0,), dtype=np.float)
        # first guess for the new component
        self.addcomp = True
        self.ncomp = 0
        while self.addcomp:
            self.ncomp += 1
            # where in the spectrum does it look like things are weirdest
            v0 = self.vel_kmps[np.argmax(np.abs(self.target 
                                                - self.fit_prev))] \
                - (len(self.target) - len(self.spec_raw)) // 2
            # seem like reasonable values for an absorber
            logN0 = 13.
            b0 = 20.
            self.components_guess = np.append(self.components_prev,
                                              np.array([logN0, b0, v0]))
            optimizeresult = so.optimize.minimize(self._fitfunc, 
                                                  self.components_guess,
                                                  args=(resolution_kmps,))
            if not optimizeresult.succes:
                msg = (f'Fitting spectrum with {self.ncomp} '
                       'components failed. scipy.optimize message:\n')
                msg = msg + optimizeresult.message
                raise RuntimeError(msg)

            self.components_curr = np.copy(optimizeresult.x)
            logcds = self.components_curr[0::3]
            bvals = self.components_curr[1::3]
            vcens = self.components_curr[2::3]
            self.fit_cur = self.getspectrum(logcds, bvals, vcens)
            self.addcomp = self._ftest(self.fit_prev, self.fit_cur,
                                       self.target, self.ncomp - 1,
                                       self.ncomp, snr)
        self.ncomp -= 1
        self.logcds_cm2_fit = self.components_prev[0::3]
        self.bvals_kmps_fit = self.components_prev[1::3]
        self.vcens_kmps_fit = self.components_curr[2::3]
        print(f'Found {self.ncomp} components')
        