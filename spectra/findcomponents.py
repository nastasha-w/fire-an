from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel
from dataclasses import dataclass
import h5py
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
                    nuspec_Hz: np.ndarray[float],
                    logcds_cm2: np.ndarray[float],
                    bvals_kmps: np.ndarray[float],
                    centers_kmps: np.ndarray[float]) -> np.ndarray[float]:
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
        print(nuzero)
        nucenters_Hz = nuzero * (1. - centers_kmps * 1e5 / c.c)
        print(nucenters_Hz)
        cds_cm2 = 10**logcds_cm2
        sigmas_Hz = bvals_kmps * 1e5 \
                    / (self.wavelength_A * 1e-8 * 2.**0.5)
        print(sigmas_Hz)
        hwhm_cauchy_Hz = self.atrans_Hz / (4. * np.pi)
        norms = (np.pi * c.electroncharge**2 / (c.electronmass * c.c) \
                 * self.fosc * cds_cm2)
        print(norms)
        xsample = (nuspec_Hz[:, np.newaxis] - nucenters_Hz[np.newaxis, :])
        z_in = (xsample + hwhm_cauchy_Hz * 1j) \
                / (sigmas_Hz[np.newaxis, :] * 2.**0.5)   
        vps = np.real(sp.wofz(z_in)) / (sigmas_Hz * (2. * np.pi)**0.5)
        print(vps)
        tau = np.sum(vps * norms[np.newaxis, :], axis=1)
        print(tau)
        normflux = np.exp(-1. * tau)
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
    
    def save(self, h5grp):
        h5grp.attrs.create('name', np.string_(self.name))
        h5grp.attrs.create('wavelength_A', self.wavelength_A)
        h5grp.attrs.create('fosc', self.fosc)
        h5grp.attrs.create('atrans_Hz', self.atrans_Hz)

# copied from Trident file
# Verner, Verner, & Ferland (1996): 
# lambda = 770.4089 A
# A_ki = 5.79e8 Hz
# fosc = 1.03e-1
ne8_770 = Line(name='Ne VIII 770', wavelength_A=770.409000,
               fosc=1.020000e-01, atrans_Hz=5.720000e+08)

class Spectrum:
    
    def __init__(self, line: Line, filen: str):
        self.line = line
        self.filen = filen
        self._readin_txtspectrum()

    def _readin_txtspectrum(self):
        spec = pd.read_csv(self.filen, comment='#', 
                           columns=['velocity_kmps', 'tau', 
                                    'flux', 'flux_error'])
        self.spec_raw = np.array(spec['flux'])
        self.vel_kmps = np.array(spec['velocity_kmps'])

    def _convolve_gauss(self, spec: np.ndarray[float], width_kmps: float=30.):
        self.dv = np.average(np.diff(self.spec_raw))
        self.width_bins = width_kmps / self.dv
        self.kernel = Gaussian1DKernel(self.width_bins)
        spec_smoothed = convolve(spec, self.kernel, 
                                 boundary='fill', fill_value=1.)
        return spec_smoothed

    def _fitfunc(self,
                 fitpars: np.ndarray[float],
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
               prev: np.ndarray[float], 
               current: np.ndarray[float], 
               target: np.ndarray[float], 
               ncomp_prev: int,
               ncomp_cur: int,
               snr: float) -> bool:
        '''
        Based on:
        https://online.stat.psu.edu/stat501/lesson/6/6.2 
        https://en.wikipedia.org/wiki/F-test
        '''
        # snr = 1 / sigma: for signal = 1 (continuum level), 
        # noise = sigma = 1
        chisq_prev = np.sum((prev - target)**2) * snr**2
        chisq_curr = np.sum((current - target)**2) * snr**2
        # larger fitting spectra are gaurantueed free of extra info
        nobs = len(self.spec_raw)
        pars_prev = 3 * ncomp_prev # each component has N, b, v
        pars_curr = 3 * ncomp_cur
        fstat = ((chisq_prev - chisq_curr) / (pars_curr - pars_prev)) \
                / (chisq_curr / (nobs - pars_curr))
        pval = st.fdist.cdf(fstat, pars_curr - pars_prev, nobs - pars_curr)
        return pval > self.minpval_newcomp
        
    def findcomponents_1pert(self, 
                             resolution_kmps: float=30., 
                             snr: float | None =30.,
                             sigma_det: float=3.) -> (int, 
                                                      np.ndarray[float],
                                                      np.ndarray[float],
                                                      np.ndarray[float]):
        self.minpval_newcomp = st.NormalDist().cdf(sigma_det)
        self.target = self._convolve_gauss(self.spec_raw, resolution_kmps)
        if snr is not None:
            self._sigma = 1. / snr
            self._delta = np.random.normal(loc=0.0, scale=self._sigma, 
                                           size=len(self.target))
            self.target += self._delta
            self.target = np.max(self.target, 0.) # don't allow values < 0
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
        self.vcens_kmps_fit = self.components_prev[2::3]
        print(f'Found {self.ncomp} components')
        return (self.ncomp, self.logcds_cm2_fit, 
                self.bvals_kmps_fit, self.vcens_kmps_fit)
    
    def findcomponents(self, 
                       resolution_kmps: float=30., 
                       snr: float | None =30.,
                       sigma_det: float=3.,
                       npert: int=20,
                       save: tuple[str, str] | None =None):
        '''
        Find the absorption components and their parameters in the 
        input spectrum.

        Parameters:
        -----------
        resolution_kmps:
            spectral resolution at which to search [km/s]. Should match
            a target instrument resolution.
        snr:
            signal-to-noise ratio at which to search. Should match a
            target observation SNR.
        sigma_det:
            the significance at which a component counts as detected.
            Used for the F-test random chance probability.
        npert:
            number of noise realizations used to determine how many 
            components to fit.
        save:
            file name (string, including the full directory path) 
            and hdf5 group within that file
            to save the results to, or None to not save.
            Only save the fitted parameters, not the spectra that go
            with them.

        Returns:
        --------
        None. (Results are stored as object attributes.)

        Raises:
        -------
        RuntimeError if one of the component fits fails.

        Notes:
        ------
        Finds components in a number of steps.
        (1) Determine how many components to fit.
            This is done by:
            (a) Convolving the target (read-in) normalized flux
                spectrum with a Gaussian of width <resolution_kmps>
            (b) Adding a random noise realization to this spectrum
                according to the given SNR.
            (c) Iteratively fitting components to this spectrum,
                until an F-test determines the newest component is
                not statistically justified at the level <sigma_det>
            (d) Repeating (a)-(c) <npert> times, and determining the
                median number of components found.
        (2) Fit components to the unpertubed spectrum.
            Fit the number of components determined in step (1d) to
            the read-in spectrum, at the resolution in step (1a), 
            but without the noise from (1b). 
        '''

        self.noisyfits = []
        for _ in range(npert):
            _ncomp, _logcds_cm2, _bvals_kmps, _vcens_kmps = \
                self.findcomponents_1pert(resolution_kmps=resolution_kmps, 
                                          snr=snr, sigma_det=sigma_det)
            self.noisyfits.append({'ncomp': _ncomp,
                                   'logN_cm2': _logcds_cm2,
                                   'bvals_kmps': _bvals_kmps,
                                   'vcens_kmps': _vcens_kmps})
        ncomp_noisy = np.array([fit['ncomp'] for fit in self.noisyfits])
        self.ncomp_final = int(np.ceil((np.median(ncomp_noisy))))
        if np.any(ncomp_noisy == self.ncomp_final):
            guessi = np.where(ncomp_noisy == self.ncomp_final)[0][0]
            guessN = self.noisyfits[guessi]['logN_cm2']
            guessb = self.noisyfits[guessi]['bvals_kmps']
            guessv = self.noisyfits[guessi]['vcens_kmps']
            self.guess = np.empty((self.ncomp_final, 3), dtype=guessN.dtype)
            self.guess[:, 0] = guessN
            self.guess[:, 1] = guessb
            self.guess[:, 2] = guessv
            self.guess.reshape((self.ncomp_final * 3,))
        else:
            # get a spectrum with the smallest number of extra components 
            opts = np.where(ncomp_noisy > self.ncomp_final)[0]
            subi = np.argmin(ncomp_noisy[opts])
            guessi = opts[subi]
            
            guessN = self.noisyfits[guessi]['logN_cm2']
            guessb = self.noisyfits[guessi]['bvals_kmps']
            guessv = self.noisyfits[guessi]['vcens_kmps']
            # largest-logN components within that fit
            compsel = np.argsort(guessN)[::-1][:self.ncomp_final]
            self.guess = np.empty((self.ncomp_final, 3), dtype=guessN.dtype)
            self.guess[:, 0] = guessN[compsel]
            self.guess[:, 1] = guessb[compsel]
            self.guess[:, 2] = guessv[compsel]
            self.guess.reshape((self.ncomp_final * 3,))

        optimizeresult = so.optimize.minimize(self._fitfunc, 
                                              self.guess,
                                              args=(resolution_kmps,))
        if not optimizeresult.succes:
            msg = (f'Fitting spectrum with {self.ncomp} '
                    'components failed. scipy.optimize message:\n')
            msg = msg + optimizeresult.message
            raise RuntimeError(msg)

        self.components_final = np.copy(optimizeresult.x)
        self.logN_cm2_final = self.components_final[0::3]
        self.bvals_kmps_final = self.components_final[1::3]
        self.vcens_kmps_final = self.components_final[2::3]

        if save is not None:
            # just save the parameters, can get the spectra from there
            with h5py.File(save[0], 'a')as f:
                grp = f.create_group(save[1])
                _grp = grp.create_group('line')
                self.line.save(_grp)
                grp.attrs.create('spectrumfilen', np.string_(self.filen))
                grp.attrs.create('resolution_kmps', resolution_kmps)
                grp.attrs.create('snr', snr)
                grp.attrs.create('sigma_det', sigma_det)
                grp.attrs.create('npert', npert)
                _grp = grp.create_group('noisy_fits')
                for i in range(len(self.noisyfits)):
                     _sgrp = _grp.create_group(f'fit{i}')
                     dct = self.noisyfits[i]
                     for key in dct:
                         _sgrp.create_dataset(key, data=dct[key])
                _grp = grp.create_group('final_fit')
                _grp.create_dataset('ncomp', self.ncomp_final)
                _grp.create_dataset('logN_cm2', self.logN_cm2_final)
                _grp.create_dataset('bvals_kmps', self.bvals_kmps_final)
                _grp.create_dataset('vcens_kmps', self.vcens_kmps_final)
                     
        