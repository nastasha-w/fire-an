from typing import Any
from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel
from dataclasses import dataclass
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize as so
import scipy.special as sp
import scipy.stats as st


import fire_an.utils.constants_and_units as c

# Bayesian model comparison, with some methods details:
# https://www.imperial.ac.uk/media/imperial-college/research-centres-and-groups/astrophysics/public/icic/data-analysis-workshop/2018/Model_comparison_TrottaSept_2018.pdf
@dataclass
class Line:
    name: str
    wavelength_A: float
    fosc: float
    atrans_Hz: float
    
    # tested: centers and widths look ok, 
    # transmission seems ballpark right
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
        #print(nuzero)
        nucenters_Hz = nuzero * (1. - centers_kmps * 1e5 / c.c)
        #print(nucenters_Hz)
        cds_cm2 = 10**logcds_cm2
        sigmas_Hz = bvals_kmps * 1e5 \
                    / (self.wavelength_A * 1e-8 * 2.**0.5)
        #print(sigmas_Hz)
        hwhm_cauchy_Hz = self.atrans_Hz / (4. * np.pi)
        norms = (np.pi * c.electroncharge**2 / (c.electronmass * c.c) \
                 * self.fosc * cds_cm2)
        xsample = (nuspec_Hz[:, np.newaxis] - nucenters_Hz[np.newaxis, :])
        z_in = (xsample + hwhm_cauchy_Hz * 1j) \
                / (sigmas_Hz[np.newaxis, :] * 2.**0.5)   
        vps = np.real(sp.wofz(z_in)) / (sigmas_Hz * (2. * np.pi)**0.5)
        tau = np.sum(vps * norms[np.newaxis, :], axis=1)
        normflux = np.exp(-1. * tau)
        return normflux
    
    # tested: consistent with spectrum generation
    def tau_to_coldens(self, spectau, specv_kmps):
        specnu = c.c \
                 / (self.wavelength_A * 1e-8 * (1. + specv_kmps * 1e5/ c.c))
        # expecting small v range compared to c
        dnu = np.average(np.diff(specnu)) 
        norm = (c.electronmass * c.c) \
               / (np.pi * c.electroncharge**2 * self.fosc)
        coldens = norm * np.sum(spectau) * np.abs(dnu)
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

class SpectrumFitFreq:
    
    def __init__(self, line: Line, filen: str | None=None):
        self.line = line
        self.filen = filen
        if self.filen is not None:
            self.readin_txtspectrum()
            self.set_fitdefaults()
    
    def set_fitdefaults(self,
                        fitu_vcens_kmps: float = 50.,
                        fitu_bvals_kmps: float = 10.,
                        bounds_logN: tuple[float | None, float | None] 
                                     = (9., None),
                        bounds_bvals_kmps: tuple[float | None, float | None]
                                     = (0., 200.)):
        '''
        Set bounds and range normalizations for the component fitting.
        '''
        self.fitu_vcens_kmps = fitu_vcens_kmps
        self.fitu_bvals_kmps = fitu_bvals_kmps
        
        bbounds = (bounds_bvals_kmps[0] if bounds_bvals_kmps[0] is None
                   else bounds_bvals_kmps[0] / self.fitu_bvals_kmps,
                   bounds_bvals_kmps[1] if bounds_bvals_kmps[1] is None
                   else bounds_bvals_kmps[1] / self.fitu_bvals_kmps,
                   )
        self.bounds_base = [bounds_logN, bbounds,
                            ((self.vel_kmps[0] - self.dv) \
                             / self.fitu_vcens_kmps, 
                             (self.vel_kmps[-1] + self.dv)\
                             / self.fitu_vcens_kmps)]
        self.rhobeg = 1.

    def readin_txtspectrum(self):
        '''
        read in spectrum from .txt file output of Trident
        '''
        spec = pd.read_csv(self.filen, comment='#',
                           sep=' ',
                           names=['velocity_kmps', 'tau', 
                                  'flux', 'flux_error'],
                           dtype=float)
        self.spec_raw = np.array(spec['flux'])
        self.vel_kmps = np.array(spec['velocity_kmps'])
        self.tau_raw = np.array(spec['tau']) # for total N checks
        self.dv = np.average(np.diff(self.vel_kmps))
        self.set_fitdefaults()
    
    # mostly for testing
    def setrawspectrum(self,
                       velocity_kmps: np.ndarray[float],
                       flux: np.ndarray[float]):
        if velocity_kmps.shape != flux.shape:
            msg = (f'velocity_kmps {velocity_kmps.shape} and '
                   f'flux {flux.shape} should have the same shape')
            raise ValueError(msg)
        if len(velocity_kmps.shape) != 1:
            msg = (f'velocity_kmps {velocity_kmps.shape} and '
                   f'flux {flux.shape} should be arrays of dimension 1')
            raise ValueError(msg)
        self.spec_raw = flux
        self.vel_kmps = self.velocity_kmps
        self.dv = np.average(np.diff(self.vel_kmps))

    def _convolve_gauss(self,
                        spec: np.ndarray[float], 
                        width_kmps: float=30.):
        '''
        Notes:
        ------
        Assumes a spectrum with constant velocity resolution.
        '''
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
                      b value [km/s] / 10., 
                      line center [kmps] / 50.)
            of each component to fit, in order.
            (array.reshape(ncomp, 3) 
            would give a table of values for each component)
            The weird units are so that reasonable changes are
            similar for the different parameters
        specres_kmps: 
            spectral resolution in km/s; should match that of the
            target array
        '''
        self.try_logcds_cm2 = fitpars[0::3]
        self.try_bvals_kmps = fitpars[1::3] * self.fitu_bvals_kmps
        self.try_centers_kmps = fitpars[2::3] * self.fitu_vcens_kmps
        #print('fitpars try: ', fitpars)
        if np.any(self.try_bvals_kmps <= 0.):
            sumsq = 1e30 # some very large value
        else:
            self.try_normflux = self.line.getspectrum(
                self.nu_Hz, self.try_logcds_cm2,
                self.try_bvals_kmps, self.try_centers_kmps)
            self.try_convflux = self._convolve_gauss(self.try_normflux,
                                                     specres_kmps)
            sumsq = np.sum((self.try_convflux - self.target)**2)
            #print('try_notmflux: ', self.try_normflux)
            #print('try_convflux: ', self.try_convflux)
            #print('target: ', self.target)
        #print('fitfunc chi^2: ', sumsq)
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
        # F-distribution
        print('chi^2 old, new; fstat: ', chisq_prev, chisq_curr, fstat)
        pval = st.f.cdf(fstat, pars_curr - pars_prev, nobs - pars_curr)
        print('p value new component: ', pval)
        return pval > self.minpval_newcomp
        
    def findcomponents_1pert(self, 
                             resolution_kmps: float=30., 
                             snr: float | None =30.,
                             sigma_det: float=3.,
                             plotinterm: bool = True) -> tuple[int, 
                                                         np.ndarray[float],
                                                         np.ndarray[float],
                                                         np.ndarray[float]]:
        # one-sided test: no sucj thing as negatiev absorption
        # (no emission line in this model or in the mock spectra)
        self.minpval_newcomp = st.norm().cdf(sigma_det)
        self.target = self._convolve_gauss(self.spec_raw, resolution_kmps)
        if snr is not None:
            self._sigma = 1. / snr
            self._delta = np.random.normal(loc=0.0, scale=self._sigma, 
                                           size=len(self.target))
            self.target += self._delta
            self.target = np.maximum(self.target, 0.) # don't allow values < 0
        self.fit_cur = np.ones(len(self.target)) # 0 components
        self.components_cur = np.ones(shape=(0,), dtype=np.float64)
        # first guess for the new component
        self.addcomp = True
        self.ncomp = 0
        while self.addcomp:
            self.ncomp += 1
            self.components_prev = np.copy(self.components_cur)
            self.fit_prev = np.copy(self.fit_cur)
            # where in the spectrum does it look like things are weirdest
            v0 = self.vel_kmps[np.argmax(np.abs(self.target 
                                                - self.fit_prev))] \
                - self.dv * (len(self.target) - len(self.spec_raw)) // 2
            # seem like reasonable values for an absorber
            logN0 = 14.5
            b0 = 20.
            print('starting with new: ', logN0, b0, v0)
            newcomp = np.array([logN0, 
                                b0 / self.fitu_bvals_kmps, 
                                v0 / self.fitu_vcens_kmps])
            self.components_guess = np.append(self.components_prev,
                                              newcomp)
            # minimum low column density (close enough to zero),
            # only positive line widths (max: spectrum width)
            # velocity center should be within the spectrum range
            self.bounds = self.bounds_base * self.ncomp
            # reaonable parameter changes
            #self.optimizeresult_cur = so.minimize(self._fitfunc, 
            #                             self.components_guess,
            #                             args=(resolution_kmps,),
            #                             bounds=self.bounds,
            #                             method='COBYLA',
            #                             options={'rhobeg': self.rhobeg,
            #                                      'tol': 1e-1})
            #                             # 'maxiter': 100,
            self.optimizeresult_cur = so.basinhopping(
                self._fitfunc, self.components_guess, 
                T=1e-2, 
                stepsize=3.,
                minimizer_kwargs={'args': (resolution_kmps,), 
                                  'bounds': self.bounds,
                                  'method': 'COBYLA',
                                  'options': {'rhobeg': self.rhobeg,
                                              }}
                )
                # 'tol': 3e-3, 1e-1 * 1./snr**2 * len(self.spec_raw),
            # assume a failed fit means there wasn't much to be done
            if not self.optimizeresult_cur.success:
                msg = (f'Fitting spectrum with {self.ncomp} '
                       'components failed. scipy.optimize message:\n')
                msg = msg + str(self.optimizeresult_cur.message)
                break 
                #raise RuntimeError(msg)
            print('extra comp. fit: ', self.optimizeresult_cur.x)

            self.components_cur = np.copy(self.optimizeresult_cur.x)
            self.logcds_cur = self.components_cur[0::3]
            self.bvals_cur = self.components_cur[1::3] \
                             * self.fitu_bvals_kmps
            self.vcens_cur = self.components_cur[2::3] \
                             * self.fitu_vcens_kmps
            self.fit_cur = self.line.getspectrum(
                self.nu_Hz, self.logcds_cur, self.bvals_cur, self.vcens_cur)
            self.addcomp = self._ftest(self.fit_prev, self.fit_cur,
                                       self.target, self.ncomp - 1,
                                       self.ncomp, snr)
        self.ncomp -= 1
        self.logcds_cm2_fit = self.components_prev[0::3]
        self.bvals_kmps_fit = self.components_prev[1::3] \
                              * self.fitu_bvals_kmps
        self.vcens_kmps_fit = self.components_prev[2::3] \
                              * self.fitu_vcens_kmps
        print(f'Found {self.ncomp} components')

        return (self.ncomp, self.logcds_cm2_fit, 
                self.bvals_kmps_fit, self.vcens_kmps_fit)
    
    def findcomponents(self, 
                       resolution_kmps: float=30., 
                       snr: float | None =30.,
                       sigma_det: float=3.,
                       npert: int=20,
                       save: tuple[str, str] | None =None,
                       plotinterm: bool = True):
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
        self.nu_Hz = c.c / (self.line.wavelength_A * 1e-8 \
                            * (1. + self.vel_kmps * 1e5 / c.c))

        self.noisyfits = []
        self.noisyspectra = []
        for _ in range(npert):
            _ncomp, _logcds_cm2, _bvals_kmps, _vcens_kmps = \
                self.findcomponents_1pert(resolution_kmps=resolution_kmps, 
                                          snr=snr, sigma_det=sigma_det,
                                          plotinterm=plotinterm)
            self.noisyfits.append({'ncomp': _ncomp,
                                   'logN_cm2': np.copy(_logcds_cm2),
                                   'bvals_kmps': np.copy(_bvals_kmps),
                                   'vcens_kmps': np.copy(_vcens_kmps)})
            self.noisyspectra.append(self.target.copy())
        ncomp_noisy = np.array([fit['ncomp'] for fit in self.noisyfits])
        self.ncomp_final = int(np.ceil((np.median(ncomp_noisy))))
        if self.ncomp_final == 0:
            self.logN_cm2_final = np.zeros(shape=(0,), dtype=np.float64)
            self.bvals_kmps_final = np.zeros(shape=(0,), dtype=np.float64)
            self.vcens_kmps_final = np.zeros(shape=(0,), dtype=np.float64)
        else:
            if np.any(ncomp_noisy == self.ncomp_final):
                guessi = np.where(ncomp_noisy == self.ncomp_final)[0][0]
                guessN = self.noisyfits[guessi]['logN_cm2']
                guessb = self.noisyfits[guessi]['bvals_kmps']
                guessv = self.noisyfits[guessi]['vcens_kmps']
                self.guess = np.empty((self.ncomp_final, 3), 
                                      dtype=guessN.dtype)
                self.guess[:, 0] = guessN
                self.guess[:, 1] = guessb
                self.guess[:, 2] = guessv
                self.guess.shape = (self.ncomp_final * 3,)
            else:
                # get a spectrum with the smallest number of 
                # extra components
                opts = np.where(ncomp_noisy > self.ncomp_final)[0]
                subi = np.argmin(ncomp_noisy[opts])
                guessi = opts[subi]
                
                guessN = self.noisyfits[guessi]['logN_cm2']
                guessb = self.noisyfits[guessi]['bvals_kmps']
                guessv = self.noisyfits[guessi]['vcens_kmps']
                # largest-logN components within that fit
                compsel = np.argsort(guessN)[::-1][:self.ncomp_final]
                self.guess = np.empty((self.ncomp_final, 3), 
                                      dtype=guessN.dtype)
                self.guess[:, 0] = guessN[compsel]
                self.guess[:, 1] = guessb[compsel]
                self.guess[:, 2] = guessv[compsel]
                self.guess.shape = (self.ncomp_final * 3,)
            self.optimizeresult = so.basinhopping(
                self._fitfunc, self.guess, 
                T=1e-2,
                stepsize=3.,
                minimizer_kwargs={'args': (resolution_kmps,), 
                                  'bounds': self.bounds_base
                                           * self.ncomp_final,
                                  'method': 'COBYLA',
                                  'options': {'rhobeg': self.rhobeg,
                                              }})
                # 'tol': 1e-3, T=1e-1 * 1./snr**2 * len(self.spec_raw)
            if not self.optimizeresult.success:
                msg = (f'Fitting spectrum with {self.ncomp} '
                        'components failed. scipy.optimize message:\n')
                msg = msg + str(self.optimizeresult.message)
                raise RuntimeError(msg)

            self.components_final = np.copy(self.optimizeresult.x)
            self.logN_cm2_final = self.components_final[0::3]
            self.bvals_kmps_final = self.components_final[1::3] \
                                    * self.fitu_bvals_kmps
            self.vcens_kmps_final = self.components_final[2::3] \
                                    * self.fitu_vcens_kmps

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
    def plotfits(self,
                 subi: int | None =None,
                 resolution_kmps: float =30.):
        plt.plot(self.vel_kmps, self.spec_raw, label='raw spectrum', 
                 color='blue', linewidth=1.7)
        convflux = self._convolve_gauss(self.spec_raw, resolution_kmps)
        plt.plot(self.vel_kmps, convflux, label='smoothed spec.',
                 color='black', linewidth=1.7, linestyle='solid')
        plt.xlabel('l.o.s. velocity [km/s]')
        plt.ylabel('transmission')
        
        if subi is None:
            for ci in range(len(self.logN_cm2_final)):
                if ci == 0:
                    label = 'comp.'
                else:
                    label = None
                cs =  self.line.getspectrum(self.nu_Hz,
                                            self.logN_cm2_final[ci, None],
                                            self.bvals_kmps_final[ci, None],
                                            self.vcens_kmps_final[ci, None])
                plt.plot(self.vel_kmps, cs, label=label,
                         color='gray', linewidth=1., linestyle='dotted')
            cs =  self.line.getspectrum(self.nu_Hz,
                                        self.logN_cm2_final,
                                        self.bvals_kmps_final,
                                        self.vcens_kmps_final)
            convflux = self._convolve_gauss(cs, resolution_kmps)
            plt.plot(self.vel_kmps, convflux, label='fit',
                     color='gray', linewidth=1.3, linestyle='dashed')
        else:
            target = self.noisyspectra[subi]
            fit = self.noisyfits[subi]
            bvals = fit['bvals_kmps']
            vcens = fit['vcens_kmps']
            logN = fit['logN_cm2']

            plt.plot(self.vel_kmps, target, label='target',
                     color='black', linewidth=1., linestyle='solid')

            for ci in range(len(logN)):
                if ci == 0:
                    label = 'comp.'
                else:
                    label = None
                cs =  self.line.getspectrum(self.nu_Hz,
                                            logN[ci, None],
                                            bvals[ci, None],
                                            vcens[ci, None])
                plt.plot(self.vel_kmps, cs, label=label,
                         color='gray', linewidth=1., linestyle='dotted')
            cs =  self.line.getspectrum(self.nu_Hz,
                                        logN, bvals, vcens)
            convflux = self._convolve_gauss(cs, resolution_kmps)
            plt.plot(self.vel_kmps, convflux, label='fit',
                     color='gray', linewidth=1.3, linestyle='dashed')                         
        plt.show()

class SpectrumFitBayes:
    
    def __init__(self, 
                 line: Line, 
                 filen: str | None = None):
        self.line = line
        self.filen = filen
        if self.filen is not None:
            self.readin_txtspectrum(self.filen)
        self._set_default_fitpars()
    
    def _set_default_fitpars(self):
        self.fitrange_b_kmps = (1., 200.)
        self.fitrange_logN_cm2 = (10., 16.)
        self.fitrange_v_kmps = (-500., 500.)

        self.fitunit_b_kmps = 40.
        self.fitunit_logN_cm2 = 1.
        self.fitunit_v_kmps = 100.
        self._set_fitrange_fitunits()
    
    def _set_fitrange_fitunits(self):
        self._fitrange_b = (None if self.fitrange_b_kmps[0] is None else
                            self.fitrange_b_kmps[0] / self.fitunit_b_kmps,
                            None if self.fitrange_b_kmps[1] is None else
                            self.fitrange_b_kmps[1] / self.fitunit_b_kmps
                            ) 
        self._fitrange_v = (None if self.fitrange_v_kmps[0] is None else
                            self.fitrange_v_kmps[0] / self.fitunit_v_kmps,
                            None if self.fitrange_v_kmps[1] is None else
                            self.fitrange_v_kmps[1] / self.fitunit_v_kmps
                            ) 
        self._fitrange_logN = (None if self.fitrange_logN_cm2[0] is None else
                            self.fitrange_logN_cm2[0] / self.fitunit_logN_cm2,
                            None if self.fitrange_logN_cm2[1] is None else
                            self.fitrange_logN_cm2[1] / self.fitunit_logN_cm2
                            ) 

    def readin_txtspectrum(self, filen: str):
        '''
        Read in the spectrum from a Trident .txt file. (Assumes a 
        in velocity space.)
        '''
        spec = pd.read_csv(self.filen, comment='#',
                           sep=' ',
                           names=['velocity_kmps', 'tau', 
                                  'flux', 'flux_error'],
                           dtype=float)
        self.spec_raw = np.array(spec['flux'])
        self.tau_raw = np.array(spec['tau'])
        self.vel_kmps = np.array(spec['velocity_kmps'])
        self.wl_A = self.line.wavelength_A * (1. + self.vel_kmps * 1e5 / c.c)
        self.nu_Hz = c.c / (self.wl_A * 1e-8)
    
    def add_mockspectrum(self,
                         components: np.ndarray[float], 
                         vbins_kmps: np.ndarray[float],
                         sigma_specres_kmps: float = 30.):
        '''
        Parameters:
        -----------
        components:
            absorption components, in fit units. (See getspectrum.)
        vbins_kmps:
            spectral bins in velocity. Should include any cosmological
            redshift.
        sigma_specres_kmps:
            spectral resolution (Gaussian sigma).
        '''
        self.vel_kmps = vbins_kmps
        self.wl_A = self.line.wavelength_A * (1. + self.vel_kmps * 1e5 / c.c)
        self.nu_Hz = c.c / (self.wl_A * 1e-8)
        self.spec_raw = self.getspectrum(components)
        #self.spec_raw = self._convolve_gauss(self._spec, sigma_specres_kmps)

    def setfitranges(self, 
                     logNranges_cm2: tuple[float, float] | None = None,
                     brange_kmps: tuple[float, float] | None = None,
                     vrange_kmps: tuple[float, float] | None = None,
                     ):
        '''
        Bayesian prior support regions. a `None` bound means no
        upper/lower bound on the value.
        '''
        if brange_kmps is not None:
            self.fitrange_b_kmps = brange_kmps
        if logNranges_cm2 is not None:
            self.fitrange_logNcm2 = logNranges_cm2
        if vrange_kmps is not None:
            self.fitrange_v_kmps = vrange_kmps
        self._set_fitrange_fitunits()
    
    def getspectrum(self,
                    components: np.ndarray[float],
                    sigma_specres_kmps: float = 30.):
        '''
        
        '''
        logcds_cm2 = 10**components[0::3] * self.fitunit_logN_cm2
        bvals_kmps = components[1::3] * self.fitunit_b_kmps
        centers_kmps = components[2::3] * self.fitunit_v_kmps
        normflux = self.line.getspectrum(self.nu_Hz, logcds_cm2, 
                                         bvals_kmps, centers_kmps)
        convflux = self._convolve_gauss(normflux, sigma_specres_kmps)
        return convflux

    def _convolve_gauss(self, spec: np.ndarray[float], width_kmps: float=30.):
        self.dv = np.average(np.diff(self.vel_kmps))
        self.width_bins = width_kmps / self.dv
        self.kernel = Gaussian1DKernel(self.width_bins)
        spec_smoothed = convolve(spec, self.kernel, 
                                 boundary='fill', fill_value=1.)
        return spec_smoothed
    
    def _noisyspec(self, 
                   rawspec: np.ndarray[float], 
                   sigma_kmps: float = 30., 
                   snr: float = 30.):
        out = self._convolve_gauss(rawspec, width_kmps=sigma_kmps)
        self.sigma_noise = 1. / snr
        self._delta = np.random.normal(loc=0.0, scale=self.sigma_noise, 
                                       size=len(out))
        out += self._delta
        return out

    def _logp_spec(self,
                   fitspec: np.ndarray[float],
                   components: np.ndarray[float],
                   snr: float = 30.,
                   sigma_specres_kmps: float = 30.):
        '''
        Parameters:
        -----------
        fitspec: 
            the spectrum to fit
        components:
            the components to evaluate the likelihood at (see getspectrum)
        Notes:
        ------
        Assumes Gaussian noise, with sigma = 1. / snr
        '''
        self._isigma = snr
        self._modelspec = self.getspectrum(components, 
                                           sigma_specres_kmps=
                                           sigma_specres_kmps)
        # single-bin log p = log (1 / (sigma * sqrt(2 pi)) 
        #                     * exp(- 0.5 * (model - value)^2 / sigma^2)
        #                  = np.log(invsigma / (2 * np.pi)) 
        #                    -0.5 * (model - value)**2 * invsigma**2 
        _logp = len(fitspec) * np.log(self._isigma / (2. * np.pi)) \
                -0.5 * np.sum(((fitspec - self._modelspec) * self._isigma)**2)
        return _logp
    
    def _fitspec(self, 
                 fitspec: np.ndarray[float], 
                 componentsguess: np.ndarray[float], 
                 sigma_specres_kmps: float = 30.) -> tuple[np.ndarray, 
                                                           np.ndarray,
                                                           Any]:
        '''
        returns the best-fit parameters, the best-fit spectrum,
        and whatever object the fitter gives back
        '''
        # some clever sampling method for high-dimensional spaces
        pass

    def _checknewcomponent(self, fitres):
        '''
        input: fitting object from _fitspec
        output: boolean: accept new component or not
        '''
        pass
    
    def _fitspec_single(self,
                        target_cur: np.ndarray, 
                        sigma_specres_kmps: float = 30.) -> tuple[int,
                                                                  np.ndarray,
                                                                  np.ndarray,
                                                                  np.ndarray]:
        guess_logN = 13.5 / self.fitunit_logN_cm2
        guess_bval_kmps = 20. / self.fitunit_b_kmps
        ncomp_cur = 0
        fitspec_cur = np.ones(len(target_cur))
        fitcomp_cur = np.zeros((0,), dtype=guess_logN.dtype)
        trynew = True
        while trynew:
            ncomp_cur += 1
            fitspec_prev = fitspec_cur
            fitcomp_prev = fitcomp_cur
            maxdiffi = np.argmax(np.abs(fitspec_prev - target_cur))
            guess_v_kmps = self.vel_kmps[maxdiffi] / self.fitunit_v_kmps
            newcomp_guess = np.array([guess_logN, 
                                      guess_bval_kmps, 
                                      guess_v_kmps])
            componentsguess = np.append(fitcomp_prev, newcomp_guess)
            fitcomp_cur, fitspec_cur, fitres = self._fitspec(
                target_cur, componentsguess,
                sigma_specres_kmps=sigma_specres_kmps)
            trynew = self._checknewcomponent(fitres)
        ncomp = ncomp_cur - 1
        logcds_cm2 = fitcomp_prev[0::3] * self.fitunit_logN_cm2
        bvals_kmps = fitcomp_prev[1::3] * self.fitunit_b_kmps
        vcens_kmps = fitcomp_prev[2::3] * self.fitunit_v_kmps
        return ncomp, logcds_cm2, bvals_kmps, vcens_kmps
            
    def fitcomponents(self, 
                      snr: float = 30.,
                      sigma_specres_kmps: float = 30.,
                      sigma_det: float = 3.,
                      npert: int = 5):
        self.snr = snr
        self.sigma_specres_kmps = sigma_specres_kmps
        self.sigma_det = sigma_det
        self.npert = npert
        
        self.noisyfits = []
        self.noisyspectra = []
        for _ in npert:
            self.target_cur = self._noisyspec(
                self.spec_raw, sigma_kmps=self.sigma_specres_kmps,
                snr=self.snr)
            _ncomp, _logcds_cm2, _bvals_kmps, _vcens_kmps = \
                self._fitspec_single(self.target_cur, self.sigma_specres_kmps)
            self.noisyfits.append({'ncomp': _ncomp,
                                   'logN_cm2': np.copy(_logcds_cm2),
                                   'bvals_kmps': np.copy(_bvals_kmps),
                                   'vcens_kmps': np.copy(_vcens_kmps)})
            self.noisyspectra.append(self.target_cur.copy())
        ncomp_noisy = np.array([fit['ncomp'] for fit in self.noisyfits])
        self.ncomp_final = int(np.ceil((np.median(ncomp_noisy))))
        if self.ncomp_final == 0:
            self.logN_cm2_final = np.zeros(shape=(0,), dtype=np.float64)
            self.bvals_kmps_final = np.zeros(shape=(0,), dtype=np.float64)
            self.vcens_kmps_final = np.zeros(shape=(0,), dtype=np.float64)
        else:
            if np.any(ncomp_noisy == self.ncomp_final):
                guessi = np.where(ncomp_noisy == self.ncomp_final)[0][0]
                guessN = self.noisyfits[guessi]['logN_cm2'] 
                guessb = self.noisyfits[guessi]['bvals_kmps']
                guessv = self.noisyfits[guessi]['vcens_kmps']
                self.guess = np.empty((self.ncomp_final, 3), 
                                      dtype=guessN.dtype)
                self.guess[:, 0] = guessN / self.fitunit_logN_cm2
                self.guess[:, 1] = guessb / self.fitunit_b_kmps
                self.guess[:, 2] = guessv / self.fitunit_v_kmps
                self.guess.shape = (self.ncomp_final * 3,)
            else:
                # get a spectrum with the smallest number of 
                # extra components
                opts = np.where(ncomp_noisy > self.ncomp_final)[0]
                subi = np.argmin(ncomp_noisy[opts])
                guessi = opts[subi]
                
                guessN = self.noisyfits[guessi]['logN_cm2']
                guessb = self.noisyfits[guessi]['bvals_kmps']
                guessv = self.noisyfits[guessi]['vcens_kmps']
                # largest-logN components within that fit
                compsel = np.argsort(guessN)[::-1][:self.ncomp_final]
                self.guess = np.empty((self.ncomp_final, 3), 
                                      dtype=guessN.dtype)
                self.guess[:, 0] = guessN[compsel] / self.fitunit_logN_cm2
                self.guess[:, 1] = guessb[compsel] / self.fitunit_b_kmps
                self.guess[:, 2] = guessv[compsel] / self.fitunit_v_kmps
                self.guess.shape = (self.ncomp_final * 3,)
        self.noiseless_target = self._convolve_gauss(
            self.spec_raw, width_kmps=sigma_specres_kmps)
        self.components_final, self.fitspec_final, self.fitobj_final = \
            self._fitspec(self.fitspec_final, self.guess, 
                          sigma_specres_kmps=sigma_specres_kmps)
        self.logN_cm2_final = self.components_final[0::3] \
                                * self.fitunit_logN_cm2
        self.bvals_kmps_final = self.components_final[1::3] \
                                * self.fitu_bvals_kmps
        self.vcens_kmps_final = self.components_final[2::3] \
                                * self.fitu_vcens_kmps

    def savefit(self, h5filen: str, h5grp: str):
        self.logN_cm2_final = self.components_final[0::3] \
                              * self.fitunit_logN_cm2
        self.bvals_kmps_final = self.components_final[1::3] \
                                * self.fitunit_b_kmps
        self.vcens_kmps_final = self.components_final[2::3] \
                                * self.fitunit_v_kmps

        # just save the parameters, can get the spectra from there
        with h5py.File(h5filen, 'a')as f:
            grp = f.create_group(h5grp)
            _grp = grp.create_group('line')
            self.line.save(_grp)
            grp.attrs.create('spectrumfilen', np.string_(self.filen))
            grp.attrs.create('resolution_kmps', self.sigma_specres_kmps)
            grp.attrs.create('snr', self.snr)
            grp.attrs.create('sigma_det', self.sigma_det)
            grp.attrs.create('npert', self.npert)
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

