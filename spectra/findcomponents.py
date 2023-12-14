from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel
import numpy as np
import pandas as pd

def readin_txtspectrum(filen):
    spec = pd.read_csv(filen, comment='#', 
                       columns=['velocity_kmps', 'tau', 'flux', 'flux_error'])
    return spec

def convolve_gauss(spectrum, width_kmps=30.):
    dv = np.average(np.diff(spectrum['velocity_kmps']))
    width_bins = width_kmps / dv
    kernel = Gaussian1DKernel(width_bins)
    smoothedflux = convolve(spectrum['flux'], kernel, 
                            boundary='fill', fill_value=1.)
    spectrum['smoothed_flux'] = smoothedflux
    return spectrum


def findcomponents(spectrum):
    pass