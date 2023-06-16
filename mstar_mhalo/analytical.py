'''
various analytical Mstar-Mhalo relations
'''


def mstar_moster_etal_2013(mvir_msun, redshift):
    # uses M200c for the Mvir definition
    # fits as a function of redshift, with uncertainty in M*
    fp_m10 = 11.590
    fp_m11 = 1.195
    fp_n10 = 0.0351
    fp_n11 = -0.0247
    fp_beta10 = 1.376
    fp_beta11 = -0.826
    fp_gamma10 = 0.608
    fp_gamma11 = 0.329
    fitpar_m1_msun = 10**(fp_m10 + fp_m11 * redshift / (1. + redshift))
    fitpar_n = fp_n10 + fp_n11 * redshift / (1. + redshift)
    fitpar_beta = fp_beta10 + fp_beta11  * redshift / (1. + redshift)
    fitpar_gamma = fp_gamma10 + fp_gamma11 * redshift / (1. + redshift)
    out = mvir_msun * 2. * fitpar_n \
          / ((mvir_msun / fitpar_m1_msun)**-fitpar_beta \
             + (mvir_msun / fitpar_m1_msun)**fitpar_gamma)
    return out


# numerically inverted to Mstar -> Mhalo in the paper
def mstar_burchett_etal_2019(mvir_msun, redshift):
    # modified Moster et al. (2013)
    fp_m10 = 11.590
    fp_m11 = 1.195
    fp_n10 = 0.0351
    fp_n11 = -0.0247
    fp_beta10 = 1.376
    fp_beta11 = -0.826
    fitpar_m1_msun = 10**(fp_m10 + fp_m11 * redshift / (1. + redshift))
    fitpar_n = fp_n10 + fp_n11 * redshift / (1. + redshift)
    fitpar_beta = fp_beta10 + fp_beta11  * redshift / (1. + redshift)
    out = mvir_msun * 2. * fitpar_n \
          / ((mvir_msun / fitpar_m1_msun)**-fitpar_beta + 1.)
    return out