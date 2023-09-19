import astropy.io.fits as apfits
import numpy as np

## From Zhijie:
#Two notes to clarify the table:
# 1) For each measurements, there are three columns representing 
#    the value, 1sigma lower uncertainty, and upper uncertainty. 
#    In particular, a lower uncertainty of -1 means it is an upper 
#    limit in observation, while an upper uncertainty of -1 is for 
#    lower limits. These limits are all 2 sigma. If both 
#    uncertainties are -1, this value is unconstrained. 
# 2) The velocity of each ion in the fits file is relative to the 
#    column of “zfit”, instead of relative to the redshift of the 
#    associated galaxy.


cubsdf = ('/Users/nastasha/ciera/projects_lead/fire3_ionabs/'
          'CUBSVII_qu_etal_2022_draft_table2.fits')
def readin_cubsdata():
    with apfits.open(cubsdf, 'readonly') as hdul:
        data = hdul[1].data
        # TODO: check MS_tot
        mstar = data['Ms_rmin'] # log10 Msun
        mstar_2s_loerr = data['Ms_rmin0'] 
        mstar_2s_hierr = data['Ms_rmin1'] 
        impactpar = data['dproj_rmin'] # kpc
        ne8col = data['NeVIII'] # log10 cm**-2
        ne8col_2s_loerr = data['NeVIII0']
        ne8col_2s_hierr = data['NeVIII1']

        ne8col[ne8col == -1.] = np.NaN # missing values
        isul_ne8 = ne8col_2s_loerr == -1.
        isll_ne8 = ne8col_2s_hierr == -1.
        if np.any(isll_ne8): # I don't think there are any, but check.
            print('Warning: lower limits on Ne VIII col. in data.')
        



