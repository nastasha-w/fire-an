import numpy as np

import fire_an.analytic_halo.model_ionprof_pl as mip
import fire_an.ionrad.ion_utils as iu

# fire_an.makeplots.litcomp.obsdataread for calculations
oddir = '/projects/b1026/nastasha/extdata/'
q23filen = oddir + 'plotdata_q23_nsigmas_1_2.dat'
b19filen = oddir + 'plotdata_b19_nsigmas_1_2.dat'

# some estimates of SSP yields, 
# given a set of single-star yields and an IMF
# https://ui.adsabs.harvard.edu/abs/2016MNRAS.455.4183V/abstract
metal_yield = 0.02 # TODO: look up a literature value
# for stellar mass MZR at z=1.6-3, and z~0 refs:
# https://iopscience.iop.org/article/10.3847/1538-4357/ac399e/pdf
# one paper that shows z=0.8, I think:
# https://openaccess.inaf.it/bitstream/20.500.12386/31805/1/2104.08295.pdf
starZ = 0.014 # TODO: look up actual stellar metallicity


def check_consistency_cieplmodel():
    '''
    stub for actual calculation. Add other parts based on
    read-in values
    '''

    totmetals = mstar * metal_yield
    starmetals = mstar * starZ
    
    plis_vc = [0.0, -0.2, -0.5]
    plis_entropy = [0.67, 1., 1.2]
    fcgm = 1.0
    solarZ = 0.3
    cds = []
    massZ_nom = []
    for pli_vc in plis_vc:
        for pli_entropy in plis_entropy:
            for mvir_msun in mvir_opts:
                model = mip.PLmodel(mvir_msun, redshift, fcgm, solarZ, 
                                    pli_vc,
                                    pli_entropy=pli_entropy)
                cd = model.coldensprof('Ne8', np.array([impactpar_pkpc]), 
                                       loslen_rvir=4.)
                cds.append(cd)
                massZ_nom.append(mvir_msun * massZ_nom 
                                 * model._tab.solarZ
                                 * model.cosmopars['omegab']
                                 / model.cosmopars['omegam'])
    cds = np.array(cds)
    massZ_nom = np.array(massZ_nom)
    massZ_inferred = massZ_nom * coldens_meas / cds 

    print(f'for abs. sys. log10 N = {np.log10(coldens_meas)}'
          f', at {impactpar_kpc} kpc')
    print(f'pli_vc = {plis_vc}, pli_entropy = {plis_entropy}')


