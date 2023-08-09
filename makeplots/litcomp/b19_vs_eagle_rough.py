'''
Rough comparison of the Burchett et al. (2019) data
to median column density curves from EAGLE (z=0.5)
in my 2020 paper with Joop and Ben.
'''

import h5py
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as so

import fire_an.makeplots.litcomp.b19_vs_analytical as bva
import fire_an.utils.cosmo_utils as cu
import fire_an.utils.constants_and_units as c
import fire_an.utils.math_utils as mu

mdir = ''
eagledatadir = ''

cosmopars_ea_23 = {'a': 0.6652884960735025,
                   'boxsize': 67.77,
                   'h': 0.6777,
                   'omegab': 0.0482519,
                   'omegalambda': 0.693,
                   'omegam': 0.307,
                   'z': 0.5031073074342141,
                   }

def readin_eagledata():
    '''
    R200c units, 
    0 - 2.5 R200c impact parameter, 0.1 R200c bin size
    0.5 dex mass bins (M200c)
    '''
    eagledatafilen = ('rdist_coldens_ne8_L0100N1504_23_test3.4_PtAb_C2Sm'
                      '_32000pix_6.25slice_zcen-all_z-projection_T4EOS'
                      '_1slice_to-100-pkpc-or-3-R200c_M200c-0p5dex-7000'
                      '_centrals_stored_profiles.hdf5')
    percdct = {}
    percentiles = [10., 50., 90.]
    pkeys = {_p: f'perc_{_p:.1f}' for _p in percentiles}
    subpath = '/R200c_bins/binset_0/'
    with h5py.File(eagledatadir + eagledatafilen, 'r') as f:
        galsets = list(f.keys())
        for galset in galsets:
            edges = f[galset + subpath + 'bin_edges'][:]
            pvals = {_p: f[galset + subpath + pkeys[_p]][:]
                     for _p in percentiles}
            mmin = f[galset].attrs['seltag'].decode()
            mmin = mmin.split['_'][0]
            mmin = float(mmin[3:])
            percdct['mmin'] = {'edges': edges, 'pvals': pvals}
    return percdct

def mvir_to_m200c(mvir_msun):
    meandens_bn98 = cu.getmeandensity('BN98', cosmopars_ea_23)
    #meandens_200c = cu.getmeandensity('200c', cosmopars_ea_23)
    mvir_g = mvir_msun * c.solar_mass
    rbn98_cgs = (mvir_g / meandens_bn98 * 3. / (4. * np.pi))**(1./3.)
    
    rbins = np.linspace(0., 5. * rbn98_cgs, 5000.)
    rcens = 0.5 * (rbins[:-1] + rbins[1:])
    dr = np.average(np.diff(rbins))
    def minfunc(M200c):
        rhovals = cu.rho_NFW(rcens, M200c, delta=200, ref='rhocrit', 
                             z=cosmopars_ea_23['z'], 
                             cosmopars=cosmopars_ea_23, c='Schaller15')
        menc = np.cumsum(4. * np.pi / 3. * rhovals * rcens**2 * dr)
        rhomean = menc / (4. * np.pi / 3. * rbins[1:]**3)
        rhomean_rbn98 = mu.linterpsolve(rbins[1:], rhomean, rbn98_cgs)
        cost = np.abs(np.log10(rhomean_rbn98) - np.log10(meandens_bn98))
        return cost
    
    res = so.minimize(minfunc, mvir_g)
    if not res.succes:
        print(res)
        return res
    else:
        return res.x



## initial copy was from the paper 2 scripts
def plot_radprof_zev():
    '''
    ions in different panels
    colors indicate different halo masses (from a rainbow color bar)
      
    techvars: 0 = 1 slice, z=0.1, all gas, all haloes
              8 = 1 slice, z=0.1, all gas, all haloes
    '''
    fontsize = 12
    units = 'R200c'
    ytype = 'perc'
    yvals_toplot = [10., 50., 90.]
    ion = 'ne8'
    ylim = (12.0, 15.5)
    imgname = 'ne8_b19_cols_vs_eaglez0p5_w20_roughcomp_v1.pdf'
    imgname = mdir + imgname
    
    shading_alpha = 0.45 
    ylabel = '$\\log_{10} \\, \\mathrm{N} \\; [\\mathrm{cm}^{-2}]$'
    xlabel = '$r_{\\perp} \\; [\\mathrm{R}_{\\mathrm{200c}}]$'
    #clabel =( '$\\log_{10}\, \\mathrm{M}_{\\mathrm{200c}} '
    #         '\\; [\\mathrm{M}_{\\odot}]$')
    
    eagledat = readin_eagledata()
    b19dat = bva.eaddata_b19(nsigmas=(1, 2))
    # define used mass ranges
    deltaM200c = 0.5
    


    #fcovticklen = 0.035
    figwidth = numcols * panelwidth + cwidth + wspace * numcols
    figheight = 2 * numrows * panelheight + legheight
    #print('{}, {}'.format(figwidth, figheight))
    fig = plt.figure(figsize=(figwidth, figheight))
    grid = gsp.GridSpec(2 * numrows + 1, numcols + 1, hspace=0.0,\
                        wspace=wspace,\
                        width_ratios=[panelwidth] * numcols + [cwidth],\
                        height_ratios=[panelheight] * numrows * 2 + [legheight],\
                        bottom=0.0)
    axes = [[fig.add_subplot(grid[(i // numcols) * 2 + j, i % numcols]) for j in range(2)]for i in range(len(ions))]
        
    plt.savefig(imgname, format='pdf', bbox_inches='tight')

