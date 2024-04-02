import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import fire_an.analytic_halo.model_ionprof_pl as mip
import fire_an.ionrad.ion_utils as iu
import fire_an.makeplots.tol_colors as tc
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu

# fire_an.makeplots.litcomp.obsdataread for calculations
oddir = '/projects/b1026/nastasha/extdata/' # quest
#oddir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/' #laptop 
q23filen = oddir + 'plotdata_q23_nsigmas_1_2.dat'
b19filen = oddir + 'plotdata_b19_nsigmas_1_2.dat'
g05filen = oddir + 'gallazzi_etal_2005_table2_stellarMZ.txt'
g05_Zsun = 0.0127 # I couldn't actually find a value in the paper

mdir = '/projects/b1026/nastasha/imgs/analytical/ion_origins/'

# some estimates of SSP yields, 
# given a set of single-star yields and an IMF
# https://ui.adsabs.harvard.edu/abs/2016MNRAS.455.4183V/abstract
# pretty wide ranges here: ~0.01 -- 0.04 (at t=infinity) depending
# on IMF and stellar yield data 
metal_yield = 0.02 

## Notes on estimating stellar metallicities
# In general, z~0.5 -- 1.0 seems to be pretty rare for Z* data (mostly 
# Fe?). There are different simulation/analytical model predictions.
# Generally, around or a bit below solar seems to be typical
# at 10^10 -- 10^10.5 Msun, 
# which is where most of the observed conterparts to detected Ne VIII
# absorption land for CUBS, and roughly at the center of the CASBaH 
# range
#
# Kashino et al. (2022)
# for stellar mass MZR at z=1.6-3, and z~0 refs:
# https://iopscience.iop.org/article/10.3847/1538-4357/ac399e/pdf
# [Fe/H] = - (0.81 +- 0.01) + (0.32 ?+-? 0.03) * log10(M* / 10^10 Msun) 
# using solar metallicity values of  12 + log(O/H) = 8.69 
#    and Zsun = 0.0142 (Asplund et al. 2009)
#    Z_{O, sun} = 0.00561 and Z_{Fe, sun} = 0.00126
# show Z* ~ Zsun at M* = 10^10 -- 10^10.5 Msun at z=0 (fig. 11)  
#
# one paper that shows z=0.8, I think:
# https://openaccess.inaf.it/bitstream/20.500.12386/31805/1/2104.08295.pdf
# GAEA analytical model predictions. Do show z=0.7.
#
# Ma et al. (2016) FIRE-2 fits:
# https://ui.adsabs.harvard.edu/abs/2016MNRAS.456.2140M/abstract
# for Zsun = 0.02, solar Fe mass frac = 0.00173
# [Fe/H] = log10(Z*/Zsun) − 0.20
# log(Z*/Zsun) = 0.40 * [log(M*/Msun) − 10] + 0.67 * exp(−0.50 * z) − 1.04
# from fig. 2: Gallazzi et al. (2005, observations): 
#     median ~0.3 -- 1 Zsun at z=0, M* = 10^10 -- 10^9.5 Msun
#     2x 1-sigma dispersion (1-sigma range) of ~0.5 -- 1 dex
#
# De Rossi et al. (2017) EAGLE fits:
# https://ui.adsabs.harvard.edu/abs/2017MNRAS.472.3354D/abstract
# using Zsun = 0.0127, ?only central galaxies?
# M* = 10^10 -- 10^10.5 Msun: a bit below to a bit above Zsun at z=1
#     ~Zsun to a bit above Zsun at z=0.5
#     show percentile 25 -- 75 ranges, seen to be 0.1 -- 0.2 dex  
#
# Leethochawalit et al. (2018) obs around a z=0.4 cluster
# https://ui.adsabs.harvard.edu/abs/2018ApJ...856...15L/abstract
# find redshift evolution of M*-Z* between z=0 and z=0.4
# but this is specifically for quiescent galaxies
# [Fe/H] ~ -0.5 -- 0 at M* = 10^10 -- 10^10.5 Msun, z=0.4
#
# Gallazzi et al. (2005): z=0, but plenty of data
# https://ui.adsabs.harvard.edu/abs/2005MNRAS.362...41G/abstract
# 
starZ = 0.014 # TODO: look up actual stellar metallicity
g05dat = pd.read_csv(g05filen, sep='\t', comment='#')
def get_starZ_range(logmstar_msun):
    lovals = np.interp(logmstar_msun,
                       g05dat['Mstar_logMsun'],
                       g05dat['Z_logZsun_p16'])
    midvals = np.interp(logmstar_msun,
                       g05dat['Mstar_logMsun'],
                       g05dat['Z_logZsun_p50'])
    hivals = np.interp(logmstar_msun,
                       g05dat['Mstar_logMsun'],
                       g05dat['Z_logZsun_p84'])
    out = 10**np.array([lovals, midvals, hivals]).T * g05_Zsun
    return out


def check_consistency_cieplmodel(infodict: dict,
                                 plis_vc: tuple = (-0.1,),
                                 plis_entropy: tuple = (1.,)):
    '''
    Can we explain the measured column densities with a hot
    phase power law model, and an amount of metals consistent
    with what the central galaxy could have produced?
    '''
    mstar_msun_opts = 10**infodict['mstar_opts_logMsun']
    totmetals = mstar_msun_opts * metal_yield
    _starZ = get_starZ_range(infodict['mstar_opts_logMsun'])
    starmetals = mstar_msun_opts[:, np.newaxis] * _starZ

    mvir_msun_opts = 10**infodict['mvir_opts_logMsun']
    redshift = infodict['z']
    impactpar_kpc = infodict['ipar_kpc']
    coldens_meas = 10**infodict['cd_opts_logcm2'][2]
    
    fcgm = 1.0
    nomZ_solar = 0.3
    cds = []
    massZ_nom = []
    for pli_vc in plis_vc:
        for pli_entropy in plis_entropy:
            for mvir_msun in mvir_msun_opts:
                model = mip.PLmodel(mvir_msun, redshift, fcgm, nomZ_solar, 
                                    pli_vc, pli_entropy=pli_entropy)
                cd = model.coldensprof('Ne8', np.array([impactpar_kpc]), 
                                       loslen_rvir=4.)
                cds.append(cd)
                massZ_nom.append(mvir_msun * nomZ_solar  
                                 * model._tab.solarZ
                                 * model.cosmopars['omegab']
                                 / model.cosmopars['omegam'])
    cds = np.array(cds)
    massZ_nom = np.array(massZ_nom)
    massZ_inferred = massZ_nom * coldens_meas / cds 

    return massZ_inferred, starmetals, totmetals


def getinfo_1sys(df: pd.DataFrame,
                 index,
                 survey: {'cubs', 'casbah'} = 'cubs'):
    '''
    For one row of the obs. dataframe, get the info we need for the
    consistency tests.
    '''
    out = {}
    if survey == 'casbah':
        mstar_mid = df.at[index, 'log_Mstar_Msun']
        dmstar = df.at[index, 'log_Mstar_Msun_err']
        out['mstar_opts_logMsun'] = np.array([mstar_mid - 2. * dmstar,
                                              mstar_mid - dmstar,
                                              mstar_mid,
                                              mstar_mid + dmstar,
                                              mstar_mid + 2. * dmstar])
        cd_mid = df.at[index, 'log_N_Ne8_pcm2']
        dcd = 0.5 * df.at[index, 'log_N_Ne8_pcm2_err']
        out['cd_opts_logcm2'] = np.array([cd_mid - 2. * dcd,
                                          cd_mid - dcd,
                                          cd_mid,
                                          cd_mid + dcd,
                                          cd_mid + 2. * dcd])
        out['ipar_kpc'] = df.at[index, 'impact_parameter_kpc']
        out['z'] = df.at[index, 'zgal']

    elif survey == 'cubs':
        calscatter = 0.2
        mstar_mid = df.at[index, 'mstar_log10Msun']
        dmstar_lo = df.at[index, 'mstar_2s_loerr_dex']
        dmstar_hi = df.at[index, 'mstar_2s_hierr_dex']
        # get total, 1-sigma errors
        dmstar_lo = np.sqrt(0.25 * dmstar_lo**2 + calscatter**2)
        dmstar_hi = np.sqrt(0.25 * dmstar_hi**2 + calscatter**2)
        out['mstar_opts_logMsun'] = np.array([mstar_mid - 2. * dmstar_lo,
                                              mstar_mid - dmstar_lo,
                                              mstar_mid,
                                              mstar_mid + dmstar_hi,
                                              mstar_mid + 2. * dmstar_hi])
        cd_mid = df.at[index, 'ne8col_logcm2']
        dcd_lo = 0.5 * df.at[index, 'ne8col_2s_loerr_dex']
        dcd_hi = 0.5 * df.at[index, 'ne8col_2s_hierr_dex']
        out['cd_opts_logcm2'] = np.array([cd_mid - 2. * dcd_lo,
                                          cd_mid - dcd_lo,
                                          cd_mid,
                                          cd_mid + dcd_hi,
                                          cd_mid + 2. * dcd_hi])
        out['ipar_kpc'] = df.at[index, 'impactpar_kpc']
        out['z'] = df.at[index, 'z_gal']
    
    out['mvir_opts_logMsun'] = \
        np.array([df.at[index, 'logmvir_msun_loer'],
                  df.at[index, 'logmvir_msun_lo'],
                  df.at[index, 'logmvir_msun_bestest'],
                  df.at[index, 'logmvir_msun_hi'],
                  df.at[index, 'logmvir_msun_hier']
                  ])
    return out


def runcheck_enough_metals_hotphase():
    # print everything at the end to get a nice overview
    # (model calculation print a lot of stuff)
    res = ''
    res = res + 'CASBaH data:\n'
    res = res + '------------\n'
    df = pd.read_csv(b19filen, sep='\t')
    df = df[np.logical_not(df['log_N_Ne8_isUL'])]
    for i in df.index:
        info = getinfo_1sys(df, i, survey='casbah')
        res = res + '\n'
        massZ_inferred, starmetals, totmetals = \
            check_consistency_cieplmodel(info)
        res = res + ('for abs. sys. log10 N = '
                     f'{info["cd_opts_logcm2"][2]}'
                     f', at {info["ipar_kpc"]} kpc\n'
                     f'inferred total Z mass: {massZ_inferred} Msun\n'
                     f'\trange: {np.min(massZ_inferred):.2e} -- '
                     f'{np.max(massZ_inferred):.2e} Msun,\n'
                     f'\tmedian: {np.median(massZ_inferred):.2e} Msun\n'
                     f'inferred stellar Z mass: {starmetals} Msun\n'
                     f'inferred produced Z mass: {totmetals} Msun\n')

    res = res + '\n\n'
    res = res + 'CUBS data:\n'
    res = res + '----------\n'
    df = pd.read_csv(q23filen, sep='\t')
    df = df[np.logical_not(df['isul_ne8'])]
    for i in df.index:
        info = getinfo_1sys(df, i, survey='cubs')
        res = res + '\n'
        massZ_inferred, starmetals, totmetals = \
            check_consistency_cieplmodel(info)
        res = res + ('for abs. sys. log10 N = '
                     f'{info["cd_opts_logcm2"][2]}'
                     f', at {info["ipar_kpc"]} kpc\n'
                     f'inferred total Z mass: {massZ_inferred} Msun\n'
                     f'\trange: {np.min(massZ_inferred):.2e} -- '
                     f'{np.max(massZ_inferred):.2e} Msun,\n'
                     f'\tmedian: {np.median(massZ_inferred):.2e} Msun\n'
                     f'inferred stellar Z mass: {starmetals} Msun\n'
                     f'inferred produced Z mass: {totmetals} Msun\n')
    print(res)    


def plotcheck_enough_metals_hotphase(info: dict,
                                     plis_vc: tuple = (0.0, -0.15, -0.3),
                                     plis_entropy: tuple = (0.67, 1., 1.2),
                                     survey='CASBaH'):
    massZ_inferred, starmetals, totmetals = \
            check_consistency_cieplmodel(info, plis_vc=plis_vc,
                                         plis_entropy=plis_entropy)
    fig = plt.figure(figsize=(5.5, 5.))
    ax = fig.add_subplot()
    colors = tc.tol_cset('vibrant')

    lMZ_cgm = np.log10(massZ_inferred.flatten())
    lMZ_stars = np.log10(starmetals.flatten())
    lMZ_prod = np.log10(totmetals.flatten())
    xmin = min([np.min(lMZ_cgm), np.min(lMZ_stars), np.min(lMZ_prod)])
    xmax = max([np.max(lMZ_cgm), np.max(lMZ_stars), np.max(lMZ_prod)])
    bins = np.linspace(0.95 * xmin, 1.05 * xmax, 25)
    ax.hist(lMZ_cgm, bins=bins, density=True, color=colors[0],
            label='CGM', weights=[1. / len(lMZ_cgm)] * len(lMZ_cgm),
            histtype='step', linewidth=2, linestyle='solid')
    ax.hist(lMZ_stars, bins=bins, density=True, color=colors[1],
            label='stars',  weights=[1. / len(lMZ_stars)] * len(lMZ_stars),
            histtype='step', linewidth=2, linestyle='dashed')
    ax.hist(lMZ_prod, bins=bins, density=True, color=colors[2],
            label='prod.',  weights=[1. / len(lMZ_prod)] * len(lMZ_prod),
            histtype='step', linewidth=2, linestyle='dotted')
    ax.legend(fontsize=10)
    ax.set_xlabel('$\\log_{10} \\, \\mathrm{M}_{\\mathrm{Z}} '
                  '\\; [\\mathrm{M}_{\\odot}]$',
                  fontsize=12)
    ax.set_ylabel('pdf', fontsize=12)
    ax.tick_params(which='both', labelsize=11, direction='in',
                   top=True, right=True)
    title = ('$\\log_{10} \\, \\mathrm{N} = '
             f'{info["cd_opts_logcm2"][2]:.1f}'
             ', \\mathrm{r}_{\\perp} = '
             f'{info["ipar_kpc"]:.0f}$ kpc')
    title = survey + ', ' + title
    fig.suptitle(title, fontsize=12)

    outname = mdir + (f'hotphase_Zbudget_{survey}_logN_'
                      f'{info["cd_opts_logcm2"][2]:.1f}'
                      f'ipar_{info["ipar_kpc"]:.0f}.pdf')
    plt.savefig(outname, bbox_inches='tight')

def runset_enough_metals_hotphase():
    # print everything at the end to get a nice overview
    # (model calculation print a lot of stuff)
    df = pd.read_csv(b19filen, sep='\t')
    df = df[np.logical_not(df['log_N_Ne8_isUL'])]
    for i in df.index:
        info = getinfo_1sys(df, i, survey='casbah')
        plotcheck_enough_metals_hotphase(info, survey='CASBaH')
    df = pd.read_csv(q23filen, sep='\t')
    df = df[np.logical_not(df['isul_ne8'])]
    for i in df.index:
        info = getinfo_1sys(df, i, survey='cubs')
        plotcheck_enough_metals_hotphase(info, survey='CUBS')
 
    
def check_pathlengths_pie(info: dict,
                          logTKrange: tuple[float, float] = (3.5, 4.5),
                          survey='cubs'):
    z = info['z']
    ion = 'Ne8'
    assumedZ_solar = 0.3
    if survey == 'cubs':
        cosmopars = {'h': 0.7, 'omegam': 0.3, 'omegalambda': 0.7}
    elif survey == 'casbah':
        cosmopars = {'h': 0.677, 'omegam': 0.31, 'omegalambda': 0.69}
    cosmopars.update({'z': z, 'a': 1. / (1. + z)})

    itab = iu.Linetable_PS20(ion, z, lintable=True)
    itab.findiontable()
    solarZ = itab.solarZ
    tableZ = itab.logZsol
    iZ = np.argmin(np.abs(tableZ - assumedZ_solar))
    ionf = itab.iontable_T_Z_nH[:, iZ, :]
    tablenH = itab.lognHcm3
    tableT = itab.logTK
    selT = tableT <= logTKrange[1]
    selT &= tableT >= logTKrange[0]
    dctZ = {'logZ': np.array([np.log10(assumedZ_solar * solarZ)])}
    nitonh = itab.find_assumedabundance(dctZ, log=False)
    iondens = ionf[selT, :] * 10**tablenH[np.newaxis, :] * nitonh
    maxid = np.max(iondens)

    pathlengths_kpc = 10**info['cd_opts_logcm2'] / maxid \
                      / (c.cm_per_mpc * 1e-3)
    rvir_cm = cu.rvir_from_mvir(10**info['mvir_opts_logMsun'] * c.solar_mass,
                                cosmopars=cosmopars, meandef='BN98')
    rvir_kpc = rvir_cm / (c.cm_per_mpc * 1e-3)

    return pathlengths_kpc, rvir_kpc

def runcheck_pathlengths_pie():
    # print everything at the end to get a nice overview
    # (model calculation print a lot of stuff)
    res = ''
    res = res + 'CASBaH data:\n'
    res = res + '------------\n'
    df = pd.read_csv(b19filen, sep='\t')
    df = df[np.logical_not(df['log_N_Ne8_isUL'])]
    for i in df.index:
        info = getinfo_1sys(df, i, survey='casbah')
        res = res + '\n'
        pathlengths_kpc, rvir_kpc = \
            check_pathlengths_pie(info, survey='casbah')
        res = res + ('for abs. sys. log10 N = '
                     f'{info["cd_opts_logcm2"][2]}'
                     f', at {info["ipar_kpc"]} kpc\n'
                     f'inferred min. pathlength: {pathlengths_kpc} kpc\n'
                     f'inferred virial radius: {rvir_kpc} kpc\n')
    res = res + '\n\n'
    res = res + 'CUBS data:\n'
    res = res + '----------\n'
    df = pd.read_csv(q23filen, sep='\t')
    df = df[np.logical_not(df['isul_ne8'])]
    for i in df.index:
        info = getinfo_1sys(df, i, survey='cubs')
        res = res + '\n'
        pathlengths_kpc, rvir_kpc = \
            check_pathlengths_pie(info, survey='cubs')
        res = res + ('for abs. sys. log10 N = '
                     f'{info["cd_opts_logcm2"][2]}'
                     f', at {info["ipar_kpc"]} kpc\n'
                     f'inferred min. pathlength: {pathlengths_kpc} kpc\n'
                     f'inferred virial radius: {rvir_kpc} kpc\n')
    print(res)    