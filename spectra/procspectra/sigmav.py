import glob

import h5py
import numpy as np

import fire_an.simlists as sl
import fire_an.spectra.findcomponents as fc
import fire_an.spectra.genspectra as gs
import fire_an.spectra.procspectra.rough_detmasks as rdm
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu

## CUBS (Qu et al., 2024): FWHM spectral resolution is ~20 km/s
sigma_lsf_cubs_kmps = 20. / (2. * np.sqrt(2. * np.log(2.)))
## from Prochaska, Burchett, et al. (2019, CASBaH galaxy survey):
# Integrating the redshift measurements from Hectospec, DEIMOS, 
# and the SDSS database, we performed internal comparisons between 
# the ∼175 sources common to two or more of the sub-surveys. 
# Ignoring catostrophic failures (described below), the measured 
# rms values between Hectosec/SDSS and DEIMOS/Hectospec are 
# ~=35 and ~=36 km/s respectively. 
# Therefore, we advise adopting a minimum redshift uncertainty of
# 35km/s for galaxies drawn from the CASBaH database.

def calc_vcen_sigmav_logN(specfilen: str, 
                          vcengal_kmps: float, 
                          line: fc.Line = fc.ne8_770):
    '''
    Returns:
    --------
    vcen: float
        optical-depth-weighted line of sight velocity 
        (km/s, trident frame)
    sigma: float
        optical-depth-weighted second moment of the line of sight
        velocity (km/s)
    logN: float
        total column density along the line of sight (log10 cm**-2)
    '''
    spec = fc.SpectrumFitBayes(line, filen=specfilen)
    vel = spec.vel_kmps
    weight = spec.tau_raw
    vcen = np.sum(weight * vel) / np.sum(weight)
    dv = vel - vcen
    sigma = np.sqrt(np.sum(weight * dv**2) / np.sum(weight))
    dvgal = vel - vcengal_kmps
    sigma_gal = np.sqrt(np.sum(weight * dvgal**2) / np.sum(weight))
    logN = np.log10(spec.line.tau_to_coldens(weight, vel))
    return vcen, sigma, sigma_gal, logN

def calc_vcen_sigmav_logN_roughcomponents(
    specfilen: str, 
    vcengal_kmps: float, 
    line: fc.Line = fc.ne8_770,
    mindetcomp_logN_cm2: float = 14.,
    sigma_lsf_kmps: float = 30.,
    maxbpar_kmps: float = 100.,
    mincutfrac: float = 0.05,
) -> tuple[float, float, float, float, float]:
    '''
    Returns:
    --------
    vcen: float
        optical-depth-weighted line of sight velocity 
        (km/s, trident frame)
    sigma: float
        optical-depth-weighted second moment of the line of sight
        velocity (km/s)
    logN: float
        total column density along the line of sight (log10 cm**-2)
    logN_det: float
        total column density along the line of sight in the detected
        components (log10 cm**-2)
    '''
    spec = fc.SpectrumFitBayes(line, filen=specfilen)
    mask, enclinds_det, cds = \
        rdm.get_roughdetmask(spec, lsf_sigma_kmps=sigma_lsf_kmps,
                             coldetlim_logN_cm2=mindetcomp_logN_cm2,
                             maxbpar_kmps=maxbpar_kmps,
                             mincutfrac=mincutfrac)
    logN_det = np.log10(np.sum(cds))
    vel = spec.vel_kmps
    weight = spec.tau_raw
    vcen = np.sum(weight[mask] * vel[mask]) / np.sum(weight[mask])
    dv = vel - vcen
    sigma = np.sqrt(np.sum(weight[mask] * dv[mask]**2) / np.sum(weight[mask]))
    dvgal = vel - vcengal_kmps
    sigma_gal = np.sqrt(np.sum(weight[mask] * dvgal[mask]**2)
                        / np.sum(weight[mask]))
    logN = np.log10(spec.line.tau_to_coldens(weight, vel))
    return vcen, sigma, sigma_gal, logN, logN_det
    
def getdata_sigmav(infofilens: list[str], 
                   filepatterns: list[str],
                   simnames: list[str],
                   snapnums: list[str],
                   outname: str,
                   detectedonly: bool = True,
                   line: fc.Line = fc.ne8_770,
                   mindetcomp_logN_cm2: float = 14.,
                   sigma_lsf_kmps: float = 30.,
                   maxbpar_kmps: float = 100.,
                   mincutfrac: float = 0.05) -> None:
    '''
    print info listed in `lineelts` for the lines of sight in
    `infofilens` and `filepatterns` 
    '''
    
    lineelts = ['simname', 'snapnum', 'los_axis',
                'impactpar_kpc', 'Ntot_logcm2', 'Ndet_logcm2',
                'vcen_Nwtd_kmps',
                'sigmav_kmps', 'sigmav_galv_kmps', 'vgal_kmps',
                'Mstar_logMsun', 'Mvir_logMsun',
                'Rvir_kpc']
    eltfmt = {'simname': '{simname}',
               'snapnum': '{snapnum}',
               'los_axis': '{los_axis}',
               'impactpar_kpc': '{impactpar_kpc:.2f}',
               'Ntot_logcm2': '{Ntot_logcm2:.3f}',
               'Ndet_logcm2': '{Ndet_logcm2:.3f}',
               'vcen_Nwtd_kmps': '{vcen_Nwtd_kmps:.2f}',
               'sigmav_kmps': '{sigmav_kmps:.2f}',
               'sigmav_galv_kmps': '{sigmav_galv_kmps:.2f}',
               'vgal_kmps': '{vgal_kmps:.2f}',
               'Mstar_logMsun': '{Mstar_logMsun:.2f}',
               'Mvir_logMsun': '{Mvir_logMsun:.2f}',
               'Rvir_kpc': '{Rvir_kpc:.2f}'}
    with open(outname, 'a') as fo:
        header = '\t'.join(lineelts) + '\n'
        fo.write(header)
    linefmt = '\t'.join(eltfmt[elt] for elt in eltfmt) + '\n'

    for infofilen, filepattern, simname, snapnum in \
            zip(infofilens, filepatterns, simnames, snapnums):
        cgpath = 'Header/cengal/halodata_doc_dict'
        with h5py.File(infofilen, 'r') as f:
            cosmopars = {key: val for key, val in 
                         f[cgpath]['cosmopars_dict'].attrs.items()}
            Mvir_Msun = f['Header/halo_data'].attrs['Mvir_g'] \
                        / c.solar_mass
            Rvir_kpc = f['Header/halo_data'].attrs['Rvir_cm'] \
                        / (c.cm_per_mpc * 1e-3)
            pgal_cm = f['Header/cengal'].attrs['pcen_cm']
            vcom_cmps = f['Header/cengal'].attrs['vcom_cmps']
            Mstar_Msun = f['Header/cengal'].attrs['mstar_gal_g'] \
                         / c.solar_mass
            hpar = cu.Hubble(cosmopars)
            axis = (f['Header/sample'].attrs['axis']).decode()
            xaxi, yaxi, losaxi = gs.getinds_ax(axis)
            starts = f['startpos_cm'][:]
            ends = f['endpos_cm'][:]
            if not (np.allclose(starts[:, losaxi], starts[0, losaxi]) 
                    and np.allclose(ends[:, losaxi], ends[0, losaxi])):
                print('Sightlines start and end at different positions,'
                      ' so they will not be on a common velocity grid.')
                vgal_kmps = 0.
            p0 = starts[0, losaxi]
            poff = pgal_cm[losaxi] - p0
            zgal = (-1. * vcom_cmps[losaxi] / c.c + 1.) \
                   * (1. - poff * hpar / c.c) - 1.
            # for vbins, snapshot redshift is factored out
                   #* (1. + cosmopars['z']) - 1
            vgal_kmps = zgal * c.c * 1e-5
            print(zgal)
            print(vgal_kmps)

            ipars = (starts[:, xaxi] - pgal_cm[xaxi])**2 \
                    + (starts[:, yaxi] - pgal_cm[yaxi])**2
            ipars_kpc = np.sqrt(ipars) / (c.cm_per_mpc * 1e-3)
        l1dct = {'simname': simname,
                 'snapnum': snapnum,
                 'los_axis': axis,
                 'vgal_kmps': vgal_kmps,
                 'Mvir_logMsun': np.log10(Mvir_Msun),
                 'Rvir_kpc': Rvir_kpc,
                 'Mstar_logMsun': np.log10(Mstar_Msun)}
        filens = glob.glob(filepattern)
        with open(outname, 'a') as fo:
            for filen in filens:
                #print(filen)
                sli = filen.split('.')[-2]
                sli = sli.split('_')[-1]
                sli = int(sli)
                ipar_kpc = ipars_kpc[sli]
                if detectedonly:
                    vcen_kmps, sigmav_kmps, sigma_gal, logN_cm2, logNdet_cm2 = \
                        calc_vcen_sigmav_logN_roughcomponents(
                            filen,
                            vgal_kmps,
                            line=line,
                            mindetcomp_logN_cm2=mindetcomp_logN_cm2,
                            sigma_lsf_kmps=sigma_lsf_kmps,
                            maxbpar_kmps=maxbpar_kmps,
                            mincutfrac=mincutfrac,
                            )
                else:
                    vcen_kmps, sigmav_kmps, sigma_gal, logN_cm2 = \
                        calc_vcen_sigmav_logN(filen, vgal_kmps, line=line)
                    logNdet_cm2 = logN_cm2

                l2dct = l1dct.copy()
                l2dct.update({'impactpar_kpc': ipar_kpc,
                              'Ntot_logcm2': logN_cm2,
                              'Ndet_logcm2': logNdet_cm2,
                              'vcen_Nwtd_kmps': vcen_kmps,
                              'sigmav_kmps': sigmav_kmps,
                              'sigmav_galv_kmps': sigma_gal})

                pline = linefmt.format(**l2dct)
                fo.write(pline)

def getsigmav_set(testset=4):
    if testset == 4:
        simnames = sl.m12_f2md
        snapshots = [sl.snaps_f2md[0], sl.snaps_f2md[1]]
        filedir = '/projects/b1026/nastasha/spectra/test4/'
        outdir =  '/projects/b1026/nastasha/imgs/spectra/test4/'
        # the one that worked the first time
        #simnames = ['crheatfix_m12i_r7100', 'crheatfix_m12f_r7100']
        #snapshots = [294, 277]
        outname = outdir + f'sigmav_testset4.dat'
        filepatterns = [filedir + f'/tridentray_{simname}_{snapnum}_*.txt'
                       for simname in simnames for snapnum in snapshots]
        infofiles = [filedir + f'/tridentray_{simname}_{snapnum}_info.hdf5'
                       for simname in simnames for snapnum in snapshots]
        simnames_arg = [simname 
                        for simname in simnames for snapnum in snapshots]
        snapnums_arg = [snapnum 
                        for simname in simnames for snapnum in snapshots]
        line = fc.ne8_770
        getdata_sigmav(infofiles, filepatterns, simnames_arg, snapnums_arg,
                       outname, line=line, detectedonly=False)
    elif testset == '4a': 
        ## vary some of the detection parameters, see effect
        simnames = sl.m12_f2md
        snapshots = [sl.snaps_f2md[0], sl.snaps_f2md[1]]
        filedir = '/projects/b1026/nastasha/spectra/test4/'
        outdir =  '/projects/b1026/nastasha/imgs/spectra/test4/'
        # the one that worked the first time
        #simnames = ['crheatfix_m12i_r7100', 'crheatfix_m12f_r7100']
        #snapshots = [294, 277]
        filepatterns = [filedir + f'/tridentray_{simname}_{snapnum}_*.txt'
                       for simname in simnames for snapnum in snapshots]
        infofiles = [filedir + f'/tridentray_{simname}_{snapnum}_info.hdf5'
                       for simname in simnames for snapnum in snapshots]
        simnames_arg = [simname 
                        for simname in simnames for snapnum in snapshots]
        snapnums_arg = [snapnum 
                        for simname in simnames for snapnum in snapshots]
        line = fc.ne8_770
        # sigma LSF: 20 km/s FWHM for the CUBS spectra 
        # -> ~10 km/s sigma (8.5 km/s)
        detkw = [{'mindetcomp_logN_cm2': 14., 'sigma_lsf_kmps': 10.,
                  'maxbpar_kmps': 100., 'mincutfrac': 0.05},
                 {'mindetcomp_logN_cm2': 13.75, 'sigma_lsf_kmps': 10.,
                  'maxbpar_kmps': 100., 'mincutfrac': 0.05},
                 {'mindetcomp_logN_cm2': 13.5, 'sigma_lsf_kmps': 10.,
                  'maxbpar_kmps': 100., 'mincutfrac': 0.05},
                 {'mindetcomp_logN_cm2': 14., 'sigma_lsf_kmps': 10.,
                  'maxbpar_kmps': 100., 'mincutfrac': 0.01},
                 {'mindetcomp_logN_cm2': 14., 'sigma_lsf_kmps': 10.,
                  'maxbpar_kmps': 100., 'mincutfrac': 0.25},
                 {'mindetcomp_logN_cm2': 13.5, 'sigma_lsf_kmps': 10.,
                  'maxbpar_kmps': 100., 'mincutfrac': 0.01},
                 {'mindetcomp_logN_cm2': 13.5, 'sigma_lsf_kmps': 10.,
                  'maxbpar_kmps': 100., 'mincutfrac': 0.25},
                  ]
        for kwargs in detkw:
            outname = (outdir + 'sigmav_testset4_'
                       f'minN{kwargs["mindetcomp_logN_cm2"]:.2f}_'
                       f'lsf{kwargs["sigma_lsf_kmps"]:.1f}_'
                       f'maxb{kwargs["maxbpar_kmps"]:.0f}_'
                       f'cutf{kwargs["mincutfrac"]:.2f}.dat')
            getdata_sigmav(infofiles, filepatterns, simnames_arg,
                           snapnums_arg, outname, detectedonly=True, 
                           line=line, **kwargs)
        


