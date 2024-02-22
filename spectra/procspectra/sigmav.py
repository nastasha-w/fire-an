import glob

import h5py
import numpy as np

import fire_an.simlists as sl
import fire_an.spectra.findcomponents as fc
import fire_an.spectra.genspectra as gs
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu

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

def getdata_sigmav(infofilens: list[str], 
                   filepatterns: list[str],
                   simnames: list[str],
                   snapnums: list[str],
                   outname: str,
                   line: fc.Line = fc.ne8_770):
    '''
    print info listed in `lineelts` for the lines of sight in
    `infofilens` and `filepatterns` 
    '''
    
    lineelts = ['simname', 'snapnum', 'los_axis',
                'impactpar_kpc', 'Ntot_logcm2', 'vcen_Nwtd_kmps',
                'sigmav_kmps', 'sigmav_galv_kmps', 'vgal_kmps',
                'Mstar_logMsun', 'Mvir_logMsun',
                'Rvir_kpc']
    eltfmt = {'simname': '{simname}',
               'snapnum': '{snapnum}',
               'los_axis': '{los_axis}',
               'impactpar_kpc': '{impactpar_kpc:.2f}',
               'Ntot_logcm2': '{Ntot_logcm2:.3f}',
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

                vcen_kmps, sigmav_kmps, sigma_gal, logN_cm2 = \
                    calc_vcen_sigmav_logN(filen, vgal_kmps, line=line)

                l2dct = l1dct.copy()
                l2dct.update({'impactpar_kpc': ipar_kpc,
                              'Ntot_logcm2': logN_cm2,
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
                   outname, line=line)
        


