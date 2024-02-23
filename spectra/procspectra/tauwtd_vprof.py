import glob

import h5py
import numpy as np
import scipy.interpolate as spi

import fire_an.simlists as sl
import fire_an.spectra.findcomponents as fc
import fire_an.spectra.genspectra as gs
import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu


def tauwtd_vprof(infofiles: list[str], 
                 filepatterns: list[str],
                 mincoldens: tuple[float] = (13.25, 13.5, 13.75, 14.0),
                 line: fc.Line = fc.ne8_770,
                 outfilen: str | None = None):
    vgrid = np.arange(-1000, 1000.5, 1.)
    minlogN = np.min(mincoldens)
    outdata = {minN: {'tau_av': np.zeros(len(vgrid)),
                      'Nspec': 0,
                      'filesincl.': []}
               for minN in mincoldens}

    for infofile, filepattern in zip(infofiles, filepatterns):
        cgpath = 'Header/cengal/halodata_doc_dict'
        with h5py.File(infofile, 'r') as f:
            cosmopars = {key: val for key, val in 
                         f[cgpath]['cosmopars_dict'].attrs.items()}
            pgal_cm = f['Header/cengal'].attrs['pcen_cm']
            vcom_cmps = f['Header/cengal'].attrs['vcom_cmps']
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
            
        filens = glob.glob(filepattern)
        for filen in filens:
            #sli = filen.split('.')[-2]
            #sli = sli.split('_')[-1]
            #sli = int(sli)
            #ipar = ipars[sli]
            spec = fc.SpectrumFitFreq(line, filen=filen)
            #print(spec.tau_raw)
            #print(spec.vel_kmps)
            logN = np.log10(spec.line.tau_to_coldens(spec.tau_raw, 
                                                     spec.vel_kmps))
            if logN >= minlogN:
                tau_interp = spi.griddata(spec.vel_kmps - vgal_kmps, 
                                          spec.tau_raw, 
                                          vgrid, method='linear', 
                                          fill_value=0., rescale=False)
                for mN in mincoldens:
                    if logN < mN:
                        continue
                    outdata[mN]['tau_av'] += tau_interp
                    outdata[mN]['Nspec'] += 1
                    outdata[mN]['filesincl.'].append(filen)
    for mN in mincoldens:
        outdata[mN]['tau_av'] = outdata[mN]['tau_av'] / outdata[mN]['Nspec']
    with h5py.File(outfilen, 'a') as f:
        for minN in mincoldens:
            grp = f.create_group(f'logN_ge_{minN:.2f}')
            grp.attrs.create('num_spectra', outdata[minN]['Nspec'])
            grp.create_dataset('v_gal_kmps', data=vgrid)
            grp.create_dataset('tau_av', data=outdata[minN]['tau_av'])
            _filens = np.array([np.string_(fn) 
                                for fn in outdata[minN]['filesincl.']])
            grp.create_dataset('included_spectra', data=_filens)

def run_tauwtd_vprof(testset=4):
    if testset == 4:
        simnames = sl.m12_f2md
        snapshots = [sl.snaps_f2md[0], sl.snaps_f2md[1]]
        filedir = '/projects/b1026/nastasha/spectra/test4/'
        outdir =  '/projects/b1026/nastasha/imgs/spectra/test4/'
        # the one that worked the first time
        #simnames = ['crheatfix_m12i_r7100', 'crheatfix_m12f_r7100']
        #snapshots = [294, 277]
        outname = outdir + f'tau_av_ne8_770_minlogN_testset4.pdf'
        filepatterns = [filedir + f'/tridentray_{simname}_{snapnum}_*.txt'
                       for simname in simnames for snapnum in snapshots]
        infofiles = [filedir + f'/tridentray_{simname}_{snapnum}_info.hdf5'
                       for simname in simnames for snapnum in snapshots]
        line = fc.ne8_770

    tauwtd_vprof(infofiles, filepatterns,
                 mincoldens=(13.25, 13.5, 13.75, 14.0),
                 line=line, outfilen=outname)