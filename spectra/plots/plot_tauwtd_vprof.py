import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fire_an.makeplots.tol_colors as tc
import fire_an.simlists as sl
import fire_an.spectra.plots.select_obsdata as so
import fire_an.spectra.findcomponents as fc
import fire_an.utils.constants_and_units as c

#TODO: add CASBaH data
#TODO: get sim mass range from actual sim data
def plot_tauwtd_vprof(simdatafilen: str,
                      line: fc.Line = fc.ne8_770,
                      outname: str | None = None):
    simcolors = tc.tol_cset('vibrant')

    simdata = {}
    with h5py.File(simdatafilen, 'r') as f:
        fkeys = list(f.keys())
        for fkey in fkeys:
            if fkey.startswith('logN_ge_'):
                key = float(fkey.split('_')[-1])
                label = f'FIRE-2, $\\log N \\geq {key:.2f}$'
            else:
                key = fkey
                label = fkey
            tau = f[fkey]['tau_av'][:]
            vel = f[fkey]['v_gal_kmps'][:]
            nsys = len(f[fkey]['included_spectra'])
            label = label + f' ({nsys})'
            simdata[key] = {'label': label,
                            'tau': tau,
                            'vel': vel,
                            'nsys': nsys}
    
    logMh_sim = np.array([11.5, 12.3]) # just set for now
    df_sys_cubs = so.selectsystems(logMh_sim,
                                   np.array([0.5, 1.]),
                                   'cubs',
                                   Mhmargin_dex=0.2,
                                   zmargin=0.05,
                                   method='bestest_inrange')
    df_comp_cubs = pd.read_csv(so.filen_cubs_components, 
                               sep='\t', comment='#',
                               converters={
                                    'sigmav [km/s]': np.float64,
                                    'QSO': str,
                                    'system': str,
                                    'impact_parameter [kpc]': np.float64,
                                    'bline [km/s]': np.float64})
    sysids = np.unique(df_sys_cubs['system'])
    vgrid_kmps = np.arange(-1000., 1000.5, 1.)
    nuspec_Hz = c.c / (line.wavelength_A * 1e-8 \
                       * (1. + vgrid_kmps * 1e5 / c.c))
    nspec = 0
    tau_cubs = np.zeros(len(vgrid_kmps), dtype=np.float64)
    for sid in sysids:
        df_comps_this = df_comp_cubs[df_comp_cubs['system'] == sid]
        logcds_cm2 = np.array(df_comps_this['Nline [logcm2]'])
        bvals_kmps = np.array(df_comps_this['bline [km/s]'])
        centers_kmps = np.array(df_comps_this['vline [km/s]'])
        # deal with line width ULs
        bline_ul = np.isnan(bvals_kmps)
        bvals_kmps[bline_ul] = \
            np.array(df_comps_this['err_bline_plus [km/s]'])[bline_ul]
        _tau = line.getspectrum(nuspec_Hz, logcds_cm2, bvals_kmps, 
                                centers_kmps)
        print(_tau)
        _tau = -1. * np.log(_tau) # optical depth from transmission
        print(_tau) 
        tau_cubs += _tau
        nspec += 1
    tau_cubs *= 1. / nspec
    label_cubs = f'CUBS: $\\log N \\gtrsim 14.0$ ({nspec})'

    fig = plt.figure(figsize=(5.5, 5.))
    ax = fig.add_subplot(1, 1, 1)
    fontsize = 12
    ax.tick_params(which='both', direction='in', labelsize=fontsize - 1.,
                   right=True, top=True)
    ax.set_xlabel('$v - v_{\\mathrm{gal}} \\; [\\mathrm{km}/\\mathrm{s}]$',
                  fontsize=fontsize)
    ax.set_ylabel('$\\langle \\tau \\rangle$')
    
    skeys = list(simdata.keys())
    skeys.sort()
    for color, skey in zip(simcolors, skeys):
        _data = simdata[skey]
        ax.plot(_data['vel'], _data['tau'], label=_data['label'],
                color=color, linewidth=1.2, linestyle='solid')
    ax.plot(vgrid_kmps, tau_cubs, label=label_cubs,
            linestyle='solid', color='black', linewidth=1.7)
    
    ax.set_xlim(-750., 750.)
    ax.legend(fontsize=fontsize - 1)

    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

    

def plotset_tauwtd_prof(testset=4):
    if testset == 4:
        #simnames = sl.m12_f2md
        #snapshots = [sl.snaps_f2md[0], sl.snaps_f2md[1]]
        #filedir = '/projects/b1026/nastasha/spectra/test4/'
        outdir =  '/projects/b1026/nastasha/imgs/spectra/test4/'
        # the one that worked the first time
        #simnames = ['crheatfix_m12i_r7100', 'crheatfix_m12f_r7100']
        #snapshots = [294, 277]
        simdatafilen = outdir + f'tau_av_ne8_770_minlogN_testset4.hdf5'
        
        # filepatterns = [filedir + f'/tridentray_{simname}_{snapnum}_*.txt'
        #                for simname in simnames for snapnum in snapshots]
        # infofiles = [filedir + f'/tridentray_{simname}_{snapnum}_info.hdf5'
        #                for simname in simnames for snapnum in snapshots]
        line = fc.ne8_770
        outname = outdir + 'datacomp_tau_av_ne8_770_minlogN_testset4.pdf'
    
    plot_tauwtd_vprof(simdatafilen, line=line,
                      outname=outname)