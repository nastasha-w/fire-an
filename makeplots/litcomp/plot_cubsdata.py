import matplotlib.pyplot as plt
import numpy as np

import makeplots.litcomp.obsdataread as cubsdr

def plot_msmh_cubs7():
    outname = ('/projects/b1026/nastasha/imgs/datacomp/smhm/'
               'cubs7_mstar_to_UM_mhalo_2s_err.pdf')
    plotdata = cubsdr.getplotdata_cubs()

    fig = plt.figure(figsize=(5.5, 5.))
    ax = fig.add_subplot(1, 1, 1)
    fontsize = 12

    ms = plotdata['mstar_log10Msun']
    ms_dlo = plotdata['mstar_2s_loerr_dex']
    ms_dhi = plotdata['mstar_2s_hierr_dex']
    mh = plotdata['logmvir_msun_bestest']
    mh_dlo = mh - plotdata['logmvir_msun_loer']
    mh_dhi = plotdata['logmvir_msun_hier'] - mh
    
    f1 = plotdata['isul_ne8']
    ax.errorbar(ms[f1], mh[f1], xerr=(ms_dlo[f1], ms_dhi[f1]), 
                yerr=(mh_dlo[f1], mh_dhi[f1]),
                color='black', alpha=0.3, linestyle='none',
                label='Ne VIII UL, Qu+23, $2\sigma$')
    f2 = np.logical_not(f1)
    ax.errorbar(ms[f2], mh[f2], xerr=(ms_dlo[f2], ms_dhi[f2]), 
                yerr=(mh_dlo[f2], mh_dhi[f2]),
                color='blue', alpha=1., linestyle='none',
                label='Ne VIII det., Qu+23, $2\sigma$')
    ax.tick_params(which='both', direction='in', right=True, top=True,
                   labelsize=fontsize - 1)
    ax.set_xlabel('$\\log_{10} \\, \\mathrm{M}_{\\star, \\mathrm{cen}}'
                  ' \\; [\\mathrm{M_{\\odot}}]$', fontsize=fontsize)
    ax.set_ylabel('$\\log_{10} \\, \\mathrm{M}_{\\mathrm{vir}}'
                  ' \\; [\\mathrm{M_{\\odot}}]$', fontsize=fontsize)
    ax.legend(fontsize=fontsize - 1, loc='upper left')
    plt.savefig(outname, bbox_inches='tight')

