'''
make it easy to find/re-run the plots for the Ne VIII column
density paper
'''
import numpy as np

import fire_an.makeplots.litcomp.plot_obscomp as poc
import fire_an.makeplots.plotmaps_overview as pmo
import fire_an.makeplots.plot_clumpiness as pcl
import fire_an.makeplots.comp3d_3model.plot_r3d_perc_indiv as pri
import fire_an.makeplots.cgm_fgas_zne.plotprop as pcp
import fire_an.makeplots.litcomp.obs_vs_analytical as ova
import fire_an.mstar_mhalo.plot_mstar_mh as psh
import fire_an.makeplots.litcomp.smhm_issues as shi
import fire_an.makeplots.makesimtable as mst

def runplots():
    ## data comparisons
    ## output: /projects/b1026/nastasha/imgs/datacomp/
    # m12 and m13 haloes vs. B+19 and Q+23 data
    # also runs z=0.5 -- 0.7 sims comp to Q+23 data
    # and runs halo mass vs. redshift for sims and B+19 and Q+23
    # figs 6, 7, 8 (5 pdfs)
    poc.runplots_obscomp()
    # compares m12 FIRE-3 NoBH main sample to the higher halo and 
    # stellar mass sample
    # also runs halo mass vs. redshift for sims and B+19 and Q+23
    # fig 14 (2 pdfs)
    poc.runplots_appendix()

    ## images (selection)
    ## output: /projects/b1026/nastasha/imgs/summary_plots/ne8_thumbnails/
    # figs. 1, 2
    # m12f, m13h113 at z=1.0 actually used in the paper
    pmo.plotoverview_ne8(['m12f', 'm13h113', 'm13h206'], [0, 5])

    ## clumpiness plot
    ## output: /projects/b1026/nastasha/imgs/clumpiness/
    # fig. 3
    pcl.plot_clumpiness_prof(target='Ne8')

    ## radial profiles nH/T/Z
    ## output: /projects/b1026/nastasha/imgs/3dprof/
    # fig. 4
    pri.plotmain()

    ## cgm property hists
    ## output: /projects/b1026/nastasha/imgs/cgmprop/
    # fig 5, fig. for higher stellar mass m12 FIRE-3 NoBH set
    pcp.plotpanels_main(massset='m12')
    pcp.plotpanels_general('m12', rrange_rvir=(0.1, 1.0),
                           trange_logk=(5.0, np.inf),
                           panels=('meanNe8col', 'fgas', 'ZNe', 'Ne8frac',
                                   'Mvir', 'Mstarcen', 'ZoverMstarcen', 
                                   'ZoverMstarhalo', 'Ztot', 'Ne8'),
                           xlabels=None, inclm12plus=True)
    pcp.plotpanels_general('m12', rrange_rvir=(0.1, 1.0),
                           trange_logk=(-np.inf, np.inf),
                           panels=('meanNe8col', 'fgas', 'ZNe', 'Ne8frac',
                                   'Mvir', 'Mstarcen', 'Ztot', 'Ne8'),
                           xlabels=None, inclm12plus=True)
    pcp.plotpanels_main(massset='m13')
    pcp.plotpanels_general('m13', rrange_rvir=(0.1, 1.0),
                           trange_logk=(5.0, np.inf),
                           panels=('meanNe8col', 'fgas', 'ZNe', 'Ne8frac',
                                   'Mvir', 'Mstarcen', 'ZoverMstarcen', 
                                   'ZoverMstarhalo', 'Ztot', 'Ne8'),
                           xlabels=None, inclm12plus=True)
    pcp.plotpanels_general('m13', rrange_rvir=(0.1, 1.0),
                           trange_logk=(-np.inf, np.inf),
                           panels=('meanNe8col', 'fgas', 'ZNe', 'Ne8frac',
                                   'Mvir', 'Mstarcen', 'Ztot', 'Ne8'),
                           xlabels=None, inclm12plus=True)
    
    ## data and analytical models
    ## output: /projects/b1026/nastasha/imgs/analytical/
    # figs 9, 10
    ova.plot_plmodel_datacomp_Kvar(obsdata=('Q+23', 'B+19'))
    ova.plot_plmodel_datacomp_parvar(obsdata=('Q+23', 'B+19'))

    ## halo mass calc. appendix
    ## output: /projects/b1026/nastasha/imgs/datacomp/smhm/
    # figs 11, 13
    psh.plot_mstar_mh([0.0, 0.5, 1.0], variation='main')
    psh.plot_mhdist_for_mstar_example_scatterdecomp(0.74)
    # fig 12
    shi.plot_hm_uncertainty_um_b19()

    ## tables 1, 2
    ## output: printed to screen
    print('\n' * 6)
    mst.maketable_main(simset='all')
    print('\n' * 6)
    mst.maketable_main(simset='m12plus')

if __name__ == '__main__':
    runplots()

