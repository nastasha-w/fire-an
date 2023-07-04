'''
convenient lists of simulations to avoid copying lists too often
'''

import fire_an.makeplots.tol_colors as tc
import fire_an.utils.opts_locs as ol

# clean: ICs for each noBH/AGN-noCR/AGN-CR run, snapshots down to z=0.5
# 15 IC + phys
m13_nobh_clean1 = [
    'm13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
]
m13_agnnocr_clean1 = [
    ('m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31'
     '_fa0.5'),
    ('m13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31'
     '_fa0.5'),
]
m13_agncr_clean1 = [
    ('m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
     '_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
     '_gacc31_fa0.5_fcr1e-3_vw3000'),
]


m12_nobh_clean1 = [
    'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm12i_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm12m_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
]
m12_agnnocr_clean1 = [
    'm12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
    'm12i_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
    'm12m_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
]
m12_agncr_clean1 = [
    ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
     '_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m12i_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
     '_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m12m_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
     '_gacc31_fa0.5_fcr1e-3_vw3000'),
]

# adding m12q
m12_nobh_clean2 = m12_nobh_clean1 + [
    'm12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
] # len 4
m12_agnnocr_clean2 = m12_agnnocr_clean1 + [
    'm12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
] # len 4
m12_agncr_clean2 = m12_agncr_clean1 + [
    ('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
     '_gacc31_fa0.5_fcr1e-3_vw3000'),
] # len 4

# rest: all other sims in the noBH/AGN-noCR/AGN-CR suites. 
# ICs are not complete sets for physics models.
# must have run to redshift 0.5
# 22 IC + phys
m13_nobh_rest1 = [
    'm13h002_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm13h007_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm13h029_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm13h223_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm13h236_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
]
m13_agnnocr_rest1 = [
]
m13_agncr_rest1 = [
    ('m13h002_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m13h007_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m13h009_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'), 
    ('m13h029_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m13h037_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'), 
    ('m13h236_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
]

m12_nobh_rest1 = [
    'm12b_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm12c_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm12r_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm12w_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm12z_m4e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
]
m12_agnnocr_rest1 = [
    'm12b_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
    'm12c_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
    'm12r_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
    'm12w_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
    'm12z_m4e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
]
m12_agncr_rest1 = [
    ('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
     '_gacc31_fa0.5_fcr1e-3_vw3000'),
]

# m12q -> clean2
m12_nobh_rest2 = m12_nobh_rest1.copy() # len 5
m12_agnnocr_rest2 = m12_agnnocr_rest1.copy() # len 5
m12_agncr_rest2 = [] # len 0

m12_f2md = ['m12z_r4200',
           'm12w_r7100',
           'm12r_r7100',
           'm12i_r7100',
           'm12c_r7100',
           'm12b_r7100',
           'm12m_r7100',
           'm12f_r7100',
           'crheatfix_m12f_r7100',
           'crheatfix_m12i_r7100'
           ]

## experimental m11 selection for C ion series
# phys variations selected by match to m12 series names
# mass res is ~ m12-sr res for all; hr/sr is for snapshot cadence
# z=0: snap 60
m11_sr_nobh_set1 = [
    'm11a_m2e3_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm11b_m2e3_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm11c_m2e3_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm11d_m7e3_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm11e_m7e3_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm11f_m1e4_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm11g_m1e4_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm11h_m7e3_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm11i_m7e3_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm11q_m7e3_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm11v_m7e3_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
] # len 11
# z=0: snap 60
m11_sr_agnnocr_set1 = [
    'm11a_m2e3_MHD_fire3_fireBH_Sep052021_crdiffc690_sdp3e-7_gacc31_fa0.5',
    'm11b_m2e3_MHD_fire3_fireBH_Sep052021_crdiffc690_sdp1e-6_gacc31_fa0.5',
    'm11c_m2e3_MHD_fire3_fireBH_Sep052021_crdiffc690_sdp1e-5_gacc31_fa0.5',
    'm11d_m7e3_MHD_fire3_fireBH_Sep052021_crdiffc690_sdp3e-5_gacc31_fa0.5',
    'm11e_m7e3_MHD_fire3_fireBH_Sep052021_crdiffc690_sdp1e-5_gacc31_fa0.5',
    'm11f_m1e4_MHD_fire3_fireBH_Sep052021_crdiffc690_sdp1e-4_gacc31_fa0.5',
    'm11g_m1e4_MHD_fire3_fireBH_Sep052021_crdiffc690_sdp1e-4_gacc31_fa0.5',
    'm11h_m7e3_MHD_fire3_fireBH_Sep052021_crdiffc690_sdp1e-5_gacc31_fa0.5',
    'm11i_m7e3_MHD_fire3_fireBH_Sep052021_crdiffc690_sdp1e-5_gacc31_fa0.5',
    'm11q_m7e3_MHD_fire3_fireBH_Sep052021_crdiffc690_sdp1e-5_gacc31_fa0.5',
    'm11v_m7e3_MHD_fire3_fireBH_Sep052021_crdiffc690_sdp3e-5_gacc31_fa0.5',
] # len 11
# z=0: snap 500
m11_hr_agncr_set1 = [
    ('m11a_m2e3_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_'
     'sdp3e-7_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m11b_m2e3_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-6_gacc31_fa0.5_fcr1e-3_vw3000'),
    #('m11c_m2e3_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
    # '_sdp1e-5_gacc31_fa0.5_fcr1e-3_vw3000'), # last snap 152 (2023-04-12)
    ('m11d_m7e3_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp3e-5_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m11e_m7e3_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-5_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m11f_m1e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m11g_m1e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m11h_m7e3_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-5_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m11q_m7e3_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-5_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m11v_m7e3_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp3e-5_gacc31_fa0.5_fcr1e-3_vw3000'),
] # len 9; left out sim that didn't reach z=1.0.
# z=0: snap 60
m11_sr_agncr_set1 = [
    ('m11i_m7e3_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-5_gacc31_fa0.5_fcr3e-3_vw3000'),
] # len 1
m11_hr_agnnocr_set1 =[]
m11_hr_nobh_set1 = []

# len 9
m11_hr_set1 = m11_hr_agncr_set1.copy()
# len 23
m11_sr_set1 = m11_sr_agncr_set1 + m11_sr_agnnocr_set1 + m11_sr_nobh_set1 


## run sets:
# set 1
m13_sr_clean1 = m13_nobh_clean1 + m13_agncr_clean1 # len 4
m13_hr_clean1 = m13_agnnocr_clean1.copy() # len 2
m12_sr_clean1 = m12_agncr_clean1.copy() # len 3
m12_hr_clean1 = m12_nobh_clean1 + m12_agnnocr_clean1 # len 6

m13_sr_rest1 = m13_nobh_rest1 + m13_agncr_rest1 # len 11
m13_hr_rest1 = m13_agnnocr_rest1.copy() # len 0
m12_sr_rest1 = m12_agncr_rest1.copy() # len 1
m12_hr_rest1 = m12_nobh_rest1 + m12_agnnocr_rest1 # len 10

m13_sr_all1 = m13_sr_clean1 + m13_sr_rest1 # len 15
m13_hr_all1 = m13_hr_clean1 + m13_hr_rest1 # len 2
m12_sr_all1 = m12_sr_clean1 + m12_sr_rest1 # len 4
m12_hr_all1 = m12_hr_clean1 + m12_hr_rest1 # len 16

# set 2: added new m12q haloes (only affects m12_sr, really)
m13_nobh_clean2 = m13_nobh_clean1.copy()
m13_agnnocr_clean2 = m13_agnnocr_clean1.copy()
m13_agncr_clean2 = m13_agncr_clean1.copy()

m13_nobh_rest2 = m13_nobh_rest1.copy()
m13_agnnocr_rest2 = m13_agnnocr_rest1.copy()
m13_agncr_rest2 = m13_agncr_rest1.copy()

m13_sr_clean2 = m13_nobh_clean2 + m13_agncr_clean2 # len 4
m13_hr_clean2 = m13_agnnocr_clean2.copy() # len 2
m12_sr_clean2 = m12_agncr_clean2.copy() # len 4
m12_hr_clean2 = m12_nobh_clean2 + m12_agnnocr_clean2 # len 8

m13_sr_rest2 = m13_nobh_rest2 + m13_agncr_rest2 # len 11
m13_hr_rest2 = m13_agnnocr_rest2.copy() # len 0
m12_sr_rest2 = m12_agncr_rest2.copy() # len 0
m12_hr_rest2 = m12_nobh_rest2 + m12_agnnocr_rest2 # len 10

m13_sr_all2 = m13_sr_clean2 + m13_sr_rest2 # len 15
m13_hr_all2 = m13_hr_clean2 + m13_hr_rest2 # len 2
m12_sr_all2 = m12_sr_clean2 + m12_sr_rest2 # len 4
m12_hr_all2 = m12_hr_clean2 + m12_hr_rest2 # len 18

## subsets that ran to z=0 (2023-04-12)
m13_sr_all2_z0 =[
    ('m13h002_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m13h007_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690'
     '_sdp1e10_gacc31_fa0.5'),
    ('m13h007_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m13h009_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m13h029_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m13h037_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690'
     '_sdp1e10_gacc31_fa0.5'),
    ('m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m13h236_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000')
] # length 10; mostly AGN-CR, a 2x noBH
m13_hr_all2_z0 = [] # none
# all m12s ran to z=0
m12_sr_all2_z0 = [
    ('m12f_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m12i_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m12m_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m12q_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1'
     '_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000') 
] # length 4; AGN-CR
m12_hr_all2_z0 = [
    ('m12b_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp1e10_gacc31_fa0.5'),
    ('m12b_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp2e-4_gacc31_fa0.5'),
    ('m12c_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp1e10_gacc31_fa0.5'),
    ('m12c_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp1e-4_gacc31_fa0.5'),
    ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp1e10_gacc31_fa0.5'),
    ('m12f_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp2e-4_gacc31_fa0.5'),
    ('m12i_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp1e10_gacc31_fa0.5'),
    ('m12i_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp2e-4_gacc31_fa0.5'),
    ('m12m_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp1e10_gacc31_fa0.5'),
    ('m12m_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp2e-4_gacc31_fa0.5'),
    ('m12q_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp1e10_gacc31_fa0.5'),
    ('m12r_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp1e10_gacc31_fa0.5'),
    ('m12r_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp1e-4_gacc31_fa0.5'),
    ('m12w_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp1e10_gacc31_fa0.5'),
    ('m12w_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp1e-4_gacc31_fa0.5'),
    ('m12z_m4e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp1e10_gacc31_fa0.5'),
    ('m12z_m4e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690'
     '_sdp1e-4_gacc31_fa0.5'),
] # len 17; noBH and AGN-noCR


# to find matching snaps
snapmatch = {
    'm13_sr': set(m13_sr_all2),
    'm13_hr': set(m13_hr_all2),
    'm12_sr': set(m12_sr_all2),
    'm12_hr': set(m12_hr_all2),
    'm11_sr': set(m11_sr_set1),
    'm11_hr': set(m11_hr_set1),
    }

snaps_sr = [45, 46, 47, 48, 49, 50]
snaps_hr = [186, 197, 210, 224, 240, 258]
snaps_f2md = [277, 294, 312, 332, 356, 382]

# m12q AGN-noCR: no snapshot 500 (z=0)
# m12q noBH: does have snapshot 500 (z=0) 
# m13_sr_all1: 
#     snap 60 (z=0) for m13h206_m3e5, 
snaps_sr_051 = [45, 50, 60]
snaps_hr_051 = [186, 258, 500]


snaplists = {
    'm13_sr': snaps_sr,
    'm13_hr': snaps_hr,
    'm12_sr': snaps_sr,
    'm12_hr': snaps_hr,
}

# simulations which may be affected by bugs (marked in Lindsey's list)
# mostly for plotting selection
buglist1 = [
    'm12i_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm12i_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp2e-4_gacc31_fa0.5',
    'm12z_m4e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31_fa0.5',
    ('m12i_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
     '_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m12m_m6e4_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
     '_gacc31_fa0.5_fcr1e-3_vw3000'),
]


## plotting sync/convenience
_iccolors = tc.tol_cset('muted')
ics_m12 = ['m12f', 'm12i', 'm12m', 'm12q', 'm12b', 'm12c',
           'm12r', 'm12w', 'm12z']
m12_iccolors = {ic: _iccolors[i] for i, ic in enumerate(ics_m12)}
ics_m13 = ['m13h002', 'm13h007', 'm13h009', 'm13h029',
           'm13h037', 'm13h113', 'm13h206', 'm13h223', 'm13h236']
m13_iccolors = {ic: _iccolors[i] for i, ic in enumerate(ics_m13)}
m13_iccolors['m13h02'] = m13_iccolors['m13h002']

_physcolors = tc.tol_cset('bright')
physcolors = {'AGN-CR': _physcolors.green,
              'AGN-noCR': _physcolors.red,
              'noBH': _physcolors.blue,
              'FIRE-2': _physcolors.yellow}

physlinestyles = {'AGN-CR': 'dotted',
                  'AGN-noCR': 'dashed',
                  'noBH': 'solid',
                  'FIRE-2': 'dashdot'}

## because the FIRE-2 directory = phys. model structure broke my
# 'store by single directory level' setup
def dirpath_from_simname(simname):
    if simname in m12_f2md:
        if simname.startswith('crheatfix'):
            _sn = simname.split('_')[1:]
            dirpath = ol.simdir_fire + 'cr_heating_fix/' + '_'.join(_sn)
        else:
            dirpath = ol.simdir_fire + simname
    else:
        dp2 = '_'.join(simname.split('_')[:2])
        if dp2.startswith('m13h02_'):
            dp2 = dp2.replace('m13h02', 'm13h002')
        dirpath = '/'.join([ol.simdir_fire, dp2, simname])
    return dirpath

def simname_from_dirpath(dirpath):
    # remove 'output' subdir, trailing slashes
    if dirpath.endswith('output'):
        simpart = dirpath.split('/')[:-2]
    elif dirpath.endswith('output/'):
        simparts = dirpath.split('/')[:-3]
    elif dirpath.endswith('/'):
        simparts = dirpath.split('/')[:-2]
    else:
        simparts = dirpath.split('/')[:-1]
    if simparts[-2] == 'cr_heating_fix':
        simname = 'crheatfix_' + simparts[-1] + '/'
    else:
        simname = simparts[-1] + '/'
    return simname
    
