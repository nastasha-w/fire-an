'''
convenient lists of simulations to avoid copying lists too often
'''

import makeplots.tol_colors as tc

# clean: ICs for each noBH/AGN-noCR/AGN-CR run, snapshots down to z=0.5
m13_nobh_clean1 = [
    'm13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
]
m13_agnnocr_clean1 = [
    ('m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
     '_gacc31_fa0.5_fcr1e-3_vw3000'),
    ('m13h113_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4'
     '_gacc31_fa0.5_fcr1e-3_vw3000'),
]
m13_agncr_clean1 = [
    ('m13h113_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp1e-4_gacc31'
     '_fa0.5'),
    ('m13h206_m3e4_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_sdp3e-4_gacc31'
     '_fa0.5'),
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

# rest: all other sims in the noBH/AGN-noCR/AGN-CR suites. 
# ICs are not complete sets for physics models.
# must have run to redshift 0.5
m13_nobh_rest1 = [
    'm13h002_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm13h007_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm13h029_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
    'm13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
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

## run sets:
m13_sr_clean1 = m13_nobh_clean1 + m13_agnnocr_clean1 # len 4
m13_hr_clean1 = m13_agncr_clean1.copy() # len 2
m12_sr_clean1 = m12_agncr_clean1.copy() # len 3
m12_hr_clean1 = m12_nobh_clean1 + m12_agnnocr_clean1 # len 6

m13_sr_rest1 = m13_nobh_rest1 + m13_agnnocr_rest1 # len 7
m13_hr_rest1 = m13_agncr_rest1.copy() # len 6
m12_sr_rest1 = m12_agncr_rest1.copy() # len 1
m12_hr_rest1 = m12_nobh_rest1 + m12_agnnocr_rest1 # len 10

m13_sr_all1 = m13_sr_clean1 + m13_sr_rest1 # len 11
m13_hr_all1 = m13_hr_clean1 + m13_hr_rest1 # len 8
m12_sr_all1 = m12_sr_clean1 + m12_sr_rest1 # len 4
m12_hr_all1 = m12_hr_clean1 + m12_hr_rest1 # len 16


# to find matching snaps
snapmatch = {
    'm13_sr': set(m13_sr_clean1),
    'm13_hr': set(m13_hr_clean1),
    'm12_sr': set(m12_sr_clean1),
    'm12_hr': set(m12_hr_clean1),
    }

snaps_sr = [45, 46, 47, 48, 49, 50]
snaps_hr = [186, 197, 210, 224, 240, 258]

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
ics_m12 = ['m12f', 'm12i', 'm12m', 'm12b', 'm12c', 'm12r', 'm12q', 'm12w',
           'm12z']
m12_iccolors = {ic: _iccolors[i] for i, ic in enumerate(ics_m12)}
ics_m13 = ['m13h113', 'm13h206', 'm13h002', 'm13h007', 'm13h009',
           'm13h029', 'm13h037', 'm13h236']
m13_iccolors = {ic: _iccolors[i] for i, ic in enumerate(ics_m13)}
m13_iccolors['h02'] = m13_iccolors['h002']

_physcolors = tc.tol_cset('bright')
physcolors = {'AGN-CR': _physcolors.green,
              'AGN-noCR': _physcolors.red,
              'noBH': _physcolors.blue,
              }

physlinestyles = {'AGN-CR': 'dotted',
                  'AGN-noCR': 'dashed',
                  'noBH': 'solid',
                  }