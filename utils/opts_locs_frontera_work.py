# -*- coding: utf-8 -*-



############################
#    setup and functions   #
############################
dir_halodata = '/work2/08466/tg877653/frontera/halodata/'
filen_halocenrvir = dir_halodata + 'cen_rvir.hdf5'

frontera_work = '/work2/08466/tg877653/frontera/'
frontera_scratch = '/scratch1/08466/tg877653/'
pre = frontera_work

c_interpfile  = pre + 'code/proj-an-c/interp2d/interp.so'
c_gridcoarser = '' #'/home/wijers/gridcoarser/gridcoarser.so'
c_halomask    = '' #'/home/wijers/halomaps/gethalomap.so'
hsml_dir = pre + 'code/proj-an-c/HsmlAndProject_OMP/'
emtab_sylvia_ssh = pre + 'iontab/PS20/UVB_dust1_CR1_G1_shield1_lines.hdf5' 
iontab_sylvia_ssh = pre + 'iontab/PS20/UVB_dust1_CR1_G1_shield1.hdf5'
simdir_fire3 = '/scratch3/01799/phopkins/fire3_suite_done/'
simdir_fire2_md = ''
simdir_fire3_m12plus = ''
simdir_fire3x_tests = '/scratch3/01799/phopkins/fire3_suite_done/'

kernel_list = ['C2','gadget']
# desngb = 58 read out from sample hdf5 file (RunTimePars)
# needed for C projection routine, but only used if smoothing lengths
# need to be calculated
desngb = 58

# copied from Parameters/ChemicalElements in EAGLE simulation hdf5 
# files.
# matches Wiersema, Schaye, Theuns et al. 2009 table 1 values
solar_abunds_ea = {'calcium':  6.435500108636916E-5,
                'carbon':   0.002066543558612466,
                'helium':   0.2805553376674652,
                'hydrogen': 0.706497848033905,
                'iron':     0.0011032151523977518,
                'magnesium':5.907064187340438E-4,
                'neon':     0.0014144604792818427,
                'nitrogen': 8.356256294064224E-4,
                'oxygen':   0.00549262436106801,
                'silicon':  6.825873861089349E-4,
                'sulfur':  4.0898521547205746E-4}

Zsun_ea = 0.0127 # Wiersma, Schaye,Theuns et al. (2009), consistent with Serena's tables

