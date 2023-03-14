# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 16:46:28 2016

@author: wijers

Helper file for make_maps: 
containts lists, directories, etc. used in calculations
Note: does not actually work on my laptop! 
"""

############################
#    setup and functions   #
############################

c_interpfile  = '/Users/Nastasha/code/proj-an-c/interp2d/interp.so'
c_gridcoarser = '' #'/home/wijers/gridcoarser/gridcoarser.so'
c_halomask    = '' #'/home/wijers/halomaps/gethalomap.so'
hsml_dir = '/Users/Nastasha/code/proj-an-c/HsmlAndProject_OMP/'
emtab_sylvia_ssh = ('/Users/Nastasha/phd/tables/lineem/lines_sp20/'
                    'UVB_dust1_CR1_G1_shield1_lines.hdf5') 
iontab_sylvia_ssh = ('/Users/Nastasha/phd/tables/ionbal/lines_sp20/'
                     'UVB_dust1_CR1_G1_shield1.hdf5')

simdir_fire = '/dummyvalue/'
dir_halodata = '/Users/Nastasha/ciera/halodata_fire/'
filen_halocenrvir = dir_halodata + 'cen_rvir.hdf5'

kernel_list = ['C2','gadget']
# desngb = 58 read out from EAGLE sample hdf5 file (RunTimePars)
# needed for projection routine, 
# but only for (not used) smoothing length calculation
desngb = 58

# copied from Parameters/ChemicalElements in EAGLE simulation hdf5 files.
# Seems to be the same accross simulations 
# (and it should be, if the same cooling tabels are used) 
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
# Wiersma, Schaye,Theuns et al. (2009), consistent with Serena's tables
Zsun_ea = 0.0127 



