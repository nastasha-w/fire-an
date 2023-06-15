'''
dtype copied from https://halos.as.arizona.edu/UniverseMachine/DR1/SMDPL_SFR/loader.py
'''
import numpy as np
import glob

dtype = np.dtype(dtype=[('id', 'i8'),('descid','i8'),('upid','i8'),
                        ('flags', 'i4'), ('uparent_dist', 'f4'),
                        ('pos', 'f4', (6)), ('vmp', 'f4'), ('lvmp', 'f4'),
                        ('mp', 'f4'), ('m', 'f4'), ('v', 'f4'), ('r', 'f4'),
                        ('rank1', 'f4'), ('rank2', 'f4'), ('ra', 'f4'),
                        ('rarank', 'f4'), ('A_UV', 'f4'), ('sm', 'f4'), 
                        ('icl', 'f4'), ('sfr', 'f4'), ('obs_sm', 'f4'), 
                        ('obs_sfr', 'f4'), ('obs_uv', 'f4'), ('empty', 'f4')],
                 align=True)

#Field explanations:
#**Note that halo masses are in Msun/h and stellar masses/SFRs are in Msun.
#ID: Unique halo ID
#DescID: ID of descendant halo (or -1 at z=0).
#UPID: -1 for central halos, otherwise, ID of largest parent halo
#Flags: Ignore
#Uparent_Dist: Ignore
#pos[6]: (X,Y,Z,VX,VY,VZ)
#X Y Z: halo position (comoving Mpc/h)
#VX VY VZ: halo velocity (physical peculiar km/s)
#M: Halo mass (Bryan & Norman 1998 virial mass, Msun/h)
#V: Halo vmax (physical km/s)
#MP: Halo peak historical mass (BN98 vir, Msun/h)
#VMP: Halo vmax at the time when peak mass was reached.
#R: Halo radius (BN98 vir, comoving kpc/h)
#Rank1: halo rank in Delta_vmax (see UniverseMachine paper)
#Rank2, RA, RARank: Ignore
#A_UV: UV attenuation (mag)
#SM: True stellar mass (Msun)
#ICL: True intracluster stellar mass (Msun)
#SFR: True star formation rate (Msun/yr)
#Obs_SM: observed stellar mass, including random & systematic errors (Msun)
#Obs_SFR: observed SFR, including random & systematic errors (Msun/yr)
#Obs_UV: Observed UV Magnitude (M_1500 AB)

#ddir = '/Users/nastasha/ciera/data_smdpl_sfr/'
ddir = '/projects/b1026/nastasha/extdata/data_smdpl_sfr/'
def loaddata(aexp):
    catopts = glob.glob(ddir + 'sfr_catalog_*.bin')
    filetrunks = [(catopt.split('/')[-1])[:-4]
                  for catopt in catopts]
    aexps = [float((trunk).split('_')[-1])
             for trunk in filetrunks]
    seli = np.argmin(np.abs(np.array(np.array(aexps) - aexp)))
    aexp_used = aexps[seli]
    print(f'selected file at aexp={aexp_used}, for target {aexp}')
    filen = catopts[seli]
    halos = np.fromfile(filen, dtype=dtype)
    return halos, aexp_used
