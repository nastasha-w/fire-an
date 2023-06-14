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

ddir = '/Users/nastasha/ciera/data_smdpl_sfr/'
def loaddata(aexp):
    catopts = glob.glob(ddir + 'sfr_catalog_*.bin')
    aexps = [float(((catopt.split('/')[-1]).split['.'][0]).split('_')[-1])
             for catopt in catopts]
    seli = np.argmin(np.abs(np.array(aexps - aexp)))
    aexp = aexps[seli]
    print(f'selected file at aexp={aexp}')
    filen = catopts[seli]
    halos = np.fromfile(filen, dtype=dtype)
    return halos, aexp
