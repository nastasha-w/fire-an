'''
dtype copied from 
https://halos.as.arizona.edu/UniverseMachine/DR1/SMDPL_SFR/loader.py
catalogues downloaded from 
https://halos.as.arizona.edu/UniverseMachine/DR1/SMDPL_SFR/
using 
wget -r --no-parent -A 'sfr_catalog_0.[4-9]*.bin' https://halos.as.arizona.edu/UniverseMachine/DR1/SMDPL_SFR/ 
The pattern matching example is for catalogs with expansion factors in
the half-open interval [0.4, 1.0). I downloaded the expansion factor 
1.0 catalog separately.
'''
import glob
import h5py
import numpy as np
import os

import fire_an.utils.h5utils as h5u
import fire_an.utils.math_utils as mu

## directory where the catalog files live
#ddir_smdpl = '/Users/nastasha/ciera/data_smdpl_sfr/'
ddir_smdpl = ('/projects/b1026/nastasha/extdata/data_smdpl_sfr/'
             'halos.as.arizona.edu/UniverseMachine/DR1/SMDPL_SFR/')
# from the website
cosmopars_smdpl = {'omegalambda': 0.693,
                   'omegam': 0.307, 
                   'omegab': 0.048,
                   'h': 0.678}
## copied from 
## https://halos.as.arizona.edu/UniverseMachine/DR1/SMDPL_SFR/loader.py
## --------------------------------------------------------------------
dtype_smdpl = np.dtype(dtype=[
    ('id', 'i8'),('descid','i8'),('upid','i8'),
    ('flags', 'i4'), ('uparent_dist', 'f4'),
    ('pos', 'f4', (6)), ('vmp', 'f4'), ('lvmp', 'f4'),
    ('mp', 'f4'), ('m', 'f4'), ('v', 'f4'), ('r', 'f4'),
    ('rank1', 'f4'), ('rank2', 'f4'), ('ra', 'f4'),
    ('rarank', 'f4'), ('A_UV', 'f4'), ('sm', 'f4'), 
    ('icl', 'f4'), ('sfr', 'f4'), ('obs_sm', 'f4'), 
    ('obs_sfr', 'f4'), ('obs_uv', 'f4'), ('empty', 'f4')],
    align=True)

def loadsmpdl(filen):
    return np.fromfile(filen, dtype=dtype_smdpl)
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
## --------------------------------------------------------------------
class SMHMhists:
    def __init__(self, ztargets, histdata='UM-smdpl', binsize=0.1):
        '''
        Parameters:
        -----------
        Class reads in halo catalogue files, calculates histograms, and
        saves the histograms. On later calls, the stored histograms are
        retrieved. 
        Methods allow for extraction of percentiles and calculation of
        uncertainty ranges for a halo/stellar mass given a probability
        distribution for the stellar/halo mass.

        ztargets: array-like of floats
            redshifts to load. closest options will be used; these are
            then the redshift options for which the closest option is 
            chosen in e.g., median retrieval
        histdata: {'UM-smdpl'}
            which halo catalog to use.
        binsize: float
            the size of the halo and stellar mass bins, in log10 Msun.
            The default is 0.1.
        '''
        self.ztargets = np.array(ztargets)
        self.atargets = 1. / (1. + self.ztargets)
        self.histdata_type = histdata
        self.binsize = binsize

        self.hists = {}
        self.mhbins = {}
        self.msbins = {}
        self.msbins_nozero = {}
        if self.histdata_type == 'UM-smdpl':
            self._finddata_smdpl(self.atargets)
        self.zs_used = np.sort(np.array(list(self.hists.keys())))
    
    def getperc_msmh(self, z_target, mode='mstomh', percvals=np.array([0.5])):
        '''
        retrieve the halo mass percentiles at a given stellar mass, 
        or stellar mass percentiles at a given halo mass.

        Parameters:
        -----------
        z_target: float
            redshift for which to get the Ms-Mh relation.
            Just use the closest loaded redshift.
        mode: {'mstomh', 'mhtoms'}
            'mstomh': percentiles of halo mass at fixed stellar mass
            'mhtoms': percentiles of stellar mass at fixed halo mass
        percvals: array-like of floats in the 0 -- 1 range
            which percentiles to calculate
        
        Returns:
        --------
        msvals: array of floats
            stellar mass values. if getting ms percentiles, 
            shape is (number of percentiles, number of mhvals).
            Units are log10 Msun.
        mhvals:array of floats
            halo mass values. if getting mh percentiles, 
            shape is (number of percentiles, number of msvals).
            Units are log10 Msun.
        '''
        self._usez = self.zs_used[np.argmin(np.abs(z_target - self.zs_used))]
        print(f'Using redshift {self._usez} for target {z_target}')
        self._msbins = self.msbins[self._usez]
        self._msbins_nozero = self.msbins_nozero[self._usez]
        self._mhbins = self.mhbins[self._usez]
        self._hist = self.hists[self._usez]
        if mode == 'mstomh':
            mhv = mu.percentiles_from_histogram(self._hist, self._mhbins,
                                                axis=0, percentiles=percvals)
            msv = 0.5 * (self._msbins_nozero[:-1] + self._msbins_nozero[1:])
        if mode == 'mhtoms':
            msv = mu.percentiles_from_histogram(self._hist, self._msbins,
                                                axis=1, percentiles=percvals)
            mhv = 0.5 * (self._mhbins[:-1] + self._mhbins[1:])
        del self._msbins, self._msbins_nozero, self._mhbins, self._hist
        del self._usez
        return msv, mhv
    
    def getbins(self, z_target, mode='ms'):
        '''
        useful to calculate the probability distribution for an input
        halo/stellar mass over the right ms/mh bins

        Parameters:
        -----------
        z_target: float
            redshift for which to get the bins.
            Just uses the closest loaded redshift.
        mode: {'mh', 'ms'}
            'mh': get halo mass bins
            'ms': get stellar mass bins 
                  (-np.inf bin set to a finite value)

        Returns:
        --------
        bins: array of floats
            bin edges, in log10 Msun
        '''
        self._usez = self.zs_used[np.argmin(np.abs(z_target - self.zs_used))]
        if mode == 'mh':
            out = self.mhbins[self._usez]
        elif mode == 'ms':
            out = self.msbins_nozero[self._usez]
        del self._usez
        return out
    
    def gethist(self, z_target):
        '''
        Parameters:
        -----------
        z_target: float
            redshift for which to get the bins.
            Just uses the closest loaded redshift.

        Returns:
        --------
        histogram: 2d-array of floats
            halo/central galaxy counts as a function of halo mass
            (axis 0) and stellar mass (axis 1) 
        '''
        self._usez = self.zs_used[np.argmin(np.abs(z_target - self.zs_used))]
        out = self.hists[self._usez]
        del self._usez
        return out

    def matrixconv(self, pdist, z_target, mode='mstomh'):
        '''
        retrieve the halo mass percentiles at a given stellar mass, 
        or stellar mass percentiles at a given halo mass.

        Parameters:
        -----------
        pdist: array of floats
            probabilities for halo mass or stellar mass 
            (regular, not nonzero) bins
        z_target: float
            redshift for which to use the Ms-Mh relation.
            Just uses the closest loaded redshift.
        mode: {'mstomh', 'mhtoms'}
            'mstomh': pdist is in ms bins, calculate mh probabilities
            'mhtoms': pdist is in mh bins, calculate ms probabilities
        
        Returns:
        --------
        pvals: array of floats
            probabilties in each bin
        bins: array of floats
            bin edges corresponding to the probabilities. 
            Units log10 Msun.
        '''        
        self._usez = self.zs_used[np.argmin(np.abs(z_target - self.zs_used))]
        self._tempmat = np.copy(self.hists[self._usez])
        if mode == 'mstomh':
            if not len(pdist) == len(self.msbins[self._usez]) - 1:
                msg = ('P distribution should be in the stellar mass bins for'
                      'mode "mstomh".')
                raise ValueError(msg)
            self._normfac = np.sum(self._tempmat, axis=0)
            self._tempmat = self._tempmat / self._normfac[np.newaxis, :]
            self._tempmat[:, self._normfac == 0.] = 0.
            pvals_out = np.sum(self._tempmat * pdist[np.newaxis, :], axis=1)
            bins_out = self.mhbins[self._usez]
        elif mode == 'mhtoms':
            if not len(pdist) == len(self.mhbins[self._usez]) - 1:
                msg = ('P distribution should be in the halo mass bins for'
                      'mode "mhtoms".')
                raise ValueError(msg)
            self._normfac = np.sum(self._tempmat, axis=1)
            self._tempmat = self._tempmat / self._normfac[:, np.newaxis]
            self._tempmat[self._normfac == 0., :] = 0.
            pvals_out = np.sum(self._tempmat * pdist[:, np.newaxis], axis=0)
            bins_out = self.msbins_nozero[self._usez]
        del self._usez, self._tempmat, self._normfac
        return pvals_out, bins_out

    def _finddata_smdpl(self, aexps):
        self._catopts = glob.glob(ddir_smdpl + 'sfr_catalog_*.bin')
        self._filetrunks = [(catopt.split('/')[-1])[:-4]
                            for catopt in self._catopts]
        self._aexp_opts = np.array([float((trunk).split('_')[-1])
                                    for trunk in self._filetrunks])
        self._selis = np.argmin(np.abs(self._aexp_opts[:, np.newaxis]
                                       - aexps[np.newaxis, :]), axis=0)
        self._selis = np.unique(self._selis)
        self.aexps_used = self._aexp_opts[self._selis]
        print((f'selected files at aexp={self.aexps_used},'
               f' for target {aexps}'))
        self.filens_used = np.array(self._catopts)[self._selis]
        for (self._filen, self._aexp) in zip(self.filens_used, 
                                             self.aexps_used):
            self._get_smhmdist_smdpl(self._filen, self._aexp)
        del self._filen, self._aexp
    
    def _get_h5name(self, catfilen):
        return catfilen[:-4] \
               + f'_savedhist_smhm_binsize{self.binsize:.3f}_aexp.hdf5'
    
    def _write_hdf5(self, h5filen, maindict, headdoc):
        with h5py.File(h5filen, 'w') as f:
            for key in maindict:
                f.create_dataset(key, data=maindict[key])
            hed = f.create_group('Header')
            h5u.savedict_hdf5(hed, headdoc)

    def _get_smhmdist_smdpl(self, filen, a_used):
        self._filen_stored = self._get_h5name(filen)
        if os.path.isfile(self._filen_stored):
            self._get_smhmdist_smdpl_fromhdf5(self._filen_stored)
        else:
            self._get_smhmdist_smdpl_fromcat(filen, a_used)
    
    def _get_smhmdist_smdpl_fromcat(self, filen, a_used):
        self._halos = loadsmpdl(filen)
        self._z_used = 1. / a_used - 1.
        self._cosmopars =  cosmopars_smdpl.copy()
        self._cosmopars.update({'z': self._z_used, 'a': a_used})
        self._hsel = self._halos['upid'] == -1 # centrals
        self._mhalo_msun = self._halos['m'][self._hsel]  \
                           / self._cosmopars['h'] #BN98
        # for observed masses, would need to check the measurement method
        self._mstar_true_msun = self._halos['sm'][self._hsel] 
        self._logmh = np.log10(self._mhalo_msun)
        self._logms = np.log10(self._mstar_true_msun)
        del self._mstar_true_msun, self._mhalo_msun, self._hsel
        del self._halos

        self._minmh = np.min(self._logmh)
        self._maxmh = np.max(self._logmh)
        self._b0 = np.floor(self._minmh / self.binsize) * self.binsize
        self._b1 = np.ceil(self._maxmh / self.binsize) * self.binsize
        self._mhbins = np.arange(self._b0, self._b1 + 0.5 * self.binsize,
                                 self.binsize)
        self._minms = np.min(self._logms[np.isfinite(self._logms)])
        self._maxms = np.max(self._logms)
        self._b0 = np.floor(self._minms / self.binsize) * self.binsize
        self._b1 = np.ceil(self._maxms / self.binsize) * self.binsize
        self._msbins = np.arange(self._b0, self._b1 + 0.5 * self.binsize, 
                                 self.binsize)
        self._msbins = np.append([-np.inf], self._msbins)
        self._msbins_nozero = np.copy(self._msbins)
        self._msbins_nozero[0] = self._msbins_nozero[1] - self.binsize
        del self._minmh, self._maxmh, self._b0, self._b1
        del self._minms, self._maxms

        self._chunksize = 1_000_000
        self._arlen = len(self._logmh)
        self._nchunks = (self._arlen - 1) // self._chunksize + 1
        # bus error without feeding smaller bits of array into histogramdd
        for i in range(self._nchunks):
            sel = slice(self._chunksize * i, self._chunksize * (i + 1), None)
            self._shist, _ = np.histogramdd([self._logmh[sel], 
                                             self._logms[sel]],
                                             bins=[self._mhbins, 
                                                   self._msbins])
            if i == 0:
                self._hist = self._shist
            else:
                self._hist += self._shist
        del self._chunksize, self._arlen, self._nchunks, self._shist
        self.mhbins[self._z_used] = self._mhbins
        self.msbins[self._z_used] = self._msbins
        self.msbins_nozero[self._z_used] = self._msbins_nozero
        self.hists[self._z_used] = self._hist
        # save hist (considerably smaller dataset than catalog)
        self._h5filen = self._get_h5name(filen)
        self._maindict = {'mhbins': self._mhbins,
                          'msbins': self._msbins,
                          'msbins_nozero': self._msbins_nozero,
                          'hist_mh_ms': self._hist}
        self._headdoc = {'cosmopars': self._cosmopars,
                         'catalog_file': filen,
                         'halo_selection': 'upid == -1 (central galaxies)',
                         'units_msbins': 'log10 Msun',
                         'units_mhbins': 'log10 Msun',
                         'units_hist': 'galaxy count',
                         'parent_sim': 'smdpl',
                         'info_mh': 'BN98 definition',
                         'info_ms': '"true" stellar mass',
                         'hist_axis0': 'mhbins',
                         'hist_axis1': 'msbins'}
        self._write_hdf5(self._h5filen, self._maindict, self._headdoc)
        del self._z_used, self._cosmopars, self._mhbins
        del self._msbins, self._msbins_nozero, self._hist
        del self._h5filen, self._maindict, self._headdoc

    def _get_smhmdist_smdpl_fromhdf5(self, filen_stored):
        with h5py.File(filen_stored, 'r') as f:
            self._mhbins = f['mhbins'][:]
            self._msbins = f['msbins'][:]
            self._msbins_nozero = f['msbins_nozero'][:]
            self._hist = f['hist_mh_ms'][:] 
            self._z_used = f['Header/cosmopars_dict'].attrs['z']
        self.mhbins[self._z_used] = self._mhbins
        self.msbins[self._z_used] = self._msbins
        self.msbins_nozero[self._z_used] = self._msbins_nozero
        self.hists[self._z_used] = self._hist
        del self._z_used, self._mhbins, self._msbins, 
        del self._msbins_nozero, self._hist


# modified from the SMHMhists class to get the SFR distribution 
# at a given halo mass
class SFRHMhists:
    def __init__(self, ztargets, histdata='UM-smdpl', binsize=0.1):
        '''
        Parameters:
        -----------
        Class reads in halo catalogue files, calculates histograms, and
        saves the histograms. On later calls, the stored histograms are
        retrieved. 
        Methods allow for extraction of percentiles of SFR at a given
        halo mass.

        ztargets: array-like of floats
            redshifts to load. closest options will be used; these are
            then the redshift options for which the closest option is 
            chosen in e.g., median retrieval
        histdata: {'UM-smdpl'}
            which halo catalog to use.
        binsize: float
            the size of the halo and star formation rate bins, in 
            log10 Msun and log10 Msun/yr, respectively.
            The default is 0.1.
        '''
        self.ztargets = np.array(ztargets)
        self.atargets = 1. / (1. + self.ztargets)
        self.histdata_type = histdata
        self.binsize = binsize

        self.hists = {}
        self.mhbins = {}
        self.sfrbins = {}
        self.sfrbins_nozero = {}
        if self.histdata_type == 'UM-smdpl':
            self._finddata_smdpl(self.atargets)
        self.zs_used = np.sort(np.array(list(self.hists.keys())))
    
    def getperc_sfrmh(self, z_target, mode='mhtosfr', percvals=np.array([0.5])):
        '''
        retrieve the halo mass percentiles at a given star formation 
        rate, or star formation rate percentiles at a given halo mass.

        Parameters:
        -----------
        z_target: float
            redshift for which to get the Ms-Mh relation.
            Just use the closest loaded redshift.
        mode: {'sfrtomh', 'mhtosfr'}
            'sfrtomh': percentiles of halo mass at fixed SFR
            'mhtosfr': percentiles of SFR at fixed halo mass
        percvals: array-like of floats in the 0 -- 1 range
            which percentiles to calculate
        
        Returns:
        --------
        sfrvals: array of floats
            SFR values. if getting SFR percentiles, 
            shape is (number of percentiles, number of mhvals).
            Units are log10 Msun / yr.
        mhvals:array of floats
            halo mass values. if getting mh percentiles, 
            shape is (number of percentiles, number of sfrvals).
            Units are log10 Msun.
        '''
        self._usez = self.zs_used[np.argmin(np.abs(z_target - self.zs_used))]
        print(f'Using redshift {self._usez} for target {z_target}')
        self._sfrbins = self.sfrbins[self._usez]
        self._sfrbins_nozero = self.sfrbins_nozero[self._usez]
        self._mhbins = self.mhbins[self._usez]
        self._hist = self.hists[self._usez]
        if mode == 'sfrtomh':
            mhv = mu.percentiles_from_histogram(self._hist, self._mhbins,
                                                axis=0, percentiles=percvals)
            msv = 0.5 * (self._sfrbins_nozero[:-1] + self._sfrbins_nozero[1:])
        if mode == 'mhtosfr':
            msv = mu.percentiles_from_histogram(self._hist, self._sfrbins,
                                                axis=1, percentiles=percvals)
            mhv = 0.5 * (self._mhbins[:-1] + self._mhbins[1:])
        del self._sfrbins, self._sfrbins_nozero, self._mhbins, self._hist
        del self._usez
        return msv, mhv
    
    def getbins(self, z_target, mode='sfr'):
        '''
        useful to calculate the probability distribution for an input
        halo/stellar mass over the right ms/mh bins

        Parameters:
        -----------
        z_target: float
            redshift for which to get the bins.
            Just uses the closest loaded redshift.
        mode: {'mh', 'sfr'}
            'mh': get halo mass bins
            'sfr': get star formation rate mass bins 
                  (-np.inf bin set to a finite value)

        Returns:
        --------
        bins: array of floats
            bin edges, in log10 Msun or log10 Msun / yr
        '''
        self._usez = self.zs_used[np.argmin(np.abs(z_target - self.zs_used))]
        if mode == 'mh':
            out = self.mhbins[self._usez]
        elif mode == 'sfr':
            out = self.sfrbins_nozero[self._usez]
        del self._usez
        return out
    
    def gethist(self, z_target):
        '''
        Parameters:
        -----------
        z_target: float
            redshift for which to get the bins.
            Just uses the closest loaded redshift.

        Returns:
        --------
        histogram: 2d-array of floats
            halo/central galaxy counts as a function of halo mass
            (axis 0) and star formation rate (axis 1) 
        '''
        self._usez = self.zs_used[np.argmin(np.abs(z_target - self.zs_used))]
        out = self.hists[self._usez]
        del self._usez
        return out

    def _finddata_smdpl(self, aexps):
        self._catopts = glob.glob(ddir_smdpl + 'sfr_catalog_*.bin')
        self._filetrunks = [(catopt.split('/')[-1])[:-4]
                            for catopt in self._catopts]
        self._aexp_opts = np.array([float((trunk).split('_')[-1])
                                    for trunk in self._filetrunks])
        self._selis = np.argmin(np.abs(self._aexp_opts[:, np.newaxis]
                                       - aexps[np.newaxis, :]), axis=0)
        self._selis = np.unique(self._selis)
        self.aexps_used = self._aexp_opts[self._selis]
        print((f'selected files at aexp={self.aexps_used},'
               f' for target {aexps}'))
        self.filens_used = np.array(self._catopts)[self._selis]
        for (self._filen, self._aexp) in zip(self.filens_used, 
                                             self.aexps_used):
            self._get_sfrhmdist_smdpl(self._filen, self._aexp)
        del self._filen, self._aexp
    
    def _get_h5name(self, catfilen):
        return catfilen[:-4] \
               + f'_savedhist_sfrhm_binsize{self.binsize:.3f}_aexp.hdf5'
    
    def _write_hdf5(self, h5filen, maindict, headdoc):
        with h5py.File(h5filen, 'w') as f:
            for key in maindict:
                f.create_dataset(key, data=maindict[key])
            hed = f.create_group('Header')
            h5u.savedict_hdf5(hed, headdoc)

    def _get_sfrhmdist_smdpl(self, filen, a_used):
        self._filen_stored = self._get_h5name(filen)
        if os.path.isfile(self._filen_stored):
            self._get_sfrhmdist_smdpl_fromhdf5(self._filen_stored)
        else:
            self._get_sfrhmdist_smdpl_fromcat(filen, a_used)
    
    def _get_sfrhmdist_smdpl_fromcat(self, filen, a_used):
        '''
        get SFR - Mh histogram from the UniverseMachine binaries
        '''
        self._halos = loadsmpdl(filen)
        self._z_used = 1. / a_used - 1.
        self._cosmopars =  cosmopars_smdpl.copy()
        self._cosmopars.update({'z': self._z_used, 'a': a_used})
        self._hsel = self._halos['upid'] == -1 # centrals
        self._mhalo_msun = self._halos['m'][self._hsel]  \
                           / self._cosmopars['h'] #BN98, stored Msun/h
        # for observed masses, would need to check the measurement method
        # stored Msun/yr
        self._sfr_true_msunpyr = self._halos['sfr'][self._hsel] 
        self._logmh = np.log10(self._mhalo_msun)
        self._logsfr = np.log10(self._sfr_true_msunpyr)
        del self._sfr_true_msunpyr, self._mhalo_msun, self._hsel
        del self._halos

        self._minmh = np.min(self._logmh)
        self._maxmh = np.max(self._logmh)
        self._b0 = np.floor(self._minmh / self.binsize) * self.binsize
        self._b1 = np.ceil(self._maxmh / self.binsize) * self.binsize
        self._mhbins = np.arange(self._b0, self._b1 + 0.5 * self.binsize,
                                 self.binsize)
        self._minsfr = np.min(self._logsfr[np.isfinite(self._logsfr)])
        self._maxsfr = np.max(self._logsfr)
        self._b0 = np.floor(self._minsfr / self.binsize) * self.binsize
        self._b1 = np.ceil(self._maxsfr / self.binsize) * self.binsize
        self._sfrbins = np.arange(self._b0, self._b1 + 0.5 * self.binsize,
                                  self.binsize)
        self._sfrbins = np.append([-np.inf], self._sfrbins)
        self._sfrbins_nozero = np.copy(self._sfrbins)
        self._sfrbins_nozero[0] = self._sfrbins_nozero[1] - self.binsize
        del self._minmh, self._maxmh, self._b0, self._b1
        del self._minsfr, self._maxsfr

        self._chunksize = 1_000_000
        self._arlen = len(self._logmh)
        self._nchunks = (self._arlen - 1) // self._chunksize + 1
        # bus error without feeding smaller bits of array into histogramdd
        for i in range(self._nchunks):
            sel = slice(self._chunksize * i, self._chunksize * (i + 1), None)
            self._shist, _ = np.histogramdd([self._logmh[sel], 
                                             self._logsfr[sel]],
                                             bins=[self._mhbins, 
                                                   self._sfrbins])
            if i == 0:
                self._hist = self._shist
            else:
                self._hist += self._shist
        del self._chunksize, self._arlen, self._nchunks, self._shist
        self.mhbins[self._z_used] = self._mhbins
        self.sfrbins[self._z_used] = self._sfrbins
        self.sfrbins_nozero[self._z_used] = self._sfrbins_nozero
        self.hists[self._z_used] = self._hist
        # save hist (considerably smaller dataset than catalog)
        self._h5filen = self._get_h5name(filen)
        self._maindict = {'mhbins': self._mhbins,
                          'sfrbins': self._sfrbins,
                          'sfrbins_nozero': self._sfrbins_nozero,
                          'hist_mh_sfr': self._hist}
        self._headdoc = {'cosmopars': self._cosmopars,
                         'catalog_file': filen,
                         'halo_selection': 'upid == -1 (central galaxies)',
                         'units_sfrbins': 'log10 Msun / yr',
                         'units_mhbins': 'log10 Msun',
                         'units_hist': 'galaxy count',
                         'parent_sim': 'smdpl',
                         'info_mh': 'BN98 definition',
                         'info_sfr': '"true" star formation rate',
                         'hist_axis0': 'mhbins',
                         'hist_axis1': 'sfrbins'}
        self._write_hdf5(self._h5filen, self._maindict, self._headdoc)
        del self._z_used, self._cosmopars, self._mhbins
        del self._sfrbins, self._sfrbins_nozero, self._hist
        del self._h5filen, self._maindict, self._headdoc

    def _get_sfrhmdist_smdpl_fromhdf5(self, filen_stored):
        with h5py.File(filen_stored, 'r') as f:
            self._mhbins = f['mhbins'][:]
            self._sfrbins = f['sfrbins'][:]
            self._sfrbins_nozero = f['sfrbins_nozero'][:]
            self._hist = f['hist_mh_sfr'][:] 
            self._z_used = f['Header/cosmopars_dict'].attrs['z']
        self.mhbins[self._z_used] = self._mhbins
        self.sfrbins[self._z_used] = self._sfrbins
        self.sfrbins_nozero[self._z_used] = self._sfrbins_nozero
        self.hists[self._z_used] = self._hist
        del self._z_used, self._mhbins, self._sfrbins, 
        del self._sfrbins_nozero, self._hist

    
    