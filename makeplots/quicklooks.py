#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
plot or do quick checks on results from fire_maps
'''

import numpy as np
import h5py
import pandas as pd
import os

import matplotlib.collections as mcol
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsp

import fire_an.makeplots.plot_utils as pu
import fire_an.utils.constants_and_units as c


def quicklook_massmap(filen, savename=None, mincol=None):
    '''
    quick plot of the mass map in the file
    '''

    with h5py.File(filen, 'r') as f:
        map = f['map'][:]
        vmin = f['map'].attrs['minfinite']
        vmax = f['map'].attrs['max']

        box_cm = f['Header/inputpars'].attrs['diameter_used_cm']
        cosmopars = {key: val for key, val in \
                     f['Header/inputpars/cosmopars'].attrs.items()}
        #print(cosmopars)
        if 'Rvir_ckpcoverh' in f['Header/inputpars/halodata'].attrs:
            rvir_ckpcoverh = f['Header/inputpars/halodata'].attrs['Rvir_ckpcoverh']
            rvir_pkpc = rvir_ckpcoverh * cosmopars['a'] / cosmopars['h']
        elif 'Rvir_cm' in f['Header/inputpars/halodata'].attrs:
            rvir_cm = f['Header/inputpars/halodata'].attrs['Rvir_cm']
            rvir_pkpc = rvir_cm / (c.cm_per_mpc * 1e-3)
        xax = f['Header/inputpars'].attrs['Axis1']
        yax = f['Header/inputpars'].attrs['Axis2']
        box_pkpc = box_cm / (1e-3 * c.cm_per_mpc)
        extent = (-0.5 * box_pkpc[xax], 0.5 * box_pkpc[xax],
                  -0.5 * box_pkpc[yax], 0.5 * box_pkpc[yax])

    if mincol is None:
        cmap = 'viridis'
    else:
        cmap = pu.paste_cmaps(['gist_yarg', 'viridis'], [vmin, mincol, vmax])
    extend = 'neither' if np.min(map) == vmin else 'min'

    fig = plt.figure(figsize=(5.5, 5.))
    grid = gsp.GridSpec(nrows=1, ncols=2, hspace=0.1, wspace=0.0, 
                        width_ratios=[10., 1.])
    ax = fig.add_subplot(grid[0, 0]) 
    cax = fig.add_subplot(grid[0, 1])
    fontsize = 12
    
    xlabel = ['X', 'Y', 'Z'][xax] + ' [pkpc]'
    ylabel = ['X', 'Y', 'Z'][yax] + ' [pkpc]'
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_title('Mass map centered on halo')

    img = ax.imshow(map.T, origin='lower', interpolation='nearest', vmin=vmin,
                    vmax=vmax, cmap=cmap, extent=extent)
    ax.tick_params(axis='both', labelsize=fontsize-1)
    plt.colorbar(img, cax=cax, extend=extend, orientation='vertical') 

    patches = [mpatch.Circle((0., 0.), rvir_pkpc)]
    collection = mcol.PatchCollection(patches)
    collection.set(edgecolor=['red'], facecolor='none', linewidth=1.5)
    ax.add_collection(collection)
    ax.text(1.05 * 2**-0.5 * rvir_pkpc, 1.05 * 2**-0.5 * rvir_pkpc, 
            '$R_{\\mathrm{vir}}$',
            color='red', fontsize=fontsize)
    
    if savename is not None:
        plt.savefig(savename, bbox_inches='tight')

def readfile_outputtimes(path):
    if path.endswith('output'):
        path = path[:-6]
    if not path.endswith('/'):
        path = path + '/'
    targetfile = path + 'snapshot_scale-factors.txt'
    if not os.path.isfile(targetfile):
        raise RuntimeError('No file {} found'.format(targetfile))
    with open(targetfile, 'r') as f:
        aopts = f.read()
    aopts = (aopts.strip()).split('\n')
    aopts = np.array([float(aopt) for aopt in aopts])
    zopts = 1. / aopts - 1.
    return zopts
        
def plotsnaps_m13noBH():
    #basedir = '/scratch3/01799/phopkins/fire3_suite_done/'
    basedir = '/Users/nastasha/ciera/fire_data/'
    
    # noBH m13s from Lindsey's spreadsheet 
    checkpaths = ['m13h002_m3e5/m13h002_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                  'm13h007_m3e5/m13h007_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                  'm13h029_m3e5/m13h029_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                  'm13h113_m3e5/m13h113_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                  'm13h206_m3e5/m13h206_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                  'm13h206_m3e5/m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct142021_crdiffc1_sdp1e10_gacc31_fa0.5',
                  'm13h206_m3e5/m13h206_m3e5_MHD_fire3_fireBH_Sep052021_crdiffc690_sdp1e10_gacc31_fa0.5',
                  'm13h206_m3e5/m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_sdp1e10_gacc31_fa0.5',
                  'm13h206_m3e5/m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct142021_crdiffc1_sdp1e10_gacc31_fa0.5',
                  'm13h206_m3e5/m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e10_gacc31_fa0.5_fcr3e-3_vw3000',
                  'm13h223_m3e5/m13h223_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                  'm13h236_m3e5/m13h236_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5',
                  ]
    data = {}
    vals = []
    title = 'snapshot redshifts: m13, m3e5, fire3_fireBH, sdp1e10_gacc31_fa0.5'
    for path in checkpaths:
        val = readfile_outputtimes(basedir + path)
        key = path.split('/')[-1]
        # shorten
        key = key.replace('m13', '')
        key = key.replace('_m3e5', '')
        key = key.replace('_fire3_fireBH', '')
        key = key.replace('_sdp1e10_gacc31_fa0.5', '')
        if len(key) > 30:
            tot = len(key)
            opts = np.where([char == '_' for char in key])[0]
            splitopt = np.argmin(np.abs(opts - 0.5 * tot))
            splitpoint = opts[splitopt]
            key = key[:splitpoint] + '\n' + key[splitpoint:]
        data.update({key: val})
        vals.append(val)
    
    commonvals = [val if np.all([val in vals[i] for i in range(len(vals))])\
                  else None for val in vals[0]]
    normalvals = [val if np.sum([val in vals[i] for i in range(len(vals))]) \
                         > 0.5 * len(vals)\
                  else None for val in vals[3]]
    while None in commonvals:
        commonvals.remove(None)
    while None in normalvals:
        normalvals.remove(None)
    
    fontsize = 12
    fig = plt.figure(figsize=(10., 5))
    ax = fig.add_axes([0.33, 0.1, 0.62, 0.8])
    fig.suptitle(title, fontsize=fontsize)
    keys = list(data.keys())
    keys.sort()
    yvals = np.arange(len(keys)) + 0.5
    xzeroval = 0.01
    for key, yval in zip(keys, yvals):
        xvals = data[key]
        if xvals[0] <= 0.:
            xvals[0] = xzeroval
        xall = [val in normalvals for val in xvals]
        colors = np.ones((len(xvals), 4)) * np.array([0.0, 0.0, 0.0, 1.0])
        colors[np.logical_not(xall)] = np.array([0.0, 1.0, 0.0, 0.4])
        ax.scatter(xvals, [yval] * len(xvals), c=colors, marker='o',
                   s=15)
        oddvals = np.where(np.logical_not(xall))[0]
        for ind in oddvals:
            xv = xvals[ind]
            if ind > 0:
                if xv == xvals[ind - 1]:
                    st = st + '={}'.format(int)
                else:
                    st = '{}'.format(ind)
            else:   
                st = '{}'.format(ind)
            ax.text(xv, yval, st,
                    horizontalalignment='left' if ind % 4 > 1 else 'right',
                    verticalalignment='bottom' if ind % 2 else 'top')
        numsnaps = len(xvals)
        ax.text(11., yval, '({})'.format(numsnaps),
                horizontalalignment='left', verticalalignment='center')
    for xv in commonvals:
        ax.axvline(xv, color='gray', alpha=0.3)
    ax.set_yticks(yvals, labels=keys)
    ax.tick_params(left=True, bottom=True, labelsize=fontsize - 1,
                   direction='in', which='both')
    ax.set_xlabel('redshift', fontsize=fontsize)
    ax.text(xzeroval, yvals[0], '$z=0$', fontsize=fontsize - 1,
            horizontalalignment='left', verticalalignment='center')
    ax.set_xlim(0.8 * xzeroval, 17.)
    ax.set_xscale('log')

    plt.savefig('/Users/nastasha/ciera/tests/fire_start/' + \
                'snaptimes_m13_noBH_sample_Lindsey.pdf')
