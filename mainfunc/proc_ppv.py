'''
process ppv data cubes into maps
'''
import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c

def smoothmax_ppv(filen, vax=2, p1ax=0, p2ax=1, 
                  smoothsigmas=(10.e5, 20.e5, 30.e5, 40.e5, 
                                50.e5, 60.e5, 70.e5, 80.e5, 90.e5, 100.e5)):
    outfilen = filen[:-5] + '_smoothed_vmaxcols.hdf5'
    with h5py.File(filen, 'r') as fi:
        hist = fi['histogram/histogram'][:]
        islog = bool(fi['histogram/histogram'].attrs['log'])
        if islog:
            hist = 10**hist
        keepaxes = {p1ax, p2ax, vax}
        allaxes = set(np.arange(len(hist.shape)))
        sumaxes = tuple(allaxes - keepaxes)
        hist = np.sum(hist, sumaxes)
        # same order as the originals, but any other axes were removed
        _p1ax, _p2ax, _vax  = np.argsort([p1ax, p2ax, vax])
        pbins1 = fi[f'Header/axis_{_p1ax}/bins'][:]
        pbins2 = fi[f'Header/axis_{_p2ax}/bins'][:]
        vbins = fi[f'Header/axis_{_vax}/bins'][:]
    
        ctot = np.sum(hist, ax=vax)
        vcens = 0.5 * (vbins[:-1] + vbins[1:])
        dv = np.average(np.diff(vbins))
        if not np.allclose(np.diff(vbins), dv):
            raise RuntimeError(f'v bins are not evenly spaced: {vbins}')
        maxi_nosmooth = np.argmax(hist, ax=vax)
        vcomp_nosmooth = vcens[maxi_nosmooth]

        with h5py.File(outfilen, 'a') as fo:
            fi.copy(fi['Header'], fo, 'Header_ppv')
            hed = fo.create_group('Header')
            hed.attrs.create('vax', vax)
            hed.attrs.create('p1ax', p1ax)
            hed.attrs.create('p2ax', p2ax)
            gtot = fo.create_dataset('map_pp', data=ctot)
            gtot.attrs.create('log', False)
            fo.create_dataset('v_maxcol_nosmooth', data=vcomp_nosmooth)

            for smoothsigma in smoothsigmas:
                nsig_calc = 3.
                vmin = np.floor(- nsig_calc * smoothsigma / dv) * dv
                vmax = np.ceil(nsig_calc * smoothsigma / dv) * dv
                vbins_conv = np.arange(vmin - 0.5 * dv, vmax + dv, dv)
                gaussv = sp.special.erf(vbins_conv[1:] / smoothsigma) \
                         - sp.special.erf(vbins_conv[:-1] / smoothsigma) 
                gaussv *= 1. / np.sum(gaussv)
                sel = [None for _ in range(3)]
                sel[_vax] = slice(None, None, None)
                sel = tuple(sel)
                smoothedv = sp.signal.convolve(hist, gaussv[sel], mode='full')
                nextra = len(gaussv) // 2 
                vextlo = np.arange(vcens[0] - nextra * dv, 
                                   vcens[0] - 0.5 * dv, dv)
                vexthi = np.arange(vcens[-1] + dv, 
                                   vcens[-1] + (nextra + 0.5) * dv, dv)
                vcens_smoothedv = np.append(vextlo, vcens)
                vcens_smoothedv = np.append(vcens_smoothedv, vexthi)
                _maxi = np.argmax(smoothedv)
                _vcomp = vcens_smoothedv[_maxi]
                _ds = fo.create_dataset(f'v_maxcol_smooth_{smoothsigma:.0f}',
                                        data=_vcomp)
                _ds.attrs.create('sigma_v_smooth_cmps', smoothsigma)

                # initial testing
                dp2 = np.average(np.diff(pbins1)) * np.average(np.diff(pbins2))
                dp2 *= (c.cm_per_mpc * 1e-3)**2
                coldens = ctot / dp2
                selpoints = np.where(coldens >= 10**12.5)
                inds = np.random.choice(len(selpoints[0]), size=10, 
                                        replace=False)
                for ind in inds:
                    # assuming axes are p p v here
                    xp = vcens_smoothedv / 1e5
                    yp = smoothedv[selpoints[0][ind], selpoints[1][ind], :]
                    mi = _maxi[selpoints[0][ind], selpoints[1][ind]]
                    _xp = vcens / 1e5
                    _yp = hist[selpoints[0][ind], selpoints[1][ind], :]
                    _mi = maxi_nosmooth[selpoints[0][ind], selpoints[1][ind]]
                    plt.plot(xp, yp, color='black', linestyle='solid',
                             label=f'smoothed {smoothsigma/1e5:.0f} km/s')
                    plt.plot(_xp, _yp, color='gray', linestyle='dashed',
                             label=f'unsmoothed, {dv/1e5:.0f} km/s')
                    plt.scatter([xp[mi]], [yp[mi]], marker='o', color='black',
                                markersize=5)
                    plt.scatter([_xp[_mi]], [_yp[_mi]], marker='o', 
                                color='gray', markersize=5)
                    plt.set_yscale('log')
                    plt.set_xlabel('doppler v [km/s]')
                    plt.legend()
                plt.show()

def run_ppv_proc(opt):
    simnames = sl.m12_sr_all2 + sl.m12_hr_all2 +\
               sl.m13_sr_all2 + sl.m13_hr_all2
    sims_sr = sl.m12_sr_all2 +  sl.m13_sr_all2
    sims_hr = sl.m12_hr_all2 +  sl.m13_hr_all2
    snaps_sr = sl.snaps_sr
    snaps_hr = sl.snaps_hr
    numsnaps = 6
    for sn in sl.buglist1:
        if sn in simnames:
            simnames.remove(sn)
    if opt == -1:
        return simnames
    paxes = ['x', 'y', 'z']
    
    wtstr = 'Ne8'
    
    smi = opt // (numsnaps * len(paxes))
    sni = (opt % (numsnaps * len(paxes))) // len(paxes) 
    pai = opt % len(paxes) 

    simname = simnames[smi]
    snaps = snaps_sr if simname in sims_sr \
            else snaps_hr if simname in sims_hr \
            else None
    snapnum = snaps[sni]
    pax = paxes[pai]
    
    ddir = ''
    filen = ddir + (f'hist_ppv_{pax}ax_by_{wtstr}_{simname}_snap{snapnum}'
                    '_bins1_v1_hvcen.hdf5')
    smooths = 1e5 * np.arange(10., 105., 10.)
    smoothmax_ppv(filen, vax=2, p1ax=0, p2ax=1, 
                  smoothsigmas=smooths)
    
    









