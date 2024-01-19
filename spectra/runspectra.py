'''
should make a separate repo eventually, but convenient to have it with
the halo property finders for now

TODO: 
x check how ion fractions in ray cells are calculated:
    - interpolated to grid first or ion fractions calculated first
      -> no grid interpolation anywhere
    - interpolation method
      -> ion density * mass / density * kernel integral over l.o.s.
      (ion density is projected onto the l.o.s. from SPH kernel 
       integration. In practice, the integration is done for a 
       limited number of normalized impact parameters, and 
       interpolated.)
    - ion-weighted temperature, velocity or something else used to
      deposit lines onto the spectrum
      -> values per particle are used directly
o check how the velocities in the spectra are centered:
  -> total redshift (for default observing redshift of 0):
     snapshot redshift, hubble flow, peculiar velocity
     from trident/light_ray.py:
     sub_data[('gas', 'redshift')] = my_segment['redshift'] - \
              (sub_data[('gas', 'l')] / ray_length) * \
              (my_segment['redshift'] - next_redshift)
     I think next_redshift comes from 
     yt_astro_analysis/cosmological_observation/cosmology_splice.py
     _deltz_forward: 'Calculate deltaz corresponding to moving a 
          comoving distance starting from some redshift'
     for simple rays.
     redshift_eff includes velocity effect if not turned off.
     velocity/lambda calculations convert lambda to velocity using
     c * (lambda_obs - lambda_0) / lambda_0, i.e., doppler velocity
     that would give the same redshift.
'''

import h5py
import numpy as np
import trident
import unyt
import yt

import fire_an.mainfunc.haloprop as hp
import fire_an.readfire.readin_fire_data as rfd
import fire_an.simlists as sl

savedir_spectra = '/projects/b1026/nastasha/spectra/'

def getinds_ax(axis):
    if axis == 'z':
        axes = (0, 1, 2)
    elif axis == 'x':
        axes = (1, 2, 0)
    elif axis == 'y':
        axes = (2, 0, 1)
    else:
        raise ValueError(f'axis should be "x", "y", or "z"; was {axis}')
    return axes

def getytds(simname, snapnum):
    simpath = sl.dirpath_from_simname(simname)
    snap = rfd.get_Firesnap(simpath, snapnum)
    # for some reason, YTarrays fail as start/end points for rays
    codelength_cm = snap.units.getunits('PartType0/Coordinates')
    simfile = snap.firstfilen
    ds = yt.load(simfile)
    return ds, simpath, codelength_cm

def getsightlines_grid(cen, totlength=unyt.unyt_quantity(500., 'kpc'),
                       gridside=unyt.unyt_quantity(500., 'kpc'),
                       gridpoints_side=10, axis='z'):
    a1, a2, a3 = getinds_ax(axis) # 'rotated' x, y, z axes
    losmin = cen[a3] - 0.5 * totlength
    losmax = cen[a3] + 0.5 * totlength
    gridx = cen[a1] - 0.5 * gridside \
            + gridside * np.linspace(0., 1., gridpoints_side)
    gridy = cen[a2] - 0.5 * gridside \
            + gridside * np.linspace(0., 1., gridpoints_side)
    start_positions = unyt.unyt_array(np.zeros((gridpoints_side**2, 3),
                                               dtype=cen.dtype),
                                      losmin.units)
    start_positions[:, a3] = losmin
    start_positions.reshape(gridpoints_side, gridpoints_side, 3)[:, :, a1] \
        = gridx[:, np.newaxis]
    start_positions.reshape(gridpoints_side, gridpoints_side, 3)[:, :, a2] \
        = gridy[np.newaxis, :]
    end_positions = unyt.unyt_array(np.zeros((gridpoints_side**2, 3),
                                             dtype=cen.dtype),
                                    losmin.units)
    end_positions[:, a3] = losmax
    end_positions.reshape(gridpoints_side, gridpoints_side, 3)[:, :, a1] \
        = gridx[:, np.newaxis]
    end_positions.reshape(gridpoints_side, gridpoints_side, 3)[:, :, a2] \
        = gridy[np.newaxis, :]
    return start_positions, end_positions
    
def runsightlines(simname, snapnum, outname_base=None,
                  settype='grid', **setargs):
    '''
    seems to work somewhat

    Parameters:
    -----------
    setargs: dict
        arguments passed to the getsightlines_<settype> function,
        with some modifications:
        for settype 'grid':
            'totlength' is muliplied by the halo virial radius (BN98)
            'gridside' is multiplied by the same
    
    '''
    ds, simpath, codelength_cm = getytds(simname, snapnum)
    halodat, todoc_halo = hp.readhalodata_shrinkingsphere(simpath, snapnum,
                                                          meandef='BN98')
    cen = unyt.unyt_array([halodat['Xc_cm'], 
                           halodat['Yc_cm'],
                           halodat['Zc_cm']], 'cm')
    rvir = unyt.unyt_quantity(halodat['Rvir_cm'], 'cm')

    if settype == 'grid':
        setargs['totlength'] = setargs['totlength'] * rvir
        setargs['gridside'] = setargs['gridside'] * rvir
        start_positions, end_positions = getsightlines_grid(cen, **setargs)
    
    if outname_base is None:
        outname_base = f'tridentray_{simname}_{snapnum}'
    _outname = savedir_spectra + outname_base

    trident.add_ion_fields(ds, ions=['Ne VIII', 'O VI', 'H I'])
    
    for ri, (_spos, _epos) in enumerate(zip(start_positions, end_positions)):
        outname = _outname + f'_{ri}.h5'
        # YTarray values don't work in make_simple_ray for some reason
        spos = np.array(_spos.to('cm')) / codelength_cm
        epos = np.array(_epos.to('cm')) / codelength_cm
        ray = trident.make_simple_ray(ds, start_position=spos, 
                                      end_position=epos,
                                      data_filename=outname,
                                      lines=['Ne VIII', 'O VI', 'H I'],
                                      fields=[('gas', 'temperature'), 
                                              ('gas', 'metallicity'),
                                              ('gas', 'density')])
        sg = trident.spectrum_generator.SpectrumGenerator(
            lambda_min='auto', lambda_max='auto', dlambda=2., 
            bin_space='velocity')
        sg.make_spectrum(ray, lines='Ne VIII 770')
        sg.save_spectrum(outname[:-3] + '.txt')


