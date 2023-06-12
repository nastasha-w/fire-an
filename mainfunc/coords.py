
import numpy as np

import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu

class CoordinateWranger:
    def __init__(self, snapobj, center_cm, rotmatrix=None,
                 parttype=0, periodic=True, vcen_cmps=None,
                 filterdct=None):
        '''
        class to get position and velocity info in different coordinate
        bases

        Parameters:
        -----------
        snapobj: Firesnap or similar
            object that allows access to cosmological parameters and 
            has a method to read in simulation arrays
        center_cm: float array, shape (3,)
            center coordinates in cm (physical, not comoving)
        rotmatrix: float array, shape (3, 3) or None
            matrix by which to multiply coordinates to get the 
            coordinates in the desired basis
            None means no rotation
        parttype: int
            which particles to get coordinates for. Matches the 
            PartType<number> groups in the simulation outputs.
        periodic: bool
            Do we need to care about coordinates wrapping around the 
            simulation volume?
        vcen_cmps: float array, shape (3,)
            bulk velocity (subtracted from simulation velocities before
            any rotations, etc.), in cm/s
        filterdct: dict with a 'filter' key or None
            the 'filter' key value should be a boolean array with 
            size matching the number of particles of the requested 
            type. Coordinate values are only returned for the 
            particles with True values.
            The default is None. (Use all particles of specified type.)
        Note:
        -----
        These objects store the arrays they used, so it's best to 
        delete them once you've got the data you want.
        '''
        self.snapobj = snapobj
        self.cen_cm = center_cm
        self.vcen_cmps = vcen_cmps
        self.rotmatrix = rotmatrix
        self.pt = parttype
        self.periodic = periodic
        self.coordaxis = 1 
        self.cosmopars = self.snapobj.cosmopars.getdct()
        self.pcalcstarted = False
        self.vcalcstarted = False
        self.__check_rotmatrix()
        self.filter = slice(None, None, None)
        if filterdct is not None:
            if 'filter' in filterdct:
                self.filter = filterdct['filter']
    
    def __check_rotmatrix(self):
        if self.rotmatrix is not None:
            if self.rotmatrix.shape != (3, 3):
                msg = ('Rotation matrix should have shape (3, 3) not input'
                       f'{self.rotmatrix.shape} for matrix\n{self.rotmatrix}')
                raise ValueError(msg)
            if not (np.allclose(np.matmul(self.rotmatrix.T, self.rotmatrix),
                                np.diag(np.ones((3,)))) 
                    and np.isclose(np.linalg.det(self.rotmatrix), 1.)):
                msg = ('input was not a valid rotation matrix.\n'
                       'transpose (should be inverse):\n'
                       f'{self.rotmatrix.T}, {self.rotmatrix}\n'
                       'determinant (should be 1.): '
                       f'{np.linalg.det(self.rotmatrix)}')
                raise ValueError(msg)
        
    def __startcalc_pos(self, subindex=None):
        '''
        read in, center, rotate coordinates
        '''
        if self.pcalcstarted:
            return None
        h5path = f'PartType{self.pt}/Coordinates'
        if self.rotmatrix is None:
            self._subindex = subindex
        else:
            self._subindex = None
        self.sel = (self.filter, slice(None, None, None)) \
                   if self._subindex is None else \
                   self.filter
        self.coords_simxyz = self.snapobj.readarray(
            h5path, subindex=self._subindex)[self.sel]
        del self.sel
        self.toCGS_coords_simxyz = self.snapobj.toCGS
        self.__center_pos()
        if self.rotmatrix is not None:
            self.__rotate_pos()
        else:
            self.coords_rotxyz = self.coords_simxyz
        self.toCGS_coords_rotxyz = self.toCGS_coords_simxyz
        if subindex is not None and self.rotmatrix is not None:
            self.coords_rotxyz = np.copy(self.coords_rotxyz[:, subindex])
        del self.coords_simxyz
        del self.toCGS_coords_simxyz
        self.pcalcstarted = True
        return None
    
    def __startcalc_vel(self, subindex=None):
        '''
        read in, center, rotate velocities
        '''
        if self.vcalcstarted:
            return None
        h5path = f'PartType{self.pt}/Velocities'
        if self.rotmatrix is None:
            self._subindex = subindex
        else:
            self._subindex = None
        self.sel = (self.filter, slice(None, None, None)) \
              if self._subindex is None else \
              self.filter
        self.vel_simxyz = self.snapobj.readarray(
            h5path, subindex=self._subindex)[self.sel]
        del self.sel
        self.toCGS_vel_simxyz = self.snapobj.toCGS
        self.__center_vel()
        if self.rotmatrix is not None:
            self.__rotate_vel()
        else:
            self.vel_rotxyz = self.vel_simxyz
        self.toCGS_vel_rotxyz = self.toCGS_vel_simxyz
        if subindex is not None and self.rotmatrix is not None:
            self.vel_rotxyz = np.copy(self.vel_rotxyz[:, subindex])
        del self.vel_simxyz
        del self.toCGS_vel_simxyz
        self.vcalcstarted = True
        return None

    def __rotate_pos(self):
        self.rotmatrix = np.asarray(self.rotmatrix, 
                                    dtype=self.coords_simxyz.dtype)
        self.coords_rotxyz = np.einsum('kj,ij->ik', self.rotmatrix,
                                       self.coords_simxyz)
        self.toCGS_coords_rotxyz = self.toCGS_coords_simxyz
    
    def __rotate_vel(self):
        self.rotmatrix = np.asarray(self.rotmatrix, 
                                    dtype=self.vel_simxyz.dtype)
        self.vel_rotxyz = np.einsum('kj,ij->ik', self.rotmatrix, 
                                    self.vel_simxyz)
        self.toCGS_vel_rotxyz = self.toCGS_vel_simxyz
    
    def __center_pos(self):
        self.center_simu = (self.cen_cm / self.toCGS_coords_simxyz).astype(
            self.coords_simxyz.dtype)
        self._sel = (np.newaxis, slice(None, None, None)) \
                    if self._subindex is None else \
                    self._subindex
        self.coords_simxyz -= self.center_simu[self._sel]
        del self._sel
        if self.periodic:
            self.boxsize_simu = self.snapobj.cosmopars.boxsize \
                                * self.snapobj.cosmopars.a \
                                / self.snapobj.cosmopars.h \
                                * c.cm_per_mpc / self.toCGS_coords_simxyz
            self.coords_simxyz += 0.5 * self.boxsize_simu
            self.coords_simxyz %= self.boxsize_simu
            self.coords_simxyz -= 0.5 * self.boxsize_simu
    
    def __center_vel(self):
        self.vcen_simu = (self.vcen_cmps / self.toCGS_vel_simxyz).astype(
            self.vel_simxyz.dtype)
        self._sel = (np.newaxis, slice(None, None, None)) \
                    if self._subindex is None else \
                    self._subindex
        self.vel_simxyz -= self.vcen_simu[self._sel]
        del self._sel
        if self.periodic:
            self.vboxsize_simu = self.snapobj.cosmopars.boxsize \
                                 * self.snapobj.cosmopars.a \
                                 / self.snapobj.cosmopars.h \
                                 * c.cm_per_mpc / self.toCGS_vel_simxyz \
                                 * cu.Hubble(cosmopars=self.cosmopars)
            self.vel_simxyz += 0.5 * self.vboxsize_simu
            self.vel_simxyz %= self.vboxsize_simu
            self.vel_simxyz -= 0.5 * self.vboxsize_simu
        
    def calccoords(self, coordspecs):
        '''
        calculate various coordinate values. Doing this all in one go
        should save some time from reading in large arrays multiple
        times.

        Parameters:
        -----------
        coordspecs: dict or list-like of dicts
            list: different coordinates to calculate
            dict or dicts in list: specify what to calculate
            dict keys and possible values:
                'pos': [0, 1, 2, 'allcart', 'rcen', 'phi', 'azimuth']
                    0, 1, 2: position along the axis with this index
                    'allcart': for all three of these cartesian axes
                    'rcen': distance to the center
                    'phi': angle w.r.t. the positive x-axis in the
                        x-y plane (after coordinate rotation)
                    'azimuth': angle w.r.t. the postitive z-axis
                        (after coordinate rotation)
                'vel': [0, 1, 2, 'allcart', 'vrad', 'vtot', 'phi',
                        'azimuth']
                     0, 1, 2: velocity along the axis with this index
                    'allcart': for all three of these cartesian axes
                    'vrad': radial velocity (relative to coordinate
                        center)
                    'vtot': total velocity (rms coordinate velocties)
                    'phi': rotational velocity about the z-axis
                        (after coordinate rotation)
                    'azimuth': in spherical coordinates, velocity in 
                        the azimuthal angle/theta direction (after 
                        coordinate rotation)
                    'alldop', 'dop0', 'dop1', 'dop2': doppler velocity,
                        i.e., peculiar velocity + hubble flow. Note 
                        that zero points for velocity and hubble flow 
                        are at the indicated position and velocity 
                        centers, and do not include the snapshot 
                        redshift, i.e., 
                        1 + z_obs = (1 + z_snap) * (1 + v_dop / c).
                        dopall ris for all three coordinates, dop<num>
                        is for the doppler velocity along the rotated 
                        axis with index <num>.
                    note: you must specifiy vcen_cmps when initializing
                    this object to calculate this. 
                indices etc. are all for the rotated coordinates, after
                centering 
        
        Returns:
        --------
        The desired coordinates in the listed order. Always returns a 
        list of 3-tuples: (coordinate [array], CGS conversion [float], 
                           doc_dictionary [e.g., used center]) 
        note that if for some reason a coordspec is requested twice, 
        the tuples will include the same object twice
        '''
        self.coordspecs_in = [[(key, dct[key]) for key in dct][0] 
                              if len(dct) == 1 else
                              ('coordspecs dicts should have length 1', dct)
                              for dct in coordspecs]
        ## this priority setting can get messy very fast if I try to
        ## implement too much here.
        # which (groups of) properties to calculate, and in what order
        # a group is calculated in a single function named 
        # __calc_<group key>
        self.calcorder = [('poscart', [('pos', 'allcart'), ('pos', 0),
                                       ('pos', 1), ('pos', 2)]),
                          ('poscen', [('pos', 'rcen')]),
                          ('posazi', [('pos', 'azimuth')]),
                          ('posphi', [('pos', 'phi')]),
                          ('velcart', [('vel', 'allcart'), ('vel', 0),
                                      ('vel', 1), ('vel', 2)]),
                          ('veldop', [('vel', 'alldop'), ('vel', 'dop0'),
                                      ('vel', 'dop1'), ('vel', 'dop2')]),
                          ('veltot', [('vel', 'vtot')]), 
                          ('velcen', [('vel', 'vrad')]),
                          ('velazi', [('vel', 'azimuth')]),
                          ('velphi', [('vel', 'phi')])
                         ]
        # what to get just because it's needed later
        # note: should include dependencies of dependencies
        self.dependencies = {('pos', 'rcen'): [('pos', 'allcart')],
                             ('pos', 'azimuth') : [('pos', 'allcart'),
                                                   ('pos', 'rcen')],
                             ('pos', 'phi') : [('pos', 'allcart')],
                             ('vel', 'vtot'): [('vel', 'allcart')],
                             ('vel', 'vrad'): [('vel', 'allcart'),
                                               ('pos', 'allcart'),
                                               ('pos', 'rcen')],
                             ('vel', 'alldop'): [('vel', 'allcart'),
                                                  ('pos', 'allcart')],
                             ('vel', 'dop0'): [('vel', 'allcart'),
                                                ('pos', 'allcart')],
                             ('vel', 'dop1'): [('vel', 'allcart'),
                                                ('pos', 'allcart')],
                             ('vel', 'dop2'): [('vel', 'allcart'),
                                                ('pos', 'allcart')],
                             ('vel', 'azimuth'): [('vel', 'allcart'),
                                                  ('pos', 'azimuth'),
                                                  ('pos', 'phi')],
                             ('vel', 'phi'): [('vel', 'allcart'),
                                              ('pos', 'phi')],
                            }
        # set up to-do list of everything that's needed (no duplicates)
        self._coords_todo = set(self.coordspecs_in.copy())
        for _coordspec in self.coordspecs_in:
            if _coordspec in self.dependencies:
                self._coords_todo |= set(self.dependencies[_coordspec])
        self.coords_todo = [[specs[i] for i in range(len(specs)) 
                             if specs[i] in self._coords_todo] 
                            for key, specs in self.calcorder]
        while [] in self.coords_todo:
            self.coords_todo.remove([])
        # holds arrays calculated for output 
        self.coords_outlist = [None for i in range(len(self.coordspecs_in))]
        # holds all calculated arrays, including those only needed as
        # dependencies. (keys are coordspecs tuples)
        self.coords_stored = {}
        print(f'calculations todo: {self.coords_todo}')
        for self.gicur, self.gcur in enumerate(self.coords_todo):
            self.gkeymatch = [key for key, specs in self.calcorder 
                              if set(self.gcur).issubset(
                                  set(specs))]
            if len(self.gkeymatch) == 0:
                msg = ('Invalid calccoords coordspec entry '
                      f'(tuple from dict):\n{self.gcur}')
                raise ValueError(msg)
            self.gkeymatch = self.gkeymatch[0]
            print(f'calculating {self.gcur}')
            if self.gkeymatch == 'poscart':
                self.__calc_poscart(self.gcur)
            elif self.gkeymatch == 'poscen':
                self.__calc_poscen(self.gcur)
            elif self.gkeymatch == 'posazi':
                self.__calc_posazimuth(self.gcur)
            elif self.gkeymatch == 'posphi':
                self.__calc_posphi(self.gcur)
            elif self.gkeymatch == 'velcart':
                self.__calc_velcart(self.gcur)
            elif self.gkeymatch == 'veldop':
                self.__calc_veldop(self.gcur)
            elif self.gkeymatch == 'veltot':
                self.__calc_veltot(self.gcur)
            elif self.gkeymatch == 'velcen':
                self.__calc_velrad(self.gcur)
            elif self.gkeymatch == 'velazi':
                self.__calc_velazimuth(self.gcur)
            elif self.gkeymatch == 'velphi':
                self.__calc_velphi(self.gcur)
            self.__update_out_todo(self.gcur)
            print(f'still todo: {self.still_todo}')
        del self.gcur, self.gkeymatch, self.coords_todo, 
        del self.coords_stored
        return self.coords_outlist
    
    def __calc_poscart(self, specs):
        if ('pos', 'allcart') in specs or len(specs) > 1:
            self.__startcalc_pos(subindex=None)
            for self.scur in specs:
                if self.scur == ('pos', 'allcart'):
                    self._todoc_cur = {'cen_cm': self.cen_cm,
                                       'rotmatrix': self.rotmatrix,
                                       'rotcoord_index': [0, 1, 2],
                                       'units': 'cm'}
                    self.coords_stored[self.scur] = (self.coords_rotxyz, 
                                                     self.toCGS_coords_rotxyz,
                                                     self._todoc_cur)
                else:
                    self._todoc_cur = {'cen_cm': self.cen_cm,
                                       'rotmatrix': self.rotmatrix,
                                       'rotcoord_index': self.scur[1],
                                       'units': 'cm'}
                    # storing a view of an array could cause unexpected
                    # side-effects
                    self._out = np.copy(self.coords_rotxyz[:, self.scur[1]])
                    self.coords_stored[self.scur] = (self._out, 
                                                     self.toCGS_coords_rotxyz,
                                                     self._todoc_cur)
                    del self._out
            del self.scur, self._todoc_cur
        else:
            self.__startcalc_pos(subindex=specs[0][1])
            self._todoc_cur = {'cen_cm': self.cen_cm,
                               'rotmatrix': self.rotmatrix,
                               'rotcoord_index': specs[0][1],
                               'units': 'cm'}
            self.coords_stored[specs[0]] = (self.coords_rotxyz, 
                                            self.toCGS_coords_rotxyz,
                                            self._todoc_cur)
            del self._todoc_cur

    def __calc_poscen(self, specs):
        self.scur = specs[0]
        self._in = self.coords_stored[('pos', 'allcart')]
        self._todoc_cur = self._in[2].copy()
        del self._todoc_cur['rotmatrix']
        del self._todoc_cur['rotcoord_index']
        self._out = np.sqrt(np.sum(self._in[0]**2, axis=self.coordaxis))
        self.coords_stored[self.scur] = (self._out, self._in[1], 
                                         self._todoc_cur)
        del self.scur, self._out, self._todoc_cur, self._in
    
    def __calc_posazimuth(self, specs):
        '''
        spherical coordinate theta: angle with the rotated z axis
        (axis index 2)
        range: -pi/2 -- pi/2 
        units: radians
        '''
        self.scur = specs[0]
        self._in = self.coords_stored[('pos', 'allcart')]
        self._todoc_cur = self._in[2].copy()
        del self._todoc_cur['rotcoord_index']
        self._todoc_cur['units'] = 'radians'
        self._todoc_cur['info'] = ('angle with rot. coord. z axis,'
                                   ' [-pi / 2, pi / 2]')
        self._r = self.coords_stored[('pos', 'rcen')]
        self._rat = self._in[0][:, 2] / self._r[0]
        # adjust if any unit differences
        if not np.isclose(self._in[1], self._r[1]):
            self._rat *= (self._in[1] / self._r[1])
        del self._r, self._in
        self._out = np.arccos(self._rat)
        self.coords_stored[self.scur] = (self._out, 1., self._todoc_cur)
        del self.scur, self._out, self._todoc_cur, self._rat

    def __calc_posphi(self, specs):
        '''
        spherical coordinate phi: angle with the rotated positive 
        x axis (axis index 0), in the rotated x-y plane (axis 
        indices 0 , 1)
        range: -pi -- pi 
        units: radians
        '''
        self.scur = specs[0]
        self._in = self.coords_stored[('pos', 'allcart')]
        self._todoc_cur = self._in[2].copy()
        del self._todoc_cur['rotcoord_index']
        self._todoc_cur['units'] = 'radians'
        self._todoc_cur['info'] = ('angle with rot. coord. +x axis,'
                                   ' in the rot. x-y plane [-pi, pi]')
        # arctan2 arguments are y, x!
        # 0., 0., -> 0.
        self._out = np.arctan2(self._in[0][:, 1], self._in[0][:, 0])
        self.coords_stored[self.scur] = (self._out, 1., self._todoc_cur)
        del self.scur, self._out, self._todoc_cur, self._in

    def __calc_velcart(self, specs):
        if ('vel', 'allcart') in specs or len(specs) > 1:
            self.__startcalc_vel(subindex=None)
            for self.scur in specs:
                if self.scur == ('vel', 'allcart'):
                    self._todoc_cur = {'vcen_cmps': self.vcen_cmps,
                                       'rotmatrix': self.rotmatrix,
                                       'rotcoord_index': [0, 1, 2],
                                       'units': 'cm * s**-1'}
                    self.coords_stored[self.scur] = (self.vel_rotxyz, 
                                                     self.toCGS_vel_rotxyz,
                                                     self._todoc_cur)
                else:
                    self._todoc_cur = {'vcen_cmps': self.vcen_cmps,
                                       'rotmatrix': self.rotmatrix,
                                       'rotcoord_index': self.scur[1],
                                       'units': 'cm * s**-1'}
                    # storing a view of an array could cause unexpected
                    # side-effects
                    self._out = np.copy(self.vel_rotxyz[:, self.scur[1]])
                    self.coords_stored[self.scur] = (self._out, 
                                                     self.toCGS_vel_rotxyz,
                                                     self._todoc_cur)
                    del self._out
            del self.scur, self._todoc_cur
        else:
            self.__startcalc_vel(subindex=specs[0][1])
            self._todoc_cur = {'vcen_cmps': self.vcen_cmps,
                               'rotmatrix': self.rotmatrix,
                               'rotcoord_index': specs[0][1],
                               'units': 'cm * s**-1'}
            self.coords_stored[specs[0]] = (self.vel_rotxyz, 
                                            self.toCGS_vel_rotxyz,
                                            self._todoc_cur)
            del self._todoc_cur
    
    def __calc_veldop(self, specs):
        self._ptov_simu = self.toCGS_coords_rotxyz \
                         * cu.Hubble(cosmopars=self.cosmopars) \
                         / self.toCGS_vel_rotxyz
        if ('vel', 'alldop') in specs or len(specs) > 1:
            for self.scur in specs:
                if self.scur == ('vel', 'alldop'):
                    self._todoc_cur = {'vcen_cmps': self.vcen_cmps,
                                       'cen_cm': self.cen_cm,
                                       'rotmatrix': self.rotmatrix,
                                       'rotcoord_index': [0, 1, 2],
                                       'units': 'cm * s**-1'}
                    self._out = np.copy(self.vel_rotxyz)
                    self._out += self.coords_rotxyz * self._ptov_simu
                    self.coords_stored[self.scur] = (self._out, 
                                                     self.toCGS_vel_rotxyz,
                                                     self._todoc_cur)
                else:
                    self._axind = int(self.scur[1][-1])
                    self._todoc_cur = {'vcen_cmps': self.vcen_cmps,
                                       'cen_cm': self.cen_cm,
                                       'rotmatrix': self.rotmatrix,
                                       'rotcoord_index': self._axind,
                                       'units': 'cm * s**-1'}
                    # storing a view of an array could cause unexpected
                    # side-effects
                    self._out = np.copy(self.vel_rotxyz[:, self._axind])
                    self._out += self.coords_rotxyz[:, self._axind] \
                                 * self._ptov_simu
                    self.coords_stored[self.scur] = (self._out, 
                                                     self.toCGS_vel_rotxyz,
                                                     self._todoc_cur)
                    del self._axind
            del self.scur, self._todoc_cur, self._out
        else:
            self._axind = int(specs[0][1][-1])
            self._todoc_cur = {'vcen_cmps': self.vcen_cmps,
                               'cen_cm': self.cen_cm,
                               'rotmatrix': self.rotmatrix,
                               'rotcoord_index': self._axind,
                               'units': 'cm * s**-1'}
            self._out = np.copy(self.vel_rotxyz[:, self._axind])
            self._out += self.coords_rotxyz[:, self._axind] \
                         * self._ptov_simu
            self.coords_stored[specs[0]] = (self._out, 
                                            self.toCGS_vel_rotxyz,
                                            self._todoc_cur)
            del self._todoc_cur, self._axind, self._out
        del self._ptov_simu,

    def __calc_veltot(self, specs):
        self.scur = specs[0]
        self._in = self.coords_stored[('vel', 'allcart')]
        self._todoc_cur = self._in[2].copy()
        del self._todoc_cur['rotmatrix']
        del self._todoc_cur['rotcoord_index']
        self._out = np.sqrt(np.sum(self._in[0]**2, axis=self.coordaxis))
        self.coords_stored[self.scur] = (self._out, self._in[1], 
                                         self._todoc_cur)
        del self.scur, self._out, self._todoc_cur, self._in

    def __calc_velrad(self, specs):
        self.scur = specs[0]
        self._cendir = np.copy(self.coords_stored[('pos', 'allcart')][0])
        self._cendir /= self.coords_stored[('pos', 'rcen')][0][:, np.newaxis]
        self._out = np.einsum('ij,ij->i', self._cendir, 
                              self.coords_stored[('vel', 'allcart')][0])
        self._units = self.coords_stored[('vel', 'allcart')][1]
        self._todoc_cur = self.coords_stored[('vel', 'allcart')][2].copy()
        del self._todoc_cur['rotmatrix']
        del self._todoc_cur['rotcoord_index']
        self.pkey = ('pos', 'allcart')
        self._todoc_cur['cen_cm'] = self.coords_stored[self.pkey][2]['cen_cm']
        self.coords_stored[self.scur] = (self._out, self._units, 
                                         self._todoc_cur)
        del self.scur, self._out, self._todoc_cur, self._cendir, self._units
        del self.pkey

    def __calc_velphi(self, specs):
        self.scur = specs[0]
        self._in = self.coords_stored[('vel', 'allcart')]
        self._phi = self.coords_stored[('pos', 'phi')][0]
        self._phidir = np.zeros(self._in[0].shape, dtype=self._in[0].dtype)
        self._phidir[:, 0] = - np.sin(self._phi)
        self._phidir[:, 1] = np.cos(self._phi)
        del self._phi
        self._out = np.einsum('ij,ij->i', self._phidir, self._in[0])
        self._units = self.coords_stored[('vel', 'allcart')][1]
        self._todoc_cur = self._in[2].copy()
        del self._todoc_cur['rotcoord_index']
        self._todoc_cur['info'] = ('velocity component in rot. x-y plane'
                                   ' rotation')
        self.coords_stored[self.scur] = (self._out, self._units, 
                                         self._todoc_cur)
        del self.scur, self._out, self._todoc_cur, self._in, self._units
        del self._phidir
    
    def __calc_velazimuth(self, specs):
        self.scur = specs[0]
        self._in = self.coords_stored[('vel', 'allcart')]
        self._phi = self.coords_stored[('pos', 'phi')][0]
        self._azi = self.coords_stored[('pos', 'azimuth')][0]
        self._thetadir = np.zeros(self._in[0].shape, dtype=self._in[0].dtype)
        self._thetadir[:, 0] = np.cos(self._azi) * np.cos(self._phi)
        self._thetadir[:, 1] = np.cos(self._azi) * np.sin(self._phi)
        self._thetadir[:, 2] =  - np.sin(self._azi)
        del self._phi, self._azi
        self._out = np.einsum('ij,ij->i', self._thetadir, self._in[0])
        self._units = self._in[1]
        self._todoc_cur = self._in[2].copy()
        del self._todoc_cur['rotcoord_index']
        self._todoc_cur['info'] = ('velocity component in rot. x-y plane'
                                   ' rotation')
        self.coords_stored[self.scur] = (self._out, self._units, 
                                         self._todoc_cur)
        del self.scur, self._out, self._todoc_cur, self._in, self._units
        del self._thetadir

    def __update_out_todo(self, specs):
        # update output list
        for self.scur in specs:
            for self.i, self.si in enumerate(self.coordspecs_in):
                if self.scur == self.si:
                    self.coords_outlist[self.i] = self.coords_stored[self.si]
        del self.scur, self.i, self.si
        # clean up stored list
        self.still_todo = self.coords_todo[self.gicur + 1:]
        if len(self.still_todo) == 0:
            pass
        else:
            self.curstored = list(self.coords_stored)
            for self.kcur in self.curstored:
                if not (np.any([self.kcur in self.dependencies[_s] 
                                for _g in self.still_todo for _s in _g
                                if _s in self.dependencies])):
                    del self.coords_stored[self.kcur]
            del self.kcur        

def rotmatrix_from_zdir(newzdir):
    '''
    calculate a rotation matrix that will move vector newzdir to
    along the new z axis.
    The new x axis is chosen along the initial x direction, unless
    the angle with newzdir is too small. In that case, the y 
    direction is used.
    '''
    # input might be something like an ang. mom. vector; 
    # normalize first 
    znew = newzdir / np.sqrt(np.sum(newzdir**2))
    xinit = np.array([1., 0., 0.])
    xin = np.inner(znew, xinit)
    if xin <= 0.9:
        xnew = xinit - xin * znew
        xnew = xnew / np.sqrt(np.sum(xnew**2))
        ynew = np.cross(znew, xnew)
    else:
        yinit = np.array([0., 1., 0.])
        yin = np.inner(znew, yinit)
        ynew = yinit - yin * znew
        ynew = ynew / np.sqrt(np.sum(ynew**2))
        xnew = np.cross(ynew, znew)
    rotmat = np.array([xnew, ynew, znew])
    return rotmat