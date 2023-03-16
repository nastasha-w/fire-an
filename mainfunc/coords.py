
import numpy as np

import fire_an.utils.constants_and_units as c
import fire_an.utils.cosmo_utils as cu

class CoordinateWranger:
    def __init__(self, snapobj, center_cm, rotmatrix=None,
                 parttype=0, periodic=True, vcen_cmps=None):
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
        self.pcalcstarted = False
        self.__check_rotmatrix()
    
    def __check_rotmatrix(self):
        if self.rotmatrix is not None:
            if self.rotmatrix.shape != (3, 3):
                msg = ('Rotation matrix should have shape (3, 3) not input'
                       f'{self.rotmatrix.shape} for matrix\n{self.rotmatrix}')
                raise ValueError(msg)
            if not (np.allclose(self.rotmatrix.T, self.rotmatrix) 
                    and np.isclose(np.linalg.det(self.rotmatrix), 1.)):
                msg = ('input was not a valid rotation matrix.\n'
                       'transpose (should be same):\n'
                       f'{self.rotmatrix.T}, {self.rotmatrix}\n'
                       'determinant (should be 1.): '
                       f'{np.linalg.det(self.rotmatrix)}')
                raise ValueError(msg)
        
    def __startcalc_pos(self, subindex=None):
        h5path = f'PartType{self.pt}/Coordinates'
        if self.rotmatrix is None:
            self._subindex = subindex
        else:
            self._subindex = None
        self.coords_simxyz = self.snapobj.readarray(h5path, 
                                                    subindex=self._subindex)
        self.toCGS_coords_simxyz = self.snapobj.toCGS
        self.__center_pos()
        if self.rotmatrix is not None:
            self.__rotate_pos()
        else:
            self.coords_rotxyz = self.coords_simxyz
            self.toCGS_coords_rotxyz = self.toCGS_coords_simxyz
        if subindex is not None:
            self.coords_rotxyz = np.copy(self.coords_rotxyz[subindex])
        del self.coords_simxyz
        del self.toCGS_coords_simxyz
        self.pcalcstarted = True

    def __startcalc_vel(self, subindex=None):
        h5path = f'PartType{self.pt}/Velocities'
        if self.rotmatrix is None:
            self._subindex = subindex
        else:
            self._subindex = None
        self.vel_simxyz = self.snapobj.readarray(h5path, 
                                                 subindex=self._subindex)
        self.toCGS_vel_simxyz = self.snapobj.toCGS
        self.__center_vel()
        if self.rotmatrix is not None:
            self.__rotate_vel()
        else:
            self.vel_rotxyz = self.vel_simxyz
            self.toCGS_vel_rotxyz = self.toCGS_vel_simxyz
        del self.vel_simxyz
        del self.toCGS_vel_simxyz
        self.vcalcstarted = True

    def __rotate_pos(self):
        self.rotmatrix = np.asarray(self.rotmatrix, 
                                    dtype=self.coords_simxyz.dtype)
        self.coords_rotxyz = np.tensordot(self.rotmatrix, self.coords_simxyz,
                                          axes=([1], [self.coordaxis]))
        self.toCGS_coords_rotxyz = self.toCGS_coords_simxyz
    
    def __rotate_vel(self):
        self.rotmatrix = np.asarray(self.rotmatrix, 
                                    dtype=self.vel_simxyz.dtype)
        self.vel_rotxyz = np.tensordot(self.rotmatrix, self.vel_simxyz,
                                       axes=([1], [self.coordaxis]))
        self.toCGS_vel_rotxyz = self.toCGS_vel_simxyz
    
    def __center_pos(self):
        self.center_simu = np.astype(self.cen_cm / self.toCGS_coords_simxyz,
                                     dtype=self.coords_simxyz.dtype)
        self.coords_simxyz -= self.center_simu
        if self.periodic:
            self.boxsize_simu = self.snapobj.cosmopars.boxsize \
                                * self.snapobj.cosmopars.a \
                                / self.snapobj.cosmopars.h \
                                * c.cm_per_mpc / self.toCGS_coords_simxyz
            self.coords_simxyx += 0.5 * self.boxsize_simu
            self.coords_simxyz %= self.boxsize_simu
            self.coords_simxyx -= 0.5 * self.boxsize_simu
    
    def __center_vel(self):
        self.vcen_simu = np.astype(self.vcen_cmps / self.toCGS_vel_simxyz,
                                   dtype=self.vel_simxyz.dtype)
        self.vel_simxyz -= self.vcen_simu
        if self.periodic:
            self.cosmopars = self.snapobj.cosmopars.getdct()
            self.vboxsize_simu = self.snapobj.cosmopars.boxsize \
                                 * self.snapobj.cosmopars.a \
                                 / self.snapobj.cosmopars.h \
                                 * c.cm_per_mpc / self.toCGS_coords_simxyz \
                                 * cu.Hubble(self.cosmopars['z'], 
                                             cosmopars=self.cosmopars)
            self.vel_simxyx += 0.5 * self.vboxsize_simu
            self.vel_simxyz %= self.vboxsize_simu
            self.vel_simxyx -= 0.5 * self.vboxsize_simu
        
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
                'pos': [0, 1, 2, 'allcart', 'rcen']
                    0, 1, 2: position along the axis with this index
                    'allcart': for all three of these cartesian axes
                    'rcen': distance to the center
                'vel': [0, 1, 2, 'allcart', 'vrad']
                     0, 1, 2: velocity along the axis with this index
                    'allcart': for all three of these cartesian axes
                    'vrad': radial velocity (relative to coordinate
                            center)
                    'vtot': total velocity (rms coordinate velocties)
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
        self.coordspecs_in = [(key, coordspecs[key]) for key in coordspecs]
        ## this priority setting can get messy very fast if I try to
        ## implement too much here.
        # which (groups of) properties to calculate, and in what order
        # a group is calculated in a single function named 
        # __calc_<group key>
        self.calcorder = {'poscart': [('pos', 'allcart'), ('pos', 0),
                                      ('pos', 1), ('pos', 2)],
                          'poscen': [('pos', 'rcen')],
                          'velcart': [('vel', 'allcart'), ('vel', 0),
                                      ('vel', 1), ('vel', 2)],
                          'veltot': [('vel', 'vtot')], 
                          'velcen': [('vel', 'vrad')],
                          }
        # what to get just because it's needed later
        # note: should include dependencies of dependencies
        self.dependencies = {('pos', 'rcen'): [('pos', 'allcart')],
                             ('vel', 'vtot'): [('vel', 'allcart')],
                             ('vel', 'vrad'): [('vel', 'allcart'),
                                               ('pos', 'allcart'),
                                               ('pos', 'rcen')]
                            }
        # set up to-do list of everything that's needed (no duplicates)
        self._coords_todo = set(self.coordspecs_in.copy())
        for _coordspec in self.coordspecs_in:
            if _coordspec in self.dependencies:
                self._coords_todo |= set(self.dependencies[_coordspec])
        self.coords_todo = [[group[key] for key in group 
                             if group[key] in self._coords_todo] 
                            for group in self.calcorder]
        # holds arrays calculated for output 
        self.coords_outlist = [None] * len(self.coordspecs_in)
        # holds all calculated arrays, including those only needed as
        # dependencies. (keys are coordspecs tuples)
        self.coords_stored = {}
        
        for self.gicur, self.gcur in enumerate(self.coords_todo):
            self.gkeymatch = [key for key in self.calcorder 
                              if set(self.gcur).issubset(
                                  set(self.calcorder[key]))]
            self.gkeymatch = self.gkeymatch[0]
            print(f'calculating {self.gcur}')
            if self.gkeymatch == 'poscart':
                self.__calc_poscart(self.gcur)
            elif self.gkeymatch == 'poscen':
                self.__calc_poscen(self.gcur)
            elif self.gkeymatch == 'velcart':
                self.__calc_velcart(self.gcur)
            elif self.gkeymatch == 'veltot':
                self.__calc_veltot(self.gcur)
            elif self.gkeymatch == 'velcen':
                self.__calc_velrad(self.gcur)
            self.__update_out_todo()
            print(f'still todo: {self.still_todo}')
        del self.gcur, self.fcur, self.gkeymatch, self.coords_todo, 
        del self.coords_stored
        return self.coords_outlist
    
    def __calc_poscart(self, specs):
        if ('pos', 'allcart') in specs or len(specs) > 1:
            self.__startcalc_pos(subindex=None)
            for self.scur in specs:
                if self.spec == ('pos', 'allcart'):
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
                    self._out = np.copy(self.coords_rotxyz[self.scur[1]])
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

    def __calc_velcart(self, specs):
        if ('vel', 'allcart') in specs or len(specs) > 1:
            self.__startcalc_vel(subindex=None)
            for self.scur in specs:
                if self.spec == ('vel', 'allcart'):
                    self._todoc_cur = {'vcen_cmps': self.vcen,
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
                    self._out = np.copy(self.vel_rotxyz[self.scur[1]])
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
        self._cendir = self.coords_stored[('pos', 'allcart')][0]
        self._cendir /= self.coords_stored[('pos', 'rcen')][0]
        self._out = np.tensordot(self.cendir, 
                                 self.coords_stored[('vel', 'allcart')][0],
                                 axes=(self.coordaxis, self.coordaxis))
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

    def __update_out_todo(self, specs):
        # update output list
        for self.scur in specs:
            for self.i, self.si in enumerate(self.coordspecs_in):
                if self.insub == self.si:
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
                                for _g in self.still_todo for _s in _g])):
                    del self.coords_stored[self.kcur]
            del self.kcur        
