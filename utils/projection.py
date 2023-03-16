# -*- coding: utf-8 -*-
"""

"""
import numpy as np
import ctypes as ct

import fire_an.utils.opts_locs as ol

def project(NumPart, Ls, Axis1, Axis2, Axis3, box3, periodic, npix_x, npix_y,
            kernel, dct, tree, ompproj=True, projmin=None, projmax=None):
    '''
    Parameters:
    -----------
    NumPart: int
        number of SPH particles to project
    Ls: length 3 indexable of floats    
        dimensions (diameter) of the box to project (same units as 
        coordinates)           
    Axis<i>: [0, 1, 2]
        for Ls and coordinates, these variables control which axis is
        the the projection axis (Axis3), and the orientation of the 
        other two. For a z-projection with a resulting array 
        (X index, Y index), use Axis1=0, Axis2=1, Axis3=2
    box3: length 3 indexable of floats
        the dimensions of the parent box (same units as Ls)
    periodic: bool
        is the projected region (perpendicular to the line of sight) a
        full slice of a periodic simulation (gas distributions at the 
        edges should be wrapped around the box), or a smaller part of 
        the simulation (gas contributions outside the projected region 
        edges should be ignored)
    npix_x,y: int
        how many pixels to use in the Axis1 and Axis2 directions, 
        respectively. Note that the minimum smoothing length is set to
        the pixel diagonal, so very non-square pixels won't actually 
        add much resolution in the higher-resolution direction.
    kernel: ['C2', 'gadget'] 
        what shape to assume for the gas distribution respresented by a
        single SPH particle
    dct: dict containing arrays 'coords', 'lsmooth', 'qW', 'qQ' 
        (prevents copying of large arrays)
        'coords': coordinates, (Numpart, 3) float array. Coordinates 
            should be transformed so that the projected region is a 
                [-Ls[0] / 2., Ls[0] / 2.,
                 -Ls[1] / 2., Ls[1] / 2.,
                 -Ls[2] / 2., Ls[2] / 2. ]
            box (if not periodic) or
                [0., Ls[0], 0., Ls[1], 0., Ls[2]] 
            if it is.
            (The reason for this assumption in the periodic case is
            that it makes it easy to determine when something needs to
            be wrapped around the edge, and for the non-periodic case,
            it allows the code to ignore periodic conditions even 
            though the simulations are periodic and the selected region
            could therefore in principle require wrapping.)
        'lsmooth': float array
            gas smoothing lengths (same units as coords)
        'qW': float array
            the array containing the particle property to directly
            project, and to weight qQ by
        'qQ': float array
            the array to get a qW-weighted average for in each pixel
    projmin, projmax: float
        maximum coordinate values in projection direction
        (override default values in Ls; I put this in for a specific
        application)
              
    Returns:
    --------
    (ResultW, ResultQ) : tuple of (npix_x, npix_y) arrays (float32)
    ResultW: qW projected onto the grid. The array contains the sum of
        qW contributions to each pixel, not a qW surface density.
        The sums of ResultW and qW shoudl be the same to floating-point
        errors when projecting a whole simulation, but apparent mass 
        loss in the projection may occur when projecting smaller 
        regions, where some particles in the qW array are (partially) 
        outside the projected region
    - ResultQ: qW-weighted average of qQ in each pixel
    '''

    # positions [Mpc / cm/s], kernel sizes [Mpc] and input quantities
    # a quirk of HsmlAndProject is that it only works for >= 100 
    # particles. Pad with zeros if less.
    if NumPart >=100:
        pos = dct['coords'].astype(np.float32)
        Hsml = dct['lsmooth'].astype(np.float32)
        qW = dct['qW'].astype(np.float32)
        qQ = dct['qQ'].astype(np.float32)

    else:
        qQ = np.zeros((100,), dtype=np.float32)
        qQ[:NumPart] = dct['qQ'].astype(np.float32)
        qW = np.zeros((100,), dtype=np.float32)
        qW[:NumPart] = dct['qW'].astype(np.float32)
        Hsml = np.zeros((100,), dtype=np.float32)
        Hsml[:NumPart] = dct['lsmooth'].astype(np.float32)
        # should put the particles outside any EAGLE or FIRE projection
        # box
        pos = np.ones((100,3), dtype=np.float32) * 1e8 
        pos[:NumPart,:] = dct['coords'].astype(np.float32)
        NumPart = 100

    # ==============================================
    # Putting everything in right format for C routine
    # ==============================================

    print('\n--- Calling findHsmlAndProject ---\n')

    # define edges of the map wrt centre coordinates [Mpc]
    # in the periodic case, the C function expects all coordinates to
    # be in the [0, BoxSize] range 
    # (though I don't think it actually reads Xmin etc. in for this)
    # these need to be defined wrt the 'rotated' axes, 
    # e.g. Zmin, Zmax are always the min/max along the projection 
    # direction
    if not periodic: # 0-centered
        Xmin = -0.5 * Ls[Axis1]
        Xmax =  0.5 * Ls[Axis1]
        Ymin = -0.5 * Ls[Axis2]
        Ymax =  0.5 * Ls[Axis2]
        if projmin is None:
            Zmin = -0.5 * Ls[Axis3]
        else:
            Zmin = projmin
        if projmax is None:
            Zmax = 0.5 * Ls[Axis3]
        else:
            Zmax = projmax
    # half box centered (BoxSize used for x-y periodic boundary 
    # conditions)
    else: 
        Xmin, Ymin = (0.,) * 2
        Xmax, Ymax = (box3[Axis1], box3[Axis2])
        if projmin is None:
            Zmin = 0.5 * (box3[Axis3] - Ls[Axis3])
        else:
            Zmin = projmin
        if projmax is None:
            Zmax = 0.5 * (box3[Axis3] + Ls[Axis3])
        else:
            Zmax = projmax

    BoxSize = box3[Axis1]

    # maximum kernel size [Mpc] (modified from Marijke's version)
    # Axis3 might be velocity; whole different units, so just ignore
    Hmax = 0.5 * min(Ls[Axis1], Ls[Axis2]) 

    # arrays to be filled with resulting maps
    ResultW = np.zeros((npix_x, npix_y)).astype(np.float32)
    ResultQ = np.zeros((npix_x, npix_y)).astype(np.float32)

    # input arrays for C routine (change in c_pos <-> change in pos)
    c_pos = pos[:,:]
    c_Hsml = Hsml[:]
    c_QuantityW = qW[:]
    c_QuantityQ = qQ[:]
    c_ResultW = ResultW[:,:]
    c_ResultQ = ResultQ[:,:]

    # check if HsmlAndProject changes
    print('Total quantity W in: %.5e' % (np.sum(c_QuantityW)))
    print('Total quantity Q in: %.5e' % (np.sum(c_QuantityQ)))

    # path to shared library
    if ompproj:
        sompproj = '_omp'
    else:
        sompproj = ''
    if tree:
        # in v3, projection can use more particles than c_int max,
        # but the tree building cannot
        if not ct.c_int(NumPart).value == NumPart:
            msg = (' ***         Warning         ***\n\n'
                   'Number of particles %i overflows C int type.\n'
                   'This will likely cause the tree building routine'
                   ' in HsmlAndProjcet_v3 to fail.\nSee notes on v3 version.'
                   '\n\n*****************************\n')
            print(msg)
        if periodic:
            lib_path = ol.hsml_dir \
                       + f'HsmlAndProject_v3_{kernel}_perbc{sompproj}.so'
        else:
            lib_path = ol.hsml_dir \
                       + f'HsmlAndProject_v3_{kernel}{sompproj}.so'
    else:
        if periodic:
            lib_path = ol.hsml_dir \
                + f'HsmlAndProject_v3_notree_{kernel}_perbc{sompproj}.so'
        else:
            lib_path = ol.hsml_dir \
                + f'HsmlAndProject_v3_notree_{kernel}{sompproj}.so'

    print('Using projection file: %s \n' % lib_path)
    # load the library
    my_library = ct.CDLL(lib_path)

    # set the parameter types (numbers with ctypes, arrays with ndpointers)
    argt = [ct.c_long,
            np.ctypeslib.ndpointer(dtype=ct.c_float,shape=(NumPart,3)),
            np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,)),
            np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,)),
            np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,)),
            ct.c_float, ct.c_float, ct.c_float,
            ct.c_float, ct.c_float, ct.c_float,
            ct.c_int, ct.c_int, ct.c_int,
            ct.c_int, ct.c_int, ct.c_int,
            ct.c_float,
            ct.c_double,
            np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(npix_x,npix_y)),
            np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(npix_x,npix_y)),
            ]
    my_library.findHsmlAndProject.argtypes = argt

    # set the return type
    my_library.findHsmlAndProject.restype = None

    print('----------')

    # call the findHsmlAndProject C routine
    my_library.findHsmlAndProject(ct.c_long(NumPart),
                                  c_pos,
                                  c_Hsml,
                                  c_QuantityW,
                                  c_QuantityQ,
                                  ct.c_float(Xmin),
                                  ct.c_float(Xmax),
                                  ct.c_float(Ymin),
                                  ct.c_float(Ymax),
                                  ct.c_float(Zmin),
                                  ct.c_float(Zmax),
                                  ct.c_int(npix_x),
                                  ct.c_int(npix_y),
                                  ct.c_int(ol.desngb),
                                  ct.c_int(Axis1),
                                  ct.c_int(Axis2),
                                  ct.c_int(Axis3),
                                  ct.c_float(Hmax),
                                  ct.c_double(BoxSize),
                                  c_ResultW,
                                  c_ResultQ)

    print('----------')

    # check if mapped quantities conserved 
    # (but first one also counts some particles outside the map 
    # actually)
    print('Total quantity W in:  %.5e' % (np.sum(c_QuantityW)))
    print('Total quantity W out: %.5e' % (np.sum(ResultW)))
    print('Total quantity Q in:  %.5e' % (np.sum(c_QuantityQ)))
    print('Total quantity Q out: %.5e' % (np.sum(ResultQ)))

    return ResultW, ResultQ