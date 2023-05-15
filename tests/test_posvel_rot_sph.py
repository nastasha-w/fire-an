
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as m3d 
import numpy as np

import fire_an.mainfunc.coords as crd
import fire_an.mainfunc.get_qty as gq
import fire_an.readfire.readin_fire_data as rfd

zdir_tests = [np.array([1., 1., 1.]),
              np.array([-1., -1., -1.]),
              np.array([1., -1., 0.])]

torotate_tests = np.diag(np.ones(3)) * 2.

def plot_posquiver(*args, colors=None):
    fig = plt.figure()
    ax = m3d.Axes3D(fig)
    maxv = 0.
    for i, arg in enumerate(args):
        if len(arg.shape) == 2:
            xv = arg[0, :]
            yv = arg[1, :]
            zv = arg[2, :]
            x, y, z = np.zeros((3, len(xv)))
            maxv = max(maxv, np.max(np.abs(arg)))
        else:
            xv = arg[0]
            yv = arg[1]
            zv = arg[2]
            x, y, z = np.zeros((3,))
            maxv = max(maxv, np.max(np.abs(arg)))
        if colors is None:
            c = f'C{i % 10}'
        else:
            c = colors[i]
        ax.quiver(x, y, z, xv, yv, zv, color=c)
    ax.set_xlim(-maxv, maxv)
    ax.set_ylim(-maxv, maxv)
    ax.set_zlim(-maxv, maxv)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

def plot_posvelquiver(poss, vels, colors=None):
    fig = plt.figure()
    ax = m3d.Axes3D(fig)
    maxv = 0.
    for i, (pos, vel) in enumerate(zip(poss, vels)):
        if len(pos.shape) == 2:
            xv = pos[:, 0]
            yv = pos[:, 1]
            zv = pos[:, 2]
            x, y, z = np.zeros((3, len(xv)))
            maxv = max(maxv, np.max(np.abs(pos + vel)),
                       np.max(np.abs(pos)))
            x2 = pos[:, 0]
            y2 = pos[:, 1]
            z2 = pos[:, 2]
            xv2 = vel[:, 0]
            yv2 = vel[:, 1]
            zv2 = vel[:, 2]
            maxv = max(maxv, np.max(np.abs(pos + vel)),
                       np.max(np.abs(pos)))
        else:
            xv = pos[0]
            yv = pos[1]
            zv = pos[2]
            x, y, z = np.zeros((3,))
            maxv = max(maxv, np.max(np.abs(pos)),
                       np.max(np.abs(pos + vel)))
            x2 = pos[0]
            y2 = pos[1]
            z2 = pos[2]
            xv2 = vel[0]
            yv2 = vel[1]
            zv2 = vel[2]
        if colors is None:
            c = f'C{i % 10}'
        else:
            c = colors[i]
        ax.quiver(x, y, z, xv, yv, zv, color=c, linestyle='dashed')
        ax.quiver(x2, y2, z2, xv2, yv2, zv2, color=c)
    ax.set_xlim(-maxv, maxv)
    ax.set_ylim(-maxv, maxv)
    ax.set_zlim(-maxv, maxv)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

def run_rottests():
    '''
    plot and check by eye, just the rotation matrix
    '''
    zdir_tests = [np.array([1., 1., 1.]),
                  np.array([-1., -1., -1.]),
                  np.array([1., -1., 0.]),
                  np.array([1., 0., 0.])]

    xin = np.array([2., 0., 0.])
    yin = np.array([0., 2., 0.])
    zin = np.array([0., 0., 2.])
    
    for znew in zdir_tests:
        rm = crd.rotmatrix_from_zdir(znew)
        xout = np.matmul(rm, xin)
        yout = np.matmul(rm, yin)
        zout = np.matmul(rm, zin)
        zdirout = np.matmul(rm, znew)
        plot_posquiver(xout, yout, zout, zdirout, znew,
                       colors=['red', 'green', 'blue', 'gray', 'black'])
        
def check_rot_cartesian():
    '''
    check p, v rotation through get_qty, CoordinateWrangler
    '''
    zdir_tests = [np.array([1., 1., 1.]),
                  np.array([1., -1., 0.])]
    pos_in =  2. * np.diag(np.ones(3))
    vel_in = np.array([[0., 0.5, 0.],
                       [0., 2., 0.],
                       [1., 0., 0.]])
    pcen = np.array([1., 2., -1.])
    vcen = np.array([-0.2, 3., 1.2])
    pos_in += pcen
    vel_in += vcen
    znewvel = np.array([1., 0., 0.])
    for znew in zdir_tests:
        _pos_in = np.append(pos_in, (znew + pcen)[np.newaxis, :], axis=0)
        _vel_in = np.append(vel_in, (znewvel + vcen)[np.newaxis, :], axis=0)
        print(_pos_in)
        print(_vel_in)
        mockcoords = {'PartType1/Coordinates': _pos_in,
                      'PartType1/Velocities': _vel_in,
                     }
        testsnap = rfd.MockFireSpec(mockcoords)
        rm = crd.rotmatrix_from_zdir(znew)
        mtargs = {'multiple': [{'pos': 'allcart'}, {'vel': 'allcart'}],
                  'vcen_cmps': vcen, 'center_cm': pcen, 'rotmatrix': rm}
        pv, pv_units, pv_todoc = gq.get_qty(testsnap, 1, 'coords', mtargs,
                                            filterdct=None)
        
        poss_plot = [pv[0][0], pv[0][1], pv[0][2], pv[0][3],
                     znew]
        vels_plot = [pv[1][0], pv[1][1], pv[1][2], pv[1][3],
                     znewvel]
        colors = ['red', 'green', 'blue', 'gray', 'black']
        plot_posvelquiver(poss_plot, vels_plot, colors=colors)

def test_sph_coords():
    '''
    simulation axis coordinates are only used in the centering and
    for the rotation calc; a simple test should suffice to see if
    the rotation is included
    r and vr are already tested, so no need to focus on those
    '''
    phiring = np.array([[1., 0., 0.],
                        [1., 1., 0.],
                        [0., 1., 0.],
                        [-1., 1., 0.],
                        [-1., 0., 0.],
                        [-1., -1., 0.],
                        [0., -1., 0.],
                        [1., -1., 0.],
                        ])
    vel_phiring = np.array([[0., 1., 1.],
                            [-1., 1., -1.],
                            [-1., 0., 1.],
                            [-1., -1, -1.],
                            [0., -1, 1.],
                            [1., -1., -1.],
                            [1., 0., 1.],
                            [1., 1., -1.],
                            ])
    phi_exp_phiring = 0.25 * np.pi * np.arange(8)
    theta_exp_phiring = 0.5 * np.pi * np.ones(8)
    r_exp_phiring = np.array([1., 2.**0.5] * 4)
    phivel_exp_phiring = np.array([1., 2**0.5] * 4)
    thetavel_exp_phiring = np.array([-1., 1.] * 4)
    rvel_exp_phiring = np.zeros(8)

    thetaslices = np.array([[0., 0., 1.],
                            [0., 0., -1.],
                            [1., 0., 1.],
                            [1., 0., -1.],
                            [0., 1., 1.],
                            [0., 1., -1.],
                            [-1, 0., 1.],
                            [-1., 0, -1.],
                            [0., -1., 1.],
                            [0., -1., -1.],
                            ]) 
    vel_thetaring = np.array([[0., 0., 1.],
                              [0., 0., 1.],
                              [-1., 0., 1.],
                              [1., 0., 1.],
                              [0., -1., 1.],
                              [0., 1., 1.],
                              [1., 0., 1.,],
                              [-1., 0., 1.],
                              [0., 1., 1.],
                              [0., 1., -1.]])
    phi_exp_thetaslices = np.array([0., 0., 
                                    0., 0., 
                                    0.5 * np.pi, 0.5 * np.pi,
                                    np.pi, np.pi,
                                    -0.5 * np.pi, -0.5 * np.pi])
    theta_exp_thetaslices = np.array([0., np.pi] + 
                                     [0.25 * np.pi, 0.75 * np.pi] * 4)
    r_exp_thetaslices = np.array([1.] * 2 + [2**0.5] * 8)
    phivel_exp_thetaslices = np.zeros(10)
    thetavel_exp_thetaslices = np.array([0., 0.] + [-2**0.5] * 7 +[2**0.5])
    rvel_exp_thetaslices = np.array([1., -1.] + [0.] * 8)

    # just swap axes x -> y -> z -> x, see if the phiring is still ok   
    rm_axswap = np.array([[0., 1., 0.], 
                          [0., 0., 1.], 
                          [1., 0., 0.]])
    phiring_prerot = np.array([[0., 1., 0.],
                               [0., 1., 1.],
                               [0., 0., 1.],
                               [0., -1., 1.],
                               [0., -1, 0.],
                               [0., -1., -1.],
                               [0., 0., -1.],
                               [0., 1., -1.],
                               ])
    vel_phiring_prerot = np.array([vel_phiring[:, 2], 
                                   vel_phiring[:, 0], 
                                   vel_phiring[:, 1]]).T
    # non-rotated coords, ring around z-axis:
    pcen = np.array([-6., 1.2, -5.])
    vcen = np.array([3., -2.5, 1.])
    mockcoords = {'PartType0/Coordinates': phiring + pcen,
                  'PartType0/Velocities': vel_phiring + vcen} 
    testsnap = rfd.MockFireSpec(mockcoords)
    rm = None
    mtargs = {'multiple': [{'vel': 'azimuth'}, {'vel': 'phi'}, 
                           {'vel': 'vrad'}, 
                           {'pos': 'azimuth'}, {'pos': 'phi'}, 
                           {'pos': 'rcen'}],
              'vcen_cmps': vcen, 'center_cm': pcen, 'rotmatrix': rm}
    pv, pv_units, pv_todoc = gq.get_qty(testsnap, 0, 'coords', mtargs,
                                        filterdct=None)
    vtheta, vphi, vr, ptheta, pphi, pr = pv
    phiring_match = True
    _m = np.allclose(pphi % (2 * np.pi), phi_exp_phiring % (2. * np.pi))
    print('phiring, phi: ', _m)
    phiring_match &= _m
    _m = np.allclose(ptheta, theta_exp_phiring)
    print('phiring, theta: ', _m)
    phiring_match &= _m
    _m = np.allclose(pr, r_exp_phiring)
    print('phiring, r: ', _m)
    phiring_match &= _m
    _m = np.allclose(vphi, phivel_exp_phiring)
    print('phiring, vphi: ', _m)
    phiring_match &= _m
    _m = np.allclose(vtheta, thetavel_exp_phiring)
    print('phiring, vtheta: ', _m)
    print(vtheta)
    print(thetavel_exp_phiring)
    phiring_match &= _m
    _m = np.allclose(vr, rvel_exp_phiring)
    print('phiring, vr: ', _m)
    phiring_match &= _m
    
    # non-rotated coordinates, slices along x-z, y-z planes
    mockcoords = {'PartType0/Coordinates': thetaslices + pcen,
                  'PartType0/Velocities': vel_thetaring + vcen} 
    testsnap = rfd.MockFireSpec(mockcoords)
    rm = None
    mtargs = {'multiple': [{'vel': 'azimuth'}, {'vel': 'phi'}, 
                           {'vel': 'vrad'}, 
                           {'pos': 'azimuth'}, {'pos': 'phi'}, 
                           {'pos': 'rcen'}],
              'vcen_cmps': vcen, 'center_cm': pcen, 'rotmatrix': rm}
    pv, pv_units, pv_todoc = gq.get_qty(testsnap, 0, 'coords', mtargs,
                                        filterdct=None)
    vtheta, vphi, vr, ptheta, pphi, pr = pv
    thetaring_match = True
    _m = np.allclose(pphi % (2 * np.pi), phi_exp_thetaslices % (2 * np.pi))
    print('thetaring, phi: ', _m)
    thetaring_match &= _m
    _m = np.allclose(ptheta, theta_exp_thetaslices)
    print('thetaring, theta: ', _m)
    thetaring_match &= _m
    _m = np.allclose(pr, r_exp_thetaslices)
    print('thetaring, r: ', _m)
    thetaring_match &= _m
    _m = np.allclose(vphi, phivel_exp_thetaslices)
    print('thetaring, vphi: ', _m)
    thetaring_match &= _m
    _m = np.allclose(vtheta, thetavel_exp_thetaslices)
    print('thetaring, vtheta: ', _m)
    thetaring_match &= _m
    _m = np.allclose(vr, rvel_exp_thetaslices)
    print('thetaring, vr: ', _m)
    thetaring_match &= _m

    # rotated coordinates:
    pcen = np.array([-6., 1.2, -5.])
    vcen = np.array([3., -2.5, 1.])
    mockcoords = {'PartType0/Coordinates': phiring_prerot + pcen,
                  'PartType0/Velocities': vel_phiring_prerot + vcen} 
    testsnap = rfd.MockFireSpec(mockcoords)
    mtargs = {'multiple': [{'vel': 'azimuth'}, {'vel': 'phi'}, 
                           {'vel': 'vrad'}, 
                           {'pos': 'azimuth'}, {'pos': 'phi'}, 
                           {'pos': 'rcen'}],
              'vcen_cmps': vcen, 'center_cm': pcen, 'rotmatrix': rm_axswap}
    pv, pv_units, pv_todoc = gq.get_qty(testsnap, 0, 'coords', mtargs,
                                        filterdct=None)
    vtheta, vphi, vr, ptheta, pphi, pr = pv
    phiring_rot_match = True
    _m = np.allclose(pphi % (2 * np.pi), phi_exp_phiring % (2 * np.pi))
    print('phiring_rot, phi: ', _m)
    phiring_rot_match &= _m
    _m = np.allclose(ptheta, theta_exp_phiring)
    print('phiring_rot, theta: ', _m)
    phiring_rot_match &= _m
    _m = np.allclose(pr, r_exp_phiring)
    print('phiring_rot, r: ', _m)
    phiring_rot_match &= _m
    _m = np.allclose(vphi, phivel_exp_phiring)
    print('phiring_rot, vphi: ', _m)
    phiring_rot_match &= _m
    _m = np.allclose(vtheta, thetavel_exp_phiring)
    print('phiring_rot, vtheta: ', _m)
    phiring_rot_match &= _m
    _m = np.allclose(vr, rvel_exp_phiring)
    print('phiring_rot, vr: ', _m)
    phiring_rot_match &= _m

    print('phi ring match: ', phiring_match)
    print('theta slices match: ', thetaring_match)
    print('phi ring w. rotation match: ', phiring_rot_match)
    return phiring_match & thetaring_match & phiring_rot_match 


