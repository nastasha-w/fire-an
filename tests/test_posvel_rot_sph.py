
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
    plot and check by eye
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




