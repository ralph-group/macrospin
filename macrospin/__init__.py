"""
Includes helper functions for macrospin
"""
from __future__ import division
import numpy as np
from mpl_toolkits.mplot3d import proj3d

__version__ = 0.1


def brillouin(x, J):
    """ Returns the Brillouin function evaluated at x and J """
    f1 = (2*J+1)/(2*J)
    f2 = 1/(2*J)
    return f1/np.tanh(f1*x)-f2/np.tanh(f2*x)


def langevin(x):
    """ Returns the Langevin function evaluated at x and J """
    return np.coth(x)-1/x


def rotation_array(R_function, angles):
    """ Returns an array of rotation matrixes in correct order based on
    a rotation function and an array of angles in degrees

    """
    rotations = np.empty((len(angles), 3, 3))
    for i, angle in enumerate(angles):
        rotations[i] = R_function(angle)
    return rotations


def Rx(angle_degrees):
    """ Constructs a rotation matrix around x based on an angle in degrees """
    angle = angle_degrees/180*np.pi
    return np.matrix([
        [1, 0, 0],
        [0, np.cos(angle), -np.sin(angle)],
        [0, np.sin(angle), np.cos(angle)],
    ])


def Ry(angle_degrees):
    """ Constructs a rotation matrix around y based on an angle in degrees """
    angle = angle_degrees/180*np.pi
    return np.matrix([
        [np.cos(angle), 0, np.sin(angle)],
        [0, 1, 0],
        [-np.sin(angle), 0, np.cos(angle)],
    ])


def Rz(angle_degrees):
    """ Constructs a rotation matrix around z based on an angle in degrees """
    angle = angle_degrees/180*np.pi
    return np.matrix([
        [np.cos(angle), -np.sin(angle), 0],
        [np.sin(angle), np.cos(angle), 0],
        [0, 0, 1],
    ])


def normalize(array):
    """ Normalizes an array row-wise """
    if len(array.shape) == 1:
        return array/np.linalg.norm(array).T
    else:
        normals = np.linalg.norm(array, axis=1)
        normals = np.array(normals.reshape(normals.size, 1))
        return array / normals


def plot_unit_square(figure):
    ax = figure.gca(projection='3d')
    from itertools import product, combinations
    r = [-1, 1]
    for s, e in combinations(np.array(list(product(r, r, r))), 2):
        if np.sum(np.abs(s-e)) == r[1]-r[0]:
            ax.plot3D(*zip(s, e), color="b")


def plot_unit_vectors(figure, vectors, color='r'):
    origin = np.array([0, 0, 0])

    from matplotlib.patches import FancyArrowPatch

    class Arrow3D(FancyArrowPatch):
        def __init__(self, xs, ys, zs, *args, **kwargs):
            FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
            self._verts3d = xs, ys, zs

        def draw(self, renderer):
            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
            self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
            FancyArrowPatch.draw(self, renderer)

    ax = figure.gca(projection='3d')
    for vector in vectors:
        distances = (
            [origin[0], vector[0]],
            [origin[1], vector[1]],
            [origin[2], vector[2]]
        )
        arrow = Arrow3D(*distances, mutation_scale=20, lw=2,
                        arrowstyle="-|>", color=color)
        ax.add_artist(arrow)
    return figure


def plot_surface(figure, normal, point=np.array([0, 0, 0])):
    # TODO: Add functionality for other directions
    d = -point.dot(normal)

    r = np.linspace(-2, 2, num=10)
    yy, zz = np.meshgrid(r, r)
    xx = (-normal[2] * zz - normal[1] * yy - d) * 1. / normal[0]

    ax = figure.gca(projection='3d')
    ax.plot_surface(xx, yy, zz, alpha=0.1)
