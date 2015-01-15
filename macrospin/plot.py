import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.patches import FancyArrowPatch
import numpy as np


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)


def unit_square(figure):
    ax = figure.gca(projection='3d')
    from itertools import product, combinations
    r = [-1, 1]
    for s, e in combinations(np.array(list(product(r, r, r))), 2):
        if np.sum(np.abs(s-e)) == r[1]-r[0]:
            ax.plot3D(*zip(s, e), color="b")


def unit_view(figure):
    ax = figure.gca(projection='3d')
    ax.view_init(elev=25, azim=45)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_zlim(-1.5, 1.5)
    ax.set_axis_off()
    axes_3D(figure)


def axes_3D(figure, mutation_scale=20, lw=1.5, arrowstyle="<|-|>", color='k', **kwargs):
    ax = figure.gca(projection='3d')
    distances = [
        ([-1.5,1.5], [0.,0.], [0.,0.]),
        ([0.,0.], [-1.5,1.5], [0.,0.]),
        ([0.,0.], [0.,0.], [-1.5,1.5]),
    ]
    for d in distances:
        arrow = Arrow3D(*d, mutation_scale=mutation_scale, lw=lw, 
                    arrowstyle=arrowstyle, color=color, **kwargs)
        ax.add_artist(arrow)


def unit_vector(figure, vector, **kwargs):
    origin = np.array([0, 0, 0])

    ax = figure.gca(projection='3d')
    distances = (
        [origin[0], vector[0]],
        [origin[1], vector[1]],
        [origin[2], vector[2]]
    )
    arrow = Arrow3D(*distances, **kwargs)
    ax.add_artist(arrow)


def unit_vectors(figure, vectors, **kwargs):
    for vector in vectors:
        unit_vector(figure, vector, **kwargs)


def surface(figure, normal, point=np.array([0, 0, 0])):
    # TODO: Add functionality for other directions
    d = -point.dot(normal)

    r = np.linspace(-1, 1, num=10)
    yy, zz = np.meshgrid(r, r)
    xx = (-normal[2] * zz - normal[1] * yy - d) * 1. /normal[0]

    ax = figure.gca(projection='3d')
    ax.plot_surface(xx, yy, zz, alpha=0.1)


def moment_time_domain(times, moments):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(times, moments[:,0], label='x')
    ax.plot(times, moments[:,1], label='y')
    ax.plot(times, moments[:,2], label='z')
    ax.set_ylim(-1.1, 1.1)
    ax.set_ylabel("M/Ms")
    ax.set_xlabel("Time (sec)")
    ax.legend()
    return ax


def moment_3d(moments):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(moments[:,0], moments[:,1], moments[:,2])
    unit_view(fig)
    #unit_vector(fig, moments[0], color='r')
    #unit_vector(fig, moments[-1], color='r')
    return ax


def energy_surface(energies):
    """ Returns the axes of an energy surface plot
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(energies[:,0], energies[:,1], energies[:,2])
    return ax