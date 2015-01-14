import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def moment_time_domain(times, moments):
    plt.plot(times, moments[:,0], label='x')
    plt.plot(times, moments[:,1], label='y')
    plt.plot(times, moments[:,2], label='z')
    plt.ylim(-1.1, 1.1)
    plt.ylabel("M/Ms")
    plt.xlabel("Time (sec)")
    plt.legend()
    return plt.gca()


def energy_surface(energies):
    """ Returns the axes of an energy surface plot
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(energies[:,0], energies[:,1], energies[:,2])
    return ax