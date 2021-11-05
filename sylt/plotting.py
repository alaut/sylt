import matplotlib.pyplot as plt
import numpy as np

from sylt.analysis import twiss, project

def plot_projections(x1, x2, ax, bins=100, hist_height=5):
    """plot projections of x, y distribution on a given axes"""

    ax1, ax2 = ax.twinx(), ax.twiny()

    xq1, h1, f1 = project(x1, bins)
    xq2, h2, f2 = project(x2, bins)

    ax1.plot(xq1, h1)
    ax2.plot(h2, xq2)

    ax1.plot(xq1, f1, ':')
    ax2.plot(f2, xq2, ':')

    ax1.set_ylim(0, hist_height*h1.max())
    ax2.set_xlim(0, hist_height*h2.max())

    ax1.set_axis_off()
    ax2.set_axis_off()

    return ax1, ax2


def plot_phase_space(data, keys, shape=None, title=None):
    """plots phase space distributions given data and keys"""

    if shape is None:
        shape(1, len(keys))

    fig, axes = plt.subplots(*shape, constrained_layout=True, num=title)

    t = np.linspace(0, 2*np.pi, 50)

    for ax, (k1, k2) in zip(axes.flatten(), keys):

        x1, x2 = data[k1], data[k2]
        tp, (xr, yr), caption = twiss(x1, x2)

        ax.plot(x1.ravel(), x2.ravel(), ',')
        ax.plot(xr(t), yr(t))

        ax.set_xlabel(k1)
        ax.set_ylabel(k2)
        
        ax.annotate(caption, (0, 0), xycoords='axes fraction')

        plot_projections(x1, x2, ax)

    return axes
