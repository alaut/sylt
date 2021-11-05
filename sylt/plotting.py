import matplotlib.pyplot as plt
import numpy as np

from sylt.analysis import twiss, project

def plot_projections(x, y, ax, bins=100, hist_height=5):
    """plot projections of x, y distribution on a given axes"""

    axh, axv = ax.twinx(), ax.twiny()

    xc, xh, x_fit = project(x, bins)
    yc, yh, y_fit = project(y, bins)

    axh.plot(xc, xh)
    axv.plot(yh, yc)

    axh.plot(xc, x_fit, ':')
    axv.plot(y_fit, yc, ':')

    axh.set_ylim(0, hist_height*xh.max())
    axv.set_xlim(0, hist_height*yh.max())

    axh.set_axis_off()
    axv.set_axis_off()

    return axh, axv


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
