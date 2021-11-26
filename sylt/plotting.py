import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import numpy as np
import matplotlib.cm as cm

from sylt.analysis import project, twiss
from sylt.fitting import oscillator, gauss
from sylt.functions import binomial

from sylt.fitting import oscillator

from scipy.special import ellipk

from datetime import datetime


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

    for ax, (k1, k2) in zip(axes, keys):

        x1, x2 = data[k1], data[k2]
        ind = np.isfinite(x1) * np.isfinite(x2)

        tp, (xr, yr), caption = twiss(x1[ind], x2[ind])

        ax.plot(x1[ind], x2[ind], ',')
        ax.plot(xr(t), yr(t))

        ax.annotate(caption, (0, 0), xycoords='axes fraction')
        ax.annotate(k1, (0.5, 0), xycoords='axes fraction',
                    horizontalalignment='center')
        ax.annotate(k2, (0, 0.5), xycoords='axes fraction',
                    rotation=90, verticalalignment='center')

        plot_projections(x1[ind], x2[ind], ax)

    return fig, axes


def plot_tune(DATA):
    """"""
    fig, ax = plt.subplots(num='tune', constrained_layout=True)

    ax.set_ylabel(r"$\mu$")
    ax.set_xlabel(r"$\hat{\phi}$")
    ax.set_xlim(0, np.pi)
    ax.set_ylim(0.5, 1.1)

    phi_hat_q = np.linspace(0, np.pi, 300)
    ax.plot(phi_hat_q, 1-phi_hat_q**2/16, 'k:')
    ax.plot(phi_hat_q, np.pi/(2*sp.ellipk(np.sin(phi_hat_q/2)**2)), 'k-')

    for i, (key, data) in enumerate(DATA.items()):
        ax.plot(data['phi_hat'], data['mu'], '.', alpha=0.33)
        ax.plot([], [], f'C{i}.', label=key)

    ax.legend(loc='upper right')


def plot_bunch_profiles(tau, t, lam, fit, waterfall=True, contour=True, decay=True):
    """plot bunch length oscillations from longitudinal profile evolution"""

    signal = oscillator(t, **fit['oscillator'])
    env = oscillator(t, A=fit['oscillator']['A'],
                     lam=fit['oscillator']['lam'])
    N = np.sum(np.gradient(tau)*lam, -1)

    figs = {}

    if waterfall:

        series = (
            (r'$\lambda(\tau)$', lam.T, '-'),
            # ('gaussian', gauss(tau[:, np.newaxis], **fit['gaussian']), '--'),
            ('binomial', binomial(tau[:, np.newaxis], **fit['binomial']), ':'),
        )

        k = t.max()/(3*lam.max())
        ind = np.arange(0, t.size, int(t.size/20))

        figs['waterfall'], ax1 = plt.subplots(constrained_layout=True)
        for key, y, linespec in series:
            ax1.plot(tau*1e9, 1e3*(t+k*y)[:, ind], linespec)
            ax1.plot([], [], f"k{linespec}", label=key)
        ax1.set_xlabel(r"$\tau$ (ns)")
        ax1.set_ylabel(r"$t$ (ms)")
        ax1.set_xlim(1e9*np.mean(4*fit['binomial']['sig'])*np.array([-0.75, 0.75]))
        ax1.legend(loc='upper right')

    if contour:
        figs['contour'], ax2 = plt.subplots(constrained_layout=True)
        pcm = ax2.pcolormesh(1e3*t, 1e9*tau, lam.T,
                             cmap='Blues', shading='auto')
        ax2.set_ylabel(r"$\tau$ (ns)")
        ax2.set_xlabel(r"$t$ (ms)")
        figs['contour'].colorbar(
            pcm, ax=ax2, label=r'$\lambda(\tau)$ (ns$^{-1}$)')
        ax2.plot(1e3*t, 1e9*fit['gaussian']['mu'],'k-',
                 label=r"$\mu_{\sigma_\tau}(t)$")
        ax2.annotate(f"$<N>$:{np.mean(N):0.2e}",
                     xy=(0, 0), xycoords='axes fraction')

        ax2.set_ylim(
            1e9*(np.mean(fit['gaussian']['mu']+3*fit['gaussian']['var']**0.5)),
            1e9*(np.mean(fit['gaussian']['mu']-3*fit['gaussian']['var']**0.5)),
        )


    if decay:

        amp, omega, lam, phi, mu = fit['oscillator'].values()

        series = (
            # (r"$\sigma_\tau(t)$", fit['gaussian']['var']**0.5, '.'),
            (r"$L_\tau(t)/4$", fit['binomial']['sig'], '.'),
            (f"${1e9*mu:0.3g}+{1e9*amp:0.3g}\exp(-t/{1e3/lam:0.3g})\cos({1e-3*omega:0.3g}t{phi:+0.3g})$", signal, '-'),
            (f"${1e9*mu:0.3g}+{1e9*amp:0.3g}\exp(-t/{1e3/lam:0.3g})$",
             fit['oscillator']['mu']+env, ':'),
            (f"${1e9*mu:0.3g}-{1e9*amp:0.3g}\exp(-t/{1e3/lam:0.3g})$",
             fit['oscillator']['mu']-env, ':'),
        )

        figs['decay'], ax3 = plt.subplots(constrained_layout=True)
        for key, y, linespec in series:
            ax3.plot(1e3*t, 1e9*y, linespec, label=key)

        ax3.set_xlim(1e3*t.min(), 1e3*t.max())
        ax3.set_xlabel(r"$t$ (ms)")
        ax3.set_ylabel(r"$\sigma_\tau$ (ns)")
        ax3.legend(loc='upper right')

    return figs
