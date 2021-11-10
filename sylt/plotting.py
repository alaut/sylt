import matplotlib.pyplot as plt
import numpy as np

from sylt.analysis import twiss, project

from sylt.fitting import gauss, fit_gaussian, oscillator, fit_oscillator

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

    for ax, (k1, k2) in zip(axes.flatten(), keys):

        x1, x2 = data[k1], data[k2]
        ind = np.isfinite(x1) * np.isfinite(x2)

        tp, (xr, yr), caption = twiss(x1[ind], x2[ind])

        ax.plot(x1[ind], x2[ind], ',')
        ax.plot(xr(t), yr(t))

        ax.annotate(caption, (0, 0), xycoords='axes fraction')
        ax.annotate(k1, (0.5, 0), xycoords='axes fraction', horizontalalignment='center')
        ax.annotate(k2, (0, 0.5), xycoords='axes fraction', rotation=90, verticalalignment='center')
        
        plot_projections(x1[ind], x2[ind], ax)

    return axes

def plot_tune(DATA):
    """"""
    fig, ax = plt.subplots(num='tune', constrained_layout=True)
    ax.set_ylabel(r"$\mu$")
    ax.set_xlabel(r"$\hat{\phi}$")
    ax.set_xlim(0, np.pi)
    ax.set_ylim(0.5, 1.1)

    phi_hat_q = np.linspace(0, np.pi, 300)
    ax.plot(phi_hat_q, 1-phi_hat_q**2/16, 'k:')
    ax.plot(phi_hat_q, np.pi/(2*ellipk(np.sin(phi_hat_q/2)**2)), 'k-')

    for i, (key, data) in enumerate(DATA.items()):
        ax.plot(data['phi_hat'], data['mu'], '.', alpha=0.33)
def plot_bunch_profiles(tau, t, lam, show=True):
    """analyze bunch length oscillations from longitudinal profile evolution"""

    fit = {}
    fit['gaussian'] = fit_gaussian(tau, lam)
    fit['oscillator'], eqn = fit_oscillator(t, fit['gaussian']['var']**0.5)

    if show:

        fig, (ax1, ax2) = plt.subplots(2, 1, num=None,
                                       constrained_layout=True, sharex=True)

        pcm = ax1.pcolormesh(t, tau, lam.T, cmap='Blues', shading='auto')

        ax2.plot(t, fit['gaussian']['var']**0.5,
                 '.', label=r"$\sigma_\tau(t)$")

        signal = oscillator(t, **fit['oscillator'])
        env = oscillator(t, A=fit['oscillator']['A'],
                         lam=fit['oscillator']['lam'])

        ax2.plot(t, signal)  # , label=eqn)
        ax2.plot(t, fit['oscillator']['mu']+env, 'k--', lw=0.5)
        ax2.plot(t, fit['oscillator']['mu']-env, 'k--', lw=0.5)

        ax2.annotate(eqn, xy=(0, 0), xycoords='axes fraction')

        ax2.set_xlabel(r"$t$ (ms)")
        ax1.set_ylabel(r"$\tau$ (ns)")

        plt.colorbar(pcm, ax=ax1, label=r'$\lambda(\tau)$ (ns$^{-1}$)')

        ax1.plot(t, fit['gaussian']['mu'], label=r"$\mu_{\sigma_\tau}(t)$")
        ax1.legend(loc='upper right')

        ax2.legend(loc='upper right')
        ax2.set_ylabel(r"$\tau$ (ns)")

        N = np.sum(np.gradient(tau)*1e-9*lam, -1)
        ax1.annotate(f"$<N>$:{np.mean(N):0.1e}",
                     xy=(0, 0), xycoords='axes fraction')

    return fit
