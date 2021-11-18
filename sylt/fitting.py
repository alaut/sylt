from scipy.optimize import curve_fit

import numpy as np


def oscillator(t, A, omega=0, mu=0, phi=0, lam=0):
    """return signal for damped oscillator with complex amplitude amp, fundamental omega and decay lam """
    return A*np.cos(omega*t+phi)*np.exp(-lam*t)+mu


def gauss(x, mu, var, amp):
    """return pdf(x) for gaussian distribution"""
    dev = x-mu
    return amp*np.exp(-0.5*dev*dev/var)


def moments_array(x, y):
    norm = 1/np.sum(y, -1)
    mu = norm*np.sum(x*y, -1)
    dev = x-mu[:, np.newaxis]
    var = norm*np.sum(y*dev*dev, -1)
    return mu, var


def moments(x, y):
    """compute first and second moment given x-y profile"""
    norm = 1/np.sum(y)
    mu = norm*np.sum(x*y)
    dev = x-mu
    var = norm*np.sum(y*dev*dev)
    return mu, var


def fit_gaussian(x, Y):
    """given x-y profile, return fit parameters for a gaussian function"""
    MU, VAR = moments_array(x, Y)
    AMP = np.max(Y, -1)
    fit = []
    for y, mu, var, amp in zip(Y, MU, VAR, AMP):
        try:
            popt, pcov = curve_fit(gauss, x, y, p0=[mu, var, amp])
        except:
            popt = [np.nan, np.nan, np.nan]
        fit.append(popt)
    fit = np.array(fit)
    return {'mu': fit[:, 0], 'var': fit[:, 1], 'amp': fit[:, 2]}


def fit_oscillator(x, y, omega=None):
    """ given damped oscillation y along x, return signal envelope"""

    y0 = np.mean(y)

    if omega is None:
        amp = np.abs(np.fft.rfft(y-y0))
        freq = np.fft.rfftfreq(x.size, d=np.mean(np.gradient(x)))
        omega = 2*np.pi*freq[amp == amp.max()][0]

    p0 = [np.max(y-y0), omega, y0, 0, 0]
    (A, omega, mu, phi, lam), _ = curve_fit(oscillator, x, y, p0)
    
    return {'A': A, 'omega': omega, 'lam': lam, 'phi': phi, 'mu': mu}
