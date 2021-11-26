import numpy as np
from scipy.stats import norm

from sylt.fitting import fit_gaussian, fit_oscillator, fit_binomial
from sylt.tools import centers


def twiss(x, y):
    """return twiss parameters and phase-space ellipse for bivariate distribution"""

    cov = np.cov(x, y)
    eps = np.sqrt(np.linalg.det(cov))

    with np.errstate(divide='ignore', invalid='ignore'):
        twiss_parameters = {
            'emittance': eps,
            'alpha': -cov[0, 1]/eps,
            'beta': cov[0, 0]/eps,
            'gamma': cov[1, 1]/eps,
        }

    ellipse_params = covariance_ellipse(cov)

    ellipse = rotated_ellipse(
        **ellipse_params,
        dx=np.nanmean(x),
        dy=np.nanmean(y),
    )

    caption = '\n'.join([
        f"$\\epsilon$={twiss_parameters['emittance']:0.3g}",
        f"$\\alpha$={twiss_parameters['alpha']:0.3g}",
        f"$\\beta$={twiss_parameters['beta']:0.3g}",
        f"$\\gamma$={twiss_parameters['gamma']:0.3g}",
    ])

    return twiss_parameters, ellipse, caption


def covariance_ellipse(cov):
    """return ellipse parameters a, b, theta from covariance matrix
    https://cookierobotics.com/007/"""

    a, b, c = cov[0, 0], cov[0, 1], cov[1, 1]

    lam_1 = (a + c) / 2 + np.sqrt(((a - c) / 2) ** 2 + b ** 2)
    lam_2 = (a + c) / 2 - np.sqrt(((a - c) / 2) ** 2 + b ** 2)

    theta = (np.pi / 2 if a < c else 0) if b == 0 else np.arctan2(lam_1 - a, b)

    return {'a': lam_1**0.5, 'b': lam_2**0.5, 'theta': theta}


def rotated_ellipse(a=1, b=1, theta=0, verbose=False, dx=0, dy=0):
    """return x(t), y(t) for rotated ellipse(a, b, theta)
    https://en.wikipedia.org/wiki/Ellipse"""

    if verbose:
        print('\n'.join([
            f'a:{a}', f'b:{b}',
            f'phi:{theta*180/np.pi} deg',
            f'dx:{dx}', f"dy:{dy}",
        ]))

    def x(t): return a*np.cos(theta)*np.cos(t)-b*np.sin(theta)*np.sin(t)+dx
    def y(t): return a*np.sin(theta)*np.cos(t)+b*np.cos(theta)*np.sin(t)+dy

    return x, y


def project(x, bins=100):
    """return domain, hist and gaussian fit of x"""
    h, xe = np.histogram(x, bins, density=True)
    xc = centers(xe)
    with np.errstate(divide='ignore'):
        fit = norm.pdf(xc, *norm.fit(x))
    return xc, h, fit


def compute_synchrotron_frequency(tau, T):
    """given .npz data, compute fft along phase coordiant and estimate synchrotron frequency"""
    amp = np.fft.rfft(tau, axis=0)
    freq = np.fft.rfftfreq(n=tau.shape[0], d=T)
    f = freq[np.argmax(amp, axis=0)]
    tau_hat = np.max(tau, axis=0)
    return tau_hat, f


def analyze_bunch_profiles(tau, t, lam, show=False):
    """return BLO fit parameters"""
    params_g = fit_gaussian(tau, lam)
    params_b = fit_binomial(tau, lam)
    # params_o = fit_oscillator(t, params_g['var']**0.5)
    params_o = fit_oscillator(t, params_b['sig'])
    return {'gaussian': params_g, 'oscillator': params_o, 'binomial': params_b}
