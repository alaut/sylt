import numpy as np
from scipy.stats import norm

from sylt.tools import centers

def twiss(x, y):
    """return twiss parameters and phase-space ellipse for bivariate distribution"""

    cov = np.cov(x, y)
    area = np.sqrt(np.linalg.det(cov))

    twiss_parameters = {
        'emittance': area,
        'alpha': -cov[0, 1]/area,
        'beta': cov[0, 0]/area,
        'gamma': cov[1, 1]/area,
    }

    ellipse_params = covariance_ellipse(cov)

    ellipse = rotated_ellipse(
        **ellipse_params,
        dx=np.nanmean(x),
        dy=np.nanmean(y),
    )

    caption = '\n'.join([
        f"$\\epsilon$={twiss_parameters['emittance']:0.2e}",
        f"$\\alpha$={twiss_parameters['alpha']:0.2f}",
        f"$\\beta$={twiss_parameters['beta']:0.2f}",
        f"$\\gamma$={twiss_parameters['gamma']:0.2f}",
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
    return xc, h, norm.pdf(xc, *norm.fit(x))
