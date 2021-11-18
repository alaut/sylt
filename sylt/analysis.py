import numpy as np
from scipy.stats import norm

from sylt.tools import centers
from sylt.fitting import fit_gaussian, fit_oscillator
import sylt.tracking as tracking
import sylt.plotting as plotting

def twiss(x, y):
    """return twiss parameters and phase-space ellipse for bivariate distribution"""

    cov = np.cov(x, y)
    area = np.sqrt(np.linalg.det(cov))

    with np.errstate(divide='ignore', invalid='ignore'):
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
def analyze_bunch_profiles(tau, t, lam, ):
    """return BLO fit parameters"""
    params_g = fit_gaussian(tau, lam)
    params_o = fit_oscillator(t, params_g['var']**0.5)
    return {'gaussian': params_g, 'oscillator': params_o}


def benchmark_bunch_profiles(tau, t, lam, show=True, simulate=False, sig_eps=3e-6):
    """broad analysis of BLO with options to benchmark with simulation"""

    fit = analyze_bunch_profiles(tau, t, lam)

    if show:
        plotting.plot_bunch_profiles(tau, t, lam, fit)

    if simulate:

        N = np.max(fit['gaussian']['amp'] *
                   np.sqrt(2*np.pi*fit['gaussian']['var']))
        print(f"N:{N*1e-10:0.2g} x 1e+10")

        bunch = tracking.Bunch(
            E=2.938272e9, sig_w=4e6, sig_tau=fit['oscillator']['mu'],
            N=N, sig_eps=sig_eps,
            eps=None, LONG_SHAPE='binomial',
        )

        ring = tracking.Ring(
            R=100, h=8, Vg=80e3, vphi_s=0, gamma_t=6.1,
            beta=17, D=(2.66, 0), b=(73e-3, 35e-3),
        )

        tracker = tracking.Tracker(bunch, ring, FIXED_MU=True)

        Omega = fit['oscillator']['omega']/2

        ring.Vg = tracker.estimate_voltage(Omega)

        tracker.__post_init__()
        k = fit['oscillator']['A']/fit['oscillator']['mu']
        tracker.match(sig_tau=fit['oscillator']['mu'], k=k)

        turns = np.arange(12_000)
        t = turns*tracker.T
        tau = np.linspace(-1, 1, 99)*bunch.sig_tau*3
        dt = np.mean(np.diff(tau))

        LAM = []
        tracker.show('start')
        for turn in turns:
            tracker.track()

            if turn % 1_000 == 0:
                lost = tracker.clean()
                print(f"lost {lost.sum()/lost.size*100:0.2g} % of particles")

            # if turn in range(0, NUM_TURNS, 2500):
            #     tracker.show(f"{turn}")

            lam, _ = np.histogram(tracker.bunch.tau, tau)
            LAM.append(lam/dt/tracker.bunch.n*N)
        tracker.show('stop')
        LAM = np.array(LAM)

        fit_sim = analyze_bunch_profiles(centers(tau), t, LAM)

        if show:
            plotting.plot_bunch_profiles(centers(tau), t, LAM, fit_sim)

    return fit, fit_sim
