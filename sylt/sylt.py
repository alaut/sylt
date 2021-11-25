import matplotlib.pyplot as plt
import numpy as np

from sylt.plotting import plot_bunch_profiles
from sylt.analysis import analyze_bunch_profiles
from sylt.tracking import Ring, Bunch, Tracker
from sylt.tools import centers


def benchmark_bunch_profiles(tau, t, lam, sig_eps, show=True, simulate=False):
    """broad analysis of BLO with options to benchmark with simulation"""

    out = {}

    out['exp'] = {}
    out['exp']['fit'] = analyze_bunch_profiles(tau, t, lam)

    if show:
        out['exp']['figs'] = plot_bunch_profiles(tau, t, lam, out['exp']['fit'])

    if simulate:
        for mode in ['SC', 'SC+TM']:
            out[mode] = {}
            eps = None if mode == 'SC+TM' else 0

            N = np.max(out['exp']['fit']['binomial']['amp'])
            print(f"N:{N:0.3e}")

            mu_sig_tau = out['exp']['fit']['oscillator']['mu']

            bunch = Bunch(
                E=2.938272e9,
                sig_w=4e6,
                sig_tau=mu_sig_tau,
                n=40_000,     # 0.5% error
                N=N,
                sig_eps=sig_eps,
                LONG_SHAPE='binomial',
                eps=eps,
            )

            tracker = Tracker(bunch, Ring(), FIXED_MU=True)

            tracker.estimate_voltage(Omega=out['exp']['fit']['oscillator']['omega']/2)

            tracker.match(
                sig_tau=mu_sig_tau,
                k=out['exp']['fit']['oscillator']['A']/mu_sig_tau,
            )

            turns = np.arange(12_000)
            t = turns*tracker.T
            tau = np.linspace(-1, 1, 99)*bunch.sig_tau*3
            dt = np.mean(np.diff(tau))

            LAM = []
            print('tracking ...')
            for turn in turns:
                tracker.track()

                if turn % 1_000 == 0:
                    tracker.clean()

                lam, _ = np.histogram(tracker.bunch.tau, tau)
                LAM.append(lam/dt/tracker.bunch.n*N)

            LAM = np.array(LAM)

            out[mode]['fit'] = analyze_bunch_profiles(centers(tau), t, LAM)

            if show:
                out[mode]['figs'] = plot_bunch_profiles(centers(tau), t, LAM, out[mode]['fit'])

    return out
