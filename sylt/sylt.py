import matplotlib.pyplot as plt
import numpy as np

from sylt.plotting import plot_bunch_profiles
from sylt.analysis import analyze_bunch_profiles
from sylt.tracking import Ring, Bunch, Tracker
from sylt.tools import centers


def benchmark_bunch_profiles(tau, t, lam, sig_eps, show=True, NUM_TURNS=12_000):
    """broad analysis of BLO with options to benchmark with simulation"""

    out = {}

    out['exp'] = {}
    out['exp']['fit'] = analyze_bunch_profiles(tau, t, lam)

    if show:
        out['exp']['figs'] = plot_bunch_profiles(
            tau, t, lam, out['exp']['fit'])

    if NUM_TURNS > 0:
        for mode in ['SC+TM']:  # , 'SC']:
            out[mode] = {}
            eps = None if mode == 'SC+TM' else 0

            N = np.max(out['exp']['fit']['binomial']['amp'])
            print(f"N:{N:0.3e}")

            mu_sig_tau = out['exp']['fit']['oscillator']['mu']

            bunch = Bunch(
                E=2.938272e9,
                sig_w=4e6,
                sig_tau=mu_sig_tau,
                n=40_000,
                N=N,
                sig_eps=sig_eps,
                eps=eps,
                FROZEN_MEAN=True,
                # FUDGE_FACTOR=0.87,
            )

            tracker = Tracker(bunch, Ring())

            tracker.estimate_voltage(
                Omega=out['exp']['fit']['oscillator']['omega']/2)

            tracker.match(
                sig_tau=mu_sig_tau,
                k=out['exp']['fit']['oscillator']['A']/mu_sig_tau,
            )

            tracker.show(f'{mode}-start')

            turns = np.arange(NUM_TURNS)
            t = turns*tracker.T
            tau = np.linspace(-1, 1, 99)*bunch.sig_tau*3
            dt = np.mean(np.diff(tau))

            LAM = []
            for turn in range(NUM_TURNS):
                tracker.track()

                if turn % 1_000 == 0:
                    tracker.clean()
                    msg = f"{mode}-{turn} of {NUM_TURNS}"
                    print(msg)
                    # tracker.show(msg)

                lam, _ = np.histogram(tracker.bunch.tau, tau)
                LAM.append(lam/dt/tracker.bunch.n*N)
            LAM = np.array(LAM)

            tracker.show(f'{mode}-stop')

            out[mode]['fit'] = analyze_bunch_profiles(centers(tau), t, LAM)

            if show:
                out[mode]['figs'] = plot_bunch_profiles(
                    centers(tau), t, LAM, out[mode]['fit'])

    return out
