import matplotlib.pyplot as plt
import numpy as np

from sylt.plotting import plot_bunch_profiles
from sylt.analysis import analyze_bunch_profiles
from sylt.tracking import Ring, Bunch, Tracker
from sylt.tools import centers


def benchmark_bunch_profiles(tau, t, lam, show=True, simulate=False, sig_eps=3e-6):
    """broad analysis of BLO with options to benchmark with simulation"""

    fit_exp = analyze_bunch_profiles(tau, t, lam)

    if show:
        figs_exp = plot_bunch_profiles(tau, t, lam, fit_exp)

    if simulate:

        N = np.max(fit_exp['gaussian']['amp'] *
                   np.sqrt(2*np.pi*fit_exp['gaussian']['var']))
        print(f"N:{N*1e-10:0.2g} x 1e+10")

        mu_sig_tau = fit_exp['oscillator']['mu']

        bunch = Bunch(
            E=2.938272e9,
            sig_w=4e6,
            sig_tau=mu_sig_tau,
            n=40_000,     # 0.5% error
            N=N, sig_eps=sig_eps,
            eps=None, LONG_SHAPE='binomial',
        )

        tracker = Tracker(bunch, Ring(), FIXED_MU=True)

        tracker.estimate_voltage(Omega=fit_exp['oscillator']['omega']/2)

        tracker.match(
            sig_tau=mu_sig_tau,
            k=fit_exp['oscillator']['A']/mu_sig_tau,
        )

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
            figs_sim = plot_bunch_profiles(centers(tau), t, LAM, fit_sim)

    return {'exp': {'fit': fit_exp, 'figs': figs_exp}, 'sim': {'fit': fit_sim, 'figs': figs_sim}}
