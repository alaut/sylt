import matplotlib.pyplot as plt
import numpy as np
from sylt.analysis import compute_synchrotron_frequency

from sylt.tracking import Bunch, Ring, Tracker
from sylt.plotting import plot_tune

ring = Ring(
    R=100,          # machine radius (m)
    h=8,            # rf harmonic
    Vg=80e3,        # gap voltage (V)
    gamma_t=6.1,    # transition

    beta=17,
    D=(2.66, 0),
    b=(73e-3, 35e-3),
)

E = 2e9+938.272e6
n = 500
sig_w = 0e6
sig_tau = 30e-9

N = 400e10
sig_eps = 1e-6

trackers = {
    'SPM': Tracker(Bunch(E, sig_tau, sig_w, n), ring),
    'SPM+SC': Tracker(Bunch(E, sig_tau, sig_w, n, N, sig_eps), ring),
    'SPM+SC+TM': Tracker(Bunch(E, sig_tau, sig_w, n, N, sig_eps, eps=None), ring)
}

DATA = {}
for i, (key, tracker) in enumerate(trackers.items()):
    tracker.UPDATE = False

    tau = []
    for turn in range(5_000):
        tracker.track()
        tau.append(tracker.bunch.tau)

    tau_hat, f = compute_synchrotron_frequency(np.array(tau), tracker.T)

    DATA[key] = {
        'phi_hat': ring.h*tracker.omega * tau_hat,
        'mu': 2*np.pi*f/tracker.Omega.real,
    }

plot_tune(DATA)

plt.show()

print('Finished !')
