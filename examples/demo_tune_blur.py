import matplotlib.pyplot as plt
import numpy as np

from sylt.tracking import Bunch, Ring, Tracker

from sylt.analysis import compute_synchrotron_frequency
from sylt.plotting import plot_tune

from properties import ring_properties, bunch_properties

ring = Ring(**ring_properties)

N = 400e10      # bunch intensity
n = 500         # number of macroparticles

trackers = {
    'SPM': Tracker(Bunch(**bunch_properties, n=n), ring),
    'SPM+SC': Tracker(Bunch(**bunch_properties, n=n, N=N), ring),
    'SPM+SC+TM': Tracker(Bunch(**bunch_properties, n=n, N=N, eps=None), ring)
}

DATA = {}
for i, (key, tracker) in enumerate(trackers.items()):
    tracker.UPDATE = False

    tau = []
    for turn in range(500_000):
        tracker.track()
        tau.append(tracker.bunch.tau)

    tau_hat, f = compute_synchrotron_frequency(np.array(tau), tracker.T)

    DATA[key] = {
        'phi_hat': ring.h*tracker.omega * tau_hat,
        'mu': 2*np.pi*f/tracker.Omega.real,
    }

plot_tune(DATA)

for ext in ['png', 'pdf', 'svg']:
    plt.savefig(f'./figs/tune.{ext}')

plt.show()

print('Finished !')
