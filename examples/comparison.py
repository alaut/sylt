import matplotlib.pyplot as plt

from sylt.tracking import Bunch, Ring, Tracker

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
n = 50_000
sig_w = 0e6
sig_tau = 30e-9

N = 200e10
sig_eps = 1e-6

trackers = {
    'SPM': Tracker(Bunch(E, sig_tau, sig_w, n), ring),
    'SPM+SC': Tracker(Bunch(E, sig_tau, sig_w, n, N, sig_eps), ring),
    'SPM+SC+TM': Tracker(Bunch(E, sig_tau, sig_w, n, N, sig_eps, eps=None), ring)
}

fig, ax = plt.subplots(constrained_layout=True)

for i, (key, tracker) in enumerate(trackers.items()):
    tracker.UPDATE = False

    for turn in range(2_000):
        tracker.track()

    tracker.clean()

    ax.plot(tracker.bunch.tau*1e9, tracker.bunch.w*1e-6,
            f'C{i}.', label=key, alpha=0.5, zorder=2-i)

ts, ws = tracker.separatrix()
ax.plot(ts*1e9, +ws*1e-6, 'm-')
ax.plot(ts*1e9, -ws*1e-6, 'm-')
ax.set_xlabel(r"$\tau$ (ns)")
ax.set_ylabel(r"$w$ (MeV)")
ax.legend()

plt.savefig('./figs/comparison.png')
plt.show()
print('Finished !')
