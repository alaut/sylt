import matplotlib.pyplot as plt
import numpy as np

from sylt.tracking import Bunch, Ring, Tracker

ring = Ring(
    R=100,          # machine radius (m)
    h=8,            # rf harmonic
    Vg=80e3,        # gap voltage (V)
    gamma_t=6.1,    # transition

    beta=17,
    D=(2.66, 0),
    b=(73e-3, 35e-3),

    # vphi_s=np.radians(10)
)

bunch = Bunch(
    E=2.938272e9,   # particle energy (eV)
    n=50_000,       # number of macroparticles
    sig_w=2e6,      # rms energy width
    sig_tau=30e-9,  # rms bunch length

    N=200e10,
    sig_eps=1e-6,

    eps=None,
)

tracker = Tracker(bunch, ring)

tracker.show('start')
plt.show()

for turn in range(5_000):

    if turn % 100 == 0:
        tracker.clean()

    tracker.track()

tracker.show('stop')

plt.show()

print('Finished !')