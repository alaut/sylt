import matplotlib.pyplot as plt

from sylt.tracking import Bunch, Ring, Tracker

from properties import ring_properties, bunch_properties

ring = Ring(**ring_properties)

bunch = Bunch(
    **bunch_properties,
    n=50_000,   # number of macroparticles
    # N=200e10,   # include intensity effects
    eps=None,   # initialize transverse distribution
)

tracker = Tracker(bunch, ring)

tracker.show('start')

tracker.match(k=0.5)
tracker.show('mismatched')

for turn in range(1_000):
    tracker.track()

tracker.clean()
tracker.show('stop')

plt.savefig('./figs/demo.png')
