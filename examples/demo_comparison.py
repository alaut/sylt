import numpy as np

from sylt.tracking import Bunch, Ring, Tracker

from properties import ring_properties, bunch_properties

ring = Ring(**ring_properties)

NUM_TURNS = 5_000       # number of tracked turns
N = 200e10              # bunch intensity
n = 30_000              # number of macro particles
ng = n - 1_000          # number of ghost particles

# Define Trackers
trackers = {
    'SPM': Tracker(Bunch(**bunch_properties, n=n), ring),
    'SPM+SC': Tracker(Bunch(**bunch_properties, n=n, N=N), ring),
    'SPM+SC+TM': Tracker(Bunch(**bunch_properties, n=n, N=N, eps=None), ring)
}

if __name__ == "__main__":

    for key, tracker in trackers.items():

        # Define "Ghost" Particles
        tracker.UPDATE = False
        tracker.bunch.w[:ng] = np.zeros(ng)
        tracker.bunch.tau[:ng] = np.random.choice(
            np.linspace(0, 1, 300)*tracker.tau_hat, ng)

        # Track Particles
        data = np.empty((2, NUM_TURNS, n))
        for turn in range(NUM_TURNS):
            tracker.track()
            data[:, turn] = tracker.bunch.tau, tracker.bunch.w

        # Clean Distribution
        data[:, :, tracker.clean(tracker.bunch.tau, tracker.bunch.w)] = np.nan

        # Export Data
        np.savez(f"./data/{key}", tau=data[0], w=data[1],
                 N=tracker.bunch.N, eps=tracker.bunch.eps,
                 **ring_properties, **bunch_properties, ng=ng,
                 )
