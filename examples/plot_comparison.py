import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import numpy as np
from glob import glob

import os

from demo_comparison import ng, trackers

files = sorted(glob('./data/*.npz'))[::-1]

# Load Data into Memory
data = {}
for file in files:
    path, ext = os.path.splitext(file)
    key = os.path.basename(path)
    tmp = np.load(file)
    data[key] = {'tau': tmp['tau']*1e9, 'w': tmp['w']*1e-6}

NUM_TURNS, n = data['SPM']['tau'].shape

fig, ax = plt.subplots(num='comparison', constrained_layout=True)

artists = []


def init(turn=0):

    for i, (key, tracker) in enumerate(trackers.items()):
        ax.plot([], [], f"C{i}", label=key)

        tau = data[key]['tau'][turn]
        w = data[key]['w'][turn]

        artists.append([
            ax.plot(tau[:ng], w[:ng], f'C{i}.', zorder=5-i)[0],
            ax.plot(tau[ng:], w[ng:], f'C{i},', zorder=2-i)[0],
        ])

        ts, ws = tracker.separatrix()
        ax.plot(ts*1e9, +ws*1e-6, 'm-')
        ax.plot(ts*1e9, -ws*1e-6, 'm-')

    ax.set_xlabel(r"$\tau$ (ns)")
    ax.set_ylabel(r"$w$ (MeV)")
    ax.legend()


init(2_500)

fig.savefig('./figs/comparison.png')


def update(turn):
    for key, (p1, p2) in zip(trackers.keys(), artists):
        tau = data[key]['tau'][turn]
        w = data[key]['w'][turn]
        p1.set_data(tau[0:ng], w[0:ng])
        p2.set_data(tau[ng:], w[ng:])
    return sum(artists, [])


ani = FuncAnimation(fig, update, frames=range(
    0, NUM_TURNS, 20), blit=True, interval=0)

plt.show()
ani.save('./figs/comparison.gif', fps=60)
