# sylt

**SY**nchrotron **L**ongitudinal **T**racker

## Implementation

The __Bunch__, __Ring__, and __Tracker__ objects must be imported.

```python
from sylt.tracking import Bunch, Ring, Tracker
```

A basic ring can be instantiated by

```python
ring = Ring(
    R=100,          # machine radius (m)
    h=8,            # rf harmonic
    Vg=80e3,        # gap voltage (V)
    gamma_t=6.1,    # transition
    )
```

A bivariate gaussian bunch can be generated with

```python
bunch = Bunch(
    E=2.938272e9,   # particle energy (eV)
    n=50e3,         # number of macroparticles
    sig_w=4e6,      # rms energy width
    sig_tau=30e-9,  # rms bunch length
    )
```

To incorporate intensity effects, a nonzero beam intensity must be defined.

```python
bunch.N = 200e12
```

To incorporate the effects of transverse motion, the ring's transverse optics and aperture can be specified by

```python
ring.beta = 17      # beta function
ring.D = 2.66       # dispersion function
ring.b = 50e-3      # beam pipe aperture radius
```
The transverse beam size will be defined by the rms transverse emittance

```python
bunch.sig_eps = 3e-6    # rms transverse emittance (m)
```

To track, define a tracker using __Tracker__ and use the __track__ method.

```python
tracker = Tracker(bunch, ring)

for turn in range(10_000):
    tracker.track()
```

At any point, the bunch's phase space distribution can be visualized using __Tracker.show__.

```python
tracker.show()
```

## Contributors
- [Alexander Laut](https://alaut.gihub.io/sylt)