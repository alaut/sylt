import numpy as np

from sylt.functions import binomial, p_binomial


def bivariate_binomial(a, b, n, mu=1):
    """return bivariate distribution of widths a, b and density projections of binomial mu"""

    k = np.sqrt((4+mu)/(3+mu))

    u = k*p_binomial(n, sig=1, mu=mu-0.5)

    th = np.random.uniform(-np.pi, np.pi, n)
    x = a*u*np.cos(th)
    y = b*u*np.sin(th)

    return x, y
