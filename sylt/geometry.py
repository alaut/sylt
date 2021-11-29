import numpy as np

from sylt.functions import binomial, p_binomial


def bivariate_binomial(a, b, n, mu=1):

    k = np.sqrt((4+mu*2)/(3+mu*2))

    u = p_binomial(n, sig=1, mu=mu)

    th = np.random.uniform(-np.pi, np.pi, n)
    x = k*a*u*np.cos(th)
    y = k*b*u*np.sin(th)

    print(f"sig_x:{np.nanstd(x):0.3f}\tsig_y:{np.nanstd(y):0.3f}")

    return x, y
