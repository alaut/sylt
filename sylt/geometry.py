import numpy as np

from sylt.functions import binomial
def bivariate_binomial(a, b, n, x0=0, y0=0, nx=299, ny=301, mu=1):

    x, y = np.meshgrid(
        x0 + a*np.linspace(-1, 1, nx),
        y0 + b*np.linspace(-1, 1, ny),
    )

    x, y = x.flatten(), y.flatten()

    H = ((x-x0)/a)**2+((y-y0)/b)**2

    # f = 1-r**2
    f = binomial(H, sig=np.sqrt(a*b), mu=mu)

    f[f < 0] = 0

    ind = np.random.choice(f.size, n, p=f/f.sum())

    x = x[ind] + 2*a/nx*(np.random.rand(n)-0.5)
    y = y[ind] + 2*b/ny*(np.random.rand(n)-0.5)

    return x, y
