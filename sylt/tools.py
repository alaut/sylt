import numpy as np
from scipy.stats import norm

def limits(x, n=2):
    """distribute points between limits of array x"""
    return np.linspace(np.nanmin(x), np.nanmax(x), n)

def centers(x):
    """return midpoints of array x"""
    return (x[1:] + x[:-1])/2

def statlim(x, n=2, N=3):
    """distribute points between statistical limits of array x"""
    mu, sig = norm.fit(x)
    return np.linspace(mu-N*sig, mu+N*sig, n)
