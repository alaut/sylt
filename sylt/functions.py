import numpy as np

from scipy.special import gamma

def parabolic(x, L):
    y = 3/(2*L)*(1-4*x**2/L**2)
    y[np.abs(x) > L/2] = 0
    return y

def gaussian(dev, var):
    return np.exp(-1/2*dev*dev/var)/np.sqrt(2*np.pi*var)

def binomial(x, x0, L, amp, mu=1):
    dev = x-x0
    y = amp*2*gamma(3/2+mu)/(L*np.sqrt(np.pi)*gamma(1+mu))*(1-4*dev*dev/L**2)**mu
    y[np.abs(dev) > L/2] = 0
    return y
