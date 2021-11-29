import numpy as np

from scipy.special import gamma
import matplotlib.pyplot as plt

def parabolic(x, L):
    y = 3/(2*L)*(1-4*x**2/L**2)
    y[np.abs(x) > L/2] = 0
    return y


def gaussian(dev, var):
    return np.exp(-1/2*dev*dev/var)/np.sqrt(2*np.pi*var)


def binomial(x, sig=1, amp=1, x0=0, mu=1.12):
    dev = x-x0
    L = sig*2*np.sqrt(3+2*mu)
    y = amp*2*gamma(3/2+mu)/(L*np.sqrt(np.pi)*gamma(1+mu)) * \
        (1-4*dev*dev/L**2)**mu
    y[np.abs(dev) > L/2] = 0
    return y.real


def p_binomial(n, sig=1, mu=1, show=False):
    L = sig*2*np.sqrt(3+2*mu)

    x = np.linspace(L/2, 0, 999, endpoint=False)

    f = 2*gamma(3/2+mu)/(L*np.sqrt(np.pi)*gamma(1+mu))*(1-4*x*x/L**2)**mu
    p = f*x

    X = np.random.choice(x, n, p=p/p.sum())
    
    return X


def binomial_der(x, sig=1, amp=1, x0=0, mu=1.12):
    dev = x-x0
    L = sig*2*np.sqrt(3+2*mu)
    dy = -amp*16*mu*dev*gamma(mu+3/2)*(1-4*dev*dev/L **
                                       2)**(mu-1)/np.sqrt(np.pi)/L**3/gamma(mu+1)
    dy[np.abs(dev) > L/2] = 0
    return dy
