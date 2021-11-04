import numpy as np

from logging import error

from dataclasses import dataclass

from scipy.constants import value, c
e = value('elementary charge')
Z_0 = value('characteristic impedance of vacuum')


def lam_parabolic(x, L):
    y = 3/(2*L)*(1-4*x**2/L**2)
    y[np.abs(x) > L/2] = 0
    return y


@dataclass
class Bunch:
    E: float    # particle energy
    sig_tau: float  # width relative time
    sig_w: float = 0    # width relative energy
    n: int = 10_000      # macroparticle number

    N: float = 0    # bunch intensity
    sig_eps: float = None  # width emittance

    eps: float = 0       # assume particle emittance is zero

    mu_w: float = 0     # center relative energy
    mu_tau: float = 0   # center relative time

    q: int = 1      # particle charge
    E_0 = 938.272e6  # particle rest energy

    FIXED_MU: bool = False
    FIXED_SIG: bool = False

    TRAN_SHAPE: str = 'rayleigh'
    LONG_SHAPE: str = 'gaussian'

    def __post_init__(self):
        self.w = np.random.normal(self.mu_w, self.sig_w, self.n)
        self.tau = np.random.normal(self.mu_tau, self.sig_tau, self.n)

        if self.eps is None:
            self.populate_emittance()

        self.gamma = self.E/self.E_0
        self.beta = np.sqrt(1-self.gamma**-2)

    def delta(self):
        """return particle momentum spread"""
        return self.w/(self.beta**2*self.E)

    def sig_delta(self):
        """return rms momentum spread"""
        return self.sig_w/(self.beta**2*self.E)

    def derivative(self):
        """compute particle derivative assuming gaussian longitudinal profile"""
        dev = self.tau-self.mu_tau
        var = self.sig_tau**2
        lam = np.exp(-0.5*dev*dev/var)/np.sqrt(2*np.pi*var)
        return -dev*lam/var

    def populate_emittance(self):
        if self.TRAN_SHAPE == 'rayleigh':
            self.eps = np.array([
                np.random.rayleigh(self.sig_eps**0.5, self.n)**2,
                np.random.rayleigh(self.sig_eps**0.5, self.n)**2,
            ])
        elif self.TRAN_SHAPE == "exponential":
            self.eps = np.array([
                np.random.exponential(self.sig_eps*2, self.n),
                np.random.exponential(self.sig_eps*2, self.n),
            ])
        elif self.TRAN_SHAPE == 'parabolic':
            eq = np.linspace(4*self.sig_eps, 0, 1000, endpoint=False)
            pq = lam_parabolic(eq, self.sig_eps*4)  # /epsq
            self.eps = np.array([
                np.random.choice(eq, self.n, p=pq/pq.sum()),
                np.random.choice(eq, self.n, p=pq/pq.sum()),
            ])
        elif self.TRAN_SHAPE == "uniform":

            A = self.sig_eps*4

            th = np.random.uniform(0, 2*np.pi, self.n)
            r = np.random.uniform(0, A, self.n)**0.5

            x = r*np.cos(th)
            y = r*np.sin(th)

            phi = np.random.uniform(0, 2*np.pi, self.n)
            xp = np.sqrt(A-r**2)*np.cos(phi)
            yp = np.sqrt(A-r**2)*np.sin(phi)

            self.eps = (x**2+xp**2, y**2+yp**2)
        else:
            error(f'Unrecognized style <{self.TRAN_SHAPE}>')


@dataclass
class Ring:
    R: float
    h: np.ndarray
    Vg: np.ndarray
    gamma_t: float

    beta: float = None
    D: float = None
    b: float = None

    vphi_s: np.ndarray = 0

    def __post_init__(self):
        self.alpha = self.gamma_t**-2
        self.beta = np.array([self.beta]).T
        self.D = np.array([self.D]).T
        self.b = np.array([self.b]).T

    def G(self, phi):
        return np.cos(self.vphi_s) - np.cos(phi+self.vphi_s) - phi*np.sin(self.vphi_s)

    def W(self, phi):
        return self.G(phi)/np.cos(self.vphi_s)


@dataclass
class Tracker:
    bunch: object
    ring: object

    UPDATE: bool = False
    FIXED_MU: bool = False
    FIXED_SIG: bool = False

    a: float = None
    g: float = None

    def __post_init__(self):
        self.eta = -(1/self.bunch.gamma**2-1/self.ring.gamma_t**2)
        self.omega = self.bunch.beta*c/self.ring.R
        self.T = 2*np.pi/self.omega
        self.Omega = np.sqrt(0j-self.eta/(self.bunch.beta**2*self.bunch.E) *
                             self.bunch.q*self.ring.Vg/self.T*self.ring.h*self.omega*np.cos(self.ring.vphi_s))
        self.kappa = self.eta / (self.bunch.beta**2*self.bunch.E)
        self.nu = self.Omega/self.omega

    def V_RF(self):
        """compute RF voltage"""
        vphi = self.ring.h*self.omega*self.bunch.tau + self.ring.vphi_s
        return self.ring.Vg*(np.sin(vphi)-np.sin(self.ring.vphi_s))

    def V_W(self, ImZ_n):
        """compute wakefield voltage"""
        return -e*self.bunch.N/self.omega*ImZ_n*self.bunch.derivative()

    def V_SC(self):
        """compute space-charge voltage"""
        bunch = self.bunch
        ring = self.ring

        a = np.prod(self.a)**(1/self.a.size)
        b = np.prod(ring.b)**(1/ring.b.size)

        r2 = np.mean(ring.beta*bunch.eps, 0) + \
            np.sum(ring.D**2*bunch.delta()**2, 0)

        g = 0.5 + np.log(b/a) - 0.5*r2/a**2

        ImZ_n = -g*Z_0/(bunch.beta*bunch.gamma**2)
        return self.V_W(ImZ_n)

    def track(self):
        if self.UPDATE:
            self.update()

        if self.bunch.N > 0:
            V = self.V_RF() + self.V_SC()
        else:
            V = self.V_RF()

        self.bunch.w = self.bunch.w + self.bunch.q*V
        self.bunch.tau = self.bunch.tau + self.T*self.kappa*self.bunch.w

    def update(self):
        bunch = self.bunch
        ring = self.ring

        if not self.FIXED_MU:
            bunch.mu_w = np.nanmean(bunch.w)
            bunch.mu_tau = np.nanmean(bunch.tau)

        if not self.FIXED_SIG:
            bunch.sig_w = np.nanstd(bunch.w)
            bunch.sig_tau = np.nanstd(bunch.tau)

        if bunch.sig_eps is not None:
            self.a = 2*np.sqrt(ring.beta*bunch.sig_eps +
                               ring.D**2*bunch.sig_delta()**2)

    def separatrix(self, tau_hat):
        tau = np.linspace(-tau_hat, tau_hat, 999)

        phi = self.ring.h*self.omega*tau
        phi_hat = self.ring.h*self.omega*tau_hat

        phi_dot = self.Omega*np.sqrt(2*(self.ring.W(phi_hat)-self.ring.W(phi)))
        tau_dot = phi_dot/(self.ring.h*self.omega)
        w = tau_dot/self.kappa
        return tau, w

    def match(self, sig_w=None, sig_tau=None, k=0):
        ring = self.ring
        bunch = self.bunch

        if sig_tau is not None:
            sig_phi = ring.h*self.omega*sig_tau
            sig_w = np.sqrt(bunch.q*ring.Vg/(np.abs(self.kappa)*np.pi*ring.h)
                            * (1-np.cos(sig_phi)))

        bunch.sig_tau = sig_tau*(1+k)
        bunch.sig_w = sig_w*(1-k)
        bunch.__post_init__()

        text = [
            f'matching k={k:0.2f}',
            f"sig_tau: {bunch.sig_tau*1e9:0.2f} ns",
            f"sig_w: {bunch.sig_w*1e-6:0.2f} MeV",
            '',
        ]
        print('\n'.join(text))
