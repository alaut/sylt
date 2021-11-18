import numpy as np
from logging import error
from dataclasses import dataclass
from scipy.constants import value, c

import sylt.plotting

from sylt.functions import parabolic
from sylt.geometry import bivariate_binomial

e = value('elementary charge')
Z_0 = value('characteristic impedance of vacuum')


@dataclass
class Bunch:
    E: float            # particle energy
    sig_tau: float      # width relative time
    sig_w: float = 0    # width relative energy
    n: int = 40_000     # macroparticle number (~0.5% error)

    N: float = 0        # bunch intensity
    sig_eps: float = np.nan  # width emittance

    eps: float = 0       # assume particle emittance is zero

    mu_w: float = 0     # center relative energy
    mu_tau: float = 0   # center relative time

    q: int = 1      # particle charge
    E_0 = 938.272e6  # particle rest energy

    TRAN_SHAPE: str = 'rayleigh'
    LONG_SHAPE: str = 'gaussian'

    def __post_init__(self):
        self.w = np.random.normal(self.mu_w, self.sig_w, self.n)

        if self.LONG_SHAPE == "gaussian":
            self.tau = np.random.normal(self.mu_tau, self.sig_tau, self.n)
        elif self.LONG_SHAPE == "parabolic":
            L = 4*self.sig_tau
            tau = np.linspace(-L/2, L/2, 1000)
            p = parabolic(tau, L)
            self.tau = np.random.choice(tau, size=self.n, p=p/p.sum())
        elif self.LONG_SHAPE == "binomial":
            self.tau, self.w = bivariate_binomial(
                a=2*self.sig_tau, b=2*self.sig_w, n=self.n,
            )
        else:
            error(f'Unrecognized longitudinal shape {self.LONG_SHAPE}!')
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

    def update(self, FIXED_MU=False, FIXED_SIG=False):
        """update representative bunch statistics"""

        if FIXED_MU:
            self.mu_w = np.nanmean(self.w)
            self.mu_tau = np.nanmean(self.tau)

        if FIXED_SIG:
            self.sig_w = np.nanstd(self.w)
            self.sig_tau = np.nanstd(self.tau)

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
            L = 4*self.sig_eps
            eq = np.linspace(L, 0, 1000, endpoint=False)
            pq = parabolic(eq, L)
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

    beta: float = 1
    D: float = 0
    b: float = np.inf

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

    UPDATE: bool = True
    FIXED_MU: bool = False
    FIXED_SIG: bool = False

    def __post_init__(self):
        self.eta = -(1/self.bunch.gamma**2-1/self.ring.gamma_t**2)
        self.omega = self.bunch.beta*c/self.ring.R
        self.T = 2*np.pi/self.omega
        self.Omega = np.sqrt(0j-self.eta/(self.bunch.beta**2*self.bunch.E) *
                             self.bunch.q*self.ring.Vg/self.T*self.ring.h*self.omega*np.cos(self.ring.vphi_s))
        self.kappa = self.eta / (self.bunch.beta**2*self.bunch.E)
        self.nu = self.Omega/self.omega

        self.tau_hat = self.T/self.ring.h/2
        # self.tau_hat = (np.pi-self.ring.vphi_s)/(self.ring.h*self.omega)

    def H(self, tau, w):
        """return particle hamiltonion"""
        phi = self.ring.h*self.omega*tau
        return 0.5*self.kappa*w*w-self.bunch.q*self.ring.Vg/(2*np.pi*self.ring.h)*self.ring.G(phi)

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

        a = self.a()
        a = np.prod(a)**(1/a.size)
        b = np.prod(ring.b)**(1/ring.b.size)

        r2 = np.mean(ring.beta*bunch.eps, 0) + \
            np.sum(ring.D**2*bunch.delta()**2, 0)

        g = 0.5 + np.log(b/a) - 0.5*r2/a**2

        ImZ_n = -g*Z_0/(bunch.beta*bunch.gamma**2)
        return self.V_W(ImZ_n)

    def track(self):
        if self.UPDATE:
            self.bunch.update(self.FIXED_MU, self.FIXED_SIG)

        if self.bunch.N > 0:
            V = self.V_RF() + self.V_SC()
        else:
            V = self.V_RF()

        self.bunch.w = self.bunch.w + self.bunch.q*V
        self.bunch.tau = self.bunch.tau + self.T*self.kappa*self.bunch.w

    def clean(self, tau=None, w=None):
        """assign NaN to particle's who's longitudinal position is external to the separatrix"""
        tau = self.bunch.tau if tau is None else tau
        w = self.bunch.w if w is None else w

        lost = ((self.H(self.tau_hat, 0) - self.H(tau, w)) > 0) + \
            (np.abs(tau) > self.tau_hat)

        self.bunch.tau[lost] = np.nan
        self.bunch.w[lost] = np.nan

        return lost

    def a(self):
        bunch = self.bunch
        ring = self.ring
        a = 2*np.sqrt(ring.beta*bunch.sig_eps +
                      ring.D**2*bunch.sig_delta()**2)
        return a

    def separatrix(self, phi_hat=None, phi=None):
        if phi_hat is None:
            phi_hat = np.pi-self.ring.vphi_s
        if phi is None:
            phi = np.linspace(-phi_hat, phi_hat, 299)

        phi_dot = self.Omega * \
            np.sqrt(0j+2*(self.ring.W(phi_hat)-self.ring.W(phi)))
        tau_dot = phi_dot/(self.ring.h*self.omega)
        w = tau_dot/self.kappa
        return phi/(self.ring.h*self.omega), w.real

    def show(self, title='distribution', alpha=0):
        """show distribution assuming stable transverse optics"""

        ring = self.ring
        bunch = self.bunch

        mu = np.array([
            np.random.uniform(-np.pi, np.pi, bunch.n),
            np.random.uniform(-np.pi, np.pi, bunch.n),
        ])

        u = np.sqrt(ring.beta*bunch.eps)*np.cos(mu)+ring.D*bunch.delta()
        up = -np.sqrt(bunch.eps/ring.beta)*(alpha*np.cos(mu)+np.sin(mu))

        data = {
            "x (mm)": u[0]*1e3, "x' (mrad)": up[0]*1e3,
            "y (mm)": u[1]*1e3, "y' (mrad)": up[1]*1e3,
            "tau (ns)": bunch.tau*1e9, "w (MeV)": bunch.w*1e-6}

        keys = [
            ['x (mm)', "x' (mrad)"],
            ['y (mm)', "y' (mrad)"],
            ["x (mm)", "y (mm)"],
            ["tau (ns)", 'w (MeV)'],
        ]

        axes = plotting.plot_phase_space(data, keys, (2, 2), f"{title}")

        ax = axes.flatten()[-1]
        for phi_max in [ring.h*self.omega*bunch.sig_tau, np.pi-ring.vphi_s]:
            xs, ys = self.separatrix(phi_max)
            ax.plot(xs*1e9, +ys.real*1e-6, 'm-')
            ax.plot(xs*1e9, -ys.real*1e-6, 'm-')

        axv = ax.secondary_yaxis('right', functions=(self.MeV2mm, self.mm2MeV))
        ax.annotate(r'$x$ (mm)', (1, 0.5), xycoords='axes fraction',
                    rotation=90, horizontalalignment='right', verticalalignment='center')

        return axes

    def mm2MeV(self, mm):
        C = self.bunch.beta**2*self.bunch.E/self.ring.alpha/self.ring.R
        return mm*1e-3*C*1e-6

    def MeV2mm(self, MeV):
        C = self.bunch.beta**2*self.bunch.E/self.ring.alpha/self.ring.R
        return MeV*1e6/C*1e3

    def match(self, sig_w=None, sig_tau=None, k=0):
        ring = self.ring
        bunch = self.bunch

        if sig_tau is not None:
            sig_phi = ring.h*self.omega*sig_tau
            bunch.sig_w = np.sqrt(bunch.q*ring.Vg/(np.abs(self.kappa)*np.pi*ring.h)
                                  * (1-np.cos(sig_phi)))
        if sig_w is not None:
            pass

        bunch.sig_tau = bunch.sig_tau*(1+k)
        bunch.sig_w = bunch.sig_w*(1-k)
        bunch.__post_init__()

        text = [
            f'matching k={k:0.2f}',
            f"sig_tau: {bunch.sig_tau*1e9:0.2f} ns",
            f"sig_w: {bunch.sig_w*1e-6:0.2f} MeV",
            '',
        ]
        print('\n'.join(text))

    def estimate_voltage(self, Omega):
        """estimate effective voltage given observed synchrotron frequency"""
        ring = self.ring
        bunch = self.bunch
        print(f"estimating voltage from Omega:{Omega*1e-3:0.3g} kHz")
        V_eff = -(Omega*self.T)**2/(self.kappa*bunch.q *
                                    2*np.pi*ring.h*np.cos(ring.vphi_s))
        print(f"V_eff:{V_eff*1e-3:0.3g} kV")
        mu_eff = 1-(ring.h*self.omega*bunch.sig_tau)**2/16
        print(f"mu_eff:{mu_eff:0.3g}")
        Vg = V_eff/mu_eff**2
        print(f"Vg:{Vg*1e-3:0.3g} kV")
        return Vg
