import numpy as np

from dataclasses import dataclass
from scipy.constants import value, c

from sylt.plotting import plot_phase_space
from sylt.functions import binomial, binomial_der, p_binomial
from sylt.geometry import bivariate_binomial

e = value('elementary charge')
Z_0 = value('characteristic impedance of vacuum')


@dataclass
class Bunch:
    E: float                # particle energy
    sig_tau: float = 30e-9    # width relative time
    sig_w: float = 4e6      # width relative energy
    n: int = 40_000         # macroparticle number (~0.5% error)

    N: float = 0        # bunch intensity
    sig_eps: float = np.nan  # width emittance

    eps: float = 0       # assume particle emittance is zero

    mu_w: float = 0     # center relative energy
    mu_tau: float = 0   # center relative time

    q: int = 1      # particle charge
    E_0 = 938.272e6  # particle rest energy

    mu_l: int = 1.1 # quasi-parabolic
    mu_t: int = 1

    def __post_init__(self):

        self.tau, self.w = bivariate_binomial(
            a=self.sig_tau,
            b=self.sig_w,
            n=self.n,
            mu=self.mu_l,
        )

        self.gamma = self.E/self.E_0
        self.beta = np.sqrt(1-self.gamma**-2)

    def delta(self):
        """return particle momentum spread"""
        return self.w/(self.beta**2*self.E)

    def sig_delta(self):
        """return rms momentum spread"""
        return self.sig_w/(self.beta**2*self.E)

    def derivative(self):
        """compute analytic particle distribution derivative"""
        dev = self.tau-self.mu_tau
        dlam = binomial_der(dev, sig=self.sig_tau, amp=1, mu=self.mu_l)
        return dlam


    def update(self, FIXED_MU=False, FIXED_SIG=False):
        """update representative bunch statistics"""

        if not FIXED_MU:
            self.mu_w = np.nanmean(self.w)
            self.mu_tau = np.nanmean(self.tau)

        if not FIXED_SIG:
            self.sig_w = np.nanstd(self.w)
            self.sig_tau = np.nanstd(self.tau)

@dataclass
class Ring:
    R: float = 100
    h: int = 8
    Vg: float = 80e3
    gamma_t: float = 6.1

    beta = (17, 17)
    beta_p = (0, 0)

    D = (2.66, 0)
    b = (73e-3, 35e-3)

    vphi_s: np.ndarray = 0

    def __post_init__(self):

        self.beta = np.array([self.beta]).T
        self.beta_p = np.array([self.beta_p]).T
        self.D = np.array([self.D]).T
        self.b = np.array([self.b]).T

        # twiss parameters
        self.alpha = -0.5*self.beta_p
        self.gamma = (1+self.alpha**2)/self.beta

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

        if self.bunch.eps is None:
            self.populate_emittance()

    def populate_emittance(self, alpha=0):
        bunch = self.bunch
        ring = self.ring

        bunch.eps = np.array([
            np.random.rayleigh(bunch.sig_eps**0.5, bunch.n)**2,
            np.random.rayleigh(bunch.sig_eps**0.5, bunch.n)**2,
        ])

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
        a = 2*((ring.beta*bunch.sig_eps +
                ring.D**2*bunch.sig_delta()**2)*(ring.beta*bunch.sig_eps))**(1/4)
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
            np.random.uniform(0, 2*np.pi, bunch.n),
            np.random.uniform(0, 2*np.pi, bunch.n),
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

        fig, axes = plot_phase_space(data, keys, (2, 2), f"{title}")

        ax = axes.flatten()[-1]
        for phi_max in [ring.h*self.omega*bunch.sig_tau, np.pi-ring.vphi_s]:
            xs, ys = self.separatrix(phi_max)
            ax.plot(xs*1e9, +ys.real*1e-6, 'm-')
            ax.plot(xs*1e9, -ys.real*1e-6, 'm-')

        ax.secondary_yaxis('right', functions=(self.MeV2mm, self.mm2MeV))
        ax.annotate(r'$x$ (mm)', (1, 0.5), xycoords='axes fraction',
                    rotation=90, horizontalalignment='right', verticalalignment='center')

        return fig, axes

    def mm2MeV(self, mm):
        C = self.bunch.beta**2*self.bunch.E/self.ring.R*self.ring.gamma_t**2
        return mm*1e-3*C*1e-6

    def MeV2mm(self, MeV):
        C = self.bunch.beta**2*self.bunch.E/self.ring.R*self.ring.gamma_t**2
        return MeV*1e6/C*1e3

    def match(self, sig_tau=None, sig_w=None, k=0):
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

        self.ring.Vg = Vg
        self.__post_init__()
        return Vg
