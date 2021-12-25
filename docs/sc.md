---
title: Space Charge
permalink: sc
author: Alex Laut
date: December 2, 2021
layout: page
---

## Wakefields

Charged particles traversing in a slowly varying long beam will induce a wakefiled voltage given by

$$V_{W} = -\frac{i}{\omega_s}\frac{Z}{n}\frac{d\lambda}{d\tau},$$

where the longitudinal charge density is given by $\lambda(\tau)$ and the space charge impedance $Z/n$ is given by

$$\frac{Z}{n}=-j\frac{Z_o}{\beta_s\gamma_s^2}{g}$$

where $g$ is the space charge geometry factor, $Z_0$ is the impedance of free space, and $n=f/f_s$.

## Induced Voltage

The induced voltage gradient due to space charge can therfore be written as

$$V'_{SC}(\tau) = \frac{g}{\omega_s}\frac{Z_0}{\beta_s\gamma_s^2}\frac{d^2\lambda}{d\tau^2}$$

## Synchrotron Frequency Shift

The synchrotron frequency is generalized by

$$\Omega^2=-\frac{\eta}{\beta_s^2E_s}\frac{q}{T_s}V'(\tau)$$

where

$$V'(\tau) = V'_{RF}(\tau)+V'_{SC}(\tau)$$

and

$$V'_{RF}=V_g\frac{d}{d\tau}g(\phi)\approx h\omega_s V_g\cos\varphi$$

The frequency tune shift is given by

$$\mu=\Omega/\Omega_s=\sqrt{1+\frac{V'_{SC}}{V'_{RF}}}$$

For a parabolic distribution given by

$$\lambda(\tau) = \frac{3Q}{2L_\tau}(1-4\frac{\tau^2}{L_\tau^2})$$

and so

$$\lambda''=-\frac{12Q}{L_\tau^3}$$

therefore the gradient is given by

$$V'_{SC}(\tau) = -\frac{12Q}{L_\tau^3}\frac{g}{\omega_s}\frac{Z_0}{\beta_s\gamma_s^2}$$

and so the frequency spread is shifted by a constant given vy $V'_{SC}(\tau)$

## Effective Voltage

Given a distribution of particles with varying synchrotron frequency, the effective synchrotron frequency is given by

$$<\mu> = \oint\lambda(\phi)\mu(\phi)d\phi$$

For a parabolic distribution we have

$$<\mu> = \int_\infty \frac{3}{2L}(1-4\frac{\tau^2}{L^2})(1-\frac{(h\omega_s\tau)^2}{16})d\tau$$

$$<\mu> = 1-\frac{\sigma_{\hat{\phi}}^2}{16}$$

from

$$\mu^2=\Omega^2/\Omega_s^2\propto V'(\tau)\propto V_g$$

and that

$$\mu = \Omega/\Omega_s$$

The implied gap voltage $V_g$ given by the interpreted effective voltage $<V>$ per the synchrotron frequency and the shape function $<\mu>$, it is given by

$$V_g = \frac{<V>}{<\mu>^2}$$

__This needs to include space charge tho__