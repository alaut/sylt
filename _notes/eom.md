---
title: Equations of Motion
date: December 2, 2021
layout: page
permalink: eom
---

# Equations of Motion

Consider the relative time and energy coordinates

$$\tau = t-t_s\qquad w=W-W_s,$$

where $t$ represents time and $W$ represents particle energy. Subscripted variables indicate the synchronous particle within a synchrotron.

## Kick


The following 

$$\Delta W = qV_g\sin\varphi,$$

describes the discrete energy gain of a particle passing through an RF gap of voltage $V_g$ with RF phase $\varphi=h\omega_st$, defined by the harmonic $h$ and angular revolution frequency $\omega_s$.

Averaged across one synchrotron period $T_s$, the instantaneous relative energy gain can be given by

$$\dot{w}=\frac{\Delta w}{T_s} = \frac{qV_g}{T_s}(\sin\varphi-\sin\varphi_s).$$

## Drift

Using logarithmic derivatives, the relationship between revolution period $T$, path circumference $C$ and relative velocity $\beta$ can be described by

$$\frac{dT}{T} = \frac{dC}{C} - \frac{d\beta}{\beta}.$$

__Momentum compaction__ is defined as follows

$$\alpha \equiv \frac{dC/C}{dp/p} \equiv \frac{dR/R}{dp/p}.$$

Using the relationship

$$\frac{d\beta}{\beta} = \frac{1}{\gamma^2}\frac{dp}{p}$$

if we subsitute for $dC/C$ and $d\beta/\beta$ and get

$$\frac{dT}{T} = \left(\alpha - \frac{1}{\gamma^2}\right)\frac{dp}{p} = \eta\frac{dp}{p},$$

where the __slippage factor__ is defined by

$$\eta =\alpha- \frac{1}{\gamma^2}.$$

From

$$\frac{dp}{p}=\frac{1}{\beta^2}\frac{dE}{E}$$

we have

$$\frac{dT}{T} = \frac{\eta}{\beta^2}\frac{dE}{E}.$$

Our equations of motion can be succinctly written as

$$\dot{\tau} = \frac{\eta}{\beta_s^2E_s}w \qquad \dot{w} = \frac{qV_g}{T_s}(\sin\varphi-\sin\varphi_s).$$
