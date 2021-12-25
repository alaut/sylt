---
title: Transverse Motion
permalink: tm
author: Alex Laut
date: December 2, 2021
layout: page
---

## Betatronic Motion

The transverse position of a particle $u$ along longitudinal coordinate $s$ is described by

$$u(s) = \sqrt{\beta(s)\epsilon}\cos\phi+D\delta$$

where $\beta$ and $D$ are the beta and dispersion functions respectively. The transverse phase advance is given by

$$\phi = \mu(s) + \mu_0$$

where $\mu_0$ is an integration constant and

$$\mu(s) = \int_0^s\frac{ds}{\beta(s)}.$$

## Effective Geometry Factor

For a particle within a uniform transverse bunch distribution in a long round beam, the longitudinal space charge fields can be described by the following geometry factor

$$g(r) = \begin{cases}\frac{1}{2}+\ln\frac{b}{a}-\frac{1}2{}\frac{r^2}{a^2} & r < a \\ \ln \frac{b}{r} & r > a
\end{cases}.$$

If we approximate that

$$\beta(s)\approx \beta \qquad D(s) \approx D$$

we can determine that

$$\overline{r^2} = \frac{\beta}{2}(\epsilon_x +\epsilon_y)+D^2\delta^2$$

and that

$$\bar{a} = 2\sqrt[4]{(\beta\sigma_\epsilon+D^2\sigma^2_\delta)(\beta\sigma_\epsilon)}$$


Accordingly

$$\bar{g}(X,Y)=\frac{1}{2}+\ln\frac{b}{\bar{a}(Y)}-\frac{1}{2}\frac{\overline{r^2}(X)}{\bar{a}^2(Y)}$$

where 

$$X \in (\delta, \epsilon_x, \epsilon_y) \qquad Y \in (\sigma_\epsilon, \sigma_\delta)$$

are evaluated at each turn.

<!-- ## Synchrotron Frequency Blur -->
