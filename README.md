# Overview

This is a Julia source code for simulating the time-dependent BCS equation, which is formulated in terms of Anderson pseudospins $\vec{\sigma}_k = (\sigma_k^x, \sigma_k^y, \sigma_k^z)$.

The equation for the time evolution is given by

$\partial_t \vec{\sigma}_k = 2 \vec{h}_k \times \vec{\sigma}_k$,

where $\vec{h}_k=(-\Delta', -\Delta'', \epsilon_k)$, $\Delta'$ and $\Delta''$ are the real and imaginary parts of the gap function, and $\epsilon_k$ is the band dispersion.

The gap function is determined self-consistently via

$\Delta' + i\Delta'' = \frac{V}{N} \sum_k (\sigma_k^x + i\sigma_k^y)$.

The code solves the time evolution of the pseudospins and gap function using the second-order implicit Runge-Kutta method.
