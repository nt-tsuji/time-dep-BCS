# Time-dependent BCS equation

## Overview

This is a Julia source code for simulating the time-dependent BCS equation, which is formulated in terms of Anderson pseudospins $\vec{\sigma}_k = (\sigma_k^1, \sigma_k^2, \sigma_k^3)$.

The pseudospins are defined by $\sigma_k^\alpha = \frac{1}{2} \langle \psi_k^\dagger \tau^\alpha \psi_k \rangle$ with the Nambu spinor $\psi_k = (c$<sub>$k\uparrow$</sub>$, c^\dagger$<sub>$-k\downarrow$</sub>$)^T$ and Pauli matrices $\tau^\alpha$ ($\alpha=1,2,3$).

The equation for the time evolution is given by

$\partial_t \vec{\sigma}_k = 2 \vec{b}_k \times \vec{\sigma}_k$,

where $\vec{b}_k=(-\Delta', -\Delta'', \varepsilon_k)$ is an effective magnetic field acting on the pseudospins, $\Delta'$ and $\Delta''$ are the real and imaginary parts of the gap function, and $\varepsilon_k$ is the band dispersion.

The gap function is determined self-consistently via

$\Delta' + i\Delta'' = \frac{V}{N_k} \sum_k (\sigma_k^x + i\sigma_k^y)$.

The code solves the time evolution of the pseudospins and gap function using the second-order implicit Runge-Kutta method.

The model is a one-dimensional tight-binding model with the dispersion $\varepsilon_k=-2t_{\rm hop}\cos k$. The interaction parameter is quenched at $t=0$: $V=V_i \to V_f$. The initial condition is given by an equilibrium state with temperature $T$.

## Parameters

$V_i$: Initial interaction strength

$V_f$: Final interaction strength

$t_{\rm hop}$: Hopping amplitude

$T$: Initial temperature

$N_k$: Number of $k$ points

$t_{\rm max}$: Maximum time

$N_t$: Number of time steps

## References

[1] R. Shimano, N. Tsuji, _Higgs mode in superconductors_, [Annu. Rev. Condens. Matter Phys. 11, 103 (2020)](https://www.annualreviews.org/content/journals/10.1146/annurev-conmatphys-031119-050813) ([arXiv:1906.09401](https://arxiv.org/abs/1906.09401)).

[2] N. Tsuji, I. Danshita, S. Tsuchiya, _Higgs and Nambu-Goldstone modes in condensed matter physics_, [Encyclopedia of Condensed Matter Physics (2nd ed.), Vol. 1, 174 (2024)](https://www.sciencedirect.com/science/article/pii/B9780323908009002560?via%3Dihub) ([arXiv:2310.17148](https://arxiv.org/abs/2310.17148)).
