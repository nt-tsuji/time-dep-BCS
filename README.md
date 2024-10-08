# Time-dependent BCS equation

## Overview

This is a Julia source code for simulating the time-dependent BCS equation, which is formulated in terms of Anderson pseudospins $\vec{\sigma}_k = (\sigma_k^1, \sigma_k^2, \sigma_k^3)$.

The pseudospins are defined by $\sigma_k^\alpha = \frac{1}{2} \langle \psi_k^\dagger \tau^\alpha \psi_k \rangle$ for each momentum $k$ with the Nambu spinor $\psi_k = (c_{k\uparrow}, c_{-k\downarrow}^\dagger)^T$ and Pauli matrices $\tau^\alpha$ ($\alpha=1,2,3$).

The equation for the time evolution is given by

$\partial_t \vec{\sigma}_k = 2 \vec{b}_k \times \vec{\sigma}_k$,

where $\vec{b}_k=(-\Delta', -\Delta'', \varepsilon_k)$ is an effective magnetic field acting on the pseudospins, $\Delta'$ and $\Delta''$ are the real and imaginary parts of the gap function $\Delta$, and $\varepsilon_k$ is the band dispersion.

The gap function $\Delta=\frac{V}{N_k}\sum_k \langle c_{k\uparrow}^\dagger c_{-k\downarrow}^\dagger \rangle$ is determined self-consistently at each time via

$\Delta' + i\Delta'' = \frac{V}{N_k} \sum_k (\sigma_k^1 + i\sigma_k^2)$.

The code solves the time evolution of the pseudospins and gap function using the second-order ([time-dep-BCS-RK2.jl](src/time-dep-BCS-RK2.jl)) and fourth-order ([time-dep-BCS-RK4.jl](src/time-dep-BCS-RK4.jl)) implicit Runge-Kutta methods.

The model is a one-dimensional tight-binding model with the dispersion $\varepsilon_k=-2t_{\rm hop}\cos k-\mu$. The interaction parameter is quenched at $t=0$: $V=V_i \to V_f$. The initial condition is given by an equilibrium state with temperature $T$.

## Parameters

$V_i$: Initial interaction strength

$V_f$: Final interaction strength

$t_{\rm hop}$: Hopping amplitude

$\mu$: Chemical potential

$T$: Initial temperature

$N_k$: Number of $k$ points

$t_{\rm max}$: Maximum time

$N_t$: Number of time steps

## Example

Time evolution of the gap function $\Delta$ for the parameters, $V_i=2, V_f=6, t_{\rm hop}=1, \mu=0, T=0.05, N_k=200, t_{\rm max}=50, N_t=1000$:

![Delta](fig/Delta.png)

## Author

Naoto Tsuji (University of Tokyo)

## References

[1] R. Shimano, N. Tsuji, _Higgs mode in superconductors_, [Annu. Rev. Condens. Matter Phys. 11, 103 (2020)](https://www.annualreviews.org/content/journals/10.1146/annurev-conmatphys-031119-050813) ([arXiv:1906.09401](https://arxiv.org/abs/1906.09401)).

[2] N. Tsuji, I. Danshita, S. Tsuchiya, _Higgs and Nambu-Goldstone modes in condensed matter physics_, [Encyclopedia of Condensed Matter Physics (2nd ed.), Vol. 1, 174 (2024)](https://www.sciencedirect.com/science/article/pii/B9780323908009002560?via%3Dihub) ([arXiv:2310.17148](https://arxiv.org/abs/2310.17148)).
