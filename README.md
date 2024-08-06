# Overview

This is a Julia source code for simulating the time-dependent BCS equation, which is formulated in terms of Anderson pseudospins \sigma_k = (\sigma_k^x, \sigma_k^y, \sigma_k^z).

The equation for the time evolution is given by

\partial_t \sigma_k = 2 h_k \times \sigma_k,

where h_k = (-\Delta_re, -\Delta_im, \epsilon_k), \Delta_re and \Delta_im are the real and imaginary parts of the gap function, and \epsilon_k is the band dispersion.

The gap function is determined self-consistently via

\Delta_re + i*\Delta_im = V/N*\sum_k (\sigma_k^x + i*\sigma_k^y).

The code solves the time evolution of the pseudospins and gap function using the second-order implicit Runge-Kutta method.
