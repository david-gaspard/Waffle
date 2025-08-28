# Waffle

[![C++](https://img.shields.io/badge/C++-%2300599C.svg?logo=c%2B%2B&logoColor=white)](https://cplusplus.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-red.svg)](https://opensource.org/licenses/MIT)

* [PRESENTATION](#presentation)
    - [Wave equation](#wave-equation)
    - [Transmission eigenvalue distribution](#transmission-eigenvalue-distribution)
    - [Profile of transmission eigenchannels](#profile-of-transmission-eigenchannels)
    - [Acknowledgements](#acknowledgements)
* [USAGE AND OPTIONS](#usage-and-options)
* [INTERNAL IMPLEMENTATION](#internal-implementation)

## PRESENTATION

Waffle is a C++ 2017 program to solve the stationary wave equation in a two-dimensional disordered medium of arbitrary shape.
The name is an contraction of *"Wave Field From Finite Elements"*.
This program primarily focuses on the computation of the transmission and reflection matrices between two or more leads, the distribution of transmission eigenvalues, and the profile of transmission eigenchannels (aka transmission eigenstates).

### Stationary wave equation

The stationary wave equation reads:

<p>$$ ( \nabla^2 + k^2 - U(\mathbf{r}) ) \psi(\mathbf{r}) = 0 $$</p>

where $\psi(\mathbf{r})$ is the wavefunction, $k$ is the wavenumber, and $U(\mathbf{r})$ is the random potential whose fluctuations are related to the scattering [mean free path](https://en.wikipedia.org/wiki/Mean_free_path) $\ell_{\mathrm{s}}$.
However, it should be noted that the program actually solves a discretized variant of the wave equation defined on a square lattice of step $h$:

<p>$$ ( \mathsf{L} + (kh)^2 - \mathsf{u} ) \mathbf{\psi} = 0 $$</p>

In this equation, the matrix $\mathsf{L}$ represents the five-point Laplacian, defined by elements $-4$ along the diagonal and $+1$ with the four nearest neighbors (4-neighborhood), $(kh)^2$ is implicitly multiplied by the identity matrix, $\mathsf{u}$ is a diagonal matrix containing the values of $U(\mathbf{r})h^2$ on all grid points, and $\mathbf{\psi}$ is a column vector containing the values of the wavefunction on the grid points.
It is important to note that the dispersion relation of this discretized wave equation is *not* isotropic. It reads

<p>$$ 4\sin^2(p_xh/2) + 4\sin^2(p_yh/2) = (kh)^2 $$</p>

and reduces to the standard isotropic dispersion $p_x^2+p_y^2=k^2$ only in the limit $h\rightarrow 0$ (for $k$ fixed).
In order to avoid forbidden directions of propagation, the wavenumber must be limited to $kh\leq 2$.
A reasonable compromise between mitigating anisotropy and reducing the number of grid points is $kh=1$, which corresponds to about 6 points per wavelength.
With this choice, the dispersion curve deviates by less than $2.5\%$ from a perfect circle.

<span style="color:red">TODO: Write a few words on open boundary conditions, and stress that they are exact...</span>

### Acknowledgements

This program was developed by David Gaspard ([Institut Langevin](https://ror.org/00kr24y60), [ESPCI Paris](https://ror.org/03zx86w41), [PSL University](https://ror.org/013cjyk83), [CNRS](https://ror.org/02feahw73)) mainly in August 2025.
This research has been supported by the [ANR](https://ror.org/00rbzpz17) project MARS_light under reference [ANR-19-CE30-0026](https://anr.fr/Project-ANR-19-CE30-0026), by the program "Investissements d'Avenir" launched by the French Government.
It also received support from a grant of the [Simons Foundation](https://ror.org/01cmst727) (No. 1027116).

## USAGE AND OPTIONS

<span style="color:red">TODO: Write a few words on the commands...................................</span>


