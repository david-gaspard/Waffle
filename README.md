# Waffle

[![C++](https://img.shields.io/badge/C++-%2300599C.svg?logo=c%2B%2B&logoColor=white)](https://cplusplus.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-red.svg)](https://opensource.org/licenses/MIT)

* [PRESENTATION](#presentation)
    - [Wave equation and transmission eigenstates](#wave-equation-and-transmission-eigenstates)
    - [Discretization](#discretization)
    - [Acknowledgements](#acknowledgements)
* [INSTALLATION AND USAGE](#installation-and-usage)
    - [Dependencies](#dependencies)
    - [Usage](#usage)
* [REFERENCES](#references)

## PRESENTATION

Waffle is a C++ 2017 program to solve the [stationary wave equation](https://en.wikipedia.org/wiki/Helmholtz_equation) in a two-dimensional disordered medium of arbitrary shape.
The name is a contraction of *"Wave Field From Finite Elements"*.
This program primarily focuses on the computation of the [transmission and reflection matrices](https://en.wikipedia.org/wiki/S-matrix) between two or more ducts, the distribution of transmission eigenvalues, and the profile of transmission eigenstates (aka transmission eigenchannels).

### Wave equation and transmission eigenstates

The starting point of all wave simulations is the computation of the retarded modal [Green's function](https://en.wikipedia.org/wiki/Green's_function) $G^+_{n}(\mathbf{r})$ defined by

<p>$$ \left( \nabla^2_{\mathbf{r}} + k^2 + \mathrm{i}\varepsilon - U(\mathbf{r}) \right) G^+_{n}(\mathbf{r}) = \chi_{{\rm in},n}(\mathbf{y}) \delta(x-x_{\rm in}) , $$</p>

where $\mathbf{r}=(x,\mathbf{y})$ is the position, $x$ the longitudinal coordinate, $\mathbf{y}$ the transverse coordinates, $k$ the [wavenumber](https://en.wikipedia.org/wiki/Wavenumber), $\varepsilon=k/\ell_{\rm a}$ the absorption parameter, $U(\mathbf{r})$ the arbitrary potential, and $x_{\rm in}$ represents some longitudinal coordinate in the input duct (typically at the edge of the domain).
The function $\chi_{{\rm in},n}(\mathbf{y})$ is the $n$-th [transverse eigenmode](https://en.wikipedia.org/wiki/Waveguide#Propagation_modes_and_cutoff_frequencies) of the input duct normalized by $\int_{\mathcal{S}_{\rm in}} \mathrm{d}{\mathbf{y}}\, \chi^*_{{\rm in},m}(\mathbf{y}) \chi_{{\rm in},n}(\mathbf{y}) = \delta_{mn}$.
At a reflecting wall, we use [Dirichlet boundary condition](https://en.wikipedia.org/wiki/Dirichlet_boundary_condition) $\psi(\mathbf{r})=0$ which conserves wave's energy.
Moreover, at the end of each open duct, we use a free-escape boundary condition ensuring that the wave leaves the system without reflection.
This boundary condition formally imposes the cancellation of all incoming waves, i.e.,

<p>$$ \left( \mathbf{n}\cdot\nabla_{\mathbf{r}} - \mathrm{i}\sqrt{k^2 + \nabla^2_{\mathbf{y}}} \right) G^+_n(\mathbf{r}) = 0 , $$</p>

where $\mathbf{n}$ is the outward [normal](https://en.wikipedia.org/wiki/Normal_(geometry)) to the duct's open edge.
Note that the square root makes the operator [nonlocal](https://en.wikipedia.org/wiki/Nonlocal_operator) across all points on the duct's open edge.
From the Green's function, we can derive the transmission and reflection matrices according to the Fisher and Lee relations [^1],

<p>$$ t_{mn} = 2\mathrm{i}\sqrt{k_{{\rm out},m} k_{{\rm in},n}} \int_{\mathcal{S}_{\rm out}} \mathrm{d}{\mathbf{y}}\, \chi_{{\rm out},m}^*(\mathbf{y}) G^+_{n}(x_{\rm out},\mathbf{y}) , $$</p>

<p>$$ r_{mn} = -\delta_{mn} + 2\mathrm{i}\sqrt{k_{{\rm in},m} k_{{\rm in},n}} \int_{\mathcal{S}_{\rm in}} \mathrm{d}{\mathbf{y}}\, \chi_{{\rm in},m}^*(\mathbf{y}) G^+_{n}(x_{\rm in},\mathbf{y}) , $$</p>

where $x_{\rm out}$ is some longitudinal coordinate in the output duct.
Transmission eigenvalues and eigenstates are defined from the [singular value decomposition](https://en.wikipedia.org/wiki/Singular_value_decomposition) (SVD) of the transmission matrix:

<p>$$ \mathsf{t} = \sum_{e=1}^{N_{\rm e}} \sqrt{T_e} \mathbf{u}_e \mathbf{v}_e^\dagger , $$</p>

where $N_{\rm e}=\min(N_{\rm in},N_{\rm out})$ is the number of transmission eigenvalues, $N_{\rm in}$ and $N_{\rm out}$ being the number of propagating modes in the input and output ducts, respectively, and $\dagger$ denotes the conjugate transpose.
The quantities $T_1,T_2,\ldots,T_{N_{\rm e}}$ are the transmission eigenvalues bounded by $0\leq T_e\leq 1~\forall e\in\{1,\ldots,N_{\rm e}\}$ due to flux conservation.
The column vectors $\mathbf{u}_e$ (of length $N_{\rm out}$) and $\mathbf{v}_e$ (of length $N_{\rm in}$) are normalized according to $\mathbf{u}_{e}^\dagger\mathbf{u}_{e'}=\delta_{ee'}$ and $\mathbf{v}^\dagger_{e}\mathbf{v}_{e'}=\delta_{ee'}$.

The transmission eigenvalue distribution $\rho(T)$ is then formally defined by

<p>$$ \rho(T) = \frac{1}{N_{\rm e}} \sum_{e=1}^{N_{\rm e}} \langle\delta(T-T_e)\rangle $$</p>

where $\langle\cdot\rangle$ denotes the average over the disorder.
In this program, all the wavefunctions are normalized such that their average incident intensity is equal to one ($\int_{\mathcal{S}_{\rm in}} \frac{\mathrm{d}{\mathbf{y}}}{S_{\rm in}} |\psi_{\rm in}(\mathbf{r})|^2 = 1$).
As a consequence, their general expression reads

<p>$$ \psi(\mathbf{r}) = \frac{\sum_{n=1}^{N_{\rm in}} a_n G^+_n(\mathbf{r}_{\rm in})}{\sqrt{\frac{1}{S_{\rm in}} \sum_{n=1}^{N_{\rm in}} \frac{|a_n|^2}{4k_{{\rm in},n}^2}}}  $$</p>

where $a_n$ are modal amplitudes determining the wavefunction uniquely, and $S_{\rm in}$ the cross-section surface of the input duct.
In particular, for transmission eigenstates, the amplitudes are given by

<p>$$ a_n = \sqrt{k_{{\rm in},n}} v_{n,e} $$</p>

where $v_{n,e}$ is the $n$-th component of $\mathbf{v}_e$ given by the SVD above.
Plane wave illumination at normal incidence can be achieved by choosing instead

<p>$$ a_n = \frac{1 - (-1)^n}{n} 2\mathrm{i} k_{{\rm in},n} $$</p>

The intensity profile of transmission eigenstate is therefore given by

<p>$$ I_T(\mathbf{r}) = \frac{1}{N_{\rm e}\rho(T)} \left\langle \sum_{e=1}^{N_{\rm e}} |\psi_{T_e}(\mathbf{r})|^2 \delta(T-T_e) \right\rangle $$</p>

where $\psi_{T_e}(\mathbf{r})$ is the transmission eigenstate defined above.

### Discretization

In the program, the equations above are discretized over a two-dimensional square lattice of step $h$.
In particular, a [five-point stencil](https://en.wikipedia.org/wiki/Five-point_stencil) is used to approximate the Laplacian in the wave equation, which affects the free-space [dispersion relation](https://en.wikipedia.org/wiki/Dispersion_relation): 

<p>$$ 4\sin^2(\frac{1}{2}k_xh) + 4\sin^2(\frac{1}{2}k_yh) = (kh)^2 $$</p>

However, by choosing $kh\leq 1$ in our simulations, which corresponds to 6 points per wavelength, this relation deviates by less than $2.5\%$ from the ideal circle, $k_x^2+k_y^2=k^2$, which can be considered acceptable.

### Acknowledgements

This program has been written by David Gaspard ([Institut Langevin](https://ror.org/00kr24y60), [ESPCI Paris](https://ror.org/03zx86w41), [PSL University](https://ror.org/013cjyk83), [CNRS](https://ror.org/02feahw73)) mainly in August 2025.
This research has been supported by the [ANR](https://ror.org/00rbzpz17) project MARS_light under reference [ANR-19-CE30-0026](https://anr.fr/Project-ANR-19-CE30-0026), by the program "Investissements d'Avenir" launched by the French Government.
It also received support from a grant of the [Simons Foundation](https://ror.org/01cmst727) (No. 1027116).

## INSTALLATION AND USAGE

The source files can be downloaded using the [`git clone`](https://git-scm.com/docs/git-clone) command:
```shell
git clone https://github.com/david-gaspard/Waffle.git
```
To compile the program, call the [`make`](https://en.wikipedia.org/wiki/Make_(software)) utility in the root directory:
```shell
make all
```

### Dependencies

The program requires a C++ compiler, such as from the [GNU Compiler Collection](https://en.wikipedia.org/wiki/GNU_Compiler_Collection), and the libraries [OpenBLAS](https://en.wikipedia.org/wiki/OpenBLAS), [UMFPACK](https://en.wikipedia.org/wiki/UMFPACK), [MUMPS](https://en.wikipedia.org/wiki/MUMPS_(software)), and [libpng](https://en.wikipedia.org/wiki/Libpng).
It also calls [Python 3](https://en.wikipedia.org/wiki/Python_(programming_language)) with the [NumPy](https://en.wikipedia.org/wiki/NumPy), [Matplotlib](https://en.wikipedia.org/wiki/Matplotlib), and [csv](https://docs.python.org/3/library/csv.html) modules, but also the command `pdflatex` from [TeX Live](https://en.wikipedia.org/wiki/TeX_Live) with the [PGFPlots](https://ctan.org/pkg/pgfplots) package to process graphical outputs automatically.

### Usage

The recommended way to launch the program is through the `run` script:
```shell
./run
```
This script defines the environment variables and calls the executable.
The executable is built from the `src/Main.cpp` file which contains the computation instructions.
The first step is to define the mesh over which the wave equation will be solved.
This can be done by instantiating a `SquareMesh` object as follows
```cpp
SquareMesh mesh("model/image.png");
```
where `model/image.png` can be replaced by a file path to any image in the `model` folder or elsewhere on your machine.
This model image completely defines the mesh structure and uses the following color code:

* White (`#FFF`): The pixel is void and should be skipped from the mesh.
* Black (`#000`): The pixel is added to the mesh. The default boudary condition between black and white pixels are Dirichlet boundary conditions $\psi(\mathbf{r})=0$.
* Red (`#F00`): Black pixels surrounding a red pixel are flagged as belonging to an input duct for the construction of the transmission matrix.
  Red pixels are not added to the mesh.
* Blue (`#00F`): Black pixels surrounding a blue pixel are flagged as belonging to an output duct for the construction of the transmission matrix. 
  Blue pixels are not added to the mesh.
* Green (`#0F0`): Black pixels surrouding a green pixels are flagged as belonging to a duct which is neither input nor output and thus ignored from the transmission matrix.
  Green pixels are not added to the mesh.

Note that red (`#F00`), blue (`#00F`), and green (`#0F0`) pixels assume free-esacpe boundary conditions presented before.
It is worth noting that `SquareMesh` objects can also be generated procedurally using the methods found in the file `src/SquareMesh.hpp`.
This can be helpful if the boundaries are complicated or change randomly from one simulation to the other.

A `WaveSystem` object (containing the potential, the discretization of the wave equation, and the Green's function) is then instantiated by the command:
```cpp
WaveSystem sys(sysname, mesh, kh, density, holscat, holabso);
```
where

* `sysname` is a string containing the full name of the system, and which will be used for file output.
* `mesh` is the `SquareMesh` object initialized from a model image (see above).
* `kh` is the product of the wavenumber and the mesh step.
* `density` is the density of pixels which are given a random potential disorder between 0 and 1.
* `holscat` is the ratio of the mesh step over the scattering [mean free path](https://en.wikipedia.org/wiki/Mean_free_path).
* `holabso` is the ratio of the mesh step over the ballistic absorption length.

The `WaveSystem` object is the main computational object of the program.
Computations can be performed by calling the methods of `WaveSystem` (see also `src/WaveSystem.hpp`).
For instance, in order to setup a new random realization of the disorder, one can call
```cpp
uint64_t seed = 1;
sys.setDisorder(seed);
```
where `seed` is a long integer uniquely representing the realization of the disorder.
The transmission matrix `tmat` associated with the propagation from input ducts (red `#F00` pixels) to output ducts (blue `#00F` pixels) can then be computed using
```cpp
ComplexMatrix tmat(sys.getNOutputProp(), sys.getNInputProp());
sys.transmissionMatrix(tmat);
```
Frequent tasks have their own methods. For instance, the computation of transmission eigenvalues can be directly achieved with
```cpp
RealMatrix tval(std::min(sys.getNOutputProp(), sys.getNInputProp()), 1);
sys.addTSpectrum(tval);
```
without having to handle the transmission matrix.
On output, the transmission eigenvalues are stored in the column matrix `tval`.
Another useful procedure is the following one, which simultaneously compute the transmission eigenvalues and the intensity profile of transmission eigenstates in specific intervals of transmission eigenvalues:
```cpp
const int nprofile = 3;
const RealMatrix trange(nprofile, 2);
trange(0, 0) = 0.998;  trange(0, 1) = 0.002;  // Interval is T=[0.996, 1].
trange(1, 0) = 0.500;  trange(1, 1) = 0.010;  // Interval is T=[0.49, 0.51].
trange(2, 0) = 0.100;  trange(2, 1) = 0.005;  // Interval is T=[0.095, 0.105].

RealMatrix tprofile(sys.getNPoint(), nprofile):
RealMatrix nsample(nprofile, 1);
RealMatrix tval(std::min(sys.getNOutputProp(), sys.getNInputProp()), 1);

sys.addITransmission(trange, tprofile, nsample, tval);
```
where

* `trange` is an input matrix with two columns, the first column defined the center of intervals over which the transmission igenstates are desired, and the second column are the half width of intervals. The number of rows of this matrix determines the number of transmission eigenstates profiles.
* `tprofile` is a matrix containing, on output, the intensity profile of transmission eigenstates, one per column. 
If multiple eigenvalues are found in the same interval, then they are summed up (not averaged).
The method adds the profiles directly to `tprofile`, so this matrix must be initialized to zero.
* `nsample` is a column matrix containing, on output, the number of transmission eigenvalues (and thus eigenstates) found in the intervals prescribed by `trange`.
This matrix must have as many rows as `trange`, and one column.
The method increments the number of found transmission eigenvalues in `nsample`, so this matrix must be initialized to zero.
* `tval` is a column matrix containing, on output, all the transmission eigenvalues (should they belong to the intervals of `trange` or not).

The full list of methods can be found in the header file `src/WaveSystem.hpp` and the documentation in the corresponding `*.cpp` file.

## REFERENCES

[^1]: D. S. Fisher and P. A. Lee, _Relation between conductivity and transmission matrix_, [Phys. Rev. B **23**, 6851-6854 (1981)](https://doi.org/10.1103/PhysRevB.23.6851).
