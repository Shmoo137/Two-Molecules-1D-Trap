# Two Molecules in a 1D Trap

This repository accompanies two research papers exploring the physics of two interacting molecules in a 1D trap:
- A. Dawid, M. Lewenstein, M. Tomza. 2018. Two ultracold interacting molecules in a one-dimensional harmonic trap. [Phys. Rev. A 97, 063618](https://doi.org/10.1103/PhysRevA.97.063618)
- A. Dawid, M. Tomza. 2020. Magnetic properties and quench dynamics of two interacting ultracold molecules.

The contained code allows to reproduce results and figures from both papers, namely:
- spectrum as a function of any Hamiltonian parameter,
- magnetization and expected J<sup>2</sup> of the system,
- the time evolution of the system after the quench of any Hamiltonian parameter.

The Hamiltonian parameters governing the molecular system are:
- the isotropic and anisotropic parts of the intermolecular interaction,
- the external electric field,
- the external magnetic field,
- the spin-rotation coupling.
Other static physical parameters can be modified in `constants.h`.

The arguments of the main() function control both the parameters of the system as the output of the program.
They are described in `main.cpp`.

The structure of the program:
- `matrix` contains the Matrix class,
- `basis` includes functions to build the basis set depending on the input parameters and chosen quantum numbers,
- `build_hamiltonian` builds the system's Hamiltonian in a created basis set as a Matrix object, depending on the input parameters, using functions from `ham_functions`,
- `constants` contains the physical constants, the rotational constants of the molecules, their electric dipole moments, and maximum $n$ and $j$ used to build a basis set,
- `eigenproblem` includes the Intel Math Kernel library functions to solve a symmetric real matrix eigenproblem,
- `looped_eigenproblem` allows to calculate a desired property as a function of any Hamiltonian parameter,
- `parameters` contains extern variables, which are mostly declared via reading the main() function arguments,
- `quench dynamics` has functions to calculate the time evolution of the chosen property of the system after the quench of any Hamiltonian parameter, using functions from `quench_functions`,
- `utility` contains functions that save results, read main() arguments, pair and sort data, etc.

Code was written by Anna Dawid (University of Warsaw & ICFO) with help of Michał Tomza (University of Warsaw)