# Performance Estimation Problem for Bregman first order methods

This code provides the MATLAB implementation of the computer-aided analysis of Bregman methods, as described in the paper:

> [1] R.A. Dragomir, A. Taylor, A. d'Aspremont, J. Bolte, "Optimal Complexity and Certification of Bregman First-OrderMethods" arXiv, 2019.

## Getting Started

This code requires [YALMIP](https://yalmip.github.io/) along with a suitable SDP solver (e.g., Sedumi, SDPT3, Mosek).
The files in the [PESTO](PESTO/) folder require the Performance Estimation Toolbox (PESTO(https://github.com/AdrienTaylor/Performance-Estimation-Toolbox))


## List of files

- [`PEP_NoLips.m`](PEP_NoLips.m) solves the Performance Estimation Problem (PEP) for computing the worst-case convergence rate of the NoLips/Bregman Gradient algorithm as described in [1, Section 4].
- [`PEP_NoLips_lowrank_worstcase.m`](PEP_NoLips_lowrank_worstcase.m) searches for particular solutions of the PEP, which are low-rank and where the gradients $\nabla f(x_i)$ are orthogonal, as described in [1, Section 4.5.4].
- [`plot_discrete_functions.m`](plot_discrete_functions.m) defines the function used to display the worst-case instances in dimension 1 or 2.

## Authors

- Radu-Alexandru Dragomir
- [**Adrien Taylor**](http://www.di.ens.fr/~ataylor/)
- [**Alexandre d'Aspremont**](https://www.di.ens.fr/~aspremon/)
- [**Jérôme Bolte**](https://www.tse-fr.eu/fr/people/jerome-bolte)
