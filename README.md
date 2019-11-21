# Performance Estimation Problem for Bregman first order methods

This code provides the MATLAB implementation of the computer-aided analysis of Bregman methods, as described in the paper:

> [1] R.A. Dragomir, A. Taylor, A. d'Aspremont, J. Bolte, "Optimal Complexity and Certification of Bregman First-OrderMethods" arXiv preprint arXiv:1911.08510, 2019.

## Getting Started

This code requires [YALMIP](https://yalmip.github.io/) along with a suitable SDP solver (e.g., Sedumi, SDPT3, Mosek).
The files in the [PESTO](PESTO/) folder require the Performance Estimation Toolbox ([PESTO](https://github.com/AdrienTaylor/Performance-Estimation-Toolbox))

For getting started, the basic performance estimation problem (PEP) for computing the worst-case performance of the NoLips/Bregman Gradient algorithm is provided in the script [`PEP_NoLips.m`](PEP_NoLips.m).


## List of files

- [`PEP_NoLips.m`](PEP_NoLips.m) solves the Performance Estimation Problem (PEP) for computing the worst-case convergence rate of the NoLips/Bregman Gradient algorithm as described in [1, Section 4].
- [`PEP_NoLips_lowrank_worstcase.m`](PEP_NoLips_lowrank_worstcase.m) implements the procedure used to discover the general worst-case functions by searching particular solutions of the PEP, as described in [1, Section 4.5.4]. 
- [`plot_discrete_functions.m`](plot_discrete_functions.m) defines the function used to display the worst-case instances in dimension 1 or 2.
- [`PESTO/`](PESTO/) provides implementation of PEPs for several Bregman methods using the [PESTO](https://github.com/AdrienTaylor/Performance-Estimation-Toolbox) toolbox:
    - [`PESTO/NoLips.m`](PESTO/NoLips.m) for the basic NoLips algorithm
    - [`PESTO/NoLips_takeII.m`](PESTO/NoLips_takeII.m) for analyzing a different convergence criterion for NoLips as in [1, Section 4.5.2]
    - [`PESTO/BregmanProximalPoint.m`](PESTO/BregmanProximalPoint.m) for the Bregman proximal point method
    - [`PESTO/IGA.m`](PESTO/IGA.m) for the Improved Interior Gradient Algorithm (IGA)
    - [`PESTO/IGA.m`](PESTO/IGA.m) for IGA with no constraints [1, Section 4.5.3]
- [`symbolic_verifications/`](symbolic_verifications/) provides the symbolic verifications of the convergence rate proofs in [1]:
    - [`symbolic_verifications/NoLips_FunctionValues.nb`](symbolic_verifications/NoLips_FunctionValues.nb) for Theorem 1
    - [`symbolic_verifications/NoLips_Gradients.nb`](symbolic_verifications/NoLips_Gradients.nb) for Proposition 4 		

## Authors

- Radu-Alexandru Dragomir
- [**Adrien Taylor**](http://www.di.ens.fr/~ataylor/)
- [**Alexandre d'Aspremont**](https://www.di.ens.fr/~aspremon/)
- [**Jérôme Bolte**](https://www.tse-fr.eu/fr/people/jerome-bolte)
