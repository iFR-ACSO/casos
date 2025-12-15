<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/iFR-ACSO/.github/blob/main/assets/logo-casos-inverted.png">
  <source media="(prefers-color-scheme: light)" srcset="https://github.com/iFR-ACSO/.github/blob/main/assets/logo-casos-trans.png">
  <img alt="CaΣoS: A nonlinear sum-of-squares optimization suite">
</picture>

----

[![Paper](https://img.shields.io/badge/Paper-ACC63710.2025.11107794-00629B?logo=ieee&style=plastic)](https://doi.org/10.23919/ACC63710.2025.11107794)
[![License](https://img.shields.io/badge/License-GPLv3-A42E2B?logo=gnu&style=plastic)](https://github.com/iFR-ACSO/casos?tab=GPL-3.0-1-ov-file#GPL-3.0-1-ov-file)

**CaΣoS** provides tools for [symbolic polynomial expressions](https://github.com/ifr-acso/casos/wiki/polynomial-data-types), [conic optimization](https://github.com/ifr-acso/casos/wiki/conic-optimization), and parametrized (convex and nonconvex) [sum-of-squares optimization](https://github.com/ifr-acso/casos/wiki/sum-of-squares-optimization), making use of the [CasADi](https://web.casadi.org) software for symbolic expressions, automatic differentiation, and numerical optimization. CaΣoS is developed by researchers at the Institute of Flight Mechanics and Controls of the University of Stuttgart and distributed open-source under the GPL-3.0 license.

### Downloads

- All platforms: [version 1.0.0](https://github.com/ifr-acso/casos/releases/tag/v1.0.0) ([zip](https://github.com/ifr-acso/casos/archive/refs/tags/v1.0.0.zip) | [tar.gz](https://github.com/ifr-acso/casos/archive/refs/tags/v1.0.0.tar.gz))

### Quick links

- [Getting started](https://github.com/ifr-acso/casos/wiki#getting-started)
- Available [conic solvers](https://github.com/ifr-acso/casos/wiki#conic-solvers)
- Convex and nonconvex [sum-of-squares optimization](https://github.com/ifr-acso/casos/wiki/sum%E2%80%90of%E2%80%90squares-optimization)
- Supported [vector, matrix, and polynomial cones](https://github.com/ifr-acso/casos/wiki/cones)
- Some [practical tipps](https://github.com/ifr-acso/casos/wiki/practical-sos-guide) to sum-of-squares
- [Transitioning](https://github.com/ifr-acso/casos/wiki/transitioning-from-other-toolboxes) from other toolboxes
- Example [code snippets](https://github.com/ifr-acso/casos/wiki/numerical-examples)

### Publications

If you use CaΣoS, please cite us:

> T. Cunis and J. Olucak, ‘CaΣoS: A nonlinear sum-of-squares optimization suite’, in _2025 American Control Conference_, (Denver, CO), pp. 1659–1666, 2025. doi: [10.23919/ACC63710.2025.11107794](https://doi.org/10.23919/ACC63710.2025.11107794)


<details>

<summary>Bibtex entry</summary>

```bibtex
@inproceedings{Cunis2025acc,
	author = {Cunis, Torbjørn and Olucak, Jan},
	title = {{CaΣoS}: {A} nonlinear sum-of-squares optimization suite},
	booktitle = {2025 American Control Conference},
	address = {Denver, CO},
	year = {2025},
	pages = {1659--1666},
	doi = {10.23919/ACC63710.2025.11107794},
}
```

</details>

#### Applications
We provide the [source code](https://doi.org/10.18419/darus-4499) for our benchmarks and comparison to other toolboxes.

Further applications of CaΣoS include:

1. R. Loureiro and T. Cunis, ‘Nonlinear Observer Synthesis for Stochastic Polynomial Dynamical Systems’, in 2025 American Control Conference, (Denver, CO), 2025, pp. 2509–2514. [DOI](https://doi.org/10.23919/ACC63710.2025.11107965)
2. R. Loureiro and T. Cunis, ‘Estimating Robust Regions of Attraction with Uncertain Equilibrium Points’, in 2025 American Control Conference, (Denver, CO), 2025, pp. 1045–1050. [DOI](https://doi.org/10.23919/ACC63710.2025.11107427)

3. J. Olucak, A. C. B. de Oliveira, and T. Cunis, ‘Safe-by-Design: Approximate Nonlinear Model Predictive Control with Realtime Feasibility’, Preprint. [arXiv](https://arxiv.org/abs/2509.22422), [Source code](https://doi.org/10.18419/darus-5297)
4. F. Geyer, F. Tuttas, W. Fichter, and T. Cunis, ‘Sum-of-Squares Stability Verification on Manifolds with Applications in Spacecraft Attitude Control’, Preprint.



----
