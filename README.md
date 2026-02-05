# PicardLefschetz -- Integration methods for highly-oscillatory integrals based on Picard--Lefschetz theory.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://anneaux.github.io/PicardLefschetz.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://anneaux.github.io/PicardLefschetz.jl/dev/)
[![Build Status](https://github.com/anneaux/PicardLefschetz.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/anneaux/PicardLefschetz.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/anneaux/PicardLefschetz.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/anneaux/PicardLefschetz.jl)


This package implements two integration methods based on Picard-Lefschetz theory, as explained [here](https://p-lpi.github.io/).
That is, we're aiming to calculate integrals of the form 
$` \int_{R^n} ~ \mathbf{p}(\mathbf{x}) ~ \mathrm{e}^{\mathrm{i} \phi(\mathbf{x})} d \mathbf{x} `$, where $\phi(\mathbf{x})$ is analytic almost everywhere, $\mathbf{p}(\mathbf{x})$ is a slowly-varying prefactor and the integration domain is the one ($n=1$) or two-dimensional ($n=2$) real space, corresponding to the 1D and 2D version of the code.

The two implemented methods are the following. 
## Method 1:
Flowing the integration domain into the complex domain, towards a steepest-descent manifold. This manifold is saved as a set of linesegments (1D) / surface elements (2D), on which the integral can be evaluated using a quadrature of choice. 

## Method 2:
Finding points of the exponentiated phase function where the first derivative vanishes, and summing over relevant saddle points' contribution ("saddle-point method"). This only gives a good approximation of the integral as long as the relevant saddle points are well-separated (if you know you know). 

## Disclaimer
This package is developped alongside the application of these methods for solving integrals in attosecond science, as explained in this publication [https://arxiv.org/abs/2510.12545].
Theoretically, the methods should work for generic phase functions $\phi(\mathbf{x})$. We will make an effort to test and further develop this package here to realise that. 
If you have any ideas, want to contribute to, or make use of the package, please do not hesitate to get in touch!!!

## To Do
- (everything)
- write tests
- include 2d code
- make nice example graphs
- write documentation
