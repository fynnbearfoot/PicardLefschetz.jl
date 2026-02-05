# PicardLefschetz

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://anneaux.github.io/PicardLefschetz.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://anneaux.github.io/PicardLefschetz.jl/dev/)
[![Build Status](https://github.com/anneaux/PicardLefschetz.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/anneaux/PicardLefschetz.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/anneaux/PicardLefschetz.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/anneaux/PicardLefschetz.jl)


This package implements two integration methods based on Picard-Lefschetz theory, as explained [here](https://p-lpi.github.io/).
That is, we're aiming to calculate integrals of the form $\int_{R^n} \vec{p}(\vec{x}) \mathrm{e}^{\mathrm{i} \phi(\vec{x})} \mathrm{d} \vec{x} $, 
for one ($n=1$) or two-dimensional ($n=2$) real-valued integration domains (1D and 2D version of the code).


Method 1:
Flowing the integration domain into the complex domain, towards a steepest-descent manifold, and then integrating this manifold.

Method 2:
Finding points of the exponentiated phase function where the first derivative vanishes, and summing over relevant saddle points' contribution ("saddle-point method"). This only gives a good approximation of the integral as long as the relevant saddle points are well-separated (if you know you know). 

# To Do
- (everything)
- write tests
- include 2d code
- make nice example graphs
- write documentation
