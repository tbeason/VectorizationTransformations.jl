# VectorizationTransformations.jl

<!-- Tidyverse lifecycle badges, see https://www.tidyverse.org/lifecycle/ Uncomment or delete as needed. -->
![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![build](https://github.com/tbeason/VectorizationTransformations.jl/workflows/CI/badge.svg)](https://github.com/tbeason/VectorizationTransformations.jl/actions?query=workflow%3ACI)
<!-- travis-ci.com badge, uncomment or delete as needed, depending on whether you are using that service. -->
<!-- [![Build Status](https://travis-ci.com/tbeason/VectorizationTransformations.jl.svg?branch=master)](https://travis-ci.com/tbeason/VectorizationTransformations.jl) -->
<!-- Coverage badge on codecov.io, which is used by default. -->
[![codecov.io](http://codecov.io/github/tbeason/VectorizationTransformations.jl/coverage.svg?branch=main)](http://codecov.io/github/tbeason/VectorizationTransformations.jl?branch=main)
<!-- Documentation -- uncomment or delete as needed -->
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://tbeason.github.io/VectorizationTransformations.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://tbeason.github.io/VectorizationTransformations.jl/dev)
-->
<!-- Aqua badge, see test/runtests.jl -->
<!-- [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) -->

This package provides common transformation operations related to matrix vectorization. 

- `vech` provides the half-vectorization
- `duplication_matrix` and `elimination_matrix` allow transformation between `vec(A)` and `vech(A)`
- `commutation_matrix` transforms `vec(A)` into `vec(A')`
- `symmetrizer_matrix` transforms `vec(A)` into `vec((A + A')/2)`

The matrix transformations are returned as sparse arrays to improve performance.

Refer to the function docstrings for additional detail.

## Usage

```julia
using VectorizationTransformations

n = 5
A = rand(n,n)
Asym = (A + A')/2

vech(Asym)

duplication_matrix(5)
elimination_matrix(5)
commutation_matrix(5)
symmetrizer_matrix(5)

? vech

? duplication_matrix
```





### Reference

Wikipedia pages and the following paper:

*The Elimination Matrix: Some Lemmas and Applications*
Jan R. Magnus and H. Neudecker
SIAM Journal on Algebraic Discrete Methods 1980 1:4, 422-449
