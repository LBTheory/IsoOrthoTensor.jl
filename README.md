# IsoOrthoTensor.jl

Isotropic and Orthogonal Tensors — Provides functions to generate such tensors

# 1. Description

This package allows for obtaining Isotropic and Orthogonal Tensors of any  order
`n`  (conversely,  of  rank  `2n`),  from  operators   that   generalize   their
definitions, available in [1], in one- to three-dimensional Euclidean spaces.

Only tensors of `Int64` elements are of interest, and a bare-bones  tensor  type
definition is provided, i.e., the one line:

    Tensor{T} = Union{T,Array{T}}

from the `basicTensor` submodule, so for all practical purposes,  the  isotropic
and orthogonality tensors returned are of  `Array{Int64}`  type  (or  plain  old
`Int64` type for one-dimensional or for zero-order/rank tensors).

The package is not concerned with implementing tensor type functionality  —  the
reason why the package name brings "Tensor" in  the  singular,  following  Julia
language documentation recommendations.

# 2. Usage

To import the package,  one  can  use  Julia's  standard  `using`  and  `import`
statements. Care has been taken as not to flood the namespace with  symbols,  so
that only submodules and user functions are exported  to  the  target  namespace
(according to Julia's import/using rules).

Obtaining  Isotropic  Tensors  of  order  `n`  (rank  `2n`)  in  `D`-dimensional
Euclidean spaces:

```julia-repl
julia> using IsoOrthoTensor

julia> Δ(1, D = 3) # Δ: U+394
3×3 Array{Int64,2}:
 1  0  0
 0  1  0
 0  0  1

julia> Δ(2, D = 2) # Δ: U+394
2×2×2×2 Array{Int64,4}:
[:, :, 1, 1] =
 3  0
 0  1

[:, :, 2, 1] =
 0  1
 1  0

[:, :, 1, 2] =
 0  1
 1  0

[:, :, 2, 2] =
 1  0
 0  3

```

Obtaining Orthogonality Tensors of order  `n`  (rank  `2n`)  in  `D`-dimensional
Euclidean spaces:

```julia-repl
julia> Ο(1, D = 3) # Ο: U+39f
3×3 Array{Int64,2}:
 1  0  0
 0  1  0
 0  0  1

julia> Ο(2, D = 2) # Ο: U+39f
2×2×2×2 Array{Int64,4}:
[:, :, 1, 1] =
 2  0
 0  1

[:, :, 2, 1] =
 0  1
 0  0

[:, :, 1, 2] =
 0  0
 1  0

[:, :, 2, 2] =
 1  0
 0  2

```

# 3. Performance

The mathematical generalization of the Tensor definitions  (obtained  by  adding
products of the Kronecker delta with either a combination of  a  permutation  of
certain free indices across the operands)  can  be  expressed  in  a  short  and
concise way, but is definitely not the most efficient algorithm.

Therefore, although the function definitions that generate such  tensors  are  a
few lines long, the resulting algorithm is of a high  (computation)  complexity,
and might be optimized in future versions of the package.

Getting a correctly working implementation was the chief goal of this  inceptive
version.

Examples:

Notice the time increment from a 6-th order to a 7-th order Isotropy tensor  (of
ranks 12 and 14, respectively) in a 2-D Euclidean space (the default), each  one
having 2¹² and 2¹⁴ components, respectively:

```julia-repl
julia> @time Set(Δ(6)) # Δ: U+394
  1.737888 seconds (9.07 M allocations: 618.298 MiB, 4.91% gc time)
Set([945, 315, 0, 10395, 225])

julia> @time Set(Δ(7)) # Δ: U+394
  7.372711 seconds (43.81 M allocations: 3.272 GiB, 5.64% gc time)
Set([0, 2835, 135135, 10395, 1575])

```

A similar comparison is made for the Orthogonality tensor. A call to `Ο(n)` will
ensue in a summation of `n!` tensor products with `n` operands each.

```julia-repl
julia> @time Set(Ο(5)) # Ο: U+39f
  1.392121 seconds (12.24 M allocations: 812.973 MiB, 9.70% gc time)
Set([120, 0, 12, 24])

julia> @time Set(Ο(6)) # Ο: U+39f
 35.093839 seconds (301.58 M allocations: 20.094 GiB, 9.48% gc time)
Set([720, 120, 36, 0, 48])

```

There isn't any hardcoded upper limit for the value of `n` accepted by `Δ()` and
by `Ο()`, so bear in mind the overwhelmingly high algorithm complexity.

# References

[1]: K. K.  Mattila,  L.  A.  Hegele  Júnior,  P.  C.  Philippi,  “High-Accuracy
Approximation of High-Rank Derivatives: Isotropic Finite  Differences  Based  of
Lattice-Boltzmann Stencils,” The Scientific World Journal, vol. 2014, article ID
142907, 16 pages, 2014.
