#------------------------------------------------------------------------------#
#                                   Package                                    #
#------------------------------------------------------------------------------#

"""
# 1. Description

    module IsoOrthoTensor

Isotropic and orthogonal tensors [1 apud 2] (as Julia multi-dimensional  arrays)
in one to three Euclidean space dimensions.

Symbolic Hermite Polynomial Tensors — or N-dimensional Hermite polynomials [1] —
in one to three Euclidean space dimensions.

# 2. Usage

```julia-repl
julia> import HermitePolyTensor

```

# 3. References

[1]: H. Grad, “Note on N-dimensional  Hermite  polynomials,”  Communications  on
Pure and Applied Mathematics, vol. 2, no. 4, pp. 325–330, 1949.

[2]: K. K.  Mattila,  L.  A.  Hegele  Júnior,  P.  C.  Philippi,  “High-Accuracy
Approximation of High-Rank Derivatives: Isotropic Finite  Differences  Based  of
Lattice-Boltzmann Stencils,” The Scientific World Journal, vol. 2014, article ID
142907, 16 pages, 2014.
"""
module IsoOrthoTensor


#------------------------------------------------------------------------------#
#                                   Imports                                    #
#------------------------------------------------------------------------------#

using Reexport


#------------------------------------------------------------------------------#
#                                   Includes                                   #
#------------------------------------------------------------------------------#

include("version.jl")

include("IOT.jl")
@reexport using .IOT


end # module
