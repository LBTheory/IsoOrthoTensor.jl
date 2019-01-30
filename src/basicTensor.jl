#------------------------------------------------------------------------------#
#                                    Module                                    #
#------------------------------------------------------------------------------#

"""
# Description

    module basicTensor

Basic tensor module.

# Usage

This module is automatically `include`d by the `iso`, and by the `orθ` packages,
which are automatically `include`d by the `IsoOrthoTensor` package.

```julia-repl
julia> using IsoOrthoTensor

julia> typeof(IOT.basicTensor)
Module

```
"""
module basicTensor

#------------------------------------------------------------------------------#
#                                   Exports                                    #
#------------------------------------------------------------------------------#

# basic Tensor Union type
export Tensor{T}

# Kronecker delta tensor in 1 to 3 Euclidean dimensions, fast and safe versions
export fK, K


#------------------------------------------------------------------------------#
#                                 Union Types                                  #
#------------------------------------------------------------------------------#

Tensor{T} = Union{T,Array{T}}


#------------------------------------------------------------------------------#
#                                  Constants                                   #
#------------------------------------------------------------------------------#

const δ = (1, [[1 0];[0 1]], [[1 0 0];[0 1 0];[0 0 1]])


#------------------------------------------------------------------------------#
#                                  Functions                                   #
#------------------------------------------------------------------------------#

"""
    fK(D::Int64 = 2)::Tensor{Int64}

Fast (no ckecks) `Int64` Kronecker  δ  tensor  in  a  `D`-dimensional  Euclidean
space, `D ∈ [1, 3]`.

This function is not exported and must be called from outside the module with  a
fully qualified name:

```julia-repl
julia> IOT.iso.fK(2)
2×2 Array{Int64,2}:
 1  0
 0  1

```
"""
fK(D::Int64 = 2)::Tensor{Int64} = δ[D]


"""
    K(D::Int64 = 2)::Tensor{Int64}

Returns the  Kronecker  Delta  tensor  in  a  `D`-dimensional  Euclidean  space,
checking bounds on `D`, i.e., whether `D ∈ [1, 3].`

```julia-repl
julia> K(3)
3×3 Array{Int64,2}:
 1  0  0
 0  1  0
 0  0  1

```
"""
function K(D::Int64 = 2)::Tensor{Int64}
    D <= 0 || D >= 4 ?
        throw(DomainError("dimension D = $D outside the valid domain [1, 3]")) :
        fK(D)
end


end
