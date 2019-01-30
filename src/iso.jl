#------------------------------------------------------------------------------#
#                                    Module                                    #
#------------------------------------------------------------------------------#

"""
# Description

    module iso

Isotropic tensor module.

# Usage

This module is automatically `include`d by the `IsoOrthoTensor` package, and  is
not used in isolation.

```julia-repl
julia> using IsoOrthoTensor

julia> typeof(IOT.iso)
Module

```


"""
module iso

#------------------------------------------------------------------------------#
#                                   Imports                                    #
#------------------------------------------------------------------------------#

import Combinatorics


#------------------------------------------------------------------------------#
#                                   Exports                                    #
#------------------------------------------------------------------------------#

export K, i𝕔, 𝕔, Δ


#------------------------------------------------------------------------------#
#                             Auxiliary Functions                              #
#------------------------------------------------------------------------------#

"""
    K(D::Int64 = 2)::Union{Int64,Array{Int64,2}}

Returns the Kronecker Delta tensor in a `D`-dimensional Euclidean space.

"""
function K(D::Int64 = 2)::Union{Int64,Array{Int64,2}}
    # Validation
    if D <= 0 || D >= 4
        throw(DomainError("dimension D = $D outside the valid domain [1; 3]"))
    end
    # Execution
    return D == 1 ? 1 : D == 2 ? [[1 0]; [0 1]] : [[1 0 0]; [0 1 0]; [0 0 1]]
end

"""
# Description

    i𝕔(
        OPD::Tuple{Vararg{Int64,N}} where N,
        FID::Tuple{Vararg{Int64,N}} where N;
        PRE = ([()], 1)
    )::Array{Tuple{Vararg{Int64,N}} where N,1} # 𝕔: U+1d554

Returns a one-dimensional array with tensor index combinations (as NTuples) that
encode summation terms for a combinatorial type of nonstandard  tensor  product,
referred to by the '⊛': Unicode U+229b symbol in some Lattice  Boltzmann  theory
literature, as [1].

`OPD` is an NTuple containing the operand dimensions, like (2, 2),  and  has  as
many elements as there are operands in the  combinatorial  type  of  nonstandard
tensor product.

`FID` is an NTuple containing the fixed indices, i.e., indices that take a fixed
position in all terms, like (1,), indicating that the first index is kept  fixed
in all summation terms.

`PRE` is **not** a user input and should be left alone. It  serves  to  pass  on
information across recursive calls.

# Usage

Nonstandard combinatorial tensor product between two rank-2 tensors, keeping the
first index fixed: arguments `(2, 2)` and  `(1,)`  to  the  `i𝕔`  function  —  a
summation of `3 choose 1 = 3` products between the operands:

```julia-repl
julia> using IsoOrthoTensor

julia> i𝕔((2, 2), (1,))
3-element Array{Tuple{Vararg{Int64,N}} where N,1}:
 (1, 2, 3, 4)
 (1, 3, 2, 4)
 (1, 4, 2, 3)

```

Nonstandard combinatorial tensor product between three rank-2  tensors,  keeping
the first and sixth indices fixed: arguments `(2, 2, 2)` and  `(1,  6)`  to  the
`i𝕔` function — a summation of `(4 choose 1)*(3 choose 2) = 4*3 =  12`  products
between the operands:

```julia-repl
julia> i𝕔((2, 2, 2), (1, 6))
12-element Array{Tuple{Vararg{Int64,N}} where N,1}:
 (1, 2, 3, 4, 5, 6)
 (1, 2, 3, 5, 4, 6)
 (1, 2, 4, 5, 3, 6)
 (1, 3, 2, 4, 5, 6)
 (1, 3, 2, 5, 4, 6)
 (1, 3, 4, 5, 2, 6)
 (1, 4, 2, 3, 5, 6)
 (1, 4, 2, 5, 3, 6)
 (1, 4, 3, 5, 2, 6)
 (1, 5, 2, 3, 4, 6)
 (1, 5, 2, 4, 3, 6)
 (1, 5, 3, 4, 2, 6)

```

# References

[1]: K. K.  Mattila,  L.  A.  Hegele  Júnior,  P.  C.  Philippi,  “High-Accuracy
Approximation of High-Rank Derivatives: Isotropic Finite  Differences  Based  of
Lattice-Boltzmann Stencils,” The Scientific World Journal, vol. 2014, article ID
142907, 16 pages, 2014.
"""
function i𝕔(
    OPD::Tuple{Vararg{Int64,N}} where N,
    FID::Tuple{Vararg{Int64,N}} where N;
    PRE = ([()], 1)
)::Array{Tuple{Vararg{Int64,N}} where N,1} # 𝕔: U+1d554
    RED = sum(OPD) + length(PRE[1][1]) # REsult Dimensions
    VID = Tuple(i for i in 1:RED if !(i in FID)) # Variable InDices
    ARR = Array{NTuple{N,Int} where N,1}()
    BEG = PRE[2]
    RNG = BEG:BEG+OPD[1]-1
    SLO = Tuple(i for i in RNG if !(i in FID)) # free SLOts in current operand
    for i in PRE[1]
        ALW = Tuple(k for k in VID if !(k in i)) # ALoWed free indices
        for j in Combinatorics.multiset_combinations(ALW, length(SLO))
            ITM = Int[]
            l = 1
            for k in RNG
                if k in FID
                    push!(ITM, k)
                else
                    push!(ITM, j[l])
                    l += 1
                end
            end
            push!(ARR, Tuple(vcat(i..., ITM...)))
        end
    end
    if length(OPD) <= 1
        return ARR
    else
        return i𝕔(OPD[2:end], FID, PRE = (ARR, BEG+OPD[1]))
    end
end

function 𝕔()
end


#------------------------------------------------------------------------------#
#                                User Functions                                #
#------------------------------------------------------------------------------#

function Δ()
end


end
