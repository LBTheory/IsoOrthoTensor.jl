#------------------------------------------------------------------------------#
#                                    Module                                    #
#------------------------------------------------------------------------------#

"""
# Description

    module orθ  # θ: U+3b8

Orthogonality tensor module.

# Usage

This module is automatically `include`d by the `IsoOrthoTensor` package, and  is
not used in isolation.

```julia-repl
julia> using IsoOrthoTensor

julia> typeof(IOT.orθ)
Module

```
"""
module orθ


#------------------------------------------------------------------------------#
#                                   Imports                                    #
#------------------------------------------------------------------------------#

include("basicTensor.jl")
using .basicTensor

import Combinatorics


#------------------------------------------------------------------------------#
#                                   Exports                                    #
#------------------------------------------------------------------------------#

export i𝕡, 𝕡, Ο, O


#------------------------------------------------------------------------------#
#                             Auxiliary Functions                              #
#------------------------------------------------------------------------------#

"""
# Description

    i𝕡(
        OPD::Tuple{Vararg{Int64,N}} where N,
        FID::Tuple{Vararg{Int64,N}} where N;
        PRE = ([()], 1)
    )::Array{Tuple{Vararg{Int64,N}} where N,1} # 𝕡: U+1d561

Returns a one-dimensional array with tensor index permutations (as NTuples) that
encode summation terms for a permutatorial type of nonstandard  tensor  product,
referred to by a particular nesting of the '⊛': Unicode U+229b  symbol  in  some
Lattice Boltzmann theory literature, as [1].

`OPD` is an NTuple containing the operand dimensions, like (2, 2),  and  has  as
many elements as there are operands in the  permutatorial  type  of  nonstandard
tensor product.

`FID` is an NTuple containing the fixed indices, i.e., indices that take a fixed
position in all terms, like (1,), indicating that the first index is kept  fixed
in all summation terms.

`PRE` is **not** a user input and should be left alone. It  serves  to  pass  on
information across recursive calls.

# Usage

Nonstandard permutatorial tensor product between two rank-2 tensors, keeping the
first index of each one fixed: arguments `(2, 2)`  and  `(1,  3)`  to  the  `i𝕡`
function — a summation of `2! = 2` products between the operands:

```julia-repl
julia> using IsoOrthoTensor

julia> i𝕡((2, 2), (1, 3))
2-element Array{Tuple{Vararg{Int64,N}} where N,1}:
 (1, 2, 3, 4)
 (1, 4, 3, 2)

```

Nonstandard permutatorial tensor product between three rank-2  tensors,  keeping
the first index of each one fixed: arguments `(2, 2, 2)` and `(1, 3, 5)` to  the
`i𝕡` function — a summation of `3! = 6` products between the operands:

```julia-repl
julia> i𝕡((2, 2, 2), (1, 3, 5))
6-element Array{Tuple{Vararg{Int64,N}} where N,1}:
 (1, 2, 3, 4, 5, 6)
 (1, 2, 3, 6, 5, 4)
 (1, 4, 3, 2, 5, 6)
 (1, 4, 3, 6, 5, 2)
 (1, 6, 3, 2, 5, 4)
 (1, 6, 3, 4, 5, 2)

```

# References

[1]: K. K.  Mattila,  L.  A.  Hegele  Júnior,  P.  C.  Philippi,  “High-Accuracy
Approximation of High-Rank Derivatives: Isotropic Finite  Differences  Based  of
Lattice-Boltzmann Stencils,” The Scientific World Journal, vol. 2014, article ID
142907, 16 pages, 2014.
"""
function i𝕡(
    OPD::Tuple{Vararg{Int64,N}} where N,
    FID::Tuple{Vararg{Int64,N}} where N;
    PRE = ([()], 1)
)::Array{Tuple{Vararg{Int64,N}} where N,1} # 𝕡: U+1d561
    RED = sum(OPD) + length(PRE[1][1]) # REsult Dimensions
    VID = Tuple(i for i in 1:RED if !(i in FID)) # Variable InDices
    ARR = Array{NTuple{N,Int} where N,1}()
    BEG = PRE[2]
    RNG = BEG:BEG+OPD[1]-1
    SLO = Tuple(i for i in RNG if !(i in FID)) # free SLOts in current operand
    for i in PRE[1]
        ALW = Tuple(k for k in VID if !(k in i)) # ALoWed free indices
        for j in Combinatorics.multiset_permutations(ALW, length(SLO))
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
        return i𝕡(OPD[2:end], FID, PRE = (ARR, BEG+OPD[1]))
    end
end


"""
# Description

    function 𝕡(
        OPS::Tuple{Vararg{Tensor{Int64},N}} where N
        FID::Tuple{Vararg{Int64,N}} where N;
        D::Int64 = 2
    )::Tensor{Int64} # 𝕡: U+1d561

Performs  a  permutatorial  type  of  nonstandard  tensor  product  between  the
operands, which are elements of the `OPS` `Tuple`,  keeping  the  indices  `FID`
fixed in a `D`-dimensional Euclidean space, and returns  the  resulting  tensor.
The permutatorial nonstandard tensor product performed is referred to by nesting
the '⊛': Unicode U+229b symbol in some Lattice Boltzmann theory  literature,  as
[1].

`OPS` is an NTuple containing the operand tensors, like `(δ, δ)`  —  a  pair  of
Kronecker-delta tensors stored in the temporary `δ` identifier.

`FID` is an NTuple containing the fixed indices, i.e., indices that take a fixed
position in all terms, like (1,), indicating that the first index is kept  fixed
in all summation terms.

`D` is the dimensionality of the Euclidean space in which  the  tensor  operands
apply.

# Usage

    Calculate the `𝝙⁽ⁿ⁾` tensor of order `2n` that is isotropic with respect to  all
    of its `2n` indices [1] with `n=2` through the nonstandard product

        Δ⁽²⁾αβγε = δαβ*δγε + δαγ*δβε + δαε*δβγ,

    in which terms combine the three free indices  `βγε`  while  keeping  the  first
    index `α` fixed.

```julia-repl
```

# References

[1]: K. K.  Mattila,  L.  A.  Hegele  Júnior,  P.  C.  Philippi,  “High-Accuracy
Approximation of High-Rank Derivatives: Isotropic Finite  Differences  Based  of
Lattice-Boltzmann Stencils,” The Scientific World Journal, vol. 2014, article ID
142907, 16 pages, 2014.
"""
function 𝕡(
    OPS::Tuple{Vararg{Array{Int64},N}} where N,
    FID::Tuple{Vararg{Int64,N}} where N;
    D::Int64 = 2
)::Tensor{Int64} # 𝕡: U+1d561
    OPD = Tuple(ndims(A) for A in OPS) # OPerand Dimensions
    RED = sum(OPD) # REsult Dimensions
    RET = fill(zero(OPS[1][1]), Tuple(D for i in 1:RED)) # RETurn tensor
    ONE = one(OPS[1][1])
    for idx in Base.product(fill(1:D, RED)...) # Return tensor index
        for IDP in i𝕡(OPD, FID) # InDex Permutations
            PID = [idx[perm] for perm in IDP] # Permuted InDices
            prod, BEG = ONE, 1
            for fact in OPS
                SIZ = ndims(fact)
                prod *= fact[PID[BEG:BEG-1+SIZ]...]
                BEG += SIZ
            end
            RET[idx...] += prod
        end
    end
    return RET
end


#------------------------------------------------------------------------------#
#                                User Functions                                #
#------------------------------------------------------------------------------#

"""
# Description

    Ο(n::Int64; D::Int64 = 2)::Tensor{Int64} # Ο: U+39f

Computes  and  returns  an  `n`-th  order  Orthogonality  Tensor   (of   `Int64`
components) in a `D`-dimensional Euclidean space, checking bounds on `D`,  i.e.,
whether `D ∈ [1, 3]`, and on `n`, i.e., whether `n ∈ [0, ∞)`.

```julia-repl
julia> Ο(1) # Ο: U+39f
2×2 Array{Int64,2}:
 1  0
 0  1

```
"""
function Ο(n::Int64; D::Int64 = 2)::Tensor{Int64} # Ο: U+39f
    # Validation
    if D <= 0 || D >= 4
        throw(DomainError("dimension D = $D outside the valid domain [1, 3]"))
    elseif n < 0
        throw(DomainError("order n = $n outside the valid domain [0, ∞)"))
        # ∞: U+221e
    end
    # Execution
    if D == 1 || n == 0
        return 1
    else
        if n == 1
            return K(D)
        else
            return 𝕡(
                Tuple(repeat([K(D)], n)),
                Tuple([i for i in 1:2:2n]),
                D = D
            )
        end
    end
end


"""
Function O (Capital Latin Letter O) is an alias of the Ο (Capital  Greek  Letter
Omicron): U+39f. For help, see the Ο (omicron) function documentation: ?Ο
"""
O = Ο # Capital o = Capital omicron


end # module
