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
#                                 Union Types                                  #
#------------------------------------------------------------------------------#

Tensor{T} = Union{T,Array{T}}


#------------------------------------------------------------------------------#
#                                  Constants                                   #
#------------------------------------------------------------------------------#

const δ::NTuple{3,Tensor{Int64}} = (1, [[1 0];[0 1]], [[1 0 0];[0 1 0];[0 0 1]])


#------------------------------------------------------------------------------#
#                             Auxiliary Functions                              #
#------------------------------------------------------------------------------#

"""
    fK(D::Int64 = 2)::Tensor{Int64}

Fast (no ckecks)`Int64` Kronecker δ tensor in a `D`-dimensional Euclidean space,
`D ∈ [1, 3]`.
"""
fK(D::Int64 = 2)::Tensor{Int64} = δ[D]


"""
    K(D::Int64 = 2)::Tensor{Int64}

Returns the  Kronecker  Delta  tensor  in  a  `D`-dimensional  Euclidean  space,
checking bounds on `D`, i.e., whether `D ∈ [1, 3].`
"""
function K(D::Int64 = 2)::Tensor{Int64}
    D <= 0 || D >= 4 ?
        throw(DomainError("dimension D = $D outside the valid domain [1, 3]")) :
        fK(D)
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

"""
# Description

    function 𝕔(
        OPS::Tuple{Vararg{Tensor{Int64},N}} where N
        FID::Tuple{Vararg{Int64,N}} where N;
        D::Int64 = 2
    )::Tensor{Int64} # 𝕔: U+1d554

Performs  a  combinatorial  type  of  nonstandard  tensor  product  between  the
operands, which are elements of the `OPS` `Tuple`,  keeping  the  indices  `FID`
fixed in a `D`-dimensional Euclidean space, and returns  the  resulting  tensor.
The combinatorial nonstandard tensor product performed is  referred  to  by  the
'⊛': Unicode U+229b symbol in some Lattice Boltzmann theory literature, as [1].

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
julia> δ = K(2)
2×2 Array{Int64,2}:
 1  0
 0  1

julia> 𝕔((δ, δ), (1,))
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

# References

[1]: K. K.  Mattila,  L.  A.  Hegele  Júnior,  P.  C.  Philippi,  “High-Accuracy
Approximation of High-Rank Derivatives: Isotropic Finite  Differences  Based  of
Lattice-Boltzmann Stencils,” The Scientific World Journal, vol. 2014, article ID
142907, 16 pages, 2014.
"""
function 𝕔(
    OPS::Tuple{Vararg{Tensor{Int64},N}} where N,
    FID::Tuple{Vararg{Int64,N}} where N;
    D::Int64 = 2
)::Tensor{Int64} # 𝕔: U+1d554
    OPD = Tuple(ndims(A) for A in OPS) # OPerand Dimensions
    RED = sum(OPD) # REsult Dimensions
    RET = fill(zero(OPS[1][1]), Tuple(D for i in 1:RED)) # RETurn tensor
    ONE = one(OPS[1][1])
    for idx in Base.product(fill(1:D, RED)...) # Return tensor index
        for IDP in i𝕔(OPD, FID) # InDex Permutations
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

    Δ(n::Int64; D::Int64 = 2)::Tensor{Int64}

Computes and returns an `n`-th order Isotropic Tensor (of `Int64` components) in
a `D`-dimensional Euclidean space, checking bounds on `D`, i.e.,  whether  `D  ∈
[1, 3]`, and on `n`, i.e., whether `n ∈ [0, ∞)`.
"""
function Δ(n::Int64; D::Int64 = 2)::Tensor{Int64}
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
        return n == 1 ? δ[D] : 𝕔((δ[D], Δ(n-1, D = D)), (1,), D = D)
    end
end


end
