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

export K, 𝕔, Δ


#------------------------------------------------------------------------------#
#                             Auxiliary Functions                              #
#------------------------------------------------------------------------------#

"""
    K(Int::D = 2)::Union{Int,Array{Int}}

Returns the Kronecker Delta tensor in a `D`-dimensional Euclidean space.

"""
function K(Int::D = 2)::Union{Int,Array{Int}}
    # Validation
    if D <= 0 || D >= 4
        throw(DomainError("dimension D = $D outside the valid domain [1; 3]"))
    end
    # Execution
    return D == 1 ? 1 : D == 2 ? [[1 0]; [0 1]] : [[1 0 0]; [0 1 0]; [0 0 1]]
end

function 𝕔()
end


#------------------------------------------------------------------------------#
#                                User Functions                                #
#------------------------------------------------------------------------------#

function Δ()
end


end
