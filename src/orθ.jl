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

# export ...


#------------------------------------------------------------------------------#
#                             Auxiliary Functions                              #
#------------------------------------------------------------------------------#



end # module
