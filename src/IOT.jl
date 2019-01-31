#------------------------------------------------------------------------------#
#                                    Module                                    #
#------------------------------------------------------------------------------#

"""
# Description

    module IOT

High-level interface module to the `IsoOrthoTensor` package.

# Usage

This module is automatically `include`d by  the  `IsoOrthoTensor`  package,  and
provides the high-level interface to  the  Isotropic  and  Orthogonality  tensor
submodules.

```julia-repl
julia> using IsoOrthoTensor

julia> typeof(IOT)
Module

```
"""
module IOT


#------------------------------------------------------------------------------#
#                                   Imports                                    #
#------------------------------------------------------------------------------#

using Reexport


#------------------------------------------------------------------------------#
#                                   Includes                                   #
#------------------------------------------------------------------------------#

include("basicTensor.jl")
@reexport using .basicTensor

include("iso.jl")
@reexport using .iso

include("orθ.jl")   # θ: U+3b8
@reexport using .orθ


#------------------------------------------------------------------------------#
#                                   Exports                                    #
#------------------------------------------------------------------------------#

# NOTE: This module is @reexport'ed by the main package module


#------------------------------------------------------------------------------#
#                             Auxiliary Functions                              #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                             Interface Functions                              #
#------------------------------------------------------------------------------#


end # module
