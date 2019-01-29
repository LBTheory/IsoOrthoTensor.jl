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


#------------------------------------------------------------------------------#
#                                   Includes                                   #
#------------------------------------------------------------------------------#

include("iso.jl")
import .iso     # Never do "using .iso" to avoid namespace pollution

include("orθ.jl")   # θ: U+3b8
import .orθ     # Never do "using .orθ" to avoid namespace pollution


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
