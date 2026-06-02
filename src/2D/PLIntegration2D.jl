# module PLIntegration2D

### specify all packages I require
# using QuadGK
# using Sobol
# # using NLsolve

using LinearAlgebra
# using Integrals
using FiniteDiff
using StaticArrays
using Contour
# using Peaks
# using LinearAlgebra
# using FastGaussQuadrature


### specify all functions I want to make usable ###
# downwards flow
# determine if saddle point is relevant (necklace)
# get thimble
# integrate thimble
# 

include("necklace-for-SA-thimble.jl")
include("necklace-for-SD-thimble.jl")

include("utils.jl")


include("saddles-generic.jl")
include("saddles-contributing.jl")
include("saddle-point-method.jl")

include("flow-integration-path-downwards.jl")
# export get_thimble, integrate_thimble



# end