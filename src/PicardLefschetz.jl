module PicardLefschetz

### make sure they're added to the Project.toml by doing ] activate .,  add PackageName
using LinearAlgebra


### for finding critical points
using NLsolve
using Sobol

### for the contour intersection
using Contour
using StaticArrays
using GeometryBasics


### I probably want this at some point
# import Base: show

include("1D/integration-path-flow.jl")

include("1D/critical-points.jl")
include("1D/line-intersection.jl")
include("1D/saddle-point-method.jl")



export square

# Write your package code here.
function hello() 
	println("hello, World")
end

"""
    square(x)

Returns the square of `x`, for whatever type of argument.
"""
square(x) = x^2


end # PicardLefschetz
