### normalised gradient

function gradN(
    f_grad::Function,
    ti::ComplexF64, tr::ComplexF64,
    thresh::Float64 = 1.)

    g = f_grad(ti,tr)
    if norm(g) > thresh # bit lower than the gradient at the saddle point
        return LinearAlgebra.normalize(g)
    else 
        return g
    end
end;



### POINT
mutable struct PointA{T}
    x::T
    y::T
    active::Bool
end
PointA(x,y) = PointA(x,y,true);

import Base.isequal
isequal(p1::PointA, p2::PointA) = isequal(p1.x,p2.x) && isequal(p1.y,p2.y)
Base.hash(p::PointA, h::UInt) = hash([p.x,p.y], h)

xy(p::PointA) = (p.x, p.y);
Base.getindex(pa::PointA, i) = xy(pa)[i]

dist(p1::PointA, p2::PointA) = norm([p1.x-p2.x, p1.y-p2.y])


function float2complex(p::PointA)
    return PointA(complex(p.x), complex(p.y), p.active)
end
toR4(p::PointA) = toR4(p.x, p.y)

midx(p1::PointA,p2::PointA) = (p2.x + p1.x)./2
# midx(ps::Vector{PointA{T}}) where T<:Number = sum([p.x for p in ps])/length(ps)
midy(p1::PointA,p2::PointA) = (p2.y + p1.y)./2
# midy(ps::Vector{PointA{T}}) where T<:Number = sum([p.y for p in ps])/length(ps)

midpoint(p1::PointA,p2::PointA) = PointA(midx(p1,p2), midy(p1,p2))
# midpoint(Ps::Vector{PointA{T}}) where T<:Number = PointA(midx(Ps), midy(Ps));



### doesn't seem to be needed
# function point2vec(p::PointA2)
#     return [reim(p.x)...,reim(p.y)...]
# end;


struct QuadC ### a quadrilateral with coordinates
    points::Vector{PointA} # length = 4
end

struct TriangleC ### a triangle with coordinates
    points::Vector{PointA} # length = 3
end





## initialising the original integration domain

function make_init_points_rectangle(t1min::Real, t1max::Real,
    t2min::Real, t2max::Real,
    flow_bounds=[true, true, true, true],
    point_type=PointA)
    [point_type(complex(t1min), complex(t2min), flow_bounds[1]), 
    point_type(complex(t1min), complex(t2max), flow_bounds[2]), 
    point_type(complex(t1max), complex(t2max), flow_bounds[3]),
    point_type(complex(t1max), complex(t2min), flow_bounds[4])]
end

function make_init_points_parallelogram(timin::Real, timax::Real,
    ttmin::Real, ttmax::Real,
    flow_bounds=[true, true, true, true],
    point_type=PointA)
    [point_type(complex(timin), complex(timin+ttmin), flow_bounds[1]), 
    point_type(complex(timin), complex(timin+ttmax), flow_bounds[2]), 
    point_type(complex(timax), complex(timax+ttmax), flow_bounds[3]),
    point_type(complex(timax), complex(timax+ttmin), flow_bounds[4])] 
end


