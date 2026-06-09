# QUADRILATERALS

#### integrate the quads
using FastGaussQuadrature



### utils



import Base.map
function map(p::Vector{Float64}, p1::PointA, p2::PointA, p3::PointA, p4::PointA)
    return ([p1.x,p1.y] .* (1. - p[1]) * (1. - p[2]) + 
            [p2.x,p2.y] .* (1. + p[1]) * (1. - p[2]) + 
            [p3.x,p3.y] .* (1. + p[1]) * (1. + p[2]) + 
            [p4.x,p4.y] .* (1. - p[1]) * (1. + p[2])) / 4.
end


function jacobian(p::Vector{Float64}, p1::PointA, p2::PointA, p3::PointA, p4::PointA)
    A = +(p1.x - p3.x) * (p2.y - p4.y) - (p1.y - p3.y) * (p2.x - p4.x)
    B = -(p1.x - p2.x) * (p3.y - p4.y) + (p1.y - p2.y) * (p3.x - p4.x)
    C = +(p2.x - p3.x) * (p1.y - p4.y) - (p2.y - p3.y) * (p1.x - p4.x)
    return (A + B * p[1] + C * p[2]) / 8.
end;

function integrate_quadrilateral(
    f::Function,
    quad::QuadC, n::Int64=7;
    prefactor::Function=(ti,tr) -> ones(2)
    )
    
        p1, p2, p3, p4 = quad.points

        x, w = gausslegendre(n);
        y = x;
        sum = [0. + 0im, 0. + 0im]
        for i=1:n, j=1:n
            jac = -jacobian([x[i], y[j]], p1, p2, p3, p4) # this minus sign here comes from that debugging experiment in the 2024-10-20 figures spectra... NB
            
            ti,tr = map([x[i], x[j]], p1, p2, p3, p4)
            action = f(ti,tr)
    
            sum = sum + jac * prefactor(ti, tr) * exp(action) * w[i] * w[j]
        end
        
    return sum
    
end


