
# TRIANGLES




# function integrate_triangle(points, triangle, integrand; order=1,dim=2)
#     int = zeros(ComplexF64,dim)
#     area = triangle_area_new(xy.(points[triangle.indices]))
#     wi = area/3

#     v1,v2,v3 = triangle.indices
#     t_edges = [(v1,v2),(v2,v3),(v3,v1)] ### indices of the specific one I'm looking at
#     for k in t_edges
#         p1,p2 = points[[k...]]
#         ### is this fine or does this need to be done in the projection???
#         midpoint = ((p1.x+p2.x)/2, (p1.y+p2.y)/2) 
#         int .+= wi*integrand(midpoint...)
#     end    
#     return int
# end

### VERSION ONE

midx(p1::PointA,p2::PointA) = (p2.x + p1.x)./2
midy(p1::PointA,p2::PointA) = (p2.y + p1.y)./2
midpoint(p1::PointA,p2::PointA) = PointA(midx(p1,p2), midy(p1,p2));


function triangle_to_fake_quadrilateral(tri::TriangleC)
    t_edges = [(1,2),(2,3),(3,1)]
    d_edges = [dist(tri.points[k[1]],tri.points[k[2]]) for k in t_edges]
    val, idx = findmax(d_edges)
    longest_edge = t_edges[idx]
    
    extra_point = midpoint(tri.points[longest_edge[1]], tri.points[longest_edge[2]])
    a,b,opp = longest_edge..., only(setdiff([1,2,3], [longest_edge...]))
    fake_quad_points = [tri.points[a], extra_point, tri.points[b], tri.points[opp]]
    return QuadC([PointA(p.x, p.y) for p in fake_quad_points])   
end

function number_of_arguments(f::Function)
    mg = methods(f).ms[1]

    arity = length(mg.sig.parameters) - 1
    return arity
end  

function integrate_triangle_cheat(f::Function, triangle::TriangleC;
    prefactor::Function=(ti,tr) -> ones(2), order=7,dim=2)
    
    @assert number_of_arguments(f) == 2
    fake_quad = triangle_to_fake_quadrilateral(triangle)
    int_quad = integrate_quadrilateral(f, fake_quad, order; prefactor=prefactor)
    return int_quad
end

### ALTERNATIVE VERSION, DOESN'T ALWAYS GIVE THE SAME RESULT AND IDK YET WHY

using SimplexQuad

function jacobian(tri::TriangleC)
    p1,p2,p3=tri.points
    [p2.x-p1.x p3.x-p1.x;
    p2.y-p1.y p3.y-p1.y]
end

function integrate_triangle(f::Function, triangle::TriangleC;prefactor=(ti, tr)->ones(2), order=1,dim=2)
    p1,p2,p3=triangle.points

    jac = jacobian(triangle)

    g1(s,t) = p1.x - (p2.x-p1.x)*s + (p3.x-p1.x)*t
    g2(s,t) = p1.y - (p2.y-p1.y)*s + (p3.y-p1.y)*t
    
    int = zeros(ComplexF64,dim)

    X, W = simplexquad(Float64,order,2) # the 2 means it's a triangle 
    for i in 1:length(W)
        int +=  W[i] * prefactor(g1(X[i,:]...), g2(X[i,:]...)) * exp(f(g1(X[i,:]...), g2(X[i,:]...))) * det(jac)
    end
    return -int  
end;

nothing