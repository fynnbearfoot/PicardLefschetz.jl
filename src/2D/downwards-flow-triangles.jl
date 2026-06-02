#### OTHER UTILS
Base.getindex(z::Iterators.Zip, i) = (it -> getindex(it, i)).(z.is)
toR4(z1::ComplexF64, z2::ComplexF64) = [real(z1), imag(z1), real(z2), imag(z2)]
toR4(z1::Number, z2::Number) = [real(z1), imag(z1), real(z2), imag(z2)]




#### POINTS
mutable struct PointA{T}
    x::T
    y::T
    active::Bool
end
PointA(x,y) = PointA(x,y,true);

xy(p::PointA) = (p.x, p.y);
Base.getindex(pa::PointA, i) = xy(pa)[i]

dist(p1::PointA, p2::PointA) = norm([p1.x-p2.x, p1.y-p2.y])

import Base.isequal
isequal(p1::PointA, p2::PointA) = isequal(p1.x,p2.x) && isequal(p1.y,p2.y)
Base.hash(p::PointA, h::UInt) = hash([p.x,p.y], h)

function float2complex(p::PointA)
    return PointA(complex(p.x), complex(p.y), p.active)
end
toR4(p::PointA) = toR4(p.x, p.y)


function float2complex(p::PointA)
    return PointA(complex(p.x), complex(p.y), p.active)
end

xy(p::PointA) = (p.x, p.y);



#### TRIANGLES
mutable struct TriangleA
    coord::MVector{3,Int} ### this could be an MVector
    active::Bool
    TriangleA(coord::AbstractVector{T}) where T <:Integer = new(coord, true)
end

function project_onto_triangle(base::AbstractVector, points::AbstractVector) # what are these types?
    @assert length(base) == 3 "Need exactly 3 vertices"

    # Convert to R^4
    vs = [toR4(z.x, z.y) for z in base]

    # Basepoint
    v1, v2, v3 = vs
    e1 = v2 - v1
    e2 = v3 - v1

    # Gram-Schmidt to find orthonormal basis
    u1 = e1 / norm(e1)
    e2_proj = e2 - (u1 ⋅ e2) * u1
    u2 = e2_proj / norm(e2_proj)

    # Projection function
    proj(v) = [(u1 ⋅ (v - v1)), (u2 ⋅ (v - v1))]

    # Project all points
    ps = [proj(toR4(p.x, p.y)) for p in points]

    return ps
end

using MiniQhull
flags_original = "qhull d Qt Qbb Qc Qz"

prod(points::AbstractVector) = (points[2][1]-points[1][1])*(points[3][2]-points[1][2]) - (points[2][2]-points[1][2])*(points[3][1]-points[1][1])

function subdivide_triangle_new!(points, t_vertices::TriangleA; Δ::Real=0.5, delaunay_flags=flags_original)
        og_angle = angle(prod(points[t_vertices.coord]))
	    vertices_here = Vector{Int}()
	    push!(vertices_here, t_vertices.coord...)
	    v1,v2,v3 = t_vertices.coord
	    t_edges = [(v1,v2),(v2,v3),(v3,v1)] ### indices of the specific one I'm looking at
	    for k in t_edges
	        p1,p2 = points[[k...]]
	        d = dist(p1,p2)
	        if d > Δ ### subdivide
	            
	            ### how many points to insert in between?
	            ### wanna make sure that the new distance is not larger than Δ
	            n = floor(Int64, d/Δ)
	            xvals = range(p1.x, stop = p2.x, length = n+2)
	            yvals = range(p1.y, stop = p2.y, length = n+2)
	            new_points = [PointA(c...) for c in zip(xvals[2:end-1], yvals[2:end-1])] # I don't want to include the actual points themselves

	            for np in new_points
	                already_created = findall(p->isequal(round.(xy(p),digits=8), round.(xy(np),digits=8)), points) 
	                if isempty(already_created)
	                    push!(points, np)
	                    push!(vertices_here, length(points))
	       #         println("point is new. new index: $(length(points)),  length points: $(length(points))")
	                else
	       #             println("point is already created. length points: $(length(points))")
	                    np_vertex = already_created[1]
	                    push!(vertices_here, np_vertex)
	                    ### add a pointer to the existing point
	                end
	            end

	        end
	    end
	    
	    if length(vertices_here) == 3
	        return points, Vector{TriangleA}[]
	    else
	        ### project only the points that are relevant for this new triangle
	        projected_points = project_onto_triangle(points[t_vertices.coord], points[vertices_here])
	        connections_here = delaunay(projected_points, flags_original) ### they are in the index system of this triangle!

	        ### check that they are not degenerate!	        
	        connections_here_vectors = [ collect(col) for col in eachcol(connections_here)]
	        filter!(ind -> triangle_area_new(projected_points[ind]) > 1e-5,
	            connections_here_vectors)

	        # new_triangles = [TriangleA(vertices_here[inds]) for inds in connections_here_vectors]
            new_triangles=TriangleA[]
            ### the following is just to ensure that the newly created triangles have the same orientation
            for inds in connections_here_vectors 
                verts = vertices_here[inds]
                if dot(og_angle,angle(prod(points[verts]))) > 0 
                    push!(new_triangles, TriangleA(verts))
                else
                    push!(new_triangles, TriangleA(reverse(verts)))
                end
            end
	        return points, new_triangles
	    end
end


function dissect_thimbles(triangles)
    active_triangles = filter(s->s.active, triangles)
    thimbles = Vector()
    visited = falses(length(active_triangles))
   
    function dfs!(active_triangles, visited, i_start)
        stack = [i_start] # "open ends" to explore (= the four corner point indices of the quadrilateral)
        trace = Int[] # trace is more like a carpet now I guess
        while !isempty(stack)
            v = pop!(stack)

            if !visited[v]
			# #                 @show v
                visited[v] = true
                push!(trace, v)

                ### find the neighbouring quads to the vertices
                for vertex in 1:3
                    # find all quads that corner to you
                    nexts = findall(sim -> in(active_triangles[v].coord[vertex], sim.coord), active_triangles)
					#                     @show nexts
                    append!(stack, filter(!isequal(v),nexts)) # obviously ignore the quad that you're in atm
                end
                unique!(stack)
            end
        end
        return trace
    end
    
    for i in 1:length(active_triangles)
		#     i = 1
        if !visited[i]
            trace = dfs!(active_triangles, visited, i)
            push!(thimbles, active_triangles[trace])
        end
    end
    
    return thimbles
end


function triangle_area_new(points::AbstractVector)
    @assert length(points) == 3
    abs((points[2][1]-points[1][1])*(points[3][2]-points[1][2]) - (points[2][2]-points[1][2])*(points[3][1]-points[1][1])) / 2
end
# triangle_area(points::AbstractVector{PointA}) = triangle_area([xy(p) for p in points])


# function integrate_triangle(points, triangle, integrand; order=1,dim=2)
#     int = zeros(ComplexF64,dim)
#     area = triangle_area_new(xy.(points[triangle.coord]))
#     wi = area/3

#     v1,v2,v3 = triangle.coord
#     t_edges = [(v1,v2),(v2,v3),(v3,v1)] ### indices of the specific one I'm looking at
#     for k in t_edges
#         p1,p2 = points[[k...]]
#         ### is this fine or does this need to be done in the projection???
#         midpoint = ((p1.x+p2.x)/2, (p1.y+p2.y)/2) 
#         int .+= wi*integrand(midpoint...)
#     end    
#     return int
# end





#### DOWNWARDS FLOW


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



function initialise_time_grid(timin, timax, ttmin, ttmax, Δ;
        flow_bounds=[true, true, true, true])
    
    points = [
        PointA(timin, timin+ttmin, flow_bounds[1]), 
        PointA(timin, timin+ttmax, flow_bounds[2]), 
        PointA(timax, timax+ttmax, flow_bounds[3]),
        PointA(timax, timax+ttmin, flow_bounds[4])]      
    
    pts = [SVector(p.x, p.y) for p in points]
    connections = delaunay(pts) #, "qhull d Qbb Qc QJ Pp")
    #     "qhull d Qbb Qc QJ Pp"
	#         display(connections)
    triangles = [TriangleA(inds) for inds in eachcol(connections)]
	#     @show points
    points = [float2complex(p) for p in points]
	#     @show points
    
    ### filter out degenerate points!
    n_old = length(triangles)
    n_new = n_old + 1

    i_while = 0
    while (n_old != n_new)
        n_old = n_new 
        i_while += 1
        i = 1    
    
        for i_t in eachindex(triangles)
            triangle = triangles[i_t]
            if triangle.active
		#                 points_this_triangle = points[triangle.coord]
                points, new_connections = subdivide_triangle_new!(points, triangle, Δ=Δ)
                if !isempty(new_connections)
                    triangle.active = false
                    append!(triangles, new_connections)
                end

            end
            if i>30 break end
        end

        if i_while > 10 break end
        filter!(sim->sim.active, triangles)
        n_new = length(triangles)

    end         
    return points, triangles
end



function initialise_triangles_on_real_plane(xmin, xmax, ymin, ymax, Δx, Δy)
    xrange = range(real(xmin), stop = real(xmax), length = Int(ceil(real(xmax-xmin)/Δx)))
    yrange = range(real(ymin), stop = real(ymax), length = Int(ceil(real(ymax-ymin)/Δy)))

    points_ini = Vector{PointA}()
    for xr in collect(xrange)
        for yr in collect(yrange)
            push!(points_ini, PointA(xr, yr,true))
        end
    end
    pts = [SVector(p.x, p.y) for p in points_ini]
    connections = delaunay(pts) #, "qhull d Qbb Qc QJ Pp")
    #     "qhull d Qbb Qc QJ Pp"
	#         display(connections)

	### filter out degenerate points!

    #### THIS REQUIRES FURTHER DEVELOPMENT

    return points_ini, connections
end


function flow_down!(triangles, points::Vector{PointA{T}},
        f::Function,
        f_grad::Function;
        threshold::Real=0.5, # for normalisation of the gradient
        δ::Real=0.5, # flowstepfactor
        h_threshold::Real=-20.
        ) where T<:Number

    for i1 in 1:length(points)
        if points[i1].active # for the active points
            step = -δ .* gradN((ti,tr) -> conj.(complex.(f_grad(ti,tr))), points[i1].x +0im, points[i1].y +0im, threshold)
            points[i1].x += step[1]
            points[i1].y += step[2]

            ### turning them inactive when they are below the threshold. Doing this here prevents lonely relict points from flowing if their triangle has turned inactive already.
            if real(f(xy(points[i1])...)) < h_threshold 
                points[i1].active = false
            end
        end
    end

    for i2 in eachindex(triangles)
        if triangles[i2].active
            for v in triangles[i2].coord
                if !points[v].active
                # if real(f(points[v].x, points[v].y)) < h_threshold 
                    triangles[i2].active = false # am I sure that I want to turn the whole simplex inactive?
                    # points[v].active = false
                end
            end
        end
    end
end


function subdivide_triangles!(points, triangles, subdividethreshold)
    n_old = length(triangles)
    n_new = n_old + 1

    i_while = 0
    while (n_old != n_new)
        n_old = n_new 
        i_while += 1
        # i = 1

        ### in subdivide sub routine ideally
        for i_t in eachindex(triangles)
            triangle = triangles[i_t]
            if triangle.active
                if all([p.active for p in points[triangle.coord]])
                    points, new_connections = subdivide_triangle_new!(points, triangle, Δ=subdividethreshold)
					# n_actives = sum([t.active for t in triangles])
                    if !isempty(new_connections)
                        triangle.active = false
                        append!(triangles, new_connections)
                    end
                elseif !any([p.active for p in points[triangle.coord]])
                    triangle.active = false
                    # println("here i have turned a triangle inactive")
                end
            end
        end
        
		#  if i_while > 10 break end
        filter!(sim->sim.active, triangles)
        n_new = length(triangles)
    end
end



#### making SD thimbles
using Graphs, SimpleWeightedGraphs


function initialise_triangulated_necklace(ti::ComplexF64, tr::ComplexF64, f_hessian::Function;
    Ninit::Int64 = 20,
    ϵ::Float64 = 0.1)
    
    points_all = Vector{PointA{ComplexF64}}()
    
    hessian = f_hessian(ti, tr)

    # this could certainly be made more julian
    eigvec_mat = eigvecs(hessian)
    eigenvectors = [[complex(vec[1:2]...), complex(vec[3:4]...)] for vec in eachcol(eigvec_mat)] 

    pointsini = ([[ti,tr] .+ ϵ * (cos(θ) * eigenvectors[1] + sin(θ) * eigenvectors[2]) for θ in range(0, stop=2π, length=Ninit+1)])[1:end-1]

    push!(points_all, PointA(ti,tr)) ### add saddle point
    push!(points_all, [PointA(p[1], p[2]) for p in pointsini]...) ### add necklace points

    ### triangulating the initial necklace around the saddle point
    ### project the initial necklace onto its eigenvectors
        #### this is very stupid!!! I'm projecting back and forth.
    u1 = [vec[1] for vec in eachrow(eigvec_mat)]
    u2 = [vec[2] for vec in eachrow(eigvec_mat)]
    #         u1 = toR4(eigenvectors[1]...)
    #         u2 = toR4(eigenvectors[2]...)
   
    
    ### this is the same as in the other triangulation functions I hope
        v1 = toR4(ti,tr)
        # Projection function
        proj(v) = [(u1 ⋅ (v - v1)), (u2 ⋅ (v - v1))]
        # Project all points
        pts = [SVector(proj(toR4(p))...) for p in points_all]

        push!(pts, SVector{2, Float64}(zeros(2))) ### adding the centre point

        ### this is the very first triangles around the saddle point
        connections = delaunay(pts)
        triangles = [TriangleA(inds) for inds in eachcol(connections)]

    indices_necklace = collect(2:Ninit+1)
    
    return points_all, triangles, indices_necklace 
end


function flow_down!(points::Vector{PointA{T}},
        f::Function,
        f_grad::Function;
        threshold::Real=0.5, # for normalisation of the gradient
        δ::Real=0.5, # flowstepfactor
        h_threshold::Real=-20.
        ) where T<:Number

    for i1 in 1:length(points)
        if points[i1].active # for the active points
            step = -δ .* gradN((ti,tr) -> conj.(complex.(f_grad(ti,tr))), points[i1].x +0im, points[i1].y +0im, threshold)
            points[i1].x += step[1]
            points[i1].y += step[2]

            ### turning them inactive when they are below the threshold. Doing this here prevents lonely relict points from flowing if their triangle has turned inactive already.
            if real(f(xy(points[i1])...)) < h_threshold 
                points[i1].active = false
            end
        end
    end
end


function mesh_edges(tris::AbstractVector{TriangleA})
    edges = Set{Tuple{Int,Int}}()
    for tri in tris
        (a, b, c) = tri.coord
        for (i, j) in ((a,b), (b,c), (c,a))
            push!(edges, (min(i,j), max(i,j)))
        end
    end

    collect(edges)
end



function triangulate_flowed_points(points_all::Vector, first_brim::Vector, second_brim::Vector)
    N_fb = length(first_brim) ### this should be Ninit
    
    ### connecting new (flowed) and old points to triangles
    points_local = vcat(first_brim, second_brim)
    
    new_triangles_local = Vector{TriangleA}()
    extra_edges = Vector{Tuple{Int64,Int64}}()
    #### Locally within the brim(s) ###
    for i in 1:N_fb
        ### i is the local index within the brim
        fb_i = first_brim[i]
        sb_i = second_brim[i]
        
        i_next = i == N_fb ? 1 : i+1 ### index in the brim
        fb_n = first_brim[i_next]
        sb_n = second_brim[i_next]
        
        
        i2 = N_fb + i
        i3 = N_fb + i_next
        i4 = i_next            

        ### global indices
        globals = [only(findall(p-> isequal(xy(p), xy(points_local[i_local])), points_all)) for i_local in [i, i2, i3, i4] ]
     
        ### certainly I could use this "actives" here instead of the pairwise testing
        #         actives = [p.active for p in [fb_i, sb_i, fb_n, sb_n]]
        if isequal(fb_i, sb_i) && isequal(fb_n, sb_n) 
        # println("case 1, no one flows.")
            push!(extra_edges, (globals[1], globals[4]))
            continue
        elseif isequal(fb_i, sb_i) && !isequal(fb_n, sb_n)
        # println("case 2, i doesn't flow, defo take triangle choice 2.")
            push!(new_triangles_local, TriangleA([i, i3, i4]) )
        elseif !isequal(fb_i, sb_i) && isequal(fb_n, sb_n)
        # println("case 3, i flows, but i+1 doesn't, defo take triangle choice 1.")
            push!(new_triangles_local, TriangleA([i, i2, i4]) )
        else 
        # println("case 4, gotta compare the distances.")
            d13 = dist(fb_i, sb_n)
            d24 = dist(sb_i, fb_n)

            if d13 < d24
                push!(new_triangles_local, TriangleA([i, i3, i4]) )      
                push!(new_triangles_local, TriangleA([i, i2, i3]) )            
            elseif d24 <= d13
                push!(new_triangles_local, TriangleA([i, i2, i4]) )
                push!(new_triangles_local, TriangleA([i2, i3, i4]) )
            else
                println("Warning in the triangle making!")
            end
        end
    
    end

    new_triangles_global = Vector{TriangleA}()
    for tri in new_triangles_local
        coords_global = Vector{Int64}()
        for i_local in tri.coord        
            i_global = first(findall(p-> isequal(xy(p),xy(points_local[i_local])), points_all))
            push!(coords_global, i_global)            
        end
        push!(new_triangles_global, TriangleA(coords_global))
    end
    return new_triangles_global, extra_edges
end


function find_new_necklace_indices(points_all::Vector, second_brim::Vector, new_triangles::Vector, extra_edges::AbstractVector)
    
    edges_global = mesh_edges(new_triangles)
    # if two neighbouring points didn't flow, then there's no new triangle for them. Hence they won't appear in the "local points", but they are still part of the brim!
    ### points_local should also include all points from the second brim, just in case they're not part of the new triangles
    push!(edges_global, extra_edges...)
    
    Noffset = minimum(unique(collect(Iterators.flatten(edges_global)))) -1
    edges_global_unique = sort(unique(collect(Iterators.flatten(edges_global))))

    ### all the points that come up in the brim
    points_local = points_all[edges_global_unique]
    
    ### global 2 local
    d = Dict([val=>idx for (idx,val) in enumerate(edges_global_unique)])  
    edges_local = [(d[e[1]], d[e[2]]) for e in edges_global] 
    
    i_extra = 1
    for point in second_brim
        if isempty(findall(p-> isequal(xy(p), xy(point)), points_local))
            @show "I'm using this!!!!"
            i_extra += 1
            i_global = only(findall(p-> isequal(xy(p), xy(point)), points_all))
            push!(points_local, point)
            d[i_global] = length(deepcopy(points_local))
        end
    end

    
    ### within the brim indices
    target_indices_local = [only( findall(p-> isequal(xy(p), xy(point)), points_local)) for point in second_brim]

    full_path = Int[]
    total_length = 0.0

    ### within brim indices
    weights= [dist(points_local[i], points_local[j]) for (i,j) in edges_local]
    
    g = SimpleWeightedGraph(length(points_local))
    for (k, (i, j)) in enumerate(edges_local)
        add_edge!(g, i, j, weights[k])
    end

    for (x1, x2) in zip(target_indices_local[1:end-1], target_indices_local[2:end])
        pathstate = dijkstra_shortest_paths(g, x1)
        path = enumerate_paths(pathstate, x2)
        append!(full_path, path[1:end-1])   # avoid repeating vertices
        total_length += pathstate.dists[x2]
    end

    ### close the loop
    pathstate_closing = dijkstra_shortest_paths(g, target_indices_local[end])
    path_closing = enumerate_paths(pathstate_closing, target_indices_local[1])
    append!(full_path, path_closing[1:end-1])   # avoid repeating vertices
    total_length += pathstate_closing.dists[target_indices_local[1]]

    new_necklace_indices_global = Int64[]
    for i in full_path
        el = findall(p->xy(p) == xy(points_local[i]), points_all)
        if length(el) ==1
            push!(new_necklace_indices_global, only(el))
        else
            @show el
        end
    end
    
    return new_necklace_indices_global
end



function get_SD_thimble_triangles(
    f::Function,
    f_grad::Function,
    f_hessian::Function,
    ti::ComplexF64, tr::ComplexF64
    ; Ninit::Int64=20, Nflow::Int64=10,
#     accuracy::Float64=1e-4,
    eigvecfactorinit::Float64 = 0.1, # I should come up with sophisticated guesses here.
    flowstepfactor::Float64 = 6., 
    subdividethreshold::Float64 = 5.,
    gradn_threshold::Real=1.,
    h_threshold::Real=-50.) ### this makes the flow stop at some point
    
    points_all, triangles, indices_necklace = initialise_triangulated_necklace(ti, tr, f_hessian, Ninit=Ninit, ϵ=eigvecfactorinit)
    
    for i_flow in 1:Nflow 
        ### doing one flow step to create a loop brim ###
        points_necklace = deepcopy(points_all)[indices_necklace]
        first_brim = deepcopy(points_necklace)
        flow_down!(points_necklace, f, f_grad, threshold = gradn_threshold, δ = flowstepfactor, h_threshold = h_threshold)
        second_brim = deepcopy(points_necklace)
            @assert length(second_brim) == length(first_brim)
        
        # [push!(points_all, deepcopy(p)) for p in points_necklace if isempty(findall(px -> xy(px)== xy(p),points_all ))]
        union!(points_all, points_necklace)

            @assert length(points_all) == length(unique(points_all))
        ### triangulate this brim
        new_triangles, extra_edges = triangulate_flowed_points(points_all, first_brim, second_brim)
        
        ### subdivide the new triangles
        subdivide_triangles!(points_all, new_triangles, subdividethreshold)
        
        ### update the outer brim
        indices_necklace = find_new_necklace_indices(points_all, second_brim, new_triangles, extra_edges)

        push!(triangles, deepcopy(new_triangles)...)
    end     
    
    return points_all, triangles, indices_necklace
    
end










nothing