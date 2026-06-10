#### making SD thimbles
using Graphs, SimpleWeightedGraphs


function initialise_triangulated_necklace(ti::ComplexF64, tr::ComplexF64, f_hessian::Function;
    Ninit::Int64 = 20,
    ϵ::Float64 = 0.1)
    
    points_all = Vector{PointA}()
    
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
    
    ### this is the same as in the other triangulation functions I hope
        v1 = toR4(ti,tr)
        # Projection function
        proj(v) = [(u1 ⋅ (v - v1)), (u2 ⋅ (v - v1))]
        # Project all points
        pts = [SVector(proj(toR4(p))...) for p in points_all]

        push!(pts, SVector{2, Float64}(zeros(2))) ### adding the centre point

        ### this is the very first triangles around the saddle point
        connections = delaunay(pts)

        ### ensure they are all oriented in the same direction
        n_ref = triangle_normal(points_all[eachcol(connections)[1]]...) # randomly choose the first triangle as a reference triangle TODO could be a smarter choice
        triangles = TriangleA[]
        for inds in eachcol(connections)
            n = triangle_normal(points_all[inds]...)
            if real(dot(n, n_ref)) > 0 
                push!(triangles, TriangleA(inds))
            else
                push!(triangles, TriangleA(reverse(inds)))
            end
        end

    indices_necklace = collect(2:Ninit+1)
    
    return points_all, triangles, indices_necklace 
end


### this flows a list of points
function flow_down!(points::Vector{<:PointA},
        f::Function,
        f_grad::Function;
        threshold::Real=0.5, # for normalisation of the gradient
        δ::Real=0.5, # flowstepfactor
        h_threshold::Real=-20.
        )

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
        (a, b, c) = tri.indices
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
        indices_global = Vector{Int64}()
        for i_local in tri.indices        
            i_global = first(findall(p-> isequal(xy(p),xy(points_local[i_local])), points_all))
            push!(indices_global, i_global)            
        end
        push!(new_triangles_global, TriangleA(indices_global))
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
    eigvecfactorinit::Float64 = 0.1, # I should come up with sophisticated guesses here.
    flowstepfactor::Float64 = 6., 
    subdividethreshold::Float64 = 5.,
    gradn_threshold::Real=1.,
    h_threshold::Real=-50.) ### this makes the flow stop at some point
    
    points_all, triangles, indices_necklace = initialise_triangulated_necklace(ti, tr, f_hessian, Ninit=Ninit, ϵ=eigvecfactorinit)

    ### find a suitable threshold for the normalisation of the gradient
    gradient0 = [norm(conj.(f_grad(p.x, p.y))) for p in points_all]
    threshold = round(minimum(gradient0), RoundDown, sigdigits=2)
    @show threshold

    
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
        new_triangles_unoriented, extra_edges = triangulate_flowed_points(points_all, first_brim, second_brim)
        
        ### make sure new triangles are well oriented w.r.t. a previous one
        n_ref = triangle_normal(points_all[triangles[end].indices]...) # choose the last triangle as a reference TODO this isn't necessarily a good choice
        new_triangles = TriangleA[]
        for inds in new_triangles_unoriented
            n = triangle_normal(points_all[inds.indices]...)

            if real(dot(n, n_ref)) > 0 
                push!(new_triangles, TriangleA(inds.indices, inds.active))
            else
                push!(new_triangles, TriangleA(reverse(inds.indices), inds.active))
            end
        end

        ### subdivide the new triangles
        subdivide_triangles!(points_all, new_triangles, subdividethreshold)
        
        ### update the outer brim
        indices_necklace = find_new_necklace_indices(points_all, second_brim, new_triangles, extra_edges)

        push!(triangles, deepcopy(new_triangles)...)
    end     
    trianglesC = [TriangleC(points_all[tri.indices]) for tri in triangles]

    return trianglesC, points_all, triangles
    
end
