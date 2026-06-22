### QUAD (stores indices only)
mutable struct QuadA
    indices::Vector{Int} ### this could be an MVector
    active::Bool
    QuadA(indices::Vector{Int64}) = new(indices, true)
end


function maxdist(simplex::Vector{PointA{T}}) where T<:Number
        p1,p2,p3,p4 = simplex
        d12 = dist(p1,p2) 
        d23 = dist(p2,p3)
        d34 = dist(p3,p4)
        d41 = dist(p4,p1)

        arr = [d12,d23,d34,d41]
    return maximum(arr)
end


function maxdist(idx::Int64, simplices::Vector{QuadA}, points::Vector{<:PointA})
    sim = simplices[idx]

    v1, v2, v3, v4 = sim.indices
    p1, p2, p3, p4 = points[v1], points[v2], points[v3], points[v4]
    return maxdist([p1,p2,p3,p4])
end

function subdivide_quads!(points::Vector{<:PointA}, simplices::Vector{QuadA}, Δ::Float64, init::Bool=false) 

    for i3 in eachindex(simplices)
        sim = simplices[i3]
        if sim.active
            v1, v2, v3, v4 = sim.indices
            p1, p2, p3, p4 = points[v1], points[v2], points[v3], points[v4]

            if (p1.active && p2.active && p3.active && p4.active) || init # == allpointsactive
                l = length(points)

                d12 = norm([p1.x-p2.x, p1.y-p2.y])
                d23 = norm([p2.x-p3.x, p2.y-p3.y])
                d34 = norm([p3.x-p4.x, p3.y-p4.y])
                d41 = norm([p4.x-p1.x, p4.y-p1.y])
            
                arr = [d12,d23,d34,d41]
                
                printcond = false
                if maximum(arr) > Δ # subdivide
                    simplex = [p1,p2,p3,p4]
                    keys = [[1,2],[2,3],[3,4],[4,1]]
                    sorted_indices = sortperm(arr, rev=true)
                    (max1_idx, max2_idx) = sorted_indices[1:2]

                    # Get the values
                    max1 = arr[max1_idx]
                    max2 = arr[max2_idx] 

                    new1 = midpoint(simplex[keys[max1_idx]]...) # new point 1
                    new2 = midpoint(simplex[keys[max2_idx]]...) # new point 2

                    already_created_1 = findall(p->isequal(p, new1), points)
                    already_created_2 = findall(p->isequal(p, new2), points)

                    if isempty(already_created_1)
                        push!(points, new1)
                        new1_idx = l+1
                    else
                        new1_idx = already_created_1[1]
                    end
                    
                    if isempty(already_created_2)
                        push!(points, new2)
                        new2_idx = l+1+ isempty(already_created_1)
                    else
                        new2_idx = already_created_2[1]
                    end

                    simplices[i3].active = false #############

                    push!(keys, [keys[max1_idx][1],5],[5,keys[max1_idx][2]] )
                    push!(keys, [keys[max2_idx][1],6],[6,keys[max2_idx][2]] )
                    deleteat!(keys,sort([max1_idx, max2_idx]) )
                    push!(simplex, [new1, new2]...)
                    
                    ### sort the keys #######
                    keys_dict = Dict()
                    for k in keys
                        keys_dict[k[1]] = k
                    end
                    # Initialize the sorted line segments vector with the first line segment
                    sorted_keys = [keys[1]]
                    # Iterate until all line segments are sorted
                    while length(sorted_keys) != length(keys)
                        # Find the next line segment based on the end index of the last sorted segment
                        matching_ls = keys_dict[last(sorted_keys)[2]]
                        push!(sorted_keys, matching_ls)
                    end
                    #############

                    sorted_points = [k[1] for k in sorted_keys]
                    shortest_diag, shortest_diag_idx = Inf, [0,0]
                    for kdiag in [[1,4],[2,5],[3,6]]
                        diag = dist(simplex[[sorted_points[kdiag]...]]...)
                        if diag < shortest_diag
                            shortest_diag = diag
                            shortest_diag_idx = kdiag
                        end
                    end
                    new_sim1_idx = sorted_points[mod1.(shortest_diag_idx[1]:(shortest_diag_idx[1]+3), length(sorted_points))]
                    new_sim2_idx = sorted_points[mod1.(shortest_diag_idx[2]:(shortest_diag_idx[2]+3), length(sorted_points))]

                    new_sim1 = simplex[new_sim1_idx]
                    new_sim2 = simplex[new_sim2_idx]
                    
                    ### now translate back to the actual indices
                    proper_indices = [v1,v2,v3,v4, new1_idx, new2_idx]
                    append!(simplices, [
                        QuadA(proper_indices[new_sim1_idx]),
                        QuadA(proper_indices[new_sim2_idx])])
                end
            
            elseif !any([p1.active, p2.active, p3.active, p4.active])
                simplices[i3].active = false
            end
        end
    end
#     return simplices, points
end



function subdivide(points::Vector{<:PointA}, simplices::Vector{QuadA}, Δ::Float64)
    n_old = length(simplices)
    n_new = n_old + 1
    while (n_old != n_new)
        n_old = n_new
        subdivide_quads!(points, simplices, Δ)
        filter!(sim->sim.active, simplices)
        n_new = length(simplices)
    end 
end


function initialise_grid_quads(points::Vector{<:PointA}, Δ::Float64)
    simplices = [QuadA([1,2,3,4])]
    
    ### subdivide_2(points, simplices, Δ) # instead of calling this I'll do it here directly
    n_old = length(simplices)
    n_new = 0 #n_old - 1
    
    while (n_old != n_new)
        n_old = n_new
        subdivide_quads!(points, simplices, Δ, n_new == 0)
        filter!(sim->sim.active, simplices)
        n_new = length(simplices)
    end    
    return (points, simplices)
end


### flowing a whole meshed surface
function flow_down!(simplices::Vector{QuadA},points::Vector{<:PointA},
        f::Function,
        f_grad::Function;
        threshold::Float64=0.5, # for normalisation of the gradient
        δ::Float64=0.5, # flowstepfactor
        h_threshold::Float64=-20.
        )

    for i1 in 1:length(points)
        if points[i1].active # for the active points
            step = -δ .* gradN((ti,tr) -> conj.(complex.(f_grad(ti,tr))), points[i1].x +0im, points[i1].y +0im, threshold)
            points[i1].x += step[1]
            points[i1].y += step[2]
        end
    end

    for i2 in eachindex(simplices)
        if simplices[i2].active
            for v in simplices[i2].indices
                if real(f(points[v].x, points[v].y)) < h_threshold 
                    simplices[i2].active = false # am I sure that I want to turn the whole simplex inactive?
                    points[v].active = false
                end
            end
        end
    end
end

### convergence of number of simplices
function has_converged(history::Vector{T}; tol::Float64=1., window_size::Int=10) where T <:Number
    # Calculate the average of the last `window_size` vectors
    if length(history) < window_size
        return false  # Not enough data points
    end
    avg_oscill = sum(history[end-window_size+1:end]) / window_size    
    return abs(avg_oscill) < tol
end


### not sure if this is needed
# isequal(q1::QuadC,q2::QuadC) = all([isequal(q1.points[i],q2.points[i]) for i in 1:4])
# Base.hash(q::QuadC, h::UInt) = hash([(p.x, p.y) for p in q.points], h)


### get simplices
function get_flowed_quads(
    f::Function,
    f_grad::Function,
    init_points::Vector{<:PointA};
    Nflow::Int64=50,
    Δinit::Float64 = 10.,
    gradnthreshold::Float64 = 0.05, # grad normalisation threshold
    flowstepfactor::Float64 = 2., # flowstepfactor
    subdividethreshold::Float64 = 8., # subdivide threshold, wants to be 4 * δ
    h_threshold::Float64 = -150.,
    maxNsimplices::Int64=5000,
    tolNsimplices::Float64=0.05,
    flow_bounds::Vector{Bool}=[true,true,true,true]
    )

    netsimplices = Vector{Int64}()
    (points, simplices) = initialise_grid_quads(init_points, Δinit)

    overboard = false
    push!(netsimplices, length(simplices))
    for i_flow in 1:Nflow
        nsimplices = length(simplices)

        flow_down!(simplices, points, f, f_grad,
                threshold = gradnthreshold, δ=flowstepfactor, h_threshold = h_threshold)
        subdivide(points, simplices, subdividethreshold)

        ### several termination conditions below
        net = (length(simplices)-nsimplices)
        if net !== 0
            push!(netsimplices, net)
            i_flow > 1 ? overboard = true : false
        # elseif (net) == 0 && overboard
        #     println("I naturally converged at flow step $i_flow !") ; break
        end

        if has_converged(netsimplices[2:end], tol = tolNsimplices*netsimplices[1]) 
            println("I converged w.r.t. number of simplices at i_flow = $i_flow."); 
            break 
        end

        if isempty(findall(sim->sim.active, simplices))
            println("I broke because all things finished flowing after $i_flow steps."); break
        end
        if length(simplices) > maxNsimplices
            println("I broke after $i_flow steps because I have more than $maxNsimplices simplices now."); break
        end        

        if i_flow == Nflow
            println("I stopped because I reached the maximum flow steps, i.e. $Nflow.")
        end
    end

    quads =  [QuadC(points[sim.indices]) for sim in simplices]

    return quads, points, simplices
end


function integrate_flowed_quads(
    f::Function,
    f_grad::Function,
    timin::Number, timax::Number,
    ttmin::Number, ttmax::Number;
    prefactor::Function = (ti,tr) -> ones(2),
    Nflow::Int64=50,
    Δinit::Float64 = 10.,
    gradnthreshold::Float64 = 0.5, # grad normalisation threshold
    flowstepfactor::Float64 = 2., # flowstepfactor
    subdividethreshold::Float64 = 8., # subdivide threshold, wants to be 4 * δ
    h_threshold::Float64 = -150.,
    maxNsimplices::Int64=5000,
    integral_accuracy::Float64=1e-7,
    integral_rel_error::Float64=0.05,
    print_message::Bool=true
    )

    netsimplices = Vector{Int64}()
    (points, simplices) = initialise_grid_parallelogram(complex(timin), complex(timax), complex(ttmin), complex(ttmax), Δinit)
    overboard = false
    prev_integral = complex(ones(2))
    int = complex(zeros(2))


    for i_flow in 1:Nflow
        nsimplices = length(simplices)
        # @show simplices

        flow_down!(simplices, points, f, f_grad,
                threshold = gradnthreshold, δ=flowstepfactor, h_threshold = h_threshold)
        # @show simplices
        subdivide(points, simplices, subdividethreshold)
        # @show simplices
        quads =  [QuadC(points[sim.indices]) for sim in simplices]
        int = complex(zeros(2))
        for quad in quads
            int += integrate_quadrilateral(f, quad, prefactor = prefactor)
        end

        abs_diff = norm(int .- prev_integral) 
        if 0 < abs_diff < integral_accuracy # if I don't want to stop here I can just set the goal incredibly low
            if print_message
            println("I broke after $i_flow iterations because the accuracy goal of $integral_accuracy was met, using $(length(simplices)) simplices.")
            end
            break
        end
        
        rel_diff = norm(filter!(!isnan,(int .- prev_integral) ./ prev_integral))
        if 0 < rel_diff < integral_rel_error
            if print_message
                println("I broke after $i_flow iterations because the accuracy goal of $integral_rel_error relative error of the integral was met, using $(length(simplices)) simplices.")
            end
            break
        end        

        prev_integral = int

#         if has_converged(netsimplices[2:end], tol = 0.05*netsimplices[1]) 
#             println("The number of simplices converged at i_flow = $i_flow."); 
#             break 
#         end


        if isempty(findall(sim->sim.active, simplices))
            println("I broke because I ran out of simplices after $i_flow steps. This might be solved by using a lower h_threshold."); break
        end
        
        
        if length(simplices) > maxNsimplices 
            print_message ? println("I broke after $i_flow steps because I have more than $maxNsimplices simplices now.") : nothing
            break
        end        

        if i_flow == Nflow && print_message
            println("I stopped because I reached the maximum flow steps, i.e. $Nflow, with $(length(simplices)) simplices.")            
        end
    end

    return int, length(simplices)
end


function integrate_flowed_quads_fixed_Nflow(
    f::Function,
    f_grad::Function,
    init_points::Vector{<:PointA};
    prefactor::Function = (ti,tr) -> ones(2),
    Nflow::Int64=50,
    Δinit::Float64 = 10.,
    gradnthreshold::Float64 = 0.05, # grad normalisation threshold
    flowstepfactor::Float64 = 2., # flowstepfactor
    subdividethreshold::Float64 = 8., # subdivide threshold, wants to be 4 * δ
    h_threshold::Float64 = -150.,
    maxNsimplices::Int64=5000,
    print_message::Bool=true
    )

    netsimplices = Vector{Int64}()
    (points, simplices) = initialise_grid(init_points, Δinit)
    int = complex(zeros(2))

    for i_flow in 1:Nflow
        nsimplices = length(simplices)

        flow_down!(simplices, points, f, f_grad,
                threshold = gradnthreshold, δ=flowstepfactor, h_threshold = h_threshold)
        subdivide(points, simplices, subdividethreshold)
        quads =  [QuadC(points[sim.indices]) for sim in simplices]
        int = complex(zeros(2))
        for quad in quads
            int += integrate_quadrilateral(f, quad, prefactor = prefactor)
        end

        if isempty(findall(sim->sim.active, simplices))
            println("I broke because I ran out of simplices after $i_flow steps. This might be solved by using a lower h_threshold."); break
        end
        
        if length(simplices) > maxNsimplices 
            print_message ? println("I broke after $i_flow steps because I have more than $maxNsimplices simplices now.") : nothing
            break
        end        

    end

    return int, length(simplices)
end

nothing