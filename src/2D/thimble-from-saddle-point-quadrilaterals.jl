### I think this needs some functions from the necklace file

function initialise_SD!(necklace::Vector{LineSeg}, points::Vector{<:PointA},
        ti::ComplexF64, tr::ComplexF64;
        f_hessian::Function,
        Ninit::Int64 = 20,
        ϵ::Float64 = 0.01)

    # hessian = my_hessian(b,Ip,q,ti,tr)
    hessian = f_hessian(ti, tr)

    # this could certainly be made more julian    
    eigenvectors = [[complex(vec[1:2]...), complex(vec[3:4]...)] for vec in eachcol(eigvecs(hessian))] 

    pointsini = ([[ti,tr] .+ ϵ * (cos(θ) * eigenvectors[1] + sin(θ) * eigenvectors[2]) for θ in range(0, stop=2π, length=Ninit+1)])[1:end-1]

    push!(points, [PointA(p[1], p[2]) for p in pointsini]...)
    push!(necklace, [LineSeg(i, i+1 , true) for i in 1:(length(points)-1)]...)
    push!(necklace, LineSeg(length(points), 1 , true)) # closing the necklace
end

function adorn_necklace(necklace::Vector{LineSeg}, points::Vector{<:PointA})
    for i in 1:length(necklace)
        necklace[i] = LineSeg(points[necklace[i].sindex], points[necklace[i].eindex])
    end
    return necklace
end;

function flow_down!(necklace::Vector{LineSeg}, points::Vector{<:PointA},
        f::Function,
        f_grad::Function;
        δ::Float64=0.1,
        threshold::Float64=0.5,
        h_threshold::Float64=-200.
        )

    for i in 1:length(points)
        # TODO check both real and imaginary part?
        if points[i].active # for the active points
            # set them to be active (= still flowing) if they are above threshold
            points[i].active = real(f(points[i].x, points[i].y)) > h_threshold #(in Job's code that's h-function > thresh, I should clearly state which sign I'm using where etc.) 
            if points[i].active
                # step = -δ .* gradN(b, Ip, q, points[i].x, points[i].y, threshold)
                step = -δ .* gradN((ti,tr)->conj.(complex.(f_grad(ti,tr))),
                 points[i].x, points[i].y, threshold)

                points[i].x += step[1]
                points[i].y += step[2]
            end
        end
    end
end


### not sure what I will need this for ever
function get_necklace_SD_solver_with_traces(
    f::Function,
    f_grad::Function,
    f_hessian::Function,
    ti::ComplexF64, tr::ComplexF64
    ; Ninit::Int64=20, Ncounter::Int64=500,
    eigvecfactorinit::Float64 = 0.02, # I should come up with sophisticated guesses here.
    flowstepfactor::Float64 = 0.1, 
    subdividethreshold::Float64 = 0.5 )
    
    necklace = Vector{LineSeg}()
    points = Vector{<:PointA}()    
    
    initialise_SD!(necklace, points, ti, tr, f_hessian=f_hessian, Ninit = Ninit, ϵ = eigvecfactorinit)
    
    points_traces = Vector{Vector{<:PointA}}()
    for i in 1:length(points)
        push!(points_traces, Vector{<:PointA}())
    end
    
    necklaces = Vector{Vector{LineSeg}}()
    push!(necklaces, adorn_necklace(sort_linesegs(necklace), points))
    
    ### find a suitable threshold for the normalisation of the gradient
    gradient0 = [norm(conj.(f_grad(p.x, p.y))) for p in points]
    threshold = round(minimum(gradient0), RoundDown, sigdigits=2)
    counter = 0
        
    while counter < Ncounter
        counter += 1
            
        for i in 1:length(points)
            if i <= length(points_traces)
                push!(points_traces[i], deepcopy(points[i]))
            else
                vcat(points_traces, [deepcopy(points[i])])
            end
        end
        push!(necklaces, deepcopy(adorn_necklace(sort_linesegs(necklace), points)))       
        
        flow_down!(necklace, points, f, f_grad, threshold = threshold, δ = flowstepfactor)
       
        if count([p.active for p in points]) == 0
    #             @debug "I broke because the flow stopped after $counter iterations"
            println("I broke because the flow stopped after $counter iterations")
            break
        end
        for i in 1:length(necklace)
            subdivide!(necklace[i], necklace, points, Δ = subdividethreshold)
        end
        keepat!(necklace, [ls.active for ls in necklace])
        
        ##
        
    end    
    if counter == Ncounter && Ncounter > 1
        println("I broke because the counter reached its max, i.e. $Ncounter")
    end
    
    necklace = sort_linesegs(necklace)
    adorn_necklace!(necklace, points)
    
    return necklace, necklaces, points_traces
end

### l121 - 213 could be deleted once I finished this refurbishing experiment
# function make_quads(necklace::Vector{LineSeg},
#         points::Vector{<:PointA}, prev_necklace::Vector{LineSeg} )

#     quads = Vector{Tuple}()
#     new_necklace = adorn_necklace(sort_linesegs(necklace), points)  
#     for idx in 1:length(prev_necklace)
#         quad = (prev_necklace[idx], new_necklace[idx])
#         push!(quads, quad)
#     end

#     return quads
# end


# function get_SD_thimble_quadrangles(
#     f::Function,
#     f_grad::Function,
#     f_hessian::Function,
#     ti::ComplexF64, tr::ComplexF64
#     ; Ninit::Int64=20, Ncounter::Int64=500,
#     accuracy::Float64=1e-4,
#     eigvecfactorinit::Float64 = 0.02, # I should come up with sophisticated guesses here.
#     flowstepfactor::Float64 = 10., 
#     subdividethreshold::Float64 = 20.)
    
#     ### check that Ninit ganzzahlig
#     necklace = Vector{LineSeg}()
#     points = Vector{<:PointA}()
    
#     initialise_SD!(necklace, points, ti, tr, f_hessian=f_hessian, Ninit = Ninit, ϵ = eigvecfactorinit)
    
#         #     points_traces = Vector{Vector{<:PointA}}()
#         #     for i in 1:length(points)
#         #         push!(points_traces, Vector{<:PointA}())
#         #     end
    
#     necklaces = Vector{Vector{LineSeg}}()
#     push!(necklaces, adorn_necklace(sort_linesegs(necklace), points))
    
#     ### find a suitable threshold for the normalisation of the gradient
#     gradient0 = [norm(conj.(f_grad(p.x, p.y))) for p in points]
#     threshold = round(minimum(gradient0), RoundDown, sigdigits=2)
    
#     counter = 0
#     quadrangles = Vector{Tuple}()
    
#     while counter < Ncounter
#         counter += 1
#         push!(necklaces, deepcopy(adorn_necklace(sort_linesegs(necklace), points)))       
        
#         flow_down!(necklace, points, f, f_grad, threshold = threshold, δ = flowstepfactor)
        
#         if count([p.active for p in points]) == 0
#             @debug "I broke because the flow stopped after $counter iterations"
#             break
#         end
#         prev_necklace = necklaces[end]
#         new_quads = make_quads(necklace, points, prev_necklace)
#         push!(quadrangles, deepcopy(new_quads)...)
        
#         for i in 1:length(necklace)
#             subdivide!(necklace[i], necklace, points, Δ = subdividethreshold)
#         end
#         keepat!(necklace, [ls.active for ls in necklace])
       
#     end
#     if counter == Ncounter && Ncounter > 1
#         println("I broke because the counter reached its max, i.e. $Ncounter")
#     end
    
#     necklace = sort_linesegs(necklace)
#     adorn_necklace!(necklace, points)
    
#     return quadrangles
# #     return necklace, necklaces, points_traces
# end


### why the fuck do i have a quadrangle and a quadrilateral ???
# function integrate_quadrangle(
#     f::Function,
#     quadrangle::Tuple{LineSeg,LineSeg}, n::Int64=7;
#     prefactor::Function=(ti,tr) -> ones(2)
#     )
    
#     p1 = quadrangle[1].s
#     p2 = quadrangle[2].s
#     p3 = quadrangle[2].e
#     p4 = quadrangle[1].e
    
#     integrate_quadrilateral(f, QuadC([p1,p2,p3,p4]), n, prefactor=prefactor)
# end


function make_quad(ls1::LineSeg, ls2::LineSeg)

    p1 = ls1.s
    p2 = ls2.s
    p3 = ls2.e
    p4 = ls1.e

    return QuadC([p1,p2,p3,p4])
end

function make_quads(necklace::Vector{LineSeg},
        points::Vector{<:PointA}, prev_necklace::Vector{LineSeg} )

    quads = Vector{QuadC}()
    new_necklace = adorn_necklace(sort_linesegs(necklace), points)  
    for idx in 1:length(prev_necklace)
        quad = make_quad(prev_necklace[idx], new_necklace[idx])
        push!(quads, quad)
    end

    return quads
end



function get_SD_thimble_quads(
    f::Function,
    f_grad::Function,
    f_hessian::Function,
    ti::ComplexF64, tr::ComplexF64
    ; Ninit::Int64=20, Ncounter::Int64=500,
    accuracy::Float64=1e-4,
    eigvecfactorinit::Float64 = 0.02, # I should come up with sophisticated guesses here.
    flowstepfactor::Float64 = 10., 
    subdividethreshold::Float64 = 20.)
    
    ### check that Ninit ganzzahlig
    necklace = Vector{LineSeg}()
    points = Vector{<:PointA}()
    
    initialise_SD!(necklace, points, ti, tr, f_hessian=f_hessian, Ninit = Ninit, ϵ = eigvecfactorinit)
    
        #     points_traces = Vector{Vector{<:PointA}}()
        #     for i in 1:length(points)
        #         push!(points_traces, Vector{<:PointA}())
        #     end
    
    necklaces = Vector{Vector{LineSeg}}()
    push!(necklaces, adorn_necklace(sort_linesegs(necklace), points))
    
    ### find a suitable threshold for the normalisation of the gradient
    gradient0 = [norm(conj.(f_grad(p.x, p.y))) for p in points]
    threshold = round(minimum(gradient0), RoundDown, sigdigits=2)
    
    counter = 0
    quads = Vector{QuadC}()
    
    while counter < Ncounter
        counter += 1
        push!(necklaces, deepcopy(adorn_necklace(sort_linesegs(necklace), points)))       
        
        flow_down!(necklace, points, f, f_grad, threshold = threshold, δ = flowstepfactor)
        
        if count([p.active for p in points]) == 0
            @debug "I broke because the flow stopped after $counter iterations"
            break
        end
        prev_necklace = necklaces[end]
        new_quads = make_quads(necklace, points, prev_necklace)
        push!(quads, deepcopy(new_quads)...)
        
        for i in 1:length(necklace)
            subdivide!(necklace[i], necklace, points, Δ = subdividethreshold)
        end
        keepat!(necklace, [ls.active for ls in necklace])
       
    end
    if counter == Ncounter && Ncounter > 1
        println("I broke because the counter reached its max, i.e. $Ncounter")
    end
    
    necklace = sort_linesegs(necklace)
    adorn_necklace!(necklace, points)
    
    return quads
#     return necklace, necklaces, points_traces
end




function integrate_SD_thimble(
    f::Function,
    f_grad::Function,
    f_hessian::Function,
    ti::ComplexF64, tr::ComplexF64;
    prefactor::Function=(ti,tr) -> ones(2),
    Ninit::Int64=20, Ncounter::Int64=500,
    accuracy::Float64=1e-4,
    eigvecfactorinit::Float64 = 0.02, # I should come up with sophisticated guesses here.
    flowstepfactor::Float64 = 4., 
    subdividethreshold::Float64 = 1.)
    
    function integrate_annulus(necklace::Vector{LineSeg},
            points::Vector{<:PointA}, prev_necklace::Vector{LineSeg})
        
        new_necklace = adorn_necklace(sort_linesegs(necklace), points)  
        int = zeros(ComplexF64, 2)
        for idx in 1:length(prev_necklace)
            quad = make_quad(prev_necklace[idx], new_necklace[idx])

            int .+= integrate_quadrangle(f, quad, prefactor=prefactor)
        end
        
        return int
    end
        

    ### TODO: check that Ninit is even
    necklace = Vector{LineSeg}()
    points = Vector{<:PointA}()
    
    initialise_SD!(necklace, points, ti, tr, f_hessian = f_hessian, Ninit = Ninit, ϵ = eigvecfactorinit)
    
    ### TODO find a suitable threshold for the normalisation of the gradient
    gradient0 = [norm(conj.(f_grad(p.x, p.y))) for p in points]
    threshold = round(minimum(gradient0), RoundDown, sigdigits=2)

    total_integral = zeros(ComplexF64,2)
    
    counter = 0

    ### sign of the intersection number
    sign_in = sign(imag(hessian_root(f_hessian(ti, tr))))

    while counter < Ncounter
        counter += 1
            
        prev_necklace = deepcopy(adorn_necklace(sort_linesegs(necklace), points))   
        flow_down!(necklace, points, f, f_grad, threshold = threshold, δ = flowstepfactor)
        
        if count([p.active for p in points]) == 0
            @debug "I broke because the flow stopped after $counter iterations"
            break
        end
        
        
        int = sign_in * integrate_annulus(necklace, points, prev_necklace)
        total_integral .+= int
        

        ### these are the convergence criteria
        # abs_diff = norm(int .- prev_integral) 
        # if 0 < abs_diff < integral_accuracy 

        # rel_diff = norm(filter!(!isnan, (int .- prev_integral) ./ prev_integral))
        # if 0 < rel_diff < accuracy


        if 0 < norm(reim(int)./reim(total_integral)) < accuracy
    #             println("I broke after $counter iterations because the accuracy goal of $accuracy was met.")
            break
        end         
        
        for i in 1:length(necklace)
            subdivide!(necklace[i], necklace, points, Δ = subdividethreshold)
        end
        keepat!(necklace, [ls.active for ls in necklace])
       
    end
    
    if counter == Ncounter && Ncounter > 1
        println("I broke because the SD thimble counter reached its max, i.e. $Ncounter.")
    end
    
    return total_integral
end


nothing