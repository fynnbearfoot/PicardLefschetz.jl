### point and lineseg
### point is now in the downward flow file. And is updated to PointA. Index => Quad
	mutable struct LineSeg
	    s::Union{UndefInitializer,PointA{T}} where T<:Union{<:Number, Vector{<:Number}}
	    e::Union{UndefInitializer,PointA{T}} where T<:Union{<:Number, Vector{<:Number}}
	    active::Bool
	    sindex::Union{UndefInitializer, Int}
	    eindex::Union{UndefInitializer, Int}
	end
	LineSeg(s::PointA, e::PointA) = LineSeg(s,e, s.active && e.active, undef,undef);
	LineSeg(s::PointA, e::PointA, active::Bool) = LineSeg(s, e, active, undef,undef);
	LineSeg(sindex::Int, eindex::Int, active::Bool) = LineSeg(undef, undef, active, sindex, eindex);
	LineSeg(sindex::Int, eindex::Int) = LineSeg(undef, undef, false, sindex, eindex);

    import Base.imag, Base.real
    imag(p::PointA) = PointA(imag.(p.x),imag.(p.y), p.active)
    real(p::PointA) = PointA(real.(p.x),real.(p.y), p.active);

    imag(ls::LineSeg) = LineSeg(imag(ls.s),imag(ls.e), ls.active)
    real(ls::LineSeg) = LineSeg(real(ls.s),real(ls.e), ls.active);

    import LinearAlgebra.norm
    function norm(ls::LineSeg)
        norm([ls.e.x, ls.e.y] .- [ls.s.x, ls.s.y] )
    end
    dist(p1::PointA, p2::PointA) = norm(@SVector[p2.x-p1.x,p2.y-p1.y])
    import Base.length 
    length(ls::LineSeg) = dist(ls.e, ls.s)


    function get_point(ls::LineSeg, which::Symbol=:s)
        if which==:s
            return PointA(ls.s.x, ls.s.y)
        elseif which == :e
            return PointA(ls.e.x, ls.e.y)
        else
            return println("You've got a problem!")
        end
    end

    ### this sorts linesegments such that they form a connected line
    function sort_linesegs(linesegs::Vector{LineSeg})
        linesegs_dict = Dict()
        for ls in linesegs
            linesegs_dict[ls.sindex] = ls
        end
        
        # Initialize the sorted line segments vector with the first line segment
        sorted_linesegs = [linesegs[1]]
        
        # Iterate until all line segments are sorted
        while length(sorted_linesegs) != length(linesegs)
            # Find the next line segment based on the end index of the last sorted segment
            matching_ls = linesegs_dict[last(sorted_linesegs).eindex]
            push!(sorted_linesegs, matching_ls)
        end
        
        return sorted_linesegs
    end



### simple Gauss area formula to find the area enclosed by the necklace (to double-check if it's not got folded into itself)
function enclosed_area(linesegs::Vector{LineSeg}, f::Function = x -> real(x))
    # Initialize the area accumulator
    area = 0.0
    
    # Iterate over each line segment
    for i in 1:length(linesegs)
        # Get the coordinates of the endpoints of the line segment
        x1 = f(linesegs[i].s.x)
        y1 = f(linesegs[i].s.y)
        x2 = f(linesegs[i].e.x)
        y2 = f(linesegs[i].e.y)        
        
        # Update the area accumulator
        area += x1*y2 - x2*y1
    end
    
    # Divide the result by 2 to get the absolute area
    area = abs(area) / 2.0
    
    return area
end


### necklacy things
function initialise!(necklace::Vector{LineSeg}, points::Vector{PointA},
        ti::ComplexF64, tr::ComplexF64,
        f::Function;
        Ninit::Int64 = 20,
        ϵ::Float64 = 0.01)

    hessian = FiniteDiff.finite_difference_hessian(
        tvec -> imag(f(complex(tvec[1:2]...), complex(tvec[3:4]...))), 
        [reim(ti)..., reim(tr)...])
    # f_hessian(ti, tr) #my_hessian(b,Ip,q,ti,tr)
# 
    # this could certainly be made more julian    
    eigenvectors = [[complex(vec[1:2]...), complex(vec[3:4]...)] for vec in eachcol(eigvecs(hessian))]

    # eigenvectors 3 and 4 are the ones with positive sign. So if I want the steepest ascent thimble, then I should use those.
    pointsini = ([[ti,tr] .+ ϵ * (cos(θ) * eigenvectors[3] + sin(θ) * eigenvectors[4]) for θ in range(0, stop=2π, length=Ninit+1)])[1:end-1]
    # because 0 and 2π are the same and I don't want the point twice, me stupid!!!

    push!(points, [PointA(p[1], p[2]) for p in pointsini]...)
    push!(necklace, [LineSeg(i, i+1 , true) for i in 1:(length(points)-1)]...)
    push!(necklace, LineSeg(length(points), 1 , true)) # closing the necklace
end

### TODO this Δ could definitely get a more sophisticated default value
function subdivide!(lineseg::LineSeg,
    necklace::Vector{LineSeg}, points::Vector{PointA};
    Δ::Float64=1.)    

    p1 = points[lineseg.sindex]
    p2 = points[lineseg.eindex]

    active = p1.active && p2.active

    Δx(p1::PointA,p2::PointA) = norm(p2.x - p1.x)
    Δy(p1::PointA,p2::PointA) = norm(p2.y - p1.y)
    midx(p1::PointA,p2::PointA) = (p2.x + p1.x)./2
    midy(p1::PointA,p2::PointA) = (p2.y + p1.y)./2

    if active && ( max(Δx(p1,p2), Δy(p1,p2)) > Δ)    
        lineseg.active = false # lineseg gets turned inactive when being divided.       
        midpoint = PointA(midx(p1,p2), midy(p1,p2)) 
        push!(points, midpoint)

        lineseg1mid = LineSeg(lineseg.sindex,length(points),true)
        linesegmid2 = LineSeg(length(points),lineseg.eindex,true)
        push!(necklace, lineseg1mid)
        push!(necklace, linesegmid2)
    end
end

### flowing up
function flow!(necklace::Vector{LineSeg}, points::Vector{PointA},
        f::Function,
        f_grad::Function;
        δ::Float64=0.1,
        threshold::Float64=0.5
        )

    for i in 1:length(points)
        # TODO check both real and imaginary part?
        if points[i].active # for the active points
            # set them to be active (= still flowing) if they are above threshold
            # points[i].active = real(-im * S(b, Ip, points[i].x, points[i].y, q)) < 0 #(in Job's code that's h-function > thresh, I should clearly state which sign I'm using where etc.) 
            points[i].active = real(f(points[i].x, points[i].y)) < 0 #(in Job's code that's h-function > thresh, I should clearly state which sign I'm using where etc.) 

            if points[i].active
                step = δ .* gradN((ti,tr)->conj.(complex.(f_grad(ti,tr))), points[i].x, points[i].y, threshold)
                points[i].x += step[1]
                points[i].y += step[2]
            end
        end
    end
end

function adorn_necklace!(necklace::Vector{LineSeg}, points::Vector{PointA})
    for i in 1:length(necklace)
        necklace[i] = LineSeg(points[necklace[i].sindex], points[necklace[i].eindex])
    end
end;

### get necklace
function get_necklace_solver(f::Function,
    f_grad::Function,
    f_hessian::Function,
    ti::ComplexF64, tr::ComplexF64; 
    Ninit::Int64=20, Ncounter::Int64=600,
    eigvecfactorinit::Float64 = 0.04, # I should come up with sophisticated guesses here.
    flowstepfactor::Float64 = 0.4, 
    subdividethreshold::Float64 = 1.8 )
       
    necklace = Vector{LineSeg}()
    points = Vector{PointA}()

    initialise!(necklace, points, ti, tr, f, Ninit = Ninit, ϵ = eigvecfactorinit)

    ### find a suitable threshold for the normalisation of the gradient
    gradient0 = [norm(conj.(f_grad(p.x, p.y))) for p in points]
    threshold = round(minimum(gradient0), RoundDown, sigdigits=2)
    
    counter = 0
        
    while counter < Ncounter
        counter += 1

        tmp = deepcopy(necklace)
        flow!(necklace, points, f, f_grad, threshold = threshold, δ = flowstepfactor)

        if count([p.active for p in points]) == 0
            @debug "I broke because the flow stopped after $counter iterations"
        # println("I broke because the flow stopped after $counter iterations")
            break
        end
        for i in 1:length(necklace)
            subdivide!(necklace[i], necklace, points, Δ = subdividethreshold)
        end
        keepat!(necklace, [ls.active for ls in necklace])
    end
        
    if counter == Ncounter && Ncounter > 1
        println("I broke because the counter reached its max, i.e. $Ncounter.")
    end

    necklace = sort_linesegs(necklace)
    adorn_necklace!(necklace, points)
    
    return necklace
end;


function get_necklace(f::Function,
    f_grad::Function,
    f_hessian::Function,
    ti::ComplexF64, tr::ComplexF64; 
    logerrors::Bool=false,
    kwargs... # this passes on all th ekeyword arguments
    )
    
    necklace = get_necklace_solver(f, f_grad, f_hessian, ti, tr; kwargs...)

    necklace_init = get_necklace_solver(f, f_grad, f_hessian, ti, tr; kwargs..., Ncounter =1)
    enclosed_area_init = enclosed_area(necklace_init,imag) + enclosed_area(necklace_init,real)
    
    if (enclosed_area(necklace,imag) + enclosed_area(necklace,real)) > enclosed_area_init
        return necklace
    else
        if (real(f(ti, tr))) > -0.2
            return necklace
        else
            @warn ("Warning (3)! The necklace is smaller than its initialisation, real(f) = $(real(f(ti, tr)))")
            # println("Warning (3)! The necklace is smaller than its initialisation for beam $b at q $q with ti $ti and tr $tr, where h was $(real(-im * S(b, Ip, ti, tr, q)))!")
            # logerrors ? log_error("necklace-errors.txt", "Warning (3) for beam $b at q $q with ti $ti and tr $tr.") : nothing
            return nothing
        end

    end

    
    # counter = 0
    # while ((enclosed_area(necklace,imag) + enclosed_area(necklace,real)) < 0.5) && 
    #     counter < 4 && length(necklace) > Ninit+1

    #     necklace = get_necklace_solver(b, Ip, q, ti, tr; Ninit = Ninit, Ncounter=Ncounter,
    #     eigvecfactorinit = eigvecfactorinit*2, # I should come up with sophisticated guesses here.
    #     flowstepfactor = flowstepfactor, 
    #     subdividethreshold = subdividethreshold )

    #     if counter==3 && logerrors
    #         log_error("necklace-errors.txt", "Warning (1) for beam $b at q $q with ti $ti and tr $tr.")
    #     end

    #     Ncounter *= 2 # random other guess to improve the necklace finding
    #     counter += 1
    # end
    # return necklace
end;


nothing