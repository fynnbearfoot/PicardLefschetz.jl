using NLsolve

function splat_complex(v::Vector{ComplexF64})
    return vec(collect(Iterators.flatten(reim.(v))))
end;



function solve_first_drv(drv::Function, t0::Vector{T};
    digits::Int64=5) where T<:Real

    try
        
        function speqs!(F,x)
            F .= vec(collect(Iterators.flatten(
            reim.(drv(complex(x[1], x[2]),
            complex(x[3], x[4]) ))
            )))
        end
        
        result = nlsolve(speqs!, t0)
        #, method = :trust_region, factor =fac )#, ftol = 1e-13)#, method = :anderson)
        if converged(result)
            t1SP = complex(result.zero[1:2]...)
            t2SP = complex(result.zero[3:4]...)

            t1SP = round(t1SP, digits = digits)
            t2SP = round(t2SP, digits = digits)
            return t1SP, t2SP
        else
            return nothing, nothing
        end
    catch e
        println("Error in solve_first_drv(): $e")
        return nothing,nothing
    end    
end 

function solve_first_drv(drv::Function, t0::Vector{T};
    digits::Int64=5) where T<:Complex
    return solve_first_drv(drv,splat_complex(t0),digits=digits )
end


function find_saddles_sobol(drv::Function,
        t1_cd::ComplexDomain, t2_cd::ComplexDomain;
        N::Int64=200, # number of seeds generated per domain
        digits::Int64=5,
        check::Function = (t1,t2) -> !isequal(t1,t2)
         ) 

    solutions = Vector{Vector{ComplexF64}}()
    
    t1_seq = SobolSeq(reim(t1_cd.min),reim(t1_cd.max))
    t2_seq = SobolSeq(reim(t2_cd.min),reim(t2_cd.max))

    for i in 1:N
        t10 = Sobol.next!(t1_seq)
        t20 = Sobol.next!(t2_seq)
        
#         @show t10
        t0 = [t10[1]; t10[2]; t20[1]; t20[2]]
#         @show t0
        t1s, t2s = solve_first_drv(drv, t0, digits = digits) 

        ### check conditions and deposit in array
        if !isnothing(t1s) && check(t1s,t2s)
            push!(solutions, [t1s, t2s])
        end
    end

    unique!( ts -> round.(ts, digits=digits), solutions)

    sort!(solutions, by = x -> real(x[1]))
    return solutions
end;