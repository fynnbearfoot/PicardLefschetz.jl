
# using NLsolve
# using Sobol
export solve_first_derivative, find_saddles_sobol

function solve_first_derivative(drv::Function, t0::Vector{Float64},
    roundDigits::Int64=5)

    try
    ### using NLsolve
        function speqs!(F, x)
            F[1] = real(drv(x[1]+ x[2]*im))
            F[2] = imag(drv(x[1]+ x[2]*im))
        end

        result = nlsolve(speqs!, t0)    

        if converged(result)
            tiSP = result.zero[1] + im*result.zero[2]
            tiSP = round(tiSP, digits = roundDigits)
            return tiSP
        else
            return nothing
        end
    catch e
        println("Error in solve_SPEqs(): $e")
        return nothing
    end
end

function find_saddles_sobol(drv::Function,
        tmin::ComplexF64, tmax::ComplexF64,        
        N::Int64=200 # number of seeds generated per domain
    ) 
    
    roundDigits = 2 # I should certainly change this, it seems a bit excessive

    saddles = Vector{ComplexF64}()
    
    t_seq = SobolSeq(reim(tmin),reim(tmax))

    for i in 1:N
        t0 = Sobol.next!(t_seq)
        t0 = [t0[1]; t0[2]]

        ts = solve_first_derivative(drv, t0, roundDigits) 
#         ### check conditons and deposit in array
        if !isnothing(ts)
            ts_r = round(ts, digits=roundDigits)

            push!(saddles, complex(
                real(ts_r)==0 ? 0. : real(ts_r),
                imag(ts_r)==0 ? 0. : imag(ts_r)    
                ))

        end
    end

    unique!(ts -> round(ts, digits = roundDigits), saddles)
    sort!(saddles, by = x -> real(x))
    return saddles
end;



function find_saddle_similar_seed(drv::Function, ts::ComplexF64;roundDigits = 2)
    
    ts = solve_first_derivative(drv, [reim(ts)...], roundDigits) 
    
    if !isnothing(ts) 
        ts_r = round(ts, digits=roundDigits)
        return ts_r
    else
        return nothing
    end    
end


nothing