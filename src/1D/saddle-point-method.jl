# using Contour, GeometryBasics
# export is_contributing, integrate_around_saddle_point, integrate_SPM

function is_contributing(ts::ComplexF64, S::Function, tmin::ComplexF64, tmax::ComplexF64;
    Ntimes::Int64=100 )
    timags = range(imag(tmin), stop = imag(tmax), length = Ntimes)
    treals = range(real(tmin), stop = real(tmax), length = Ntimes)
    tlength = real(tmax-tmin)
    
    # Δi = timags[2]-timags[1]
    # Δr = treals[2]-treals[1]
    crossthresh = sqrt((step(timags))^2 + (step(treals))^2)

    Svals = [S(tr + im*ti) for tr in treals, ti in timags] 

    real_axis = Curve2([(real(tmin)-tlength, 0.), (real(tmax)+tlength, 0.)])    

    relevant = false
    S_saddle = S(ts)


    con_S_saddle_real = Contour.contour(collect(treals), collect(timags), real.(Svals), real.(S_saddle))
    for curve in con_S_saddle_real.lines              
        saddle_ip = crosses_point(curve, GeometryBasics.Point(reim(ts)...), crossthresh)
        ### filter for those level lines that actually intersect the saddle point
        if !(saddle_ip==false)

            segs = dissect_curve(ts, curve, saddle_ip, crossthresh)          
#             saddle_ip = closest_intersection(saddle_ip, curve, GeometryBasics.Point(reim(ts)...), crossthresh)
# #                 ### disect curve
#             c1 = Curve2(curve.vertices[1:saddle_ip+1])
#             c2 = Curve2(curve.vertices[saddle_ip+1:end])

            for c in segs
                actiondiff = S(complex(c.vertices[minimum([5, length(c.vertices)])]...)) - S(ts)
                if imag(-actiondiff) > 0 # ascent lines
                    if !isempty(intersection(c, real_axis)) 
                        relevant = true
                    end
                end
            end
        end
    end
    
    return relevant
    
end




function integrate_around_saddle_point(ts::ComplexF64,
    S::Function, drv::Function, drv2::Function
    ; prefactor::Function=t->1.,
    )

    S_ts = S(ts)
    
    int = prefactor(ts) * sqrt(-im*2π/drv2(ts)) * exp(im * S_ts)
end


function integrate_SPM(S::Function, drv::Function, drv2::Function,
    tmin::ComplexF64, tmax::ComplexF64
    ; prefactor::Function = t-> 1.)
    
    saddles = filter(ts -> real(tmin) < real(ts) < real(tmax), 
        find_saddles_sobol(drv, tmin, tmax, 300)
        )
    int_SPM = complex(0.)
    for ts in saddles      
        if is_contributing(ts, S, tmin, tmax)
            int_SPM += prefactor(ts) * integrate_around_saddle_point(ts, S, drv, drv2)  
        end
    end
    return int_SPM
end



nothing