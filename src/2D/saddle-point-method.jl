function hessian_root(h::AbstractArray) 
    ### I should definitely work this out properly and also make sure this actually gets rid of the branch cuts...

#     h = f_hessian(ti, tr)
    xd2S_dti2 = -conj(complex(h[1:2,1]...))
    xd2S_dtr2 = -conj(complex(h[3:4,3]...))
    xd2S_dtitr = -conj(complex(h[3:4,1]...))
    
    # ### RBSFA hessian_root
    sqrt1 = sqrt(2π/ (im*xd2S_dti2 ))
    sqrt2 = sqrt(2π * xd2S_dti2 / (im*(xd2S_dtr2 * xd2S_dti2 - xd2S_dtitr * xd2S_dtitr)) )
    return sqrt1 * sqrt2
    
    ### hessian_determinant < = works better for now. But maybe needs to be changed
    # hdet = xd2S_dtr2 * xd2S_dti2 - xd2S_dtitr * xd2S_dtitr
    # return im * 2*π/sqrt(hdet)

end



function saddles_gaussian_contribution(f::Function,
    f_hessian::Function,
    ti::ComplexF64, tr::ComplexF64;
    prefactor::Function = (ti,tr) -> ones(2)
    )  
        
    ### prefactor for the saddle-point method
    prefactor_spm = hessian_root(f_hessian(ti,tr))

    return prefactor(ti,tr) .* prefactor_spm .* exp(f(ti,tr))
        
end