############################################################
# Complex domain

import Base.+, Base.*, Base.in

struct ComplexDomain
    min::ComplexF64
    max::ComplexF64 #Union{ComplexF64,Nothing}
    
    ComplexDomain(rmin::Real,rmax::Real,imin::Real,imax::Real) = new(rmin+imin*im,rmax+imax*im)
    
    ComplexDomain(min::ComplexF64,max::ComplexF64) = new(min,max)

    ComplexDomain() = new(zero(ComplexF64),zero(ComplexF64))
    
end

function +(cd1::ComplexDomain,cd2::ComplexDomain)
    rmin = minimum([real(cd1.min),real(cd2.min)])
    rmax = maximum([real(cd1.max),real(cd2.max)])
    imin = minimum([imag(cd1.min),imag(cd2.min)])
    imax = maximum([imag(cd1.max),imag(cd2.max)])
    
    return ComplexDomain(rmin,rmax,imin,imax)
end

function +(cd1::ComplexDomain, z::Number) # this function shifts the whole domain by the number specified
    
    rmin = real(cd1.min) + real(z)
    rmax = real(cd1.max) + real(z)
    imin = imag(cd1.min) + imag(z)
    imax = imag(cd1.max) + imag(z)
    
    return ComplexDomain(rmin,rmax,imin,imax)
end

function *(cd1::ComplexDomain, z::Real) # this function multiplies the whole domain by the number specified. Useful when I multiply by TCycle
    
    rmin = real(cd1.min) * z
    rmax = real(cd1.max) * z
    imin = imag(cd1.min) * z
    imax = imag(cd1.max) * z
    
    return ComplexDomain(rmin,rmax,imin,imax)
end;

function in(z::Complex,cd::ComplexDomain)
    return (real(cd.min) <= real(z) < real(cd.max) ) && (imag(cd.min) <= imag(z) < imag(cd.max) )
end;

function realrange(cd::ComplexDomain, N::Int64=50)
    return range(real(cd.min), stop=real(cd.max), length=N) 
end
function imagrange(cd::ComplexDomain, N::Int64=50)
    return range(imag(cd.min), stop=imag(cd.max), length=N) 
end

##########################################


scalarproduct(a::AbstractVector{}, b::AbstractVector{}) = sum( a .* b )
scalarproduct(a::AbstractVector{}) = scalarproduct(a,a)
scalarproduct2(a::AbstractVector{}) = scalarproduct(a,a)
scalarproduct(a::Number, b::Number) = a .* b;
scalarproduct2(a::Number) = scalarproduct(a,a);


#####################
function log_error(file_path::String, error_message::String)
    # Check if the file exists
    if !isfile(file_path)
        # Create the file if it doesn't exist
        open(file_path, "w") do file
            # Just open and close the file to create it
        end
    end

    # Append the error message to the file
    open(file_path, "a") do file
        write(file, error_message * "\n")
    end
end

#############

import Base.in
function in(x::Real, xmin::Real, xmax::Real; inclusive::Bool=true)
    if inclusive
        return xmin <= x <= xmax
    else
        return xmin < x < xmax
    end
end;


#########




nothing