module Unwrap

export unwrap1d

# Unwrap two neighbouring values
function γ1d(a,b)
    d = a-b
    if abs(d) > π
        d-round(d/2π)*2π
    else
        d
    end
end

# 1d vector unwrapping, seems faster than unwrap() from DSP.jl
function unwrap1d(v::AbstractVector)
    u = Array{eltype(v)}(length(v))
    u[1] = v[1]
    for i = 2:length(v)
        u[i] = u[i-1] + γ1d(v[i],v[i-1])
    end
    u
end

function unwrap1d(M::AbstractMatrix, dim=1)
    m,n = size(M)
    U = Matrix{eltype(M)}(m, n)
    if dim == 1
        for j in 1:n
            U[:,j] = unwrap1d(M[:,j])
        end
    else
        for j in 1:n
            U[i,:] = unwrap1d(M[1,:])
        end
    end
    U
end

include("unwrap2d.jl")

end # module
