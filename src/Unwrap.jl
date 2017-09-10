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

include("unwrap2d.jl")

end # module
