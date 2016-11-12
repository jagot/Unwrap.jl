#=
Implementation of the 2d unwrapping algorithm described in:

- Herráez, M. A., Burton, D. R., Lalor, M. J., & Gdeisat,
  M. A. (2002). Fast two-dimensional phase-unwrapping algorithm based
  on sorting by reliability following a noncontinuous path. Applied
  Optics, 41(35), 7437. http://dx.doi.org/10.1364/ao.41.007437

=#

using ProgressMeter

export unwrap2d

γ2d(v) = v - trunc(v/2π)*2π
γ2d(a,b) = γ2d(b-a)

# Returns integer number of 2π necessary to bring a to b's branch.
function γi(a,b)
    d = b-a
    if abs(d) > 0.5
        round(Int,d)
    else
        0
    end
end

function diff2(φ)
    D² = zeros(φ)

    m,n = size(φ)

    for i = 1:m
        D²[i,1] = D²[i,end] = Inf
    end
    for j = 1:n
        D²[1,j] = D²[end,j] = Inf
    end

    for i = 2:m-1
        @simd for j = 2:n-1
            @inbounds D²[i,j]  = (γ2d(φ[i-1,j  ], φ[i,j]) - γ2d(φ[i,j], φ[i+1,j  ]))^2 # H
            @inbounds D²[i,j] += (γ2d(φ[i  ,j-1], φ[i,j]) - γ2d(φ[i,j], φ[i  ,j+1]))^2 # V
            @inbounds D²[i,j] += (γ2d(φ[i-1,j-1], φ[i,j]) - γ2d(φ[i,j], φ[i+1,j+1]))^2 # D₁
            @inbounds D²[i,j] += (γ2d(φ[i-1,j+1], φ[i,j]) - γ2d(φ[i,j], φ[i+1,j-1]))^2 # D₂
            @inbounds D²[i,j]  = √D²[i,j]
        end
    end
    D²
end

type Edge
    i
    j
    k
    r
end

function find_edges(R²)
    m,n = size(R²)
    # Find all edges
    V = R²[1:end-1,:] + R²[2:end,:]
    H = R²[:,1:end-1] + R²[:,2:end]
    # Construct edge vector
    edges = Any[nothing,nothing]
    for (t,d,M) in ((1,:v,V), (2,:h, H))
        edges[t] = vec(Edge[Edge(i,j,d,M[i,j]) for i = 1:size(M,1), j = 1:size(M,2)])
    end
    edges = [edges...;]
    # Sort edges according to reliability
    sort!(edges, by = e -> e.r, rev=true)
    # for (i,j) in ((1,1),(1,n),(m,1),(m-1,n-1))
    #     filter!(e -> e.i != i && e.j != j, edges)
    # end
    edges
end

function coords(e::Edge)
    if e.k == :v
        [e.i, e.i+1], e.j
    else
        e.i, [e.j, e.j+1]
    end
end

function coords(c::Tuple)
    if length(c[1]) == 1
        [[c[1];c[1]] c[2]]
    else
        [c[1] [c[2];c[2]]]
    end
end
function coords(c::Tuple,i::Integer)
    if length(c[1]) == 1
        c[1],c[2][i]
    else
        c[1][i],c[2]
    end
end

function form_group!(c, φ, groups, group_map, shifts)
    g_id = length(groups) + 1
    group_map[c...] = g_id
    cc = coords(c)
    gg = γi(φ[c...]...)
    shifts[cc[2,:]] = γi(φ[c...]...)
    append!(groups, Matrix{Int}[cc])
end

function add_to_group!(c, φ, groups, group_map, shifts)
    g = group_map[c...]
    ds = γi(φ[c...]...)
    cc = coords(c)
    if g[1] == 0
        group_map[coords(c,1)...] = g[2]
        groups[g[2]] = [groups[g[2]]; cc[1,:]']
        shifts[cc[1,:]] = shifts[cc[2,:]] + ds 
    else
        group_map[coords(c,2)...] = g[1]
        groups[g[1]] = [groups[g[1]]; cc[2,:]']
        shifts[cc[2,:]] = shifts[cc[1,:]] - ds 
    end
end

function join_groups!(c, φ, groups, group_map, shifts)
    g = group_map[c...]
    ds = γi(φ[c...]...)
    g1,g2 = groups[g[1]], groups[g[2]]
    l1,l2 = size(g1,1), size(g2,1)
    if l1 > l2
        ds = shifts[g1[1,:]...] - ds
        for t = 1:l2
            group_map[g2[t,:]...] = g[1]
            shifts[g2[t,:]...] += ds
        end
        groups[g[1]] = [groups[g[1]]; g2]
        groups[g[2]] = Matrix{Int}()
    else
        ds = shifts[g1[2,:]...] + ds
        for t = 1:l1
            group_map[g1[t,:]...] = g[2]
            shifts[g1[t,:]...] += ds
        end
        groups[g[2]] = [groups[g[2]]; g1]
        groups[g[1]] = Matrix{Int}()
    end
end

function unwrap2d(φ)
    φ = copy(φ)
    D² = diff2(φ)
    R² = 1./D²
    φ　/= 2π
    
    edges = find_edges(R²)

    groups = Matrix{Int}[]
    group_map = zeros(Int, size(φ))

    shifts = zeros(Int, size(φ))

    mg = zeros(length(edges))
    lg = zeros(length(edges))

    @showprogress "Unwrapping " for t = 1:length(edges)
        e = edges[t]
        c = coords(e)
        g = group_map[c...]
        
        if g[1] == g[2] == 0 # Neither point in any group
            form_group!(c, φ, groups, group_map, shifts)
        elseif g[1] == 0 || g[2] == 0 # One point in a group
            add_to_group!(c, φ, groups, group_map, shifts)
        elseif g[1] != g[2] # Both points in group
            join_groups!(c, φ, groups, group_map, shifts)
        end
        mg[t] = maximum(map(length, groups))
        lg[t] = length(groups)
    end
    2π*(φ+shifts)
end

function unwrap2d_debug(φ)
    φ = copy(φ)
    println("Calculating 2nd difference")
    @time D² = diff2(φ)
    @time R² = 1./D²
    φ　/= 2π
    
    println("Constructing edges")
    edges = find_edges(R²)

    println("Building groups")
    groups = Matrix{Int}[]
    group_map = zeros(Int, size(φ))

    shifts = zeros(Int, size(φ))

    mg = zeros(length(edges))
    lg = zeros(length(edges))

    @time for t = 1:length(edges)
        e = edges[t]
        c = coords(e)
        g = group_map[c...]
        
        if g[1] == g[2] == 0 # Neither point in any group
            form_group!(c, φ, groups, group_map, shifts)
        elseif g[1] == 0 || g[2] == 0 # One point in a group
            add_to_group!(c, φ, groups, group_map, shifts)
        elseif g[1] != g[2] # Both points in group
            join_groups!(c, φ, groups, group_map, shifts)
        end
        mg[t] = maximum(map(length, groups))
        lg[t] = length(groups)
    end
    println()
    2π*(φ+shifts),D²,R²,mg,lg,group_map,shifts
end
