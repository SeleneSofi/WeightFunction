struct nGum
    a::Float64
    n::Int8
    l::Int8
    r::Array{Float64, 1}
    c::Array{Float64, 1}
    d::Array{Float64, 1}
end

struct nGum_opt
    a::Bool
    r::Array{Bool, 1}
    c::Array{Bool, 1}
    d::Array{Bool, 1}
end

struct nSig
    n::Int8
    l::Int8
    a::Float64
    v::Float64
    r::Array{Float64, 1}
    c::Array{Float64, 1}
    d::Array{Float64, 1}
end

struct nSig_opt
    a::Bool
    v::Bool
    r::Array{Bool, 1}
    c::Array{Bool, 1}
    d::Array{Bool, 1}
end

struct nGuma
    n::Int8
    l::Int8
    r::Array{Float64, 1}
    a::Array{Float64, 1}
    c::Array{Float64, 1}
    d::Array{Float64, 1}
end

struct nGuma_opt
    r::Array{Bool, 1}
    a::Array{Bool, 1}
    c::Array{Bool, 1}
    d::Array{Bool, 1}
end

struct nGumqe
    a::Float64
    n::Int8
    l::Int8
    r::Array{Float64, 1}
    c::Array{Float64, 1}
    q::Array{Float64, 1}
    d::Array{Float64, 1}
end

struct nGumqe_opt
    a::Bool
    r::Array{Bool, 1}
    c::Array{Bool, 1}
    q::Array{Bool, 1}
    d::Array{Bool, 1}
end

# struct nGLF
#     n::Int8
#     l::Int8
#     r::Array{Float64, 1}
#     a::Array{Float64, 1}
#     b::Array{Float64, 1}
#     z::Array{Float64, 1}
#     u::Array{Float64, 1}
#     v::Array{Float64, 1}
# end

# struct nGLF_opt
#     r::Array{Bool, 1}
#     a::Array{Bool, 1}
#     b::Array{Bool, 1}
#     z::Array{Bool, 1}
#     u::Array{Bool, 1}
#     v::Array{Bool, 1}
# end

function qexp(x::Real, q)
    if q == 1
        return exp.(x)
    elseif q < 1
        if x < -1/(1-q)
            return 0
        else
            return (1 + (1-q)*x)^(1/(1-q))
        end
    else
        if x > -1/(1-q)
            return Inf
        else
            return (1 + (1-q)*x)^(1/(1-q))
        end
    end
end

function qexp(x::Array, q)
    return [qexp(ix, q) for ix = x]
end


function weight_funtion(x::nGum, lnα)
    gum = exp.(-exp.(-abs(x.a)*(lnα .- x.d'))) * abs.(x.c)
    gum .*= exp.(-1/2*lnα)
    for i = 1:x.n-x.l
        gum .*= (exp.(-(lnα .- x.r[i])) .- 1)
    end
    return gum
end

function weight_funtion(x::nSig, lnα)
    a = abs(x.a)*exp(-1)*(1+exp(-x.v))^(1+exp(x.v))
    d = x.d .- x.v/a
    gum = 1 ./ (1 .+ exp.(-a*(lnα .- d'))).^(exp(x.v)) * abs.(x.c)
    gum .*= exp.(-1/2*lnα)
    for i = 1:x.n-x.l
        gum .*= (exp.(-(lnα .- x.r[i])) .- 1)
    end
    return gum
end

function weight_funtion(x::nGuma, lnα)
    gum = exp.(-exp.(-abs.(x.a)'.*(lnα .- x.d'))) * abs.(x.c)
    gum .*= exp.(-1/2*lnα)
    for i = 1:x.n-x.l
        gum .*= (exp.(-(lnα .- x.r[i])) .- 1)
    end
    return gum
end

function weight_funtion(x::nGumqe, lnα)
    gum = [qexp(-qexp(-abs(x.a)*(ilnα - x.d[j]),x.q[j]),x.q[j]) for ilnα = lnα, j = eachindex(x.c)]
    gum = gum * abs.(x.c)
    gum .*= exp.(-1/2*lnα)
    for i = 1:x.n-x.l
        gum .*= (exp.(-(lnα .- x.r[i])) .- 1)
    end
    return gum
end

# function weight_funtion(x::nGLF, lnα)
#     a, z, u, v = abs.((x.a, x.z, x.u, x.v))
#     b = x.b
#     gum = 
#     gum = gum * abs.(x.c)
#     gum .*= exp.(-1/2*lnα)
#     for i = 1:x.n-x.l
#         gum .*= (exp.(-(lnα .- x.r[i])) .- 1)
#     end
#     return gum
# end

function make_array(x::Array, op::Array)
    xar = Array{Float64}(undef, 0)
    car = Array{Float64}(undef, 0)
    orbs = Array{Tuple}(undef, 0)
    for i = 1:size(x,1)
        nci = size(car,1)
        n = x[i].n
        l = x[i].l
        for name = fieldnames(typeof(x[i]))
            name in [:n, :l] ? continue :
            prop = getproperty(x[i], name)
            opprop = getproperty(op[i], name)
            if opprop isa Bool
                if opprop
                    append!(xar, prop)
                    append!(car, fill(NaN, size(prop,1)))
                else
                    append!(car, prop)
                end
            else
                for j = eachindex(opprop)
                    if opprop[j]
                        append!(xar, prop[j])
                        append!(car, NaN)
                    else
                        append!(car, prop[j])
                    end
                end
            end
        end
        push!(orbs, (size(car,1)-nci, n, l))
    end
    return xar, car, orbs
end

# function make_array(x::Array{nGum}, op::Array{nGum_opt})
#     xar = Array{Float64}(undef, 0)
#     car = Array{Float64}(undef, 0)
#     orbs = Array{Tuple}(undef, 0)
#     for i = 1:size(x,1)
#         nci = size(car,1)
#         n = x[i].n
#         l = x[i].l
#         for name = fieldnames(typeof(x[i]))
#             name in [:n, :l] ? continue :
#             prop = getproperty(x[i], name)
#             opprop = getproperty(op[i], name)
#             if opprop isa Bool
#                 if opprop
#                     append!(xar, prop)
#                     append!(car, fill(NaN, size(prop,1)))
#                 else
#                     append!(car, prop)
#                 end
#             else
#                 for j = eachindex(opprop)
#                     if opprop[j]
#                         append!(xar, prop[j])
#                         append!(car, NaN)
#                     else
#                         append!(car, prop[j])
#                     end
#                 end
#             end
#         end
#         push!(orbs, (size(car,1)-nci, n, l))
#     end
#     return xar, car, orbs
# end

# function make_array(x::Array{nSig}, op::Array{nSig_opt})
#     xar = Array{Float64}(undef, 0)
#     car = Array{Float64}(undef, 0)
#     orbs = Array{Tuple}(undef, 0)
#     for i = 1:size(x,1)
#         nci = size(car,1)
#         n = x[i].n
#         l = x[i].l
#         for name = fieldnames(typeof(x[i]))
#             name in [:n, :l] ? continue :
#             prop = getproperty(x[i], name)
#             opprop = getproperty(op[i], name)
#             if opprop isa Bool
#                 if opprop
#                     append!(xar, prop)
#                     append!(car, fill(NaN, size(prop,1)))
#                 else
#                     append!(car, prop)
#                 end
#             else
#                 for j = eachindex(opprop)
#                     if opprop[j]
#                         append!(xar, prop[j])
#                         append!(car, NaN)
#                     else
#                         append!(car, prop[j])
#                     end
#                 end
#             end
#         end
#         push!(orbs, (size(car,1)-nci, n, l))
#     end
#     return xar, car, orbs
# end

# function make_array(x::Array{nGuma}, op::Array{nGuma_opt})
#     xar = Array{Float64}(undef, 0)
#     car = Array{Float64}(undef, 0)
#     orbs = Array{Tuple}(undef, 0)
#     for i = 1:size(x,1)
#         nci = size(car,1)
#         n = x[i].n
#         l = x[i].l
#         for name = fieldnames(typeof(x[i]))
#             name in [:n, :l] ? continue :
#             prop = getproperty(x[i], name)
#             opprop = getproperty(op[i], name)
#             if opprop isa Bool
#                 if opprop
#                     append!(xar, prop)
#                     append!(car, fill(NaN, size(prop,1)))
#                 else
#                     append!(car, prop)
#                 end
#             else
#                 for j = eachindex(opprop)
#                     if opprop[j]
#                         append!(xar, prop[j])
#                         append!(car, NaN)
#                     else
#                         append!(car, prop[j])
#                     end
#                 end
#             end
#         end
#         push!(orbs, (size(car,1)-nci, n, l))
#     end
#     return xar, car, orbs
# end

# function make_array(x::Array{nGumqe}, op::Array{nGumqe_opt})
#     xar = Array{Float64}(undef, 0)
#     car = Array{Float64}(undef, 0)
#     orbs = Array{Tuple}(undef, 0)
#     for i = 1:size(x,1)
#         nci = size(car,1)
#         n = x[i].n
#         l = x[i].l
#         for name = fieldnames(typeof(x[i]))
#             name in [:n, :l] ? continue :
#             prop = getproperty(x[i], name)
#             opprop = getproperty(op[i], name)
#             if opprop isa Bool
#                 if opprop
#                     append!(xar, prop)
#                     append!(car, fill(NaN, size(prop,1)))
#                 else
#                     append!(car, prop)
#                 end
#             else
#                 for j = eachindex(opprop)
#                     if opprop[j]
#                         append!(xar, prop[j])
#                         append!(car, NaN)
#                     else
#                         append!(car, prop[j])
#                     end
#                 end
#             end
#         end
#         push!(orbs, (size(car,1)-nci, n, l))
#     end
#     return xar, car, orbs
# end

function make_wfstr_nGum(xcar::Array{Float64, 1}, orbs::Array{Tuple, 1})
    xc = Array{nGum}(undef, 0)
    for orb = orbs
        n, l = orb[2:3]
        a = popfirst!(xcar)
        nr = n - l
        r = [popfirst!(xcar) for _ = 1:nr]
        ncd = orb[1] - (nr + 1)
        c = [popfirst!(xcar) for _ = 1:ncd/2]
        d = [popfirst!(xcar) for _ = 1:ncd/2]
        push!(xc, nGum(a, n, l, r, c, d))
    end
    return xc
end

function make_wfstr_nSig(xcar::Array{Float64, 1}, orbs::Array{Tuple, 1})
    xc = Array{nSig}(undef, 0)
    for orb = orbs
        n, l = orb[2:3]
        a = popfirst!(xcar)
        v = popfirst!(xcar)
        nr = n - l
        r = [popfirst!(xcar) for _ = 1:nr]
        ncd = orb[1] - (nr + 2)
        # ncd = size(xcar,1)
        c = abs.([popfirst!(xcar) for _ = 1:ncd/2])
        c = c/sum(c)
        d = [popfirst!(xcar) for _ = 1:ncd/2]
        id = sortperm(d)
        c = c[id]
        d = d[id]
        push!(xc, nSig(n, l, a, v, r, c, d))
    end
    return xc
end

function make_wfstr_nGuma(xcar::Array{Float64, 1}, orbs::Array{Tuple, 1})
    xc = Array{nGuma}(undef, 0)
    for orb = orbs
        n, l = orb[2:3]
        nr = n - l
        r = [popfirst!(xcar) for _ = 1:nr]
        ncd = orb[1] - nr
        a = [popfirst!(xcar) for _ = 1:ncd/3]
        c = [popfirst!(xcar) for _ = 1:ncd/3]
        d = [popfirst!(xcar) for _ = 1:ncd/3]
        push!(xc, nGuma(n, l, r, a, c, d))
    end
    return xc
end

function make_wfstr_nGumqe(xcar::Array{Float64, 1}, orbs::Array{Tuple, 1})
    xc = Array{nGumqe}(undef, 0)
    for orb = orbs
        n, l = orb[2:3]
        a = popfirst!(xcar)
        nr = n - l
        r = [popfirst!(xcar) for _ = 1:nr]
        ncd = orb[1] - (nr + 1)
        c = [popfirst!(xcar) for _ = 1:ncd/3]
        q = [popfirst!(xcar) for _ = 1:ncd/3]
        d = [popfirst!(xcar) for _ = 1:ncd/3]
        push!(xc, nGumqe(a, n, l, r, c, q, d))
    end
    return xc
end