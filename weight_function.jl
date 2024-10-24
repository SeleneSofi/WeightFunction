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

function weight_funtion(x::nGum, lnα)
    gum = exp.(-exp.(-x.a*(lnα .- x.d'))) * x.c
    gum .*= exp.(-1/2*lnα)
    for i = 1:x.n-x.l
        gum .*= (exp.(-(lnα .- x.r[i])) .- 1)
    end
    return gum
end

function make_array(x::Array{nGum}, op::Array{nGum_opt})
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

function make_wfstructs(xcar::Array{Float64, 1}, orbs::Array{Tuple, 1})
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