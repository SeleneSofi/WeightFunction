using LinearAlgebra
using TensorOperations

function Sabss(a::Float64, b::Float64)
    return (π/(a+b))^(3/2)
end

function Sabpp(a::Float64, b::Float64)
    return (π/(a+b))^(3/2)/2/(a+b)
end

function Habss(a::Float64, b::Float64, ζ::Int16)
    return 3*a*b/(a+b)*Sabss(a,b) - 2*π*ζ/(a+b)
end

function Habpp(a::Float64, b::Float64, ζ::Int16)
    return 5*a*b/(a+b)*Sabpp(a,b) - 2*π*ζ/3/(a+b)^2
end

function Vabcdssss(a::Float64, b::Float64, c::Float64, d::Float64)
    return 2*sqrt(π^5)/((a+b)*(c+d)*sqrt(a+b+c+d))
end

function Vabcdpppp(a::Float64, b::Float64, c::Float64, d::Float64)
    # return 2*sqrt(π^5)/(a+b)/(c+d)/sqrt(a+b+c+d)*(1/6/(a+b)/(c+d) + 3/20/(a+b+c+d)^2)
    return sqrt(π^5)/30/((a+b)^2*(c+d)^2*sqrt((a+b+c+d)^5))*(10*(a+b+c+d)^2 + 9*(a+b)*(c+d))
end

function Vabcdsspp(a::Float64, b::Float64, c::Float64, d::Float64)
    # return 2*sqrt(π^5)/(a+b)/(c+d)/sqrt(a+b+c+d)/2/(c+d)*(1-(a+b)/3/(a+b+c+d))
    return sqrt(π^5)/3/((a+b)*(c+d)^2*sqrt((a+b+c+d)^3))*(3*(a+b+c+d)-(a+b))
end

function Vabcdspsp(a::Float64, b::Float64, c::Float64, d::Float64)
    # return 2*sqrt(π^5)/(a+b)/(c+d)/sqrt(a+b+c+d)/6/(a+b+c+d)
    return sqrt(π^5)/3/((a+b)*(c+d)*sqrt((a+b+c+d)^3))
end

function Vabcdpxpxpypy(a::Float64, b::Float64, c::Float64, d::Float64)
    # return 2*sqrt(π^5)/(a+b)/(c+d)/sqrt(a+b+c+d)*(1/6/(a+b)/(c+d) + 1/20/(a+b+c+d)^2)
    return sqrt(π^5)/30/((a+b)^2*(c+d)^2*sqrt((a+b+c+d)^5))*(10*(a+b+c+d)^2 + 3*(a+b)*(c+d))
end

function Vabcdpxpypxpy(a::Float64, b::Float64, c::Float64, d::Float64)
    # return sqrt(π^5)/(a+b)/(c+d)/sqrt(a+b+c+d)/10/(a+b+c+d)^2
    return sqrt(π^5)/10/((a+b)*(c+d)*sqrt((a+b+c+d)^5))
end

function normalized(
    als::Array{Float64, 1},
    gs::Array{Float64, 2},
)
    gsn = similar(gs)
    Ssst = [Sabss(i[1],i[2]) for i = Iterators.product(als, als)]
    for i = 1:size(gs,2)
        gsn[:,i] .= 1/sqrt(gs[:,i]'*Ssst*gs[:,i]) * gs[:,i]
    end

    return gsn
end

function normalized(
    als::Array{Float64, 1},
    alp::Array{Float64, 1},
    gs::Array{Float64, 2},
    gp::Array{Float64, 2},
)
    gsn = similar(gs)
    gpn = similar(gp)
    Ssst = [Sabss(i[1],i[2]) for i = Iterators.product(als, als)]
    for i = 1:size(gs,2)
        gsn[:,i] .= 1/sqrt(gs[:,i]'*Ssst*gs[:,i]) * gs[:,i]
    end

    Sppt = [Sabpp(i[1],i[2]) for i = Iterators.product(alp, alp)]
    for i = 1:size(gp,2)
        gpn[:,i] .= 1/sqrt(gp[:,i]'*Sppt*gp[:,i]) * gp[:,i]
    end

    return gsn, gpn
end

function orthonormalized(
    als::Array{Float64, 1},
    gs::Array{Float64, 2},
)
    gsn = similar(gs)
    Ssst = [Sabss(i[1],i[2]) for i = Iterators.product(als, als)]
    ng = size(gs,2)
    Sss = Array{Float64}(undef, ng, ng)
    for i = 1:ng, j = 1:ng
        Sss[i,j] = gs[:,i]'*Ssst*gs[:,j]
    end
    gsn .= gs*Sss^(-1/2)

    return gsn
end

function orthonormalized(
    als::Array{Float64, 1},
    alp::Array{Float64, 1},
    gs::Array{Float64, 2},
    gp::Array{Float64, 2},
)
    gsn = similar(gs)
    gpn = similar(gp)
    Ssst = [Sabss(i[1],i[2]) for i = Iterators.product(als, als)]
    ng = size(gs,2)
    Sss = Array{Float64}(undef, ng, ng)
    for i = 1:ng, j = 1:ng
        Sss[i,j] = gs[:,i]'*Ssst*gs[:,j]
    end
    gsn .= gs*Sss^(-1/2)

    Sppt = [Sabpp(i[1],i[2]) for i = Iterators.product(alp, alp)]
    ng = size(gp,2)
    Spp = Array{Float64}(undef, ng, ng)
    for i = 1:ng, j = 1:ng
        Spp[i,j] = gp[:,i]'*Sppt*gp[:,j]
    end
    gpn .= gp*Spp^(-1/2)

    return gsn, gpn
end

@fastmath function HValpha(
    als::Array{Float64, 1},
    ζ::Int16
)
    ns = size(als,1)
    Hsst = Array{Float64}(undef, ns, ns)
    @inbounds for i = Iterators.product(1:ns, 1:ns)
        @fastmath @views Hsst[i[1],i[2]] = Habss(als[i[1]],als[i[2]],ζ)
    end

    Vsssst = Array{Float64}(undef, ns, ns, ns, ns)
    @inbounds for i = Iterators.product(1:ns, 1:ns, 1:ns, 1:ns)
        @fastmath @views Vsssst[i[1],i[2],i[3],i[4]] = Vabcdssss(als[i[1]],als[i[2]],als[i[3]],als[i[4]])
    end
    return (Hsst,), (Vsssst,)
end

@fastmath function HValpha(
    als::Array{Float64, 1},
    alp::Array{Float64, 1},
    ζ::Int16
)
    ns = size(als,1)
    np = size(alp,1)
    Hsst = Array{Float64}(undef, ns, ns)
    @inbounds for i = Iterators.product(1:ns, 1:ns)
        @fastmath @views Hsst[i[1],i[2]] = Habss(als[i[1]],als[i[2]],ζ)
    end

    Hppt = Array{Float64}(undef, np, np)
    @inbounds for i = Iterators.product(1:np, 1:np)
        @fastmath @views Hppt[i[1],i[2]] = Habpp(alp[i[1]],alp[i[2]],ζ)
    end

    Vsssst = Array{Float64}(undef, ns, ns, ns, ns)
    @inbounds for i = Iterators.product(1:ns, 1:ns, 1:ns, 1:ns)
        @fastmath @views Vsssst[i[1],i[2],i[3],i[4]] = Vabcdssss(als[i[1]],als[i[2]],als[i[3]],als[i[4]])
    end

    Vssppt = Array{Float64}(undef, ns, ns, np, np)
    @inbounds for i = Iterators.product(1:ns, 1:ns, 1:np, 1:np)
        @fastmath @views Vssppt[i[1],i[2],i[3],i[4]] = Vabcdsspp(als[i[1]],als[i[2]],alp[i[3]],alp[i[4]])
    end

    Vspspt = Array{Float64}(undef, ns, np, ns, np)
    @inbounds for i = Iterators.product(1:ns, 1:np, 1:ns, 1:np)
        @fastmath @views Vspspt[i[1],i[2],i[3],i[4]] = Vabcdspsp(als[i[1]],alp[i[2]],als[i[3]],alp[i[4]])
    end

    Vppppt = Array{Float64}(undef, np, np, np, np)
    @inbounds for i = Iterators.product(1:np, 1:np, 1:np, 1:np)
        @fastmath @views Vppppt[i[1],i[2],i[3],i[4]] = Vabcdpppp(alp[i[1]],alp[i[2]],alp[i[3]],alp[i[4]])
    end

    Vpxpxpypyt = Array{Float64}(undef, np, np, np, np)
    @inbounds for i = Iterators.product(1:np, 1:np, 1:np, 1:np)
        @fastmath @views Vpxpxpypyt[i[1],i[2],i[3],i[4]] = Vabcdpxpxpypy(alp[i[1]],alp[i[2]],alp[i[3]],alp[i[4]])
    end

    Vpxpypxpyt = Array{Float64}(undef, np, np, np, np)
    @inbounds for i = Iterators.product(1:np, 1:np, 1:np, 1:np)
        @fastmath @views Vpxpypxpyt[i[1],i[2],i[3],i[4]] = Vabcdpxpypxpy(alp[i[1]],alp[i[2]],alp[i[3]],alp[i[4]])
    end

    return (Hsst, Hppt), (Vsssst, Vssppt, Vspspt, Vppppt, Vpxpxpypyt, Vpxpypxpyt)
end

@fastmath function integral2g(
    gs::Array{Float64, 2},
    Ht::Tuple,
)
    Hss = dot.(eachrow(gs'*Ht[1]),eachcol(gs))

    return (Hss,)

end

@fastmath function integral2g(
    gs::Array{Float64, 2},
    gp::Array{Float64, 2},
    Ht::Tuple,
)
    Hss = dot.(eachrow(gs'*Ht[1]),eachcol(gs))
    Hpp = dot.(eachrow(gp'*Ht[2]),eachcol(gp))

    return Hss, Hpp

end

function integral4g(
    gs::Array{Float64, 2},
    Vt::Tuple,
)
    ngs = size(gs,2)

    Vsssst = Vt[1]

    @tensor Jab[i,j,c,d] := gs[a,i]*Vsssst[a,b,c,d]*gs[b,j]
    @tensor J[i,j,k,l] := gs[c,k]*Jab[i,j,c,d]*gs[d,l]

    Jssss = Array{Float64}(undef, ngs, ngs)
    Kssss = Array{Float64}(undef, ngs, ngs)
    @fastmath @views @inbounds for i = 1:ngs, j = i:ngs
        Jssss[i,j] = J[i,i,j,j]
        Kssss[i,j] = J[i,j,i,j]
        if i ≠ j
            Jssss[j,i] = Jssss[i,j]
            Kssss[j,i] = Kssss[i,j]
        end
    end

    return (Jssss,), (Kssss,)
end

function integral4g(
    gs::Array{Float64, 2},
    gp::Array{Float64, 2},
    Vt::Tuple,
    # ne::Int16
)
    ngs = size(gs,2)

    Vsssst = Vt[1]
    @tensor Jab[i,j,c,d] := gs[a,i]*Vsssst[a,b,c,d]*gs[b,j]
    @tensor J[i,j,k,l] := gs[c,k]*Jab[i,j,c,d]*gs[d,l]

    Jssss = Array{Float64}(undef, ngs, ngs)
    Kssss = Array{Float64}(undef, ngs, ngs)
    @fastmath @views @inbounds for i = 1:ngs, j = i:ngs
        Jssss[i,j] = J[i,i,j,j]
        Kssss[i,j] = J[i,j,i,j]
        if i ≠ j
            Jssss[j,i] = Jssss[i,j]
            Kssss[j,i] = Kssss[i,j]
        end
    end

    ngp = size(gp,2)

    Vssppt = Vt[2]
    @tensor Jab[i,j,c,d] := gs[a,i]*Vssppt[a,b,c,d]*gs[b,j]
    @tensor J[i,j,k,l] := gp[c,k]*Jab[i,j,c,d]*gp[d,l]

    Jsspp = Array{Float64}(undef, ngs, ngp)
    @inbounds for i = 1:ngs, j = 1:ngp
        @views Jsspp[i,j] = J[i,i,j,j]
    end

    Vspspt = Vt[3]
    @tensor Jab[i,j,c,d] := gs[a,i]*Vspspt[a,b,c,d]*gp[b,j]
    @tensor J[i,j,k,l] := gs[c,k]*Jab[i,j,c,d]*gp[d,l]

    Kspsp = Array{Float64}(undef, ngs, ngp)
    @inbounds for i = 1:ngs, j = 1:ngp
        @views Kspsp[i,j] = J[i,j,i,j]
    end

    Vppppt = Vt[4]
    @tensor Jab[i,j,c,d] := gp[a,i]*Vppppt[a,b,c,d]*gp[b,j]
    @tensor J[i,j,k,l] := gp[c,k]*Jab[i,j,c,d]*gp[d,l]

    Jpppp = Array{Float64}(undef, ngp, ngp)
    Kpppp = Array{Float64}(undef, ngp, ngp)
    @inbounds for i = 1:ngp, j = i:ngp
        @views Jpppp[i,j] = J[i,i,j,j]
        @views Kpppp[i,j] = J[i,j,i,j]
        if i ≠ j
            @views Jpppp[j,i] = Jpppp[i,j]
            @views Kpppp[j,i] = Kpppp[i,j]
        end
    end

    # if ne > 5
        Vppppt = Vt[5]
        @tensor Jab[i,j,c,d] := gp[a,i]*Vppppt[a,b,c,d]*gp[b,j]
        @tensor J[i,j,k,l] := gp[c,k]*Jab[i,j,c,d]*gp[d,l]

        Jpxpxpypy = Array{Float64}(undef, ngp, ngp)
        @inbounds for i = 1:ngp, j = 1:ngp
            @views Jpxpxpypy[i,j] = J[i,i,j,j]
        end

        Vppppt = Vt[6]
        @tensor Jab[i,j,c,d] := gp[a,i]*Vppppt[a,b,c,d]*gp[b,j]
        @tensor J[i,j,k,l] := gp[c,k]*Jab[i,j,c,d]*gp[d,l]

        Kpxpypxpy = Array{Float64}(undef, ngp, ngp)
        @inbounds for i = 1:ngp, j = 1:ngp
            @views Kpxpypxpy[i,j] = J[i,j,i,j]
        end
    # else
    #     Jpxpxpypy = 0
    #     Kpxpypxpy = 0
    # end

    return (Jssss, Jsspp, Jpppp, Jpxpxpypy), (Kssss, Kspsp, Kpppp, Kpxpypxpy)
end

# function get_αβ(neo::Int16, ntype::Int16)
#     αβ = [(neo-1) % ntype + 1, 0]
#     αβ[2] = neo - αβ[1]
#     return sort(αβ, rev=true)
# end

function ne2distribution(ne::Int16)
    n = 1
    l = ['s','p','p','p']
    li = 1
    s = 0
    dist = Array{Tuple}(undef, 0)
    while ne > 0
        ne -= 1
        push!(dist,(n,l[li],li-1,s))
        if s == 0
            if li == 1
                s = 1
            elseif li < 4
                li += 1
            else
                s = 1
                li = 2
            end
        elseif s == 1
            if n == 1
                s = 0
                n += 1
            elseif li == 1
                s = 0
                li += 1
            elseif li < 4
                li += 1
            else
                s = 0
                n += 1
                li = 1 
            end
        end
    end
    return dist
end

function density_prep(
    als::Array{Float64, 1},
    ζ::Int16,
)
    Ht, Vt = HValpha(als,ζ)
    return Ht, Vt
end

function density_prep(
    als::Array{Float64, 1},
    alp::Array{Float64, 1},
    ζ::Int16,
)
    Ht, Vt = HValpha(als,alp,ζ)
    return Ht, Vt
end

function energy(
    als::Array{Float64, 1},
    gs::Array{Float64, 2},
    Ht::Tuple,
    Vt::Tuple,
    ne::Int16
)
    gs = normalized(als, gs)
    gs = orthonormalized(als, gs)
    H = integral2g(gs, Ht)
    J, K = integral4g(gs, Vt)

    dist = ne2distribution(ne)
    ndist = size(dist)[1]

    Hsum = 0.
    Jsum = 0.
    Ksum = 0.
    for i = 1:ndist
        idist = dist[i]
        if idist[2] == 's'
            Hsum += H[1][idist[1]]
        elseif idist[2] == 'p'
            Hsum += H[2][idist[1]-1]
        end
        for j = i+1:ndist
            jdist = dist[j]
            if idist[2] * jdist[2] == "ss"
                Jsum += J[1][idist[1],jdist[1]]
                if idist[4] == jdist[4]
                    Ksum += K[1][idist[1],jdist[1]]
                end
            elseif idist[2] * jdist[2] == "sp"
                Jsum += J[2][idist[1],jdist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[2][idist[1],jdist[1]-1]
                end
            elseif idist[2] * jdist[2] == "ps"
                Jsum += J[2][jdist[1],idist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[2][jdist[1],idist[1]-1]
                end
            elseif idist[2] * jdist[2] == "pp" && idist[3] == jdist[3]
                Jsum += J[3][idist[1]-1,jdist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[3][idist[1]-1,jdist[1]-1]
                end
            elseif idist[2] * jdist[2] == "pp" && idist[3] ≠ jdist[3]
                Jsum += J[4][idist[1]-1,jdist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[4][idist[1]-1,jdist[1]-1]
                end
            end
        end
    end

    return Hsum + Jsum - Ksum

end

function energy(
    als::Array{Float64, 1},
    gs::Array{Float64, 2},
    ζ::Int16,
    ne::Int16
)
    gs = normalized(als, gs)
    gs = orthonormalized(als, gs)
    Ht, Jt = density_prep(als, ζ)
    H = integral2g(gs, Ht)
    J, K = integral4g(gs, Vt)

    dist = ne2distribution(ne)
    ndist = size(dist)[1]

    Hsum = 0.
    Jsum = 0.
    Ksum = 0.
    for i = 1:ndist
        idist = dist[i]
        if idist[2] == 's'
            Hsum += H[1][idist[1]]
        elseif idist[2] == 'p'
            Hsum += H[2][idist[1]-1]
        end
        for j = i+1:ndist
            jdist = dist[j]
            if idist[2] * jdist[2] == "ss"
                Jsum += J[1][idist[1],jdist[1]]
                if idist[4] == jdist[4]
                    Ksum += K[1][idist[1],jdist[1]]
                end
            elseif idist[2] * jdist[2] == "sp"
                Jsum += J[2][idist[1],jdist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[2][idist[1],jdist[1]-1]
                end
            elseif idist[2] * jdist[2] == "ps"
                Jsum += J[2][jdist[1],idist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[2][jdist[1],idist[1]-1]
                end
            elseif idist[2] * jdist[2] == "pp" && idist[3] == jdist[3]
                Jsum += J[3][idist[1]-1,jdist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[3][idist[1]-1,jdist[1]-1]
                end
            elseif idist[2] * jdist[2] == "pp" && idist[3] ≠ jdist[3]
                Jsum += J[4][idist[1]-1,jdist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[4][idist[1]-1,jdist[1]-1]
                end
            end
        end
    end

    return Hsum + Jsum - Ksum

end

function energy(
    als::Array{Float64, 1},
    alp::Array{Float64, 1},
    gs::Array{Float64, 2},
    gp::Array{Float64, 2},
    Ht::Tuple,
    Vt::Tuple,
    ne::Int16
)
    gs, gp = normalized(als, alp, gs, gp)
    gs, gp = orthonormalized(als, alp, gs, gp)
    H = integral2g(gs, gp, Ht)
    J, K = integral4g(gs, gp, Vt)

    dist = ne2distribution(ne)
    ndist = size(dist)[1]

    Hsum = 0.
    Jsum = 0.
    Ksum = 0.
    for i = 1:ndist
        idist = dist[i]
        if idist[2] == 's'
            Hsum += H[1][idist[1]]
        elseif idist[2] == 'p'
            Hsum += H[2][idist[1]-1]
        end
        for j = i+1:ndist
            jdist = dist[j]
            if idist[2] * jdist[2] == "ss"
                Jsum += J[1][idist[1],jdist[1]]
                if idist[4] == jdist[4]
                    Ksum += K[1][idist[1],jdist[1]]
                end
            elseif idist[2] * jdist[2] == "sp"
                Jsum += J[2][idist[1],jdist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[2][idist[1],jdist[1]-1]
                end
            elseif idist[2] * jdist[2] == "ps"
                Jsum += J[2][jdist[1],idist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[2][jdist[1],idist[1]-1]
                end
            elseif idist[2] * jdist[2] == "pp" && idist[3] == jdist[3]
                Jsum += J[3][idist[1]-1,jdist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[3][idist[1]-1,jdist[1]-1]
                end
            elseif idist[2] * jdist[2] == "pp" && idist[3] ≠ jdist[3]
                Jsum += J[4][idist[1]-1,jdist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[4][idist[1]-1,jdist[1]-1]
                end
            end
        end
    end

    return Hsum + Jsum - Ksum

end

function energy(
    als::Array{Float64, 1},
    alp::Array{Float64, 1},
    gs::Array{Float64, 2},
    gp::Array{Float64, 2},
    ζ::Int16,
    ne::Int16
)
    gs, gp = normalized(als, alp, gs, gp)
    gs, gp = orthonormalized(als, alp, gs, gp)
    Ht, Vt = density_prep(als, alp, ζ)
    H = integral2g(gs, gp, Ht)
    J, K = integral4g(gs, gp, Vt)

    dist = ne2distribution(ne)
    ndist = size(dist)[1]

    Hsum = 0.
    Jsum = 0.
    Ksum = 0.
    for i = 1:ndist
        idist = dist[i]
        if idist[2] == 's'
            Hsum += H[1][idist[1]]
        elseif idist[2] == 'p'
            Hsum += H[2][idist[1]-1]
        end
        for j = i+1:ndist
            jdist = dist[j]
            if idist[2] * jdist[2] == "ss"
                Jsum += J[1][idist[1],jdist[1]]
                if idist[4] == jdist[4]
                    Ksum += K[1][idist[1],jdist[1]]
                end
            elseif idist[2] * jdist[2] == "sp"
                Jsum += J[2][idist[1],jdist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[2][idist[1],jdist[1]-1]
                end
            elseif idist[2] * jdist[2] == "ps"
                Jsum += J[2][jdist[1],idist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[2][jdist[1],idist[1]-1]
                end
            elseif idist[2] * jdist[2] == "pp" && idist[3] == jdist[3]
                Jsum += J[3][idist[1]-1,jdist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[3][idist[1]-1,jdist[1]-1]
                end
            elseif idist[2] * jdist[2] == "pp" && idist[3] ≠ jdist[3]
                Jsum += J[4][idist[1]-1,jdist[1]-1]
                if idist[4] == jdist[4]
                    Ksum += K[4][idist[1]-1,jdist[1]-1]
                end
            end
        end
    end

    return Hsum + Jsum - Ksum

end