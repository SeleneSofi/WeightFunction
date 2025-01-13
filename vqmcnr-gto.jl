# Program to run VQMC simulation for Ne atom and optimize boys parâmeters using
# Newton-Raphson method.
#
# v 1.4.9 (1 - features added, 2 - revision, 3 - atom added)
#
### Matrix of ψ are transposed, Slater determinant has rows for electrons and
### columns for orbitals.
#
# Loading packages.
using JLD
using Printf
using Random
using LinearAlgebra
using Plots
using DelimitedFiles
gr()
default(legend = false)
include("quadrature.jl")
include("weight_function.jl")
include("dict_neutrals.jl")
include("energyN.jl")

# Setting general auxiliar functions

function make_wf(
    qt::Array{ETabn},
    x0::Array{nSig},
)
    lnα = [calc_lnα(iqt) for iqt = qt]
    al = [exp.(ilnα) for ilnα = lnα]
    gnt = [Array{Array}(undef, 0) for _ = lnα]
    for ix0 = x0
        push!(gnt[ix0.l], weight_funtion(ix0, lnα[ix0.l]))
    end
    gn = [hcat(ignt...) for ignt = gnt]
    if size(gn)[1] == 1
        gn = normalized(al..., gn...)
        gn = [orthonormalized(al..., gn)]
    else
        gn = [normalized(al..., gn...)...]
        gn = [orthonormalized(al..., gn...)...]
    end
    return al, gn
end

function ne2nαβ(ne, al, gn, num)
    net = copy(ne)
    n = 1
    epn = [
        [0 0 0]
        [1 0 0]
        [0 1 0]
        [0 0 1]
    ]
    li = 1
    s = 0
    nα, nβ = 0, 0
    np = 0
    ep = Array{Int}(undef, (0,3))
    ζ = Array{Float64}(undef, 0)
    cαp = Array{Int}(undef, (1,0))
    cβp = Array{Int}(undef, (1,0))
    while net > 0
        net -= 1
        lit = li == 1 ? 1 : 2
        nsp = n-(lit-1)
        if s == 0
            nα += 1
            np += num[lit]
            ep = [ep; repeat(epn[li,:]', num[lit])]
            ζ = [ζ; al[lit]]
            if size(cαp,2) == 0
                cαp = [cαp gn[lit][:,nsp]']
            else
                cαp = [
                    cαp                     zeros(size(cαp,1),num[lit])
                    zeros(1,size(cαp,2))    gn[lit][:,nsp]'
                ]
            end
        else
            nβ += 1
            if size(cβp,2) == 0
                cβp = [cβp gn[lit][:,nsp]']
            else
                cβp = [
                    cβp                     zeros(size(cβp,1),num[lit])
                    zeros(1,size(cβp,2))    gn[lit][:,nsp]'
                ]
            end
        end
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
    cpdiff = size(cαp,2) - size(cβp,2)
    if cpdiff > 0
        cβp = [cβp zeros(size(cβp,1),cpdiff)]
    elseif cpdiff < 0
        cαp = [cαp zeros(size(cαp,1),cpdiff)]
    end
    return nα, nβ, np, ep, ζ, cαp, cβp
end

function dfact(n::Int)
    n ≥ -1 || error("n must be greater than -2")
    -1 ≤ n ≤ 0 && return 1
    n * dfact(n-2)
end

function sumd(A::AbstractArray, dims)
    return dropdims(sum(A, dims=dims), dims=dims)
end

function prodd(A::AbstractArray, dims)
    return dropdims(prod(A, dims=dims), dims=dims)
end

function printprg(en::Float64, σen::Float64, acc::Float64, τr::Float64, ict::Int, nrun::Int)
    pct = round(100ict/nrun, digits=1)
    pi = convert(Int,round(50ict/nrun))
    if ict == 100
        @printf("E = %.8f ± %.8f Eh\n", en, σen)
        # println("E = $en ± $σen")
        @printf("Step = %i\n", ict)
        # println("Step = $ict")
        @printf("Acceptance = %.4f\n", acc)
        # println("Acceptance = $acc")
        @printf("τ = %.4f\n", τr)
        # println("τ = $τr")
        println("[" * "="^pi * " "^(50-pi) * "] $(@sprintf("%5.1f",pct)) %")
        # println("[" * " "^50 * "] $(@sprintf("%5.1f",pct)) %")
    else
        # pi = convert(Int,round(50ict/nrun))
        print("\u1b[5F")
        @printf("E = %.8f ± %.8f Eh\n", en, σen)
        # println("E = $en ± $σen")
        @printf("Step = %i\n", ict)
        # println("Step = $ict")
        @printf("Acceptance = %.4f\n", acc)
        # println("Acceptance = $acc")
        @printf("τ = %.4f\n", τr)
        # println("τ = $τr")
        println("[" * "="^pi * " "^(50-pi) * "] $(@sprintf("%5.1f",pct)) %")
        # println("[" * "="^pi * " "^(50-pi) * "] $pct %")
        print("\u1b[0K")
    end
end

# Setting structures
struct QMCstr
    nw::Int
    seed::Int
    tnstp::Int
    nstp::Int
    τ::Float64
    acc::Float64
    nrstp::Int
    Δlim::Float64
end

struct GTOstr
    base::String
    zn::Int
    nα::Int
    nβ::Int
    ne::Int
    np::Int
    ep::Array{Int,2}
    ζ::Array{Float64,1}
    # iap::Array{Int,1}
    cαp::Array{Float64,2}
    cβp::Array{Float64,2}
end

struct Corstr
    c9::Array{Float64,1}
    c9i::BitArray{1}
end

mutable struct Ψstr
    ψ::Array{Float64,1}
    ψs::Array{Array{Float64,3},1}
    ∇ψs::Array{Array{Float64,4},1}
    ∇²ψs::Array{Array{Float64,3},1}
    ψ⁻¹s::Array{Array{Float64,3},1}
end

struct NRstr
    E::Array{Float64,1}
    σ::Array{Float64,1}
    Ek::Array{Float64,1}
    Vne::Array{Float64,1}
    Vee::Array{Float64,1}
    Vt::Array{Float64,1}
    c9::Array{Float64,2}
    fm::Array{Float64,2}
    hmn::Array{Float64,3}
end

struct NRμstr
    μ∂lψ∂m::Array{Float64,2}
    μ∂el∂m::Array{Float64,2}
    μel∂lψ∂m::Array{Float64,2}
    μ∂²lψ∂mn::Array{Float64,3}
    μel∂²lψ∂mn::Array{Float64,3}
    μ∂el∂m∂lψ∂n::Array{Float64,3}
end

# Setting especific auxiliar functions
function qmcinp(iask::Bool, icor::Bool, atom::Int, ion::Int, spin::Char, base::String)
    if iask
        println("Enter number of thermalization steps:")
        tnstp = parse(Int,readline())
        println("Enter number of VQMC steps:")
        nstp = parse(Int,readline())
        println("Enter number of NR steps:")
        nrstp = parse(Int,readline())
        println("Enter seed number:")
        seed = parse(Int,readline())
        Random.seed!(seed)
        qmcp = QMCstr(
            100, # nw
            seed, # seed
            tnstp, # tnstp
            nstp, # nstp
            1., # τ
            .5, # acc
        nrstp, # nrstp
        1. # Δlim
        )
    end
    # atom = 10
    # ion = 0
    # spin = 'α'
    # spin = 'β'
    # spin = 'n'
    _base = base
    if base == "STO-3G"
        if atom == 2 # He
            _zn = 2
            _nα = 1
            _nβ = 1
            _ne = 2
            _np = 3
            _ep = [ 0 0 0;
                    0 0 0;
                    0 0 0;]
            _ζ = [  6.362421394;
                    1.158922999;
                    0.3136497915;]
            _cαp = [ 0.1543289673  0.5353281423  0.4446345422;]
            _cβp = [ 0.1543289673  0.5353281423  0.4446345422;]
        elseif atom == 17 # Cl
            _zn = 17
            _nα = 9
            _nβ = 8
            _ne = 17
            _np = 27
            _ep = [ 0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    0 0 0;
                    1 0 0;
                    1 0 0;
                    1 0 0;
                    0 1 0;
                    0 1 0;
                    0 1 0;
                    0 0 1;
                    0 0 1;
                    0 0 1;
                    1 0 0;
                    1 0 0;
                    1 0 0;
                    0 1 0;
                    0 1 0;
                    0 1 0;
                    0 0 1;
                    0 0 1;
                    0 0 1;]
            _ζ = [ 601.3456136;
                   109.5358542;
                    29.64467686;
                    38.96041889;
                     9.053563477;
                     2.944499834;
                     2.129386495;
                     0.5940934274;
                     0.2325241410;
                    38.96041889;
                     9.053563477;
                     2.944499834;
                    38.96041889;
                     9.053563477;
                     2.944499834;
                    38.96041889;
                     9.053563477;
                     2.944499834;
                     2.129386495;
                     0.5940934274;
                     0.2325241410;
                     2.129386495;
                     0.5940934274;
                     0.2325241410;
                     2.129386495;
                     0.5940934274;
                     0.2325241410;]
            _cαp = [ 0.1543289673   0.5353281423   0.4446345422   0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000 -0.09996722919  0.3995128261   0.7001154689   0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000 -0.2196203690   0.2255954336   0.9003984260   0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.1559162750   0.6076837186   0.3919573931   0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.1559162750   0.6076837186   0.3919573931   0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.1559162750   0.6076837186   0.3919573931   0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.01058760429  0.5951670053   0.4620010120   0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.01058760429  0.5951670053   0.4620010120   0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.01058760429  0.5951670053   0.4620010120 ;]
            _cβp = [ 0.1543289673   0.5353281423   0.4446345422   0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000 -0.09996722919  0.3995128261   0.7001154689   0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000 -0.2196203690   0.2255954336   0.9003984260   0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.1559162750   0.6076837186   0.3919573931   0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.1559162750   0.6076837186   0.3919573931   0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.1559162750   0.6076837186   0.3919573931   0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.01058760429  0.5951670053   0.4620010120   0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000;
                     0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.00000000000  0.01058760429  0.5951670053   0.4620010120   0.00000000000  0.00000000000  0.00000000000;]
        end
    elseif base == "wf-et-50s50p"
        num = [50, 50]
        _zn = atom
        _ne = copy(_zn)
        if ion ≠ 0
            _ne -= ion
        end
        qt = database[_zn][_ne]["qt"]
        x0 = database[_zn][_ne]["x"]
        al, gn = make_wf(qt, x0)
        _nα, _nβ, _np, _ep, _ζ, _cαp, _cβp = ne2nαβ(_ne, al, gn, num)
    end
    if ion ≠ 0 && base ≠ "wf-et-50s50p" # 1s¹2s²2p⁴
        _ne -= 1
        if spin == 'α'
            _nα -= 1
            _cαp = _cαp[1:end .≠ ion, :]
        else
            _nβ -= 1
            _cβp = _cβp[1:end .≠ ion, :]
        end
    # end
    _cαp .*= sqrt.(2 .^(2sumd(_ep,2).+3/2) .* _ζ.^(sumd(_ep,2).+3/2) ./ (pi.^(3/2)*prodd(dfact.(2_ep .-1),2)))'
    if size(_cβp,2) ≠ 0
        _cβp .*= sqrt.(2 .^(2sumd(_ep,2).+3/2) .* _ζ.^(sumd(_ep,2).+3/2) ./ (pi.^(3/2)*prodd(dfact.(2_ep .-1),2)))'
    end
end
    gto = GTOstr(
            _base, # base
            _zn, # zn
            _nα, # nα
            _nβ, # nβ
            _ne, # ne
            _np, # np
            _ep, # ep
            _ζ, # ζ
            _cαp, # cαp
            _cβp  # cβp
        )
    if _ne == 1
        corp = Corstr(
            [ 0.0;  0.0;  0.0;  0.0;  # ee
              0.0;  0.0;  0.0; # ne
              0.0;  0.0 ], # nee e nne
            [ 0; 0; 0; 0;
              0; 0; 0;
              0; 0 ]
        )
    elseif icor
        corp = Corstr(
            [ 0.5;  0.0;  0.0;  0.0;  # ee
              0.0;  0.0;  0.0; # ne
              0.0;  0.0 ], # nee e nne
            [ 0; 1; 1; 1;
              1; 1; 1;
              1; 1 ]
        )
    end
    if iask && icor
        return qmcp, gto, corp
    elseif iask
        return qmcp, gto
    elseif icor
        return gto, corp
    else
        return gto
    end
end

@fastmath @views function gtof(iαβ::Int, rci::Array{Float64,2}, ri::Array{Float64,1}, gto::GTOstr)
    np, ep, ζ = gto.np, copy(gto.ep), gto.ζ
    nw = size(ri,1)
    if iαβ == 1
        no, csp = gto.nα, gto.cαp
    else
        no, csp = gto.nβ, gto.cβp
    end
    ri² = ri.^2
    dpw = Array{Float64}(undef, np, nw)
    gpw = Array{Float64}(undef, np, nw, 3)
    lpw = Array{Float64}(undef, np, nw)
    mep123 = maximum(ep[:,1:3])
    if mep123 ≠ 0
        rci⁻¹ = rci.^-1
        if mep123 ≠ 1
            rci⁻² = rci.^-2
            for i = 2:mep123
                rcie[:,:,i-1] .= rci .^ i
            end
        end
    end
    @inbounds @simd for p = 1:np
        dpw[p,:] .= exp.(-ζ[p] .* ri²)
        gpw[p,:,:] .= -2ζ[p] .* rci
        lpw[p,:] .= 4ζ[p].^2 .* ri² .- 6ζ[p]
        @simd for i = 1:3
            if ep[p,i] ≠ 0
                gpw[p,:,i] .+= ep[p,i] .* rci⁻¹[:,i]
                lpw[p,:] .-= 4ζ[p] .* ep[p,i]
                if ep[p,i] == 1
                    dpw[p,:] .*= rci[:,i]
                else
                    dpw[p,:] .*= rcie[:,i,ep[p,i]-1]
                    lpw[p,:] .+= ep[p,i] .* (ep[p,i] - 1) .* rci⁻²[:,i]
                end
            end
        end
    end
    gpw .*= dpw
    lpw .*= dpw
    # ψi = csp * dpw
    ψi = dpw' * csp'
    ∇ψi = Array{Float64}(undef, (size(ri,1),no,3))
    ∇ψi[:,:,1] .= gpw[:,:,1]' * csp'
    ∇ψi[:,:,2] .= gpw[:,:,2]' * csp'
    ∇ψi[:,:,3] .= gpw[:,:,3]' * csp'
    ∇²ψi = lpw' * csp'
    return ψi, ∇ψi, ∇²ψi
end

### Matrix of ψ are transposed, Slater determinant has rows for electrons and
### columns for orbitals.
function wavef(rc::Array{Float64,3}, gto::GTOstr, c9r::Array{Float64,1})
    nw = size(rc,1)
    nα, nβ, ne = gto.nα, gto.nβ, gto.ne
    ψ = ones(nw)
    ψs = [zeros(nw,nα,nα), zeros(nw,nβ,nβ)]
    ∇ψs = [zeros(nw,nα,nα,3), zeros(nw,nβ,nβ,3)]
    ∇²ψs = [zeros(nw,nα,nα), zeros(nw,nβ,nβ)]
    ψ⁻¹s = [zeros(nw,nα,nα), zeros(nw,nβ,nβ)]
    r2 = rc[:,:,1].^2 .+ rc[:,:,2].^2 .+ rc[:,:,3].^2
    r = sqrt.(r2)
    rb = r./(1 .+ r)
    rb2 = rb.^2
    rb3 = rb.^3
    rb4 = rb.^4
    isp = 1
    iαβ = 0
    for i = 1:ne
        i ≠ nα + 1 || ((isp, iαβ) = (2, nα))
        for j = i+1:ne
            # if (i ≤ nα && j ≤ nα) || (i > nα && j > nα)
            #     ees = 0.5
            # else
                ees = 1.0
            # end
            rij = sqrt.((rc[:,i,1] - rc[:,j,1]).^2 .+ (rc[:,i,2] - rc[:,j,2]).^2 .+ (rc[:,i,3] - rc[:,j,3]).^2)
            rbij = rij ./ (1 .+ rij)
            rbij2 = rbij.^2
            rbij3 = rbij.^3
            rbij4 = rbij.^4
            corr = c9r[1] * ees * rbij
            corr += c9r[2] * rbij2
            corr += c9r[3] * rbij3
            corr += c9r[4] * rbij4
            corr += c9r[5] * (rb2[:,i] + rb2[:,j])
            corr += c9r[6] * (rb3[:,i] + rb3[:,j])
            corr += c9r[7] * (rb4[:,i] + rb4[:,j])
            corr += c9r[8] * (rb2[:,i] .* rb2[:,j])
            corr += c9r[9] * ((rb2[:,i] + rb2[:,j]) .* rbij2)
            ψ .*= exp.(corr)
        end
        ψs[isp][:,i-iαβ,:], ∇ψs[isp][:,i-iαβ,:,:], ∇²ψs[isp][:,i-iαβ,:] = gtof(isp,rc[:,i,:],r[:,i],gto)
    end
    for w = 1:nw
        ψ[w] *= det(ψs[1][w,:,:]) * det(ψs[2][w,:,:])
        ψ⁻¹s[1][w,:,:] = inv(ψs[1][w,:,:])
        nβ ≠ 0 && (ψ⁻¹s[2][w,:,:] = inv(ψs[2][w,:,:]))
    end
    return Ψstr(ψ, ψs, ∇ψs, ∇²ψs, ψ⁻¹s)
end

### Matrix of ψ are transposed, Slater determinant has rows for electrons and
### columns for orbitals.
@fastmath function wavefup(i::Int, rci::Array{Float64,2}, ψ::Ψstr, rc::Array{Float64,3},
                 gto::GTOstr, c9r::Array{Float64,1})
    @views begin
        nw = size(rc,1)
        nα, ne = gto.nα, gto.ne
        nei = 1:ne .≠ i
        ψn = copy(ψ.ψ)
        r2 = rc[:,:,1].^2 .+ rc[:,:,2].^2 .+ rc[:,:,3].^2
        r = sqrt.(r2)
        rb = r./(1 .+ r)
        rb2 = rb.^2
        rb3 = rb .* rb2
        rb4 = rb2.^2
        rn2 = rci[:,1].^2 .+ rci[:,2].^2 .+ rci[:,3].^2
        rn = sqrt.(rn2)
        rbn = rn./(1 .+ rn)
        rbn2 = rbn.^2
        rbn3 = rbn .* rbn2
        rbn4 = rbn2.^2
    end
    rij = sqrt.((rc[:,i,1] .- rc[:,nei,1]).^2 .+ (rc[:,i,2] .- rc[:,nei,2]).^2 .+ (rc[:,i,3] .- rc[:,nei,3]).^2)
    @views begin
        rbij = rij ./ (1 .+ rij)
        rbij2 = rbij.^2
        rbij3 = rbij .* rbij2
        rbij4 = rbij2.^2
    end
    rinj = sqrt.((rci[:,1] .- rc[:,nei,1]).^2 .+ (rci[:,2] .- rc[:,nei,2]).^2 .+ (rci[:,3] .- rc[:,nei,3]).^2)
    @views begin
        rbinj = rinj ./ (1 .+ rinj)
        rbinj2 = rbinj.^2
        rbinj3 = rbinj .* rbinj2
        rbinj4 = rbinj2.^2
        ees = 1.0
        corr = c9r[1] .* ees .* (rbinj .- rbij)
        corr .+= c9r[2] .* (rbinj2 .- rbij2)
        corr .+= c9r[3] .* (rbinj3 .- rbij3)
        corr .+= c9r[4] .* (rbinj4 .- rbij4)
        corr .+= c9r[8] .* (rbn2 .- rb2[:,i]) .* rb2[:,nei]
        corr .+= c9r[9] .* ( (rbn2 .+ rb2[:,nei]) .* rbinj2 .- (rb2[:,i] .+ rb2[:,nei]) .* rbij2 )
        corr = sumd(corr, 2)
        corr .+= (ne-1) .* c9r[5] .* (rbn2 .- rb2[:,i])
        corr .+= (ne-1) .* c9r[6] .* (rbn3 .- rb3[:,i])
        corr .+= (ne-1) .* c9r[7] .* (rbn4 .- rb4[:,i])
        ψn .*= exp.(corr)
        if i ≤ nα
            ns = nα
            isp = 1
            iαβ = 0
        else
            ns = ne - nα
            isp = 2
            iαβ = nα
        end
        ψsi, ∇ψsi, ∇²ψsi = gtof(isp,rci,rn,gto)
        rno = zeros(nw)
        @inbounds for j = 1:ns
            rno .+= ψsi[:,j] .* ψ.ψ⁻¹s[isp][:,j,i-iαβ]
        end
        ψn .*= rno
        ψst = copy(ψ.ψs[isp])
        ψst[:,i-iαβ,:] = ψsi
        ψ⁻¹si = Array{Float64}(undef, nw, ns, ns)
        nei = 1:ns .≠ i-iαβ
        @inbounds @simd for ii = 1:ns
            if ii ≠ i-iαβ
                for jj = 1:ns
                    ψ⁻¹ψ = zeros(nw)
                    for kk = 1:ns
                        ψ⁻¹ψ .+= ψ.ψ⁻¹s[isp][:,jj,i-iαβ] .* ψsi[:,kk] .* ψ.ψ⁻¹s[isp][:,kk,ii]
                    end
                    ψ⁻¹si[:,jj,ii] .= ψ.ψ⁻¹s[isp][:,jj,ii] .- ψ⁻¹ψ ./ rno
                end
            else
                ψ⁻¹si[:,:,ii] .= ψ.ψ⁻¹s[isp][:,:,i-iαβ] ./ rno
            end
        end
    end
    return ψn, ψsi, ∇ψsi, ∇²ψsi, ψ⁻¹si
end

@fastmath @views function elocal(ψ::Ψstr, gto::GTOstr, rc::Array{Float64,3}, c9r::Array{Float64,1})
    nw = size(ψ.ψ,1)
    zn, nα, ne = gto.zn, gto.nα, gto.ne
    ek = zeros(nw)
    vne = zeros(nw)
    vee = zeros(nw)
    ∂lψ∂m = zeros(nw, 9)
    ∂el∂m = zeros(nw, 9)
    rcij = Array{Float64}(undef, nw, 3)
    rij2 = Array{Float64}(undef, nw)
    rij = Array{Float64}(undef, nw)
    rij3 = Array{Float64}(undef, nw)
    rij4 = Array{Float64}(undef, nw)
    rdij = Array{Float64}(undef, nw)
    rdij2 = Array{Float64}(undef, nw)
    rdij3 = Array{Float64}(undef, nw)
    rdij4 = Array{Float64}(undef, nw)
    rdij5 = Array{Float64}(undef, nw)
    rdij6 = Array{Float64}(undef, nw)
    ririj = Array{Float64}(undef, nw)
    ∇ψc = Array{Float64}(undef, nw, 3)
    r2 = rc[:,:,1].^2 .+ rc[:,:,2].^2 .+ rc[:,:,3].^2
    r = sqrt.(r2)
    r3 = r .* r2
    r4 = r2.^2
    rd = (1 .+ r).^-1
    rd2 = rd.^2
    rd3 = rd .* rd2
    rd4 = rd2.^2
    rd5 = rd2 .* rd3
    rd6 = rd3.^2
    isp = 1
    iαβ = 0
    no = nα
    @inbounds for i = 1:ne
        i ≠ nα + 1 || ((isp, iαβ, no) = (2, nα, ne-nα))
        # Setting up determinant and correlation gradients.
        ∇ψd = zeros(nw,3)
        # Calculating gradient and kinetic energy for correlation,
        # parameters 5, 6, 7.
        ∇ψc .= (ne-1) .* rc[:,i,:] .* ( 2 .* c9r[5] .* rd3[:,i] .+
                                       3 .* c9r[6] .* r[:,i] .* rd4[:,i] .+
                                       4 .* c9r[7] .* r2[:,i] .* rd5[:,i] )
        ek .-= 0.5 .* (ne-1) .* ( 6 .* c9r[5] .* rd4[:,i] .+
                                 12 .* c9r[6] .* r[:,i] .* rd5[:,i] .+
                                 20 .* c9r[7] .* r2[:,i] .* rd6[:,i] )
        # ek[w] -= 0.5(ne-1) * ( 6c9r[5] / rdi4 +
        #                        12c9r[6] * ri/rdi5 +
        #                        20c9r[7] * ri2/rdi6 )
        # Setting up ∂ln(ψ)/∂cₘ gradient
        ∇∂lψ∂m = zeros(nw,3,9)
        # Calculating ∂ln(ψ)/∂cₘ gradient and ∂eₗ/∂cₘ,
        # parameters 5, 6, 7.
        ∇∂lψ∂m[:,:,5] .+= (ne-1) .* 2 .* rc[:,i,:] .* rd3[:,i]
        ∇∂lψ∂m[:,:,6] .+= (ne-1) .* 3 .* rc[:,i,:] .* r[:,i] .* rd4[:,i]
        ∇∂lψ∂m[:,:,7] .+= (ne-1) .* 4 .* rc[:,i,:] .* r2[:,i] .* rd5[:,i]
        ∂el∂m[:,5] .-= 0.5 .* (ne-1) .* 6 .* rd4[:,i]
        ∂el∂m[:,6] .-= 0.5 .* (ne-1) .* 12 .* r[:,i] .* rd5[:,i]
        ∂el∂m[:,7] .-= 0.5 .* (ne-1) .* 20 .* r2[:,i] .* rd6[:,i]
        # Setting ij variables
        #r2 = sumd(rc.^2, 3)
        #r = sqrt.(r2)
        for j = 1:no
            ∇ψd .+= ψ.ψ⁻¹s[isp][:,j,i-iαβ] .* ψ.∇ψs[isp][:,i-iαβ,j,:]
            ek .-= 0.5 .* ψ.ψ⁻¹s[isp][:,j,i-iαβ] .* ψ.∇²ψs[isp][:,i-iαβ,j]
        end
        @inbounds for j = 1:ne
            if j == i
                continue
            end
            # Differentiating like and unlik spin electrons.
            # if (i ≤ nα && j ≤ nα) || (i > nα && j > nα)
            #     ees = 0.5
            # else
                ees = 1.0
            # end
            # Setting up electrons ij variables
            rcij .= rc[:,i,:] .- rc[:,j,:]
            rij2 .= rcij[:,1].^2 .+ rcij[:,2].^2 .+ rcij[:,3].^2
            rij .= sqrt.(rij2)
            rij3 .= rij .* rij2
            rij4 .= rij2.^2
            rdij .= (1 .+ rij).^-1
            rdij2 .= rdij.^2
            rdij3 .= rdij .* rdij2
            rdij4 .= rdij2.^2
            rdij5 .= rdij2 .* rdij3
            rdij6 .= rdij3.^2
            ririj .= rc[:,i,1] .* rcij[:,1] .+ rc[:,i,2] .* rcij[:,2] .+ rc[:,i,3] .* rcij[:,3]
            # Calculating gradient and kinetic energy for correlation,
            # parameters 1, 2, 3, 4, 8, 9.
            ∇ψc .+= rc[:,i,:] .* ( 2c9r[8] .* r2[:,j] .* rd2[:,j] .* rd3[:,i] .+
                                   2c9r[9] .* rij2 .* rdij2 .* rd3[:,i] )
            ∇ψc .+= rcij .* ( c9r[1] .* ees .* rdij2 ./ rij .+
                              2c9r[2] .* rdij3 .+
                              3c9r[3] .* rij .* rdij4 .+
                              4c9r[4] .* rij2 .* rdij5 .+
                              2c9r[9] .* (r2[:,i].*rd2[:,i] .+ r2[:,j].*rd2[:,j]) .* rdij3 )
            ek  .-= 0.5( 2c9r[1] .* ees .* rdij3 ./ rij .+
                         6c9r[2] .* rdij4 .+
                        12c9r[3] .* rij .* rdij5 .+
                        20c9r[4] .* rij2 .* rdij6 .+
                         6c9r[8] .* r2[:,j] .* rd2[:,j] .* rd4[:,i] .+
                         2c9r[9] .* ( 3 .* rij2 .* rdij2 .* rd4[:,i] .+
                                      4 .* ririj .* rd3[:,i] .* rdij3 .+
                                      3 .* (r2[:,i].*rd2[:,i] .+ r2[:,j].*rd2[:,j]) .* rdij4 ))
            # Calculating ∂ln(ψ)/∂cₘ
            ∂lψ∂m[:,1] .+= 0.5 .* rij .* rdij
            ∂lψ∂m[:,2] .+= 0.5 .* rij2 .* rdij2
            ∂lψ∂m[:,3] .+= 0.5 .* rij3 .* rdij3
            ∂lψ∂m[:,4] .+= 0.5 .* rij4 .* rdij4
            ∂lψ∂m[:,5] .+= 0.5 .* (r2[:,i].*rd2[:,i] .+ r2[:,j].*rd2[:,j])
            ∂lψ∂m[:,6] .+= 0.5 .* (r3[:,i].*rd3[:,i] .+ r3[:,j].*rd3[:,j])
            ∂lψ∂m[:,7] .+= 0.5 .* (r4[:,i].*rd4[:,i] .+ r4[:,j].*rd4[:,j])
            ∂lψ∂m[:,8] .+= 0.5 .* r2[:,i].*rd2[:,i] .* r2[:,j].*rd2[:,j]
            ∂lψ∂m[:,9] .+= 0.5 .* (r2[:,i].*rd2[:,i] .+ r2[:,j].*rd2[:,j]) .* rij2.*rdij2
            # Calculating ∂ln(ψ)/∂cₘ gradient and ∂eₗ/∂cₘ,
            # parameters 1, 2, 3, 4, 8, 9.
            ∇∂lψ∂m[:,:,1] .+= rcij .* rdij2 ./ rij
            ∇∂lψ∂m[:,:,2] .+= rcij .* 2 .* rdij3
            ∇∂lψ∂m[:,:,3] .+= rcij .* 3 .* rij .* rdij4
            ∇∂lψ∂m[:,:,4] .+= rcij .* 4 .* rij2 .* rdij5
            ∇∂lψ∂m[:,:,8] .+= rc[:,i,:] .* 2 .* r2[:,j] .* rd2[:,j] .* rd3[:,i]
            ∇∂lψ∂m[:,:,9] .+= rc[:,i,:] .* 2 .* rij2 .* rdij2 .* rd3[:,i] .+
                              rcij .* 2 .* (r2[:,i].*rd2[:,i] .+ r2[:,j].*rd2[:,j]) .* rdij3
            ∂el∂m[:,1] .-= 0.5 *  2 .* rdij3 ./ rij
            ∂el∂m[:,2] .-= 0.5 *  6 .* rdij4
            ∂el∂m[:,3] .-= 0.5 * 12 .* rij .* rdij5
            ∂el∂m[:,4] .-= 0.5 * 20 .* rij2 .* rdij6
            ∂el∂m[:,8] .-= 0.5 *  6 .* r2[:,j] .* rd2[:,j] .* rd4[:,i]
            ∂el∂m[:,9] .-= 0.5 *  2 .* ( 3 .* rij2 .* rdij2 .* rd4[:,i] .+
                                         4 .* ririj .* rd3[:,i] .* rdij3 .+
                                         3 .* (r2[:,i].*rd2[:,i] .+ r2[:,j].*rd2[:,j]) .* rdij4 )
            # Calculating electron-electron repulsion energy
            vee .+= 0.5 .* rij.^-1
        end
        # Calculating ∂eₗ/∂cₘ for gradients of determinant, correlation, eₗ.
        @inbounds for k = 1:3
            ∂el∂m[:,1] .-= 0.5 * 2 .* (∇ψd[:,k] .+ ∇ψc[:,k]) .* ∇∂lψ∂m[:,k,1]
            ∂el∂m[:,2] .-= 0.5 * 2 .* (∇ψd[:,k] .+ ∇ψc[:,k]) .* ∇∂lψ∂m[:,k,2]
            ∂el∂m[:,3] .-= 0.5 * 2 .* (∇ψd[:,k] .+ ∇ψc[:,k]) .* ∇∂lψ∂m[:,k,3]
            ∂el∂m[:,4] .-= 0.5 * 2 .* (∇ψd[:,k] .+ ∇ψc[:,k]) .* ∇∂lψ∂m[:,k,4]
            ∂el∂m[:,5] .-= 0.5 * 2 .* (∇ψd[:,k] .+ ∇ψc[:,k]) .* ∇∂lψ∂m[:,k,5]
            ∂el∂m[:,6] .-= 0.5 * 2 .* (∇ψd[:,k] .+ ∇ψc[:,k]) .* ∇∂lψ∂m[:,k,6]
            ∂el∂m[:,7] .-= 0.5 * 2 .* (∇ψd[:,k] .+ ∇ψc[:,k]) .* ∇∂lψ∂m[:,k,7]
            ∂el∂m[:,8] .-= 0.5 * 2 .* (∇ψd[:,k] .+ ∇ψc[:,k]) .* ∇∂lψ∂m[:,k,8]
            ∂el∂m[:,9] .-= 0.5 * 2 .* (∇ψd[:,k] .+ ∇ψc[:,k]) .* ∇∂lψ∂m[:,k,9]
            ek .-= (∇ψd[:,k] .+ 0.5 .* ∇ψc[:,k]) .* ∇ψc[:,k]
        end
        # Calculating kinetic energy for gradients of determinant, correlation.
        # Calculating nuclear-electron attraction energy
        vne .-= zn .* r[:,i].^-1
    end
    el = ek .+ vne .+ vee
    return (el, ek, vne, vee, ∂lψ∂m, ∂el∂m)
end

@fastmath @views function vqmcnr(qmcp::QMCstr, gto::GTOstr, corp::Corstr)
    # Setting initial seed
    nw    = qmcp.nw
    tnstp = qmcp.tnstp
    nstp  = qmcp.nstp
    τ     = qmcp.τ
    acc   = qmcp.acc
    nrstp = qmcp.nrstp
    Δlim  = qmcp.Δlim
    c9, c9i = corp.c9, corp.c9i

    Enrf = Array{Float64}(undef, nrstp)
    σnrf = Array{Float64}(undef, nrstp)
    Ekf  = Array{Float64}(undef, nrstp)
    Vnef = Array{Float64}(undef, nrstp)
    Veef = Array{Float64}(undef, nrstp)
    Vtf  = Array{Float64}(undef, nrstp)
    c9f  = Array{Float64}(undef, 9, nrstp)
    fmf  = Array{Float64}(undef, 9, nrstp)
    hmnf = Array{Float64}(undef, 9, 9, nrstp)
    μ∂lψ∂mf      = Array{Float64}(undef, 9, nrstp)
    μ∂el∂mf      = Array{Float64}(undef, 9, nrstp)
    μel∂lψ∂mf    = Array{Float64}(undef, 9, nrstp)
    μ∂²lψ∂mnf    = Array{Float64}(undef, 9, 9, nrstp)
    μel∂²lψ∂mnf  = Array{Float64}(undef, 9, 9, nrstp)
    μ∂el∂m∂lψ∂nf = Array{Float64}(undef, 9, 9, nrstp)

    c9   = corp.c9
    τ    = qmcp.τ
    acc = qmcp.acc
    nα, ne = gto.nα, gto.ne
    inr  = 0
    swt  = 1
    τr   = 0.0 # only to initialize as Float64
    accr = 0.0 # only to initialize as Float64
    c9r  = copy(c9)
    rc   = Array{Float64}(undef, nw, ne, 3)
    # default(legend = false)
    # Starting main loop
    while inr ≠ nrstp
        # Prepering thermalization or main calculation
        if swt == 1
            nrun = tnstp
            τr   = copy(τ)
            accr = copy(acc)
            rc   = randn(Float64, nw, ne, 3)
            println("-"^60)
            println("")
            println("Starting Thermalization of Step $(inr+1)")
            println("")
        else
            nrun = nstp
            inr += 1
            println("-"^60)
            println("")
            println("Starting VQMC of Step $inr")
            println("")
        end

        # aict = Array{Float64}(undef, div(nrun,100))
        # Enct = Array{Float64}(undef, div(nrun,100))
        # First wave and energy estimations
        ψo = wavef(rc,gto,c9r)
        el, ek, vne, vee, ∂lψ∂m, ∂el∂m = elocal(ψo, gto, rc, c9r)
        # Intial cycle VQMC variables
        μel = sum(el)/nw
        μel² = sum(el.^2)/nw
        μek = sum(ek)/nw
        μvne = sum(vne)/nw
        μvee = sum(vee)/nw
        ict = 1
        en = μel/ict
        σen = sqrt( (μel²/ict - en.^2) / (ict*nw - 1) )
        ekt = μek/ict
        vnet = μvne/ict
        veet = μvee/ict
        vt = vnet + veet
        if swt ≠ 1
            # Setting initil variables of VMQ cicles
            μ∂lψ∂m = sumd(∂lψ∂m, 1)./nw
            μ∂el∂m = sumd(∂el∂m, 1)./nw
            μel∂lψ∂m = sumd(el .* ∂lψ∂m, 1) ./ nw
            μ∂²lψ∂mn = ∂lψ∂m' * ∂lψ∂m ./ nw
            μel∂²lψ∂mn = (el.*∂lψ∂m)' * ∂lψ∂m ./ nw
            μ∂el∂m∂lψ∂n = (∂el∂m' * ∂lψ∂m .+ ∂lψ∂m' * ∂el∂m) ./ nw
        end
        while ict ≠ nrun
            ict += 1
            nacc = 0
            isp = 1
            iαβ = 0
            @inbounds for i = 1:ne
                i ≠ nα + 1 || ((isp, iαβ) = (2, nα))
                rci = rc[:,i,:] .+ τr.*(0.5 .- rand(nw,3))
                ψi, ψsi, ∇ψsi, ∇²ψsi, ψ⁻¹si = wavefup(i,rci,ψo,rc,gto,c9r)
                ψtrans = ψi.^2 ./ ψo.ψ.^2
                iacc = ψtrans .≥ rand(nw)
                nacc += sum(iacc)/ne
                rc[iacc,i,:] = rci[iacc,:]
                ψo.ψ[iacc] = ψi[iacc]
                ψo.ψs[isp][iacc,i-iαβ,:] = ψsi[iacc,:]
                ψo.∇ψs[isp][iacc,i-iαβ,:,:] = ∇ψsi[iacc,:,:]
                ψo.∇²ψs[isp][iacc,i-iαβ,:] = ∇²ψsi[iacc,:]
                ψo.ψ⁻¹s[isp][iacc,:,:] = ψ⁻¹si[iacc,:,:]
            end
            acc = nacc/nw
            el, ek, vne, vee, ∂lψ∂m, ∂el∂m = elocal(ψo, gto, rc, c9r)
            μel += sum(el)/nw
            μel² += sum(el.^2)/nw
            μek += sum(ek)/nw
            μvne += sum(vne)/nw
            μvee += sum(vee)/nw
            if swt ≠ 1
                μ∂lψ∂m .+= sumd(∂lψ∂m, 1)./nw
                μ∂el∂m .+= sumd(∂el∂m, 1)./nw
                μel∂lψ∂m .+= sumd(el .* ∂lψ∂m, 1) ./ nw
                μ∂²lψ∂mn .+= ∂lψ∂m' * ∂lψ∂m ./ nw
                μel∂²lψ∂mn .+= (el.*∂lψ∂m)' * ∂lψ∂m ./ nw
                μ∂el∂m∂lψ∂n .+= (∂el∂m' * ∂lψ∂m .+ ∂lψ∂m' *∂el∂m) ./ nw
            end
            en = μel/ict
            σen = sqrt( (μel²/ict - en.^2) / (ict*nw - 1) )
            ekt = μek/ict
            vnet = μvne/ict
            veet = μvee/ict
            vt = vnet + veet
            if ict % 100 == 0
                # ictd = div(ict,100)
                # aict[ictd] = ict/100
                # Enct[ictd] = en
                printprg(en,σen,acc,τr,ict,nrun)
                # display( plot(aict[1:ictd],Enct[1:ictd], ylims = (-2.861679995612-0.0005,-2.861679995612+0.0005) ))
                # rept = 1
                # while rept < ictd && rept <= 10
                #     Enrpt = (collect((rept+1):ictd).*Enct[(rept+1):ictd].-rept*Enct[rept])./(collect((rept+1):ictd).-rept)
                #     display( plot!(aict[rept+1:ictd],Enrpt, ylims = (-2.861679995612-0.0005,-2.861679995612+0.0005) ))
                #     rept += 1   
                # end
                # println("-------------------------------------------------------")
                # println("E = $en ± $σen")
                # println("Step = $ict")
                # println("Acceptance = $acc")
                # println("τ = $τr")
                τr = acc*τr/accr
            end
        end
        if swt == 1
            swt = 0
        else
            swt = 1
            Enrf[inr] = en
            σnrf[inr] = σen
            Ekf[inr] = ekt
            Vnef[inr] = vnet
            Veef[inr] = veet
            Vtf[inr] = vt
            c9f[:,inr] = c9r
            μ∂lψ∂mf[:,inr] = μ∂lψ∂m
            μ∂el∂mf[:,inr] = μ∂el∂m
            μel∂lψ∂mf[:,inr] = μel∂lψ∂m
            μ∂²lψ∂mnf[:,:,inr] = μ∂²lψ∂mn
            μel∂²lψ∂mnf[:,:,inr] = μel∂²lψ∂mn
            μ∂el∂m∂lψ∂nf[:,:,inr] = μ∂el∂m∂lψ∂n
            fmr = zeros(9)
            @inbounds for mk = 1:9
                fmr[mk] = 2( μel∂lψ∂m[mk] - en*μ∂lψ∂m[mk] ) / nstp
            end
            hmnr = zeros(9,9)
            @inbounds for mk = 1:9
                for nk = 1:9
                    hmnr[mk,nk] = ( 4(μel∂²lψ∂mn[mk,nk] - en*μ∂²lψ∂mn[mk,nk]) -
                                    2(μ∂lψ∂m[mk]*fmr[nk] + μ∂lψ∂m[nk]*fmr[mk]) +
                                    μ∂el∂m∂lψ∂n[mk,nk] -
                                    (μ∂el∂m[mk]*μ∂lψ∂m[nk] + μ∂el∂m[nk]*μ∂lψ∂m[mk]) / nstp
                                    ) / nstp
                end
            end
            fmf[:,inr] = fmr
            display( plot( plot(c9f[:,1:inr]', marker=(:circle,3)),
                           plot(fmf[:,1:inr]', marker=(:circle,3)) ) )
            hmnf[:,:,inr] = hmnr
            if sum(c9i) ≠ 0
                eigd, eigv = eigen(hmnr[c9i,c9i])
                hmnc⁻¹ = eigv*diagm(abs.(eigd.^-1))*eigv'
                Δc9 = - hmnc⁻¹*fmr[c9i]
                MΔc9 = maximum(abs.(Δc9))
                if MΔc9 > Δlim
                    Δc9 = Δlim*Δc9/MΔc9
                end
                # Δc9[abs.(Δc9) .> 1] .= sign.(Δc9[abs.(Δc9) .> 1])
                c9r[c9i] += Δc9
            end
        end
    end
    return Corstr(c9r,c9i), NRstr(Enrf, σnrf, Ekf, Vnef, Veef, Vtf, c9f, fmf, hmnf),
           NRμstr(μ∂lψ∂mf, μ∂el∂mf, μel∂lψ∂mf, μ∂²lψ∂mnf, μel∂²lψ∂mnf, μ∂el∂m∂lψ∂nf)
end

function main(fout::String, atom::Int, ion::Int, spin::Char, base::String)
    nra = Array{NRstr}(undef, 0)
    nrμa = Array{NRμstr}(undef, 0)
    qmcpa = Array{QMCstr}(undef, 0)
    gtoa = Array{GTOstr}(undef, 0)
    corpa = Array{Corstr}(undef, 0)
    qmcp, gto, corp = qmcinp(true, true, atom, ion, spin, base)
    push!(qmcpa, qmcp)
    push!(gtoa, gto)
    push!(corpa, corp)
    c9s, o1, o2 = vqmcnr(qmcp, gto, corp)
    push!(nra, o1)
    push!(nrμa, o2)
    save(fout, "nr", nra, "qmcp", qmcpa, "gto", gtoa, "corp", corpa, "c9s", c9s, "nrμ", nrμa)
    return c9s, nra, qmcpa, gtoa, corpa, nrμa
end

function main(fout::String, c9s::Corstr, nr::Array{NRstr}, qmcp::Array{QMCstr}, gto::Array{GTOstr}, corp::Array{Corstr}, nrμ::Array{NRμstr}, atom::Int, ion::Int, spin::Char, base::String)
    qmcpn, gton = qmcinp(true, false, atom, ion, spin, base)
    push!(qmcp, qmcpn)
    push!(gto, gton)
    push!(corp, c9s)
    c9s, o1, o2 = vqmcnr(qmcpn, gton, c9s)
    push!(nr, o1)
    push!(nrμ, o2)
    save(fout, "nr", nr, "qmcp", qmcp, "gto", gto, "corp", corp, "c9s", c9s, "nrμ", nrμ)
    return c9s
end

function main(fout::String, atom::Int, ion::Int, spin::Char, base::String, script::Int)
    nra = Array{NRstr}(undef, 0)
    nrμa = Array{NRμstr}(undef, 0)
    qmcpa = Array{QMCstr}(undef, 0)
    gtoa = Array{GTOstr}(undef, 0)
    corpa = Array{Corstr}(undef, 0)
    println("Enter seed number:")
    seed = parse(Int,readline())
    Random.seed!(seed)
    qmcp = QMCstr(
        100, # nw
        seed, # seed
        100, # tnstp
        1000, # nstp
        1., # τ
        .5, # acc
    20, # nrstp
    1. # Δlim
    )
    gto, corp = qmcinp(false, true, atom, ion, spin, base)
    push!(qmcpa, qmcp)
    push!(gtoa, gto)
    push!(corpa, corp)
    c9s, o1, o2 = vqmcnr(qmcp, gto, corp)
    push!(nra, o1)
    push!(nrμa, o2)
    save(fout, "nr", nra, "qmcp", qmcpa, "gto", gtoa, "corp", corpa, "c9s", c9s, "nrμ", nrμa)

    qmcp = QMCstr(
        100, # nw
        seed, # seed
        1000, # tnstp
        20000, # nstp
        1., # τ
        .5, # acc
    5, # nrstp
    1. # Δlim
    )
    gto = qmcinp(false, false, atom, ion, spin, base)
    push!(qmcpa, qmcp)
    push!(gtoa, gto)
    push!(corpa, c9s)
    c9s, o1, o2 = vqmcnr(qmcp, gto, c9s)
    push!(nra, o1)
    push!(nrμa, o2)
    save(fout, "nr", nra, "qmcp", qmcpa, "gto", gtoa, "corp", corpa, "c9s", c9s, "nrμ", nrμa)

    qmcp = QMCstr(
        100, # nw
        seed, # seed
        1000, # tnstp
        20000, # nstp
        1., # τ
        .5, # acc
    20, # nrstp
    1. # Δlim
    )
    gto = qmcinp(false, false, atom, ion, spin, base)
    push!(qmcpa, qmcp)
    push!(gtoa, gto)
    push!(corpa, c9s)
    c9s, o1, o2 = vqmcnr(qmcp, gto, c9s)
    push!(nra, o1)
    push!(nrμa, o2)
    save(fout, "nr", nra, "qmcp", qmcpa, "gto", gtoa, "corp", corpa, "c9s", c9s, "nrμ", nrμa)

    c9s = Corstr(sumd(nra[3].c9[:,:], 2) / size(nra[3].E[:], 1), c9s.c9i)

    qmcp = QMCstr(
        100, # nw
        seed, # seed
        5000, # tnstp
        200000, # nstp
        1., # τ
        .5, # acc
    1, # nrstp
    1. # Δlim
    )
    gto = qmcinp(false, false, atom, ion, spin, base)
    push!(qmcpa, qmcp)
    push!(gtoa, gto)
    push!(corpa, c9s)
    c9s, o1, o2 = vqmcnr(qmcp, gto, c9s)
    push!(nra, o1)
    push!(nrμa, o2)
    save(fout, "nr", nra, "qmcp", qmcpa, "gto", gtoa, "corp", corpa, "c9s", c9s, "nrμ", nrμa)

    return c9s, nra, qmcpa, gtoa, corpa, nrμa
end
