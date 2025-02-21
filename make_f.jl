include("quadrature.jl")
include("weight_function.jl")
include("energyN.jl")

function make_f(
    qt::Array{<:QtUnion},
    x0::Array{<:WfUnion},
    op::Array{<:WfUnion_opt},
    ζ::Int64,
    ne::Int64
)
    lnα = [calc_lnα(iqt) for iqt = qt]
    al = [exp.(ilnα) for ilnα = lnα]
    Ht, Vt = density_prep(al..., ζ)
    xar0, car0, orbs = make_wf_array(x0, op)
    function f_opt(
        xar::Array{Float64, 1},
        car0::Array{Float64, 1},
        orbs::Array{Tuple, 1},
        lnα::Array,
        al::Array,
        Ht::Tuple,
        Vt::Tuple,
        ne::Int64,
        wf_type::DataType
    )
        xcar = copy(car0)
        xcar[isnan.(xcar)] .= xar
        xc = make_wfstr(xcar, orbs, wf_type())
        gnt = [Array{Array}(undef, 0) for _ = lnα]
        for ixc = xc
            Δlnα = [
                lnα[1][2] - lnα[1][2],
                ((lnα[1][3:end] - lnα[1][1:end-2])/2)...,
                lnα[1][end] - lnα[1][end-1]
            ]
            push!(gnt[ixc.l], weight_funtion(ixc, lnα[ixc.l]) .* Δlnα)
        end
        gn = [hcat(ignt...) for ignt = gnt]
        return energy(al..., gn..., Ht, Vt, ne)
    end

    function get_x(
        xar::Array{Float64, 1},
        car0::Array{Float64, 1},
        orbs::Array{Tuple, 1},
    )
        xcar = copy(car0)
        xcar[isnan.(xcar)] .= xar
        return make_wfstr(xcar, orbs, wf_type())
    end

    return f_opt, xar0, (car0, orbs, lnα, al, Ht, Vt, ne, typeof(x0[1])), get_x
end

function make_f(
    qt::Array{<:QtUnion},
    qtop::Array{<:QtUnion_opt},
    x0::Array{<:WfUnion},
    op::Array{<:WfUnion_opt},
    ζ::Int64,
    ne::Int64
)
    sumbool = 0
    for iqtop = qtop
        for fld in fieldnames(typeof(iqtop))
            sumbool += getproperty(iqtop, fld)
        end
    end
    if sumbool == 0
        return make_f(qt, x0, op, ζ, ne)
    end
    
    xarwf, carwf, orbs = make_wf_array(x0, op)
    xarqt, carqt, n_items_n = make_qt_array(qt, qtop)
    xar0 = [xarqt; xarwf]
    car0 = [carqt; carwf]

    function f_opt(
        xar::Array{Float64, 1},
        car0::Array{Float64, 1},
        orbs::Array{Tuple, 1},
        n_items_n::Array{Tuple, 1},
        ne::Int64,
        wf_type::DataType,
        qt_type::DataType
    )
        xcar = copy(car0)
        xcar[isnan.(xcar)] .= xar
        qt = make_qtstr!(xcar, n_items_n, qt_type())
        lnα = [calc_lnα(iqt) for iqt = qt]
        al = [exp.(ilnα) for ilnα = lnα]
        xc = make_wfstr(xcar, orbs, wf_type())
        gnt = [Array{Array}(undef, 0) for _ = lnα]
        for ixc = xc
            Δlnα = [
                lnα[1][2] - lnα[1][2],
                ((lnα[1][3:end] - lnα[1][1:end-2])/2)...,
                lnα[1][end] - lnα[1][end-1]
            ]
            push!(gnt[ixc.l], weight_funtion(ixc, lnα[ixc.l]) .* Δlnα)
        end
        gn = [hcat(ignt...) for ignt = gnt]
        Ht, Vt = density_prep(al..., ζ)
        return energy(al..., gn..., Ht, Vt, ne)
    end

    function get_x(
        xar::Array{Float64, 1},
        car0::Array{Float64, 1},
        orbs::Array{Tuple, 1},
    )
        xcar = copy(car0)
        xcar[isnan.(xcar)] .= xar
        return make_wfstr(xcar, orbs, wf_type())
    end

    return f_opt, xar0, (car0, orbs, n_items_n, ne, typeof(x0[1]), typeof(qt[1])), get_x
end

function make_f_mult(
    qt::Array,
    x0::Array{nSig},
    op::Array{nSig_opt},
    ζ::Int64,
    ne::Int64
)
    lnα = [Array{Array}(undef, 0) for _ = qt]
    al = [Array{Array}(undef, 0) for _ = qt]
    Ht = Array{Tuple}(undef, 0)
    Vt = Array{Tuple}(undef, 0)
    for i = eachindex(qt)
        lnα[i] = [calc_lnα(iqt) for iqt = qt[i]]
        al[i] = [exp.(ilnα) for ilnα = lnα[i]]
        iHt, iVt = density_prep(al[i]..., ζ)
        push!(Ht, iHt)
        push!(Vt, iVt)
    end

    xar0, car0, orbs = wf_make_array(x0, op)
    function f_opt(
        xar::Array{Float64, 1},
        car0::Array{Float64, 1},
        orbs::Array{Tuple, 1},
        lnα::Array,
        al::Array,
        Ht::Array,
        Vt::Array,
        ne::Int64
    )
        xcar = copy(car0)
        xcar[isnan.(xcar)] .= xar
        xc = make_wfstr_nSig(xcar, orbs)
        m_en = 0
        for i = eachindex(qt)
            gnt = [Array{Array}(undef, 0) for _ = lnα[i]]
            for ixc = xc
                push!(gnt[ixc.l], weight_funtion(ixc, lnα[i][ixc.l]))
            end
            gn = [hcat(ignt...) for ignt = gnt]
            m_en += energy(al[i]..., gn..., Ht[i], Vt[i], ne)
        end
        return m_en/size(qt,1)
    end

    function get_x(
        xar::Array{Float64, 1},
        car0::Array{Float64, 1},
        orbs::Array{Tuple, 1},
    )
        xcar = copy(car0)
        xcar[isnan.(xcar)] .= xar
        return make_wfstr_nSig(xcar, orbs)
    end

    return f_opt, xar0, (car0, orbs, lnα, al, Ht, Vt, ne), get_x
end