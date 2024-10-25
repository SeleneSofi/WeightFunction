function make_f(
    qt::Array{ETabn},
    x0::Array{nGum},
    op::Array{nGum_opt},
    ζ::Int16,
    ne::Int16
)
    lnα = [calc_lnα(iqt) for iqt = qt]
    al = [exp.(ilnα) for ilnα = lnα]
    Ht, Vt = density_prep(al..., ζ)
    xar0, car0, orbs = make_array(x0, op)
    function f_opt(
        xar::Array{Float64, 1},
        car0::Array{Float64, 1},
        orbs::Array{Tuple, 1},
        lnα::Array,
        al::Array,
        Ht::Tuple,
        Vt::Tuple,
        ne::Int16
    )
        xcar = copy(car0)
        xcar[isnan.(xcar)] .= xar
        xc = make_wfstr_nGum(xcar, orbs)
        gnt = [Array{Array}(undef, 0) for _ = lnα]
        for ixc = xc
            push!(gnt[ixc.l], weight_funtion(ixc, lnα[ixc.l]))
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
        return make_wfstr_nGum(xcar, orbs)
    end

    return f_opt, xar0, (car0, orbs, lnα, al, Ht, Vt, ne), get_x
end

function make_f(
    qt::Array{ETabn},
    x0::Array{nGuma},
    op::Array{nGuma_opt},
    ζ::Int16,
    ne::Int16
)
    lnα = [calc_lnα(iqt) for iqt = qt]
    al = [exp.(ilnα) for ilnα = lnα]
    Ht, Vt = density_prep(al..., ζ)
    xar0, car0, orbs = make_array(x0, op)
    function f_opt(
        xar::Array{Float64, 1},
        car0::Array{Float64, 1},
        orbs::Array{Tuple, 1},
        lnα::Array,
        al::Array,
        Ht::Tuple,
        Vt::Tuple,
        ne::Int16
    )
        xcar = copy(car0)
        xcar[isnan.(xcar)] .= xar
        xc = make_wfstr_nGuma(xcar, orbs)
        gnt = [Array{Array}(undef, 0) for _ = lnα]
        for ixc = xc
            push!(gnt[ixc.l], weight_funtion(ixc, lnα[ixc.l]))
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
        return make_wfstr_nGuma(xcar, orbs)
    end

    return f_opt, xar0, (car0, orbs, lnα, al, Ht, Vt, ne), get_x
end

function make_f(
    qt::Array{ETabn},
    x0::Array{nGumqe},
    op::Array{nGumqe_opt},
    ζ::Int16,
    ne::Int16
)
    lnα = [calc_lnα(iqt) for iqt = qt]
    al = [exp.(ilnα) for ilnα = lnα]
    Ht, Vt = density_prep(al..., ζ)
    xar0, car0, orbs = make_array(x0, op)
    function f_opt(
        xar::Array{Float64, 1},
        car0::Array{Float64, 1},
        orbs::Array{Tuple, 1},
        lnα::Array,
        al::Array,
        Ht::Tuple,
        Vt::Tuple,
        ne::Int16
    )
        xcar = copy(car0)
        xcar[isnan.(xcar)] .= xar
        xc = make_wfstr_nGumqe(xcar, orbs)
        gnt = [Array{Array}(undef, 0) for _ = lnα]
        for ixc = xc
            push!(gnt[ixc.l], weight_funtion(ixc, lnα[ixc.l]))
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
        return make_wfstr_nGumqe(xcar, orbs)
    end

    return f_opt, xar0, (car0, orbs, lnα, al, Ht, Vt, ne), get_x
end