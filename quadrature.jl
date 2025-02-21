quad_list = []
quad_opt_list = []

# Adding Even-Tempered Quadrature that uses first point (a), last point (b)  and number of points (n)

struct ETabn
    n_items::Int64
    a::Float64
    b::Float64
    n::Int64
    ETabn(a,b,n) = new(2,a,b,n)
    ETabn() = new()
end

struct ETabn_opt
    a::Bool
    b::Bool
end

push!(quad_list, ETabn)
push!(quad_opt_list, ETabn_opt)

function calc_lnα(x::ETabn)
    return Array(LinRange(x.a,x.b,x.n))
end

function make_qtstr!(xcar::Array{Float64, 1}, n_items_n::Array{Tuple, 1}, ::ETabn)
    qtar = splice!(xcar, 1:sum(n_qt_items))
    qt = Array{ETabn}(undef, 0)
    for i = n_items_n
        a = popfirst!(qtar)
        b = popfirst!(qtar)
        n = i[2]
        push!(qt, ETabn(a, b, n))
    end
    return qt
end

# Adding Polynomial Quadrature, using 

struct QtPol
    n_items::Int64
    a::Array{Float64}
    n::Int64
    QrPol(a, n) = new(size(a,1), a, n)
    QrPol() = new()
end

struct QtPol_opt
    a::Array{Bool}
end

push!(quad_list, QtPol)
push!(quad_opt_list, QtPol_opt)

function calc_lnα(x::QtPol)
    lnα = zeros(x.n)
    i = collect(1:x.n)
    for d = 1:size(x.a)
        lnα += (i-1)^(d-1) * x.a[d]
    end
    return lnα
end

function make_qtstr!(xcar::Array{Float64, 1}, n_items_n::Array{Tuple, 1}, ::QtPol)
    qtar = splice!(xcar, 1:sum(n_qt_items))
    qt = Array{QtPol}(undef, 0)
    for i = n_items_n
        a = splice!(qtar, 1:n_items)
        n = i[2]
        push!(xc, QtPol(a, n))
    end
    return xc
end

# Defining Union of Quadrature Types

QtUnion = Union{quad_list...}
QtUnion_opt = Union{quad_opt_list...}

function make_qt_array(
    qt::Array{<:QtUnion},
    qtop::Array{<:QtUnion_opt}
)
    xar = Array{Float64}(undef, 0)
    car = Array{Float64}(undef, 0)
    n_items_n = Array{Tuple}(undef, 0)
    qt_n_items = [iqt.n_items for iqt = qt]
    for i = 1:size(qt,1)
        for name = fieldnames(typeof(qt[i]))
            name in [:n_items, :n] ? continue :
            prop = getproperty(qt[i], name)
            opprop = getproperty(qtop[i], name)
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
        push!(n_items_n, (qt[i].n_items, qt[i].n))
    end
    return xar, car, n_items_n
end