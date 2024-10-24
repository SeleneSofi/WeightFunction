struct ETabn
    a::Float64
    b::Float64
    n::Int16
end

function calc_lnα(x::ETabn)
    return Array(LinRange(x.a,x.b,x.n))
end