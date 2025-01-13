include("quadrature.jl")
include("weight_function.jl")
include("energyN.jl")

function graph_wf(
    qti::Array{ETabn},
    x0::Array{nSig}
)
    qt = [ETabn(iqti.a, iqti.b, 1000) for iqti in qti]
    lnα = [calc_lnα(iqt) for iqt = qt]
    al = [exp.(ilnα) for ilnα = lnα]
    gnt = [Array{Array}(undef, 0) for _ = lnα]
    for ixc = x0
        push!(gnt[ixc.l], weight_funtion(ixc, lnα[ixc.l]))
    end
    gn = [hcat(ignt...) for ignt = gnt]
    if size(gn)[1] == 1
        gn = normalized(al..., gn...)
        gn = orthonormalized(al..., gn)
        l = repeat([(blank=false,)],1,size(gn)[2])
        plt = [plot(lnα[1], gn[:,i]) for i = 1:size(gn)[2]]
        plot(plt..., layout=l, size=(400*size(gn)[2],300))
    else
        gn = [normalized(al..., gn...)...]
        gn = [orthonormalized(al..., gn...)...]
        l = [
            repeat([(blank=false,)],1,size(gn[1])[2]),
            [(blank=true,) repeat([(blank=false,width=1/size(gn[1])[2])],1,size(gn[2])[2]) (blank=true,)]
        ]
        plt = vcat(
            [plot(lnα[1], gn[1][:,i]) for i = 1:size(gn[1])[2]],
            [plot(lnα[2], gn[2][:,i]) for i = 1:size(gn[2])[2]]
        )
        plot(plt..., layout=l, size=(400*size(gn[1])[2],600))
    end
end