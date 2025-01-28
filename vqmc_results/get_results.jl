include("../vqmcnr-gto.jl")

ζs = collect(1:19)
res_arr = Array{String, 1}(undef, 0)
push!(res_arr, "Element,Atomic Number,Charge,Multiplicity,Energy,Standard Deviation")

for ζ = ζs
    dc = database[ζ]
    sym = dc["symbol"]
    for chg = [-1,0,1]
        ne = ζ - chg
        if haskey(dc, ne)
            S = dc[ne]["S"]
            file = sym * "_" * string(chg) * "_" * string(S) * ".jld"
            res = load(file)
            E = res["nr"][end].E
            σ = res["nr"][end].σ
            res_line = sym * "," * string(ζ) * "," * string(chg) * "," * string(S) * "," * string(E[1])  * "," * string(σ[1])
            push!(res_arr, res_line)
        end
    end
end

open("results.csv", "w") do fl
    for ln = res_arr
        write(fl, ln * '\n')
    end
end