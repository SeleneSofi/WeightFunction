include("including.jl")

nqt = 35
op0 = [
    nSig_opt(0, 0, Float64[], [0,0,0,0], [0,0,0,0])
    nSig_opt(0, 0, [0], [0,0,0,0], [0,0,0,0])
    nSig_opt(0, 0, [0, 0], [0,0,0,0], [0,0,0,0])
    nSig_opt(0, 0, Float64[], [0,0,0,0], [0,0,0,0])
    nSig_opt(0, 0, [0], [0,0,0,0], [0,0,0,0])
]

for ikey = keys(database)
    for jkey = keys(database[ikey])
        if jkey isa Number
            ζ::Int16 = ikey
            ne::Int16 = jkey
            print(string(ζ)*" "*string(ne)*"|")
            if ne == 1
                op = [nSig_opt(0, 0, Float64[], [0], [0])]
            elseif ne <= 2
                op = op0[[1]]
            elseif ne <= 4
                op = op0[[1,2]]
            elseif ne <= 10
                op = op0[[1,2,4]]
            elseif ne <= 12
                op = op0[[1,2,3,4]]
            else
                op = op0
            end
            data = database[ikey][jkey]
            Elim = data["Elim"]
            x = data["x"]
            qt50 = data["qt"]
            qt = [typeof(iqt50)(iqt50.a,iqt50.b,nqt) for iqt50 = qt50]
            f_opt, xar0, fixparams, get_x = make_f(qt,x,op,ζ,ne)
            data["Err"] = f_opt(xar0, fixparams...) - Elim
            print(string(data["Err"])*"\n")
        end
    end
end