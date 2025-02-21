z::Int64 = 5
ne::Int64 = 5
qt = database[z][ne]["qt"]
x = database[z][ne]["x"]
Elim = database[z][ne]["Elim"]
op0 = [
    nSig_opt(0, 0, Float64[], [0,0,0,0], [0,0,0,0])
    nSig_opt(0, 0, [0], [0,0,0,0], [0,0,0,0])
    nSig_opt(0, 0, [0, 0], [0,0,0,0], [0,0,0,0])
    nSig_opt(0, 0, Float64[], [0,0,0,0], [0,0,0,0])
    nSig_opt(0, 0, [0], [0,0,0,0], [0,0,0,0])
]
if ne <= 2
    op = op0[1]
elseif ne <= 4
    op = op0[[1,2]]
elseif ne <= 10
    op = op0[[1,2,4]]
elseif ne <= 12
    op = op0[[1,2,3,4]]
else
    op = op0
end
    
nums = [25,30,35,40,45,50]
qts = [[ETabn(iqt.a,iqt.b,n) for iqt in qt] for n in nums]
Es = Array{Float64, 1}(undef, 0)

for iqt in qts
    f_opt, xar0, fixparams, get_x = make_f(iqt,x,op,z,ne)
    push!(Es, f_opt(xar0, fixparams...))
end