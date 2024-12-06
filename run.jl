qt1000 = [
    ETabn(qt[1].a,qt[1].b,1000),
    ETabn(qt[2].a,qt[2].b,1000)
]
lnα = [calc_lnα(iqt) for iqt = qt1000]
gnt = [Array{Array}(undef, 0) for _ = lnα]
for ixc = x0
    push!(gnt[ixc.l], weight_funtion(ixc, lnα[ixc.l]))
end
gn = [hcat(ignt...) for ignt = gnt]

idx = 1000 - sum(gn[1][:,2] .< 0)
ab = inv([lnα[1][idx:idx+1] ones(2)])*gn[1][idx:idx+1,2]
gn[1][idx:idx+1,2]
-ab[2]/ab[1]

cs = cumsum(gn[1][:,3])
idx = argmax(cs)
ab = inv([lnα[1][idx:idx+1] ones(2)])*gn[1][idx:idx+1,2]
gn[1][idx:idx+1,3]
-ab[2]/ab[1]

idx = argmin(cs)
ab = inv([lnα[1][idx:idx+1] ones(2)])*gn[1][idx:idx+1,2]
gn[1][idx:idx+1,3]
-ab[2]/ab[1]

idx = 1000 - sum(gn[2][:,2] .< 0)
ab = inv([lnα[2][idx:idx+1] ones(2)])*gn[2][idx:idx+1,2]
gn[2][idx:idx+1,2]
-ab[2]/ab[1]

function runfull(qt, x0, ζ, ne, Elim)
    optopt = Optim.Options(g_tol=5e-14, iterations=2_000, show_trace=true, show_every=100)
    nsum = 4
    op0 = [
        nSig_opt(0, 0, Float64[], zeros(nsum), zeros(nsum))
        nSig_opt(0, 0, [0],       zeros(nsum), zeros(nsum))
        nSig_opt(0, 0, [0, 0],    zeros(nsum), zeros(nsum))
        nSig_opt(0, 0, Float64[], zeros(nsum), zeros(nsum))
        nSig_opt(0, 0, [0],       zeros(nsum), zeros(nsum))
    ]
    op01 = [
        nSig_opt(1, 0, Float64[], ones(nsum), ones(nsum))
        nSig_opt(1, 0, [1],       ones(nsum), ones(nsum))
        nSig_opt(1, 0, [1, 1],    ones(nsum), ones(nsum))
        nSig_opt(1, 0, Float64[], ones(nsum), ones(nsum))
        nSig_opt(1, 0, [1],       ones(nsum), ones(nsum))
    ]
    op1 = [
        nSig_opt(1, 1, Float64[], ones(nsum), ones(nsum))
        nSig_opt(1, 1, [1],       ones(nsum), ones(nsum))
        nSig_opt(1, 1, [1, 1],    ones(nsum), ones(nsum))
        nSig_opt(1, 1, Float64[], ones(nsum), ones(nsum))
        nSig_opt(1, 1, [1],       ones(nsum), ones(nsum))
    ]
    x0_ini = deepcopy(x0)
    x_scr = run_script("1s2s3s,2p3p,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,3s,2p,3p,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,3s,2p,3p,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,3s,2p,3p", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x0bak = deepcopy(x0)
    f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
    err = f_opt(xar0, fixparams...) - Elim
    for i = 1:5
        x0[i] = x0_ini[i]
        op = deepcopy(op0)
        op[i] = op01[i]
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op[i] = op1[i]
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
        if f_opt(xar0, fixparams...) - Elim < err
            err = f_opt(xar0, fixparams...) - Elim
        else
            x0[i] = x0bak[i]
        end
    end
    x_scr = run_script("1s,2s,3s,2p,3p", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,3s,2p,3p", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    # x_scr = run_script("1s,2s,3s,2p,3p,ref", qt, x0, ζ, ne, Elim)
    # x0 .= x_scr
    # x_scr = run_script("1s,2s,3s,2p,3p,ref", qt, x0, ζ, ne, Elim)
    # x0 .= x_scr
    f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
    err = f_opt(xar0, fixparams...) - Elim
    return x_scr, err
end

function runfull(qt, x0, ζ, ne, Elim)
    optopt = Optim.Options(g_tol=5e-14, iterations=2_000, show_trace=true, show_every=100)
    nsum = 4
    op0 = [
        nSig_opt(0, 0, Float64[], zeros(nsum), zeros(nsum))
        nSig_opt(0, 0, [0],       zeros(nsum), zeros(nsum))
        nSig_opt(0, 0, [0, 0],    zeros(nsum), zeros(nsum))
        nSig_opt(0, 0, Float64[], zeros(nsum), zeros(nsum))
    ]
    op01 = [
        nSig_opt(1, 0, Float64[], ones(nsum), ones(nsum))
        nSig_opt(1, 0, [1],       ones(nsum), ones(nsum))
        nSig_opt(1, 0, [1, 1],    ones(nsum), ones(nsum))
        nSig_opt(1, 0, Float64[], ones(nsum), ones(nsum))
    ]
    op1 = [
        nSig_opt(1, 1, Float64[], ones(nsum), ones(nsum))
        nSig_opt(1, 1, [1],       ones(nsum), ones(nsum))
        nSig_opt(1, 1, [1, 1],    ones(nsum), ones(nsum))
        nSig_opt(1, 1, Float64[], ones(nsum), ones(nsum))
    ]
    x0_ini = deepcopy(x0)
    x_scr = run_script("1s2s3s,2p,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,3s,2p,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,3s,2p,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,3s,2p", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x0bak = deepcopy(x0)
    f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
    err = f_opt(xar0, fixparams...) - Elim
    for i = 1:4
        x0[i] = x0_ini[i]
        op = deepcopy(op0)
        op[i] = op01[i]
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op[i] = op1[i]
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
        if f_opt(xar0, fixparams...) - Elim < err
            err = f_opt(xar0, fixparams...) - Elim
        else
            x0[i] = x0bak[i]
        end
    end
    x_scr = run_script("1s,2s,3s,2p", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,3s,2p", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
    err = f_opt(xar0, fixparams...) - Elim
    return x_scr, err
end

function runfull(qt, x0, ζ, ne, Elim)
    optopt = Optim.Options(g_tol=5e-14, iterations=2_000, show_trace=true, show_every=100)
    nsum = 4
    op0 = [
        nSig_opt(0, 0, Float64[], zeros(nsum), zeros(nsum))
        nSig_opt(0, 0, [0],       zeros(nsum), zeros(nsum))
        nSig_opt(0, 0, Float64[], zeros(nsum), zeros(nsum))
    ]
    op01 = [
        nSig_opt(1, 0, Float64[], ones(nsum), ones(nsum))
        nSig_opt(1, 0, [1],       ones(nsum), ones(nsum))
        nSig_opt(1, 0, Float64[], ones(nsum), ones(nsum))
    ]
    op1 = [
        nSig_opt(1, 1, Float64[], ones(nsum), ones(nsum))
        nSig_opt(1, 1, [1],       ones(nsum), ones(nsum))
        nSig_opt(1, 1, Float64[], ones(nsum), ones(nsum))
    ]
    x0_ini = deepcopy(x0)
    x_scr = run_script("1s2s,2p,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,2p,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,2p,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,2p", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x0bak = deepcopy(x0)
    f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
    err = f_opt(xar0, fixparams...) - Elim
    for i = 1:3
        x0[i] = x0_ini[i]
        op = deepcopy(op0)
        op[i] = op01[i]
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op[i] = op1[i]
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
        if f_opt(xar0, fixparams...) - Elim < err
            err = f_opt(xar0, fixparams...) - Elim
        else
            x0[i] = x0bak[i]
        end
    end
    x_scr = run_script("1s,2s,2p", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,2p", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
    err = f_opt(xar0, fixparams...) - Elim
    return x_scr, err
end

function runfull(qt, x0, ζ, ne, Elim)
    optopt = Optim.Options(g_tol=5e-14, iterations=2_000, show_trace=true, show_every=100)
    nsum = 4
    op0 = [
        nSig_opt(0, 0, Float64[], zeros(nsum), zeros(nsum))
        nSig_opt(0, 0, [0],       zeros(nsum), zeros(nsum))
    ]
    op01 = [
        nSig_opt(1, 0, Float64[], ones(nsum), ones(nsum))
        nSig_opt(1, 0, [1],       ones(nsum), ones(nsum))
    ]
    op1 = [
        nSig_opt(1, 1, Float64[], ones(nsum), ones(nsum))
        nSig_opt(1, 1, [1],       ones(nsum), ones(nsum))
    ]
    x0_ini = deepcopy(x0)
    x_scr = run_script("1s2s,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x0bak = deepcopy(x0)
    f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
    err = f_opt(xar0, fixparams...) - Elim
    for i = 1:2
        x0[i] = x0_ini[i]
        op = deepcopy(op0)
        op[i] = op01[i]
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
        if f_opt(xar0, fixparams...) - Elim < err
            err = f_opt(xar0, fixparams...) - Elim
        else
            x0[i] = x0bak[i]
        end
    end
    x_scr = run_script("1s,2s", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
    err = f_opt(xar0, fixparams...) - Elim
    return x_scr, err
end

function runfull(qt, x0, ζ, ne, Elim)
    nsum = 4
    op1 = [
        nSig_opt(1, 1, Float64[], ones(nsum), ones(nsum))
    ]
    x_scr = run_script("1s,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
    err = f_opt(xar0, fixparams...) - Elim
    return x_scr, err
end

function runfull(qt, x0, ζ, ne, Elim)
    nsum = 1
    op1 = [
        nSig_opt(1, 1, Float64[], ones(nsum), ones(nsum))
    ]
    x_scr = run_script("1s,no_v,H", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,no_v,H", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,no_v,H", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,H", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,H", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,H", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
    err = f_opt(xar0, fixparams...) - Elim
    return x_scr, err
end

function runfull(qt, x0, ζ, ne, Elim)
    optopt = Optim.Options(g_tol=5e-14, iterations=10_000, show_trace=true, show_every=100)
    nsum = 4
    op0 = [
        nSig_opt(0, 0, Float64[], zeros(nsum), zeros(nsum))
        nSig_opt(0, 0, [0],       zeros(nsum), zeros(nsum))
        nSig_opt(0, 0, Float64[], zeros(nsum), zeros(nsum))
    ]
    op1 = [
        nSig_opt(1, 1, Float64[], ones(nsum), ones(nsum))
        nSig_opt(1, 1, [1],       ones(nsum), ones(nsum))
        nSig_opt(1, 1, Float64[], ones(nsum), ones(nsum))
    ]
    x0_ini = deepcopy(x0)
    x_scr = run_script("1s2s,2p,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,2p,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,2p,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,2p", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x0bak = deepcopy(x0)
    f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
    err = f_opt(xar0, fixparams...) - Elim
    for i = 1:3
        x0[i] = x0_ini[i]
        op = deepcopy(op0)
        op[i] = op1[i]
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
        if f_opt(xar0, fixparams...) - Elim < err
            err = f_opt(xar0, fixparams...) - Elim
        else
            x0[i] = x0bak[i]
        end
    end
    x_scr = run_script("1s,2s,2p", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
    err = f_opt(xar0, fixparams...) - Elim
    return x_scr, err
end

function runfull(qt, x0, ζ, ne, Elim)
    optopt = Optim.Options(g_tol=5e-16, iterations=10_000, show_trace=true, show_every=100)
    nsum = 4
    op0 = [
        nSig_opt(0, 0, Float64[], zeros(nsum), zeros(nsum))
        nSig_opt(0, 0, [0],       zeros(nsum), zeros(nsum))
    ]
    op1 = [
        nSig_opt(1, 1, Float64[], ones(nsum), ones(nsum))
        nSig_opt(1, 1, [1],       ones(nsum), ones(nsum))
    ]
    x0_ini = deepcopy(x0)
    x_scr = run_script("1s2s,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s,no_v", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x_scr = run_script("1s,2s", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    x0bak = deepcopy(x0)
    f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
    err = f_opt(xar0, fixparams...) - Elim
    for i = 1:2
        x0[i] = x0_ini[i]
        op = deepcopy(op0)
        op[i] = op1[i]
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
        if f_opt(xar0, fixparams...) - Elim < err
            err = f_opt(xar0, fixparams...) - Elim
        else
            x0[i] = x0bak[i]
        end
    end
    x_scr = run_script("1s,2s", qt, x0, ζ, ne, Elim)
    x0 .= x_scr
    f_opt, xar0, fixparams, get_x = make_f(qt,x0,op1,ζ,ne)
    err = f_opt(xar0, fixparams...) - Elim
    return x_scr, err
end