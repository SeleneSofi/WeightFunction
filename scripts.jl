include("quadrature.jl")
include("weight_function.jl")
include("energyN.jl")

function run_script(script, qt, x0, ζ, ne, Elim)
    nsum = 4
    if script == "1s2s3s,2p3p"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0, 0],             [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-10, iterations=5_000, show_trace=true, show_every=100)
        for _ = 1:3
            op = deepcopy(op0)
            op[1] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
            op[2] = nSig_opt(1,1,[1],   ones(nsum),ones(nsum))
            op[3] = nSig_opt(1,1,[1,1],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
            op = deepcopy(op0)
            op[4] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
            op[5] = nSig_opt(1,1,[1],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
        end
    elseif script == "1s2s3s,2p"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0, 0],             [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-10, iterations=5_000, show_trace=true, show_every=100)
        for _ = 1:3
            op = deepcopy(op0)
            op[1] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
            op[2] = nSig_opt(1,1,[1],   ones(nsum),ones(nsum))
            op[3] = nSig_opt(1,1,[1,1],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
            op = deepcopy(op0)
            op[4] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
        end
    elseif script == "1s2s,2p"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-10, iterations=5_000, show_trace=true, show_every=100)
        for _ = 1:3
            op = deepcopy(op0)
            op[1] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
            op[2] = nSig_opt(1,1,[1],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
            op = deepcopy(op0)
            op[3] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
        end
    elseif script == "1s2s"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-10, iterations=5_000, show_trace=true, show_every=100)
        for _ = 1:3
            op = deepcopy(op0)
            op[1] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
            op[2] = nSig_opt(1,1,[1],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
        end
    elseif script == "1s,2s,3s,2p,3p"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0, 0],             [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-13, iterations=10_000, show_trace=true, show_every=100)
        op = deepcopy(op0)
        op[1] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[2] = nSig_opt(1,1,[1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[3] = nSig_opt(1,1,[1,1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[4] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[5] = nSig_opt(1,1,[1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
    elseif script == "1s,2s,3s,2p,3p,ref"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0, 0],             [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-16, iterations=2_000, show_trace=true, show_every=100)
        op = deepcopy(op0)
        op[1] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[2] = nSig_opt(1,1,[1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[3] = nSig_opt(1,1,[1,1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[4] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[5] = nSig_opt(1,1,[1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
    elseif script == "1s,2s,3s,2p,3p,no_v"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0, 0],             [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-13, iterations=10_000, show_trace=true, show_every=100)
        op = deepcopy(op0)
        op[1] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[2] = nSig_opt(1,0,[1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[3] = nSig_opt(1,0,[1,1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[4] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[5] = nSig_opt(1,0,[1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
    elseif script == "1s,2s,3s,2p,3p,no_r"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0, 0],             [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-13, iterations=10_000, show_trace=true, show_every=100)
        op = deepcopy(op0)
        op[1] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[2] = nSig_opt(1,1,[0],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[3] = nSig_opt(1,1,[0,0],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[4] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[5] = nSig_opt(1,1,[0],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
    elseif script == "1s,2s,3s,2p"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0, 0],             [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-13, iterations=10_000, show_trace=true, show_every=100)
        op = deepcopy(op0)
        op[1] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[2] = nSig_opt(1,1,[1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[3] = nSig_opt(1,1,[1,1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[4] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
    elseif script == "1s,2s,3s,2p,no_v"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0, 0],             [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-13, iterations=10_000, show_trace=true, show_every=100)
        op = deepcopy(op0)
        op[1] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[2] = nSig_opt(1,0,[1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[3] = nSig_opt(1,0,[1,1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[4] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
    elseif script == "1s,2s,2p"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-13, iterations=10_000, show_trace=true, show_every=100)
        op = deepcopy(op0)
        op[1] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[2] = nSig_opt(1,1,[1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[3] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
    elseif script == "1s,2s,2p,no_v"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-13, iterations=10_000, show_trace=true, show_every=100)
        op = deepcopy(op0)
        op[1] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[2] = nSig_opt(1,0,[1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[3] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
    elseif script == "1s,2s"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-13, iterations=10_000, show_trace=true, show_every=100)
        op = deepcopy(op0)
        op[1] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[2] = nSig_opt(1,1,[1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
    elseif script == "1s,2s,no_v"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-13, iterations=10_000, show_trace=true, show_every=100)
        op = deepcopy(op0)
        op[1] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
        op = deepcopy(op0)
        op[2] = nSig_opt(1,0,[1],   ones(nsum),ones(nsum))
        f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
        results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
        xar_res = results.minimizer
        x0 = get_x(results.minimizer, fixparams[1:2]...)
    elseif script == "1s"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-15, iterations=10_000, show_trace=true, show_every=100)
        for _ = 1:3
            op = deepcopy(op0)
            op[1] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
        end
    elseif script == "1s,H"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0], [0])
        ]
        optopt = Optim.Options(g_tol=5e-15, iterations=10_000, show_trace=true, show_every=100)
        for _ = 1:3
            op = deepcopy(op0)
            op[1] = nSig_opt(1,1,[],   ones(1),ones(1))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
        end
    elseif script == "1s2s3s,2p3p,no_v"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0, 0],             [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-10, iterations=5_000, show_trace=true, show_every=100)
        for _ = 1:3
            op = deepcopy(op0)
            op[1] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
            op[2] = nSig_opt(1,0,[1],   ones(nsum),ones(nsum))
            op[3] = nSig_opt(1,0,[1,1],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
            op = deepcopy(op0)
            op[4] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
            op[5] = nSig_opt(1,0,[1],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
        end
    elseif script == "1s2s3s,2p3p,no_r"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0, 0],             [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-10, iterations=5_000, show_trace=true, show_every=100)
        for _ = 1:3
            op = deepcopy(op0)
            op[1] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
            op[2] = nSig_opt(1,1,[0],   ones(nsum),ones(nsum))
            op[3] = nSig_opt(1,1,[0,0],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
            op = deepcopy(op0)
            op[4] = nSig_opt(1,1,[],   ones(nsum),ones(nsum))
            op[5] = nSig_opt(1,1,[0],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
        end
    elseif script == "1s2s3s,2p,no_v"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0, 0],             [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-10, iterations=5_000, show_trace=true, show_every=100)
        for _ = 1:3
            op = deepcopy(op0)
            op[1] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
            op[2] = nSig_opt(1,0,[1],   ones(nsum),ones(nsum))
            op[3] = nSig_opt(1,0,[1,1],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
            op = deepcopy(op0)
            op[4] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
        end
    elseif script == "1s2s,2p,no_v"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-10, iterations=5_000, show_trace=true, show_every=100)
        for _ = 1:3
            op = deepcopy(op0)
            op[1] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
            op[2] = nSig_opt(1,0,[1],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
            op = deepcopy(op0)
            op[3] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
        end
    elseif script == "1s2s,no_v"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
            nSig_opt(0, 0, [0],                [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-10, iterations=5_000, show_trace=true, show_every=100)
        for _ = 1:3
            op = deepcopy(op0)
            op[1] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
            op[2] = nSig_opt(1,0,[1],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
        end
    elseif script == "1s,no_v"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0,0,0,0], [0,0,0,0])
        ]
        optopt = Optim.Options(g_tol=5e-15, iterations=10_000, show_trace=true, show_every=100)
        for _ = 1:3
            op = deepcopy(op0)
            op[1] = nSig_opt(1,0,[],   ones(nsum),ones(nsum))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
        end
    elseif script == "1s,no_v,H"
        op0 = [
            nSig_opt(0, 0, Float64[],          [0], [0])
        ]
        optopt = Optim.Options(g_tol=5e-15, iterations=10_000, show_trace=true, show_every=100)
        for _ = 1:3
            op = deepcopy(op0)
            op[1] = nSig_opt(1,0,[],   ones(1),ones(1))
            f_opt, xar0, fixparams, get_x = make_f(qt,x0,op,ζ,ne)
            results = optimize(x -> f_opt(x, fixparams...) - Elim, xar0, NelderMead(), optopt)
            x0 = get_x(results.minimizer, fixparams[1:2]...)
        end
    end
    return x0
end