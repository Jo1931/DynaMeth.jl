function Yield(sol::T, cstr::Reactor) where {T<:Matrix{Float64}}
    (sol[1, :] .* sol[9, :]) / (cstr.properties.ndotin * sum(cstr.properties.yin[[2, 3]]))
end
function Yield(sol::T, cstr::Reactor) where {T<:Vector{Float64}}
    (sol[1] .* cstr.properties.yin[6]) / (sol[6] * sum(cstr.properties.yin[[2, 3]]))
end
function calcRHSequi(T, y, reac::Reactor)
    # unpack parameter

    @unpack kinetic, properties = reac
    @unpack cp_mix, k, H, yin, P, cp_cat, reflux, Trec, Tc = properties
    @unpack parameter, rates, stoich, catalyst = kinetic

    r = rates(y, T, P / 1e5, parameter.beta, parameter)

    RHS16 = (Matrix(I, 6, 6) - y[1:6] * ones(1, 6)) * stoich * r
    RHS = RHS16

    return RHS
end
ode_equi(u, cstr, t) = calcRHSequi(cstr.properties.T, u, cstr)
function equilibrium(yin, cstr; Tvec=300:1:700)
    u0 = yin
    Y = []
    tspan = (0.0, 50000000.0)
    for i in Tvec
        cstr.properties.T = i
        prob = ODEProblem(ode_equi, u0, tspan)
        sol = solve(prob, Rodas5(), reltol=1e-8, abstol=1e-8, p=cstr)
        Y = push!(Y, Yield(sol[:, end], cstr))
    end
    return Y, Tvec
end


function yield_dyn(sol, yin, ndotin)
    y1out = sol[end-8, :]
    nout = sol[end, :]
    out1 = maximum(y1out .* nout / (yin[2] * ndotin + yin[3] * ndotin))
    out2 = minimum(y1out .* nout / (yin[2] * ndotin + yin[3] * ndotin))
    return out1, out2
end

function dae_bif(residual, dy, y, reac, t)
    x = y

    #reac.properties.yin[2] = 0.2 * (1 + fac * 0.8 * sin(t * 2 * pi / 3600 / tau))
    #reac.properties.yin[4] = 0.8 * (1 - fac * 0.2 * sin(t * 2 * pi / 3600 / tau))
    #reac.properties.ndotin = 0.06 * (1 + fac * sin(t * 2 * pi / 3600 + 1 * pi / tau))
    reac.properties.Tin = per_T(t)
    reac.properties.Tc = reac.properties.Tin
    RHS = calcRHScascade(x, reac)

    for i = 1:reac.properties.num_casc
        v = (9*i-8):9*i
        LHS = calcLHSbif(x[v], reac)
        #print(LHS)
        residual[v] = -LHS * dy[v] + RHS[v]
    end
end



function limitCycle(sol_Bussche, cstr, i_start, i_end)
    YU = []
    YD = []
    TY = []
    Tmax = []
    Tmin = []
    states = []
    for idx = i_start:i_end

        res = sim_one_limit(idx, sol_Bussche, cstr)

        push!(Tmax, res.Tmax)
        push!(Tmin, res.Tmax)
        push!(YU, res.YU)
        push!(YD, res.YD)
        push!(TY, res.TY)
        push!(states, res.states)
        print(idx)
    end

    return SolutionCycle(YU, YD, Tmax, Tmin, TY, states)
end


mutable struct SolutionBranches
    sol::Solution_2d
    stable
    instable
    eig
    limit_cycle
    cstr
    function SolutionBranches(sol::Solution_2d, cstr::Reactor)
        stable, instable, out = YieldBranches(sol::Solution_2d)
        new(sol, stable, instable, out, [], cstr.properties)
    end
end
function (sol::SolutionBranches)(limit::SolutionCycle)
    sol.limit_cycle = limit
end
function StableSet!(sol::SolutionBranches, idx::Int=1)
    sol.stable[idx] = hcat(sol.stable[idx], sol.instable[idx][:, 1])
    sol.instable[idx] = sol.instable[idx][:, 2:end]
end
function StableRem!(sol::SolutionBranches, idx::Int=1)
    sol.instable[idx] = hcat(sol.stable[idx][:, end], sol.instable[idx])
    sol.stable[idx] = sol.stable[idx][:, 1:end-1]
end
function InstableSet!(sol::SolutionBranches, idx::Int=1)
    sol.instable[idx] = hcat(sol.instable[idx], sol.stable[idx+1][:, 1])
    sol.stable[idx+1] = sol.stable[idx+1][:, 2:end]
end
function InstableRem!(sol::SolutionBranches, idx::Int=1)
    sol.stable[idx+1] = hcat(sol.instable[idx][:, end], sol.stable[idx+1])
    sol.instable[idx] = sol.instable[idx][:, 1:end-1]
end
function Plots.plot(sol::SolutionBranches; kwargs...)
    st = sol.stable
    in = sol.instable

    p1 = plot(xlim=(450, 520), lw=2, xlabel=L"T_\mathrm{in} \  \mathrm{in \ K}", ylabel="Yield in %", framestyle=:box; kwargs...)
    for i = 1:length(st)
        plot!(st[i][2, :], st[i][1, :] * 100, lw=2, c=get_color_palette(:auto, 1)[1], label=false; kwargs...)
    end
    for i = 1:length(in)
        plot!(in[1][2, :], in[1][1, :] * 100, lw=2, c=get_color_palette(:auto, 1)[1], ls=:dash, label=false, alpha=0.5; kwargs...)
    end
    if typeof(sol.limit_cycle) == SolutionCycle
        scatter!(sol.limit_cycle.Tcycle, sol.limit_cycle.Ymax * 100, c=get_color_palette(:auto, 1)[1], marker=:circle, label=false, ; kwargs)
        scatter!(sol.limit_cycle.Tcycle, sol.limit_cycle.Ymin * 100, c=get_color_palette(:auto, 1)[1], marker=:circle, label=false, ; kwargs)
    end
    plot(p1)
end
function Plots.plot!(sol::SolutionBranches; kwargs...)
    st = sol.stable
    in = sol.instable

    p1 = plot!(xlim=(450, 520), lw=2, xlabel=L"T_\mathrm{in} \  \mathrm{in \ K}", ylabel="Yield in %", framestyle=:box; kwargs...)
    for i = 1:length(st)
        plot!(st[i][2, :], st[i][1, :] * 100, lw=2, c=get_color_palette(:auto, 1)[1], label=false; kwargs...)
    end
    for i = 1:length(in)
        plot!(in[1][2, :], in[1][1, :] * 100, lw=2, c=get_color_palette(:auto, 1)[1], ls=:dash, label=false, alpha=0.5; kwargs...)
    end
    if typeof(sol.limit_cycle) == SolutionCycle
        plot!(sol.limit_cycle.Tcycle, sol.limit_cycle.Ymax * 100, c=get_color_palette(:auto, 1)[1], lw=2, label=false; kwargs)
        plot!(sol.limit_cycle.Tcycle, sol.limit_cycle.Ymin * 100, c=get_color_palette(:auto, 1)[1], lw=2, label=false; kwargs)
    end
    plot!(p1)
end


function save_combined(p12, p34, p56, p78, name)
    savefig(p12, "p1.pdf")
    savefig(p34, "p2.pdf")
    savefig(p56, "p3.pdf")
    savefig(p78, "p4.pdf")

    combine_pdfs_standalone("p1.pdf", "p2.pdf", "p1234.pdf", width="600pt", height="600pt")
    combine_pdfs_standalone("p3.pdf", "p4.pdf", "p5678.pdf", width="600pt", height="600pt")

    combine_pdfs_standalone("p1234.pdf", "p5678.pdf", name, width="600pt", height="1201pt")
    rm("p1.pdf"; force=true)
    rm("p2.pdf"; force=true)
    rm("p3.pdf"; force=true)
    rm("p4.pdf"; force=true)
    rm("p1234.pdf"; force=true)
    rm("p5678.pdf"; force=true)
end

function save_combined(p1, p2, name; width="600pt", height="600pt")
    savefig(p1, "p1.pdf")
    savefig(p2, "p2.pdf")

    combine_pdfs_standalone("p1.pdf", "p2.pdf", name, width=width, height=height)
    rm("p1.pdf"; force=true)
    rm("p2.pdf"; force=true)
end
function combine_pdfs_standalone(pdf1::String, pdf2::String, output::String="combined.pdf"; width="600pt", height="600pt")
    tex = """
    \\documentclass[varwidth]{article}
    \\usepackage[paperwidth=$(width), paperheight=$(height), left=0cm, right=0cm, top=0cm, bottom=0cm]{geometry}
    \\usepackage{graphicx}
    \\begin{document}
    \\pagestyle{empty}
    \\noindent
    \\includegraphics{$pdf1}\\\\
    \\noindent
    \\centering
    \\includegraphics{$pdf2}
    \\end{document}
    """

    write("combine.tex", tex)

    run(`pdflatex -interaction=nonstopmode combine.tex`)
    mv("combine.pdf", output; force=true)

    # Aufräumen: Lösche Hilfsdateien
    for ext in [".tex", ".aux", ".log"]
        try
            rm("combine" * ext; force=true)
        catch e
            @warn "Konnte Datei combine$ext nicht löschen" exception = (e, catch_backtrace())
        end
    end

    println("✅ PDF wurde erfolgreich erstellt: $output (untereinander!)")
end





function shiraz()
    pN = 1.01325e5
    VN = 3.95e-6
    TN = 273.15
    R = 8.314472
    hcat1 = 1.12
    hcat1 = 7.022
    dint = 1.3e-2
    dint = 3.8e-2
    dext = 1.6e-2
    dext = 5.3e-2
    Aint = pi / 4 * dint^2
    V = Aint * hcat1
    eps = 0.39
    rho_bulk = 1134 # bezogen auf gesamte reactor volumen rho_bulk = mkat/V
    GHSV = 12000
    Vgas = V * eps
    A = dext * pi * hcat1
    k = 200.0
    mcat = V * rho_bulk
    ndotin = pN * (V) * GHSV / 60 / 60 / R / TN
    return mcat, ndotin, A, Vgas, k
end
function freiburg()
    pN = 1.01325e5
    VN = 3.95e-6
    TN = 273.15
    R = 8.314472
    hcat1 = 1.12
    dint = 1.3e-2
    dext = 1.6e-2
    Aint = pi / 4 * dint^2
    V = Aint * hcat1
    eps = 0.39
    rho_bulk = 1134 # bezogen auf gesamte reactor volumen rho_bulk = mkat/V
    GHSV = 12000
    Vgas = V * eps
    A = dext * pi * hcat1
    k = 88.5
    mcat = V * rho_bulk
    ndotin = pN * (V) * GHSV / 60 / 60 / R / TN
    return mcat, ndotin, A, Vgas, k
end







function sim_one_limit(idx, sol_Bussche, cstr)
    cstr.properties.Tin = sol_Bussche.Tv[idx]
    cstr.properties.Tc = sol_Bussche.Tv[idx]
    y0 = copy(sol_Bussche.y[idx]) * 1 #repeat(vcat(cstr.properties.yin, 0.5, 499.0674839614466, cstr.properties.ndotin), cstr.properties.num_casc)
    y0[end-1] = y0[end-1] * 1.03
    dy0 = zero.(y0)  # Anfangswerte für die Ableitungen
    tspan = (0.0, 20000.0)
    differential_vars = one.(y0)
    differential_vars[9:9:end] .= false

    prob = DAEProblem(dae_bif, dy0, y0, tspan, cstr, differential_vars=differential_vars)
    sol = solve(prob, IDA())

    prob = DAEProblem(dae_bif, dy0, sol[:, end], (0.0, 100000.0), cstr, differential_vars=differential_vars)
    sol = solve(prob, IDA())

    Yup, Ydp = yield_dyn(sol, cstr_Bussche.properties.yin, cstr_Bussche.properties.ndotin)

    prob = DAEProblem(dae_bif, dy0, sol[:, end], (0.0, 20000.0), cstr, differential_vars=differential_vars)
    sol = solve(prob, IDA())

    Yu, Yd = yield_dyn(sol, cstr_Bussche.properties.yin, cstr_Bussche.properties.ndotin)
    # if abs(Yu - Yup) > 1e-2
    #     prob = DAEProblem(dae_bif, dy0, sol[:, end], (0.0, 200000.0), cstr, differential_vars=differential_vars)
    #     sol = solve(prob, IDA())

    #     prob = DAEProblem(dae_bif, dy0, sol[:, end], (0.0, 20000.0), cstr, differential_vars=differential_vars)
    #     sol = solve(prob, IDA())
    # end
    Ty = sol_Bussche.Tv[idx]
    return (Tmax=maximum(sol[end-1, :]), Tmin=minimum(sol[end-1, :]), YU=Yu, YD=Yd, TY=Ty, states=(t=sol.t, u=sol[:, :]))
end



function limitCyclePara(sol_Bussche, cstr, i_start, i_end, step=1)


    function compute(idx)
        return sim_one_limit(idx, sol_Bussche, cstr)
    end

    res = pmap(compute, i_start:step:i_end)


    return res
end


function return_sol_cycle(res)
    YU = []
    YD = []
    TY = []
    Tmax = []
    Tmin = []
    states = []

    for idx = 1:length(res)

        push!(Tmax, res[idx].Tmax)
        push!(Tmin, res[idx].Tmin)
        push!(YU, res[idx].YU)
        push!(YD, res[idx].YD)
        push!(TY, res[idx].TY)
        push!(states, res[idx].states)
    end

    return SolutionCycle(YU, YD, Tmax, Tmin, TY, states)
end
