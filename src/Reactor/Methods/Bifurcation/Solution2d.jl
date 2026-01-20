Base.@kwdef mutable struct Solution_2d
    Tv = []
    y = []
    obj1 = []
    obj2 = []
    eig = []
    nin = []
    limit_cycle = []
end
mutable struct SolutionCycle
    Ymax
    Ymin
    Tmax
    Tmin
    Tcycle
    states
    function SolutionCycle(Ymax, Ymin, Tmax, Tmin, Tcycle, states)
        Ymaxn = []
        Yminn = []
        Tmaxn = []
        Tminn = []
        Tvn = []
        statesn = []
        if length(Ymax) < 3
            for i in 1:length(Ymax)
                if abs(Ymax[i] - Ymin[i]) > 1e-3 #|| abs(Ymax[i+1] - Ymin[i+1]) > 1e-6 || abs(Ymax[i-1] - Ymin[i-1]) > 1e-6
                    push!(Ymaxn, Ymax[i])
                    push!(Yminn, Ymin[i])
                    push!(Tmaxn, Tmax[i])
                    push!(Tminn, Tmin[i])
                    push!(Tvn, Tcycle[i])
                    push!(statesn, states[i])
                end
            end
        else
            for i in 2:length(Ymax)-1
                if abs(Ymax[i] - Ymin[i]) > 1e-3 #|| abs(Ymax[i+1] - Ymin[i+1]) > 1e-6 || abs(Ymax[i-1] - Ymin[i-1]) > 1e-6
                    push!(Ymaxn, Ymax[i])
                    push!(Yminn, Ymin[i])
                    push!(Tmaxn, Tmax[i])
                    push!(Tminn, Tmin[i])
                    push!(Tvn, Tcycle[i])
                    push!(statesn, states[i])
                end
            end
        end
        new(Ymaxn, Yminn, Tmaxn, Tminn, Tvn, statesn)
    end
end
function Plots.plot(sol::SolutionCycle; kwargs...)
    plot(sol.Tcycle, sol.Ymax, c=:black; kwargs...)
    plot!(sol.Tcycle, sol.Ymin, c=:black; kwargs...)

end
function (sol::Solution_2d)(limit::SolutionCycle)
    sol.limit_cycle = limit
end
function plot2d(p::Solution_2d)
    plot(p.Tv, p.obj2)
end

function Plots.plot(sol::Solution_2d, set=1; kwargs...)
    if set == 1
        fac = 100
        ylabel = "Yield in %"
    else
        fac = 1
        ylabel = L"T_{out} \ \mathrm{in} \ K"
    end
    st, inst, out = YieldBranches(sol)
    begin
        p4 = plot(lw=2, xlabel=L"T_\mathrm{in} \  \mathrm{in \ K}", ylabel=ylabel, framestyle=:box; kwargs...)
        for i = 1:length(st)
            plot!(st[i][2, :], st[i][set, :] * fac, c=:black, label=false; kwargs...)
        end
        for i = 1:length(inst)
            plot!(inst[i][2, :], inst[i][set, :] * fac, c=:black, ls=:dash, label=false, alpha=0.5; kwargs...)
        end
    end
    if typeof(sol.limit_cycle) == SolutionCycle
        if set == 1
            scatter!(sol.limit_cycle.Tcycle, sol.limit_cycle.Ymax * 100, c=:black, marker=:circle, markersize=1, label=false, ; kwargs)
            scatter!(sol.limit_cycle.Tcycle, sol.limit_cycle.Ymin * 100, c=:black, marker=:circle, markersize=1, label=false, ; kwargs)
        else
            scatter!(sol.limit_cycle.Tcycle, sol.limit_cycle.Tmax, c=:black, marker=:circle, markersize=1, label=false, ; kwargs)
            scatter!(sol.limit_cycle.Tcycle, sol.limit_cycle.Tmin, c=:black, marker=:circle, markersize=1, label=false, ; kwargs)
        end
    end
    p4
end
function Plots.plot!(sol::Solution_2d, set=1; kwargs...)
    if set == 1
        fac = 100
    else
        fac = 1
    end
    st, inst, out = YieldBranches(sol)
    begin
        plot!(lw=2, xlabel=L"T_\mathrm{in} \  \mathrm{in \ K}", ylabel="Yield in %", framestyle=:box; kwargs...)
        for i = 1:length(st)
            plot!(st[i][2, :], st[i][set, :] * fac, c=:black, label=false; kwargs...)
        end
        for i = 1:length(inst)
            plot!(inst[i][2, :], inst[i][set, :] * fac, c=:black, ls=:dash, label=false, alpha=0.5; kwargs...)
        end
    end
    plot!()
end
function plotEig8(sol::Solution_2d)
    eig8 = []
    for i in eachindex(sol.eig)
        push!(eig8, sol.eig[i][8])
    end
    print(length(eig8))
    t = sol.Tv
    plot(t, real.(eig8) * 1, label=false)
    plot!(t, imag.(eig8), label=false)
end
function plot_heat(hp1, hr1, Tvec, Tc0, cstr)
    p1 = plot(ylabel="Heat in W",
        xlabel="Temperature in K",
        fontfamily="Computer Modern")

    if cstr.properties.num_casc == 1
        cc = 1
    else
        cc = LinRange(1, 0.1, cstr.properties.num_casc)
    end
    for i = 1:cstr.properties.num_casc

        plot!(Tvec, hp1[:, i], ylims=(minimum([hp1]) * 0.3, maximum([hp1]) * 1.1), color=RGB(1 - cc[i], 1 - cc[i], 1 - cc[i]), linewidth=2, label=label = false)
        #plot!(Tvec, hr1[:, i], color=RGB(1, 1 - cc[i], 1 - cc[i]), linewidth=2, linestyle=:dash, label=false)

    end
    i = 1
    plot!(Tvec, hr1[:, i], color=RGB(1, 1 - cc[i], 1 - cc[i]), linewidth=2, linestyle=:dash, label=false, ylims=(minimum(hp1) * 0.95, maximum(hp1) * 1.05))
    plot!([Tc0, Tc0], [minimum(hp1) - (maximum(hp1) - minimum(hp1)) * 0.05, maximum(hp1) + (maximum(hp1) - minimum(hp1)) * 0.05], color=RGB(1, 1 - cc[i], 1 - cc[i]), linewidth=2, linestyle=:dash, label=false, ylims=(minimum(hp1) - (maximum(hp1) - minimum(hp1)) * 0.05, maximum(hp1) + (maximum(hp1) - minimum(hp1)) * 0.05), framestyle=:box)
end

function YieldBranches(sol::Solution_2d)

    eig = sol.eig
    Tv = sol.Tv
    Y = sol.obj2
    x = sol.y
    eigv = []
    out = []
    for i in eachindex(eig[1, :])
        push!(out, maximum(real.(eig[:, i])))
        push!(eigv, (eig[:, i]))
    end
    # out = [-1]
    # set = true
    # for i = 2:length(Tv)
    #     if Tv[i] > Tv[i-1] && set
    #         push!(out, -1)
    #     elseif Tv[i] <= Tv[i-1] && set
    #         push!(out, 1)
    #         set = false
    #     elseif Tv[i] < Tv[i-1] && set == false
    #         push!(out, 1)
    #     elseif Tv[i] > Tv[i-1] && set == false
    #         push!(out, -1)
    #         set = true
    #     end
    # end

    stable = []
    instable = []
    help = [Y[1], Tv[1], x[1][end-1]]
    for i = 2:length(out)
        if out[i] < 0 && out[i-1] >= 0
            help = hcat(help, [Y[i], Tv[i], x[i][end-1]])
            push!(instable, help)
            help = [Y[i], Tv[i], x[i][end-1]]
        elseif out[i] < 0
            help = hcat(help, [Y[i], Tv[i], x[i][end-1]])
        elseif out[i] >= 0 && out[i-1] < 0
            push!(stable, help)
            help = [Y[i-1], Tv[i-1], x[i-1][end-1]]
            help = hcat(help, [Y[i], Tv[i], x[i][end-1]])
        elseif out[i] >= 0
            help = hcat(help, [Y[i], Tv[i], x[i][end-1]])
        else
            error()
        end
    end
    if out[end] < 0
        push!(stable, help)
    else
        out[end] >= 0
        push!(stable, help)
    end

    return stable, instable, eigv
end