@kwdef mutable struct SolutionMPC
    t = nothing
    tu = nothing
    x = nothing
    n_co = nothing
    n_co2 = nothing
    n_h2 = nothing
    obj1 = nothing
    obj2 = nothing
    inputs = []
    para = nothing
    setpoint = nothing
end
deadt = 0
function (sol::SolutionMPC)(t, x, Y, cstr)
    nin = cstr.properties.ndotin
    yin = cstr.properties.yin
    tu = t[1]:1:t[end]
    if isnothing(sol.t)
        sol.t = t
        sol.tu = tu
        sol.x = x
        sol.n_co = cstr.methods.input.co(tu .+ deadt)
        sol.n_co2 = cstr.methods.input.co2(tu .+ deadt)
        sol.n_h2 = nin * yin[4] * (1 .+ cstr.methods.disturb(tu))
        sol.obj1 = x[1, :] .* x[8, :]
        sol.obj2 = Y
        sol.setpoint = [cstr.methods.setpoint, cstr.methods.setpoint]
    else
        sol.t = vcat(sol.t, t[2:end])
        sol.tu = vcat(sol.tu, tu[2:end])
        sol.x = hcat(sol.x, x[:, 2:end])
        sol.n_co = vcat(sol.n_co, cstr.methods.input.co(tu[2:end] .+ deadt))
        sol.n_co2 = vcat(sol.n_co2, cstr.methods.input.co2(tu[2:end] .+ deadt))
        sol.n_h2 = vcat(sol.n_h2, nin * yin[4] * (1 .+ cstr.methods.disturb(tu[2:end])))
        sol.obj1 = vcat(sol.obj1, x[1, 2:end] .* x[8, 2:end])
        sol.obj2 = vcat(sol.obj2, Y[2:end])
        sol.setpoint = vcat(sol.setpoint, cstr.methods.setpoint)
    end
end

function (sol::SolutionMPC)(model::T) where {T<:AbstractModel}
    push!(sol.inputs, model)
end
function (sol::SolutionMPC)(para)
    if isnothing(sol.para)
        sol.para = para
    else
        sol.para = hcat(sol.para, para)
    end
end

function (sol::SolutionMPC)(model::T, contr) where {T<:AbstractModel}
    model[:u] = value.(model[:u]) * (1 .+ contr)
    push!(sol.inputs, model)
end

function plot_solMPC(solMPC::SolutionMPC)

    t1 = value.(solMPC.inputs[1][:ts])
    x1 = value.(solMPC.inputs[1][:x])
    ob11 = value.(solMPC.inputs[1][:obj1])
    ob21 = value.(solMPC.inputs[1][:obj2])

    tc = t1 .< solMPC.t[2]
    #tc = 1:length(t1)
    tjum = t1[tc]
    xjum = x1[:, tc]
    ob1jum = ob11[tc]
    ob2jum = ob21[tc]

    for i = 1:length(solMPC.inputs)
        t1 = value.(solMPC.inputs[i][:ts])
        x1 = value.(solMPC.inputs[i][:x])
        ob11 = value.(solMPC.inputs[i][:obj1])
        ob21 = value.(solMPC.inputs[i][:obj2])

        tc = t1 .< solMPC.t[i+1]
        #tc = 1:length(t1)
        tjum = vcat(tjum, t1[tc])
        xjum = hcat(xjum, x1[:, tc])
        ob1jum = vcat(ob1jum, ob11[tc])
        ob2jum = vcat(ob2jum, ob21[tc])
    end




    p1 = plot(solMPC.t / 60, solMPC.setpoint * 100, ls=:dash, c=:black, label=false, xlims=(minimum(solMPC.t / 60), maximum(solMPC.t / 60)), ylims=(0, 4))# ylims=(minimum(solMPC.obj2 * 100 * 1.1), maximum(solMPC.obj2 * 100 * 1.1))
    #plot!(tjum / 60, ob2jum * 100, label=false)
    scatter!(solMPC.t / 60, solMPC.obj2 * 100, c=:blue, marker=:x, ylabel=L"Yield \ in \  \%", label=false, ytickfontcolor=:blue, yguidefontcolor=:blue)
    ax2 = twinx()
    scatter!(ax2, solMPC.t / 60, solMPC.obj1 * 6e4, c=:red, marker=:x, ylabel=L"\dot{n}_{CH_3OH} \ in \ \frac{mmol}{min}", label=false, ytickfontcolor=:red, guidefontcolor=:red, xlims=(minimum(solMPC.t / 60), maximum(solMPC.t / 60)))
    #plot!(ax2, tjum / 60, ob1jum, label=false)

    p2 = plot(solMPC.tu / 60, solMPC.n_co * 6e4, label=L"CO", c=:blue, left_margin=30Plots.px, bottom_margin=30Plots.px, xlabel=L"t \ in \ min", right_margin=20Plots.mm)
    plot!(solMPC.tu / 60, solMPC.n_co2 * 6e4, label=L"CO_2", ylabel=L"\dot{n}_{in} \ in \ \frac{mmol}{min}", c=:green, xlims=(minimum(solMPC.t / 60), maximum(solMPC.t / 60)))
    p3 = plot(solMPC.tu / 60, solMPC.n_h2 * 6e4, label=false, ylabel=L"\dot{n}_{in,H_2} \ in \ \frac{mmol}{min}", c=:red, left_margin=30Plots.px, xlims=(minimum(solMPC.t / 60), maximum(solMPC.t / 60)))
    p4 = plot(solMPC.t / 60, solMPC.para[1, :], c=:blue, label=false, ylabel=L"P", left_margin=30Plots.px, bottom_margin=30Plots.px, xlabel=L"t \ in \ min")
    plot!(solMPC.t / 60, solMPC.para[2, :], c=:blue, label=false, ylabel=L"P", left_margin=30Plots.px, bottom_margin=30Plots.px)
    plot!(solMPC.t / 60, solMPC.para[3, :], c=:blue, label=false, ylabel=L"P", left_margin=30Plots.px, bottom_margin=30Plots.px)

    p23 = plot(p3, p2, layout=(2, 1))



    p4 = plot(p1, p23, layout=(2, 1), size=(600, 500))

    ms = 2
    mg = 0.0Plots.mm
    mg1 = 15.0Plots.mm
    p5 = scatter(solMPC.t / 60, solMPC.x[1, :] * 100, m=:x, label=false, ms=ms, c=:blue, ylabel=L"y_{CH_3OH} \ in \ \%", left_margin=40Plots.px)
    plot!(tjum / 60, xjum[1, :] * 100, c=:blue, label=false, right_margin=mg, left_margin=mg1, xlims=(minimum(solMPC.t / 60), maximum(solMPC.t / 60)))
    p6 = scatter(solMPC.t / 60, solMPC.x[2, :] * 100, m=:x, label=false, ms=ms, c=:blue, ylabel=L"y_{CO_2} \ in \  \%")
    plot!(tjum / 60, xjum[2, :] * 100, c=:blue, label=false, left_margin=mg, ymirror=true, right_margin=mg1, xlims=(minimum(solMPC.t / 60), maximum(solMPC.t / 60)))
    p7 = scatter(solMPC.t / 60, solMPC.x[3, :] * 100, m=:x, label=false, ms=ms, c=:blue, ylabel=L"y_{CO} \ in \  \%", left_margin=40Plots.px)
    plot!(tjum / 60, xjum[3, :] * 100, c=:blue, label=false, right_margin=mg, left_margin=mg1, xlims=(minimum(solMPC.t / 60), maximum(solMPC.t / 60)))
    p8 = scatter(solMPC.t / 60, solMPC.x[4, :] * 100, m=:x, label=false, ms=ms, c=:blue, ylabel=L"y_{H_2} \ in \  \%")
    plot!(tjum / 60, xjum[4, :] * 100, c=:blue, label=false, left_margin=mg, ymirror=true, right_margin=mg1, xlims=(minimum(solMPC.t / 60), maximum(solMPC.t / 60)))
    p9 = scatter(solMPC.t / 60, solMPC.x[5, :] * 100, m=:x, label=false, ms=ms, c=:blue, ylabel=L"y_{H_2O} \ in \  \%", left_margin=40Plots.px)
    plot!(tjum / 60, xjum[5, :] * 100, c=:blue, label=false, right_margin=mg, left_margin=mg1, xlims=(minimum(solMPC.t / 60), maximum(solMPC.t / 60)))
    p10 = scatter(solMPC.t / 60, solMPC.x[6, :] * 100, m=:x, label=false, ms=ms, c=:blue, ylabel=L"y_{N_2} \ in \  \%")
    plot!(tjum / 60, xjum[6, :] * 100, c=:blue, label=false, left_margin=mg, ymirror=true, right_margin=mg1, xlims=(minimum(solMPC.t / 60), maximum(solMPC.t / 60)))
    p11 = scatter(solMPC.t / 60, solMPC.x[7, :], m=:x, label=false, ms=ms, xlabel=L"t \ in  \ min", c=:blue, ylabel=L"\phi", left_margin=30Plots.px, bottom_margin=40Plots.px)
    plot!(tjum / 60, xjum[7, :], c=:blue, label=false, right_margin=mg, left_margin=mg1, xlims=(minimum(solMPC.t / 60), maximum(solMPC.t / 60)))
    p12 = scatter(solMPC.t / 60, solMPC.x[8, :] * 6e4, m=:x, label=false, ms=ms, xlabel=L"t \ in \ min", c=:blue, ylabel=L"\dot{n} \ in \ \frac{mmol}{min}", bottom_margin=30Plots.px)
    plot!(tjum / 60, xjum[8, :] * 6e4, c=:blue, label=false, left_margin=mg, ymirror=true, right_margin=mg1, xlims=(minimum(solMPC.t / 60), maximum(solMPC.t / 60)))

    p512 = plot(p5, p6, p7, p8, p9, p10, p11, p12, layout=(4, 2))

    plot(p4, p512, size=(1400, 800), framestyle=:box, guidefontsize=16, tickfont=12, legendfontsize=12, legend=:bottomleft)
end


