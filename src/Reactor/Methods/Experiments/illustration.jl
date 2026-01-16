function plotReactorSteady(cstr::Reactor)
    xsim, xmess = solveSS(cstr::Reactor)


    sum(((xsim[1:3, :] - xmess[1:3, :])) .^ 2)

    p1 = scatter(xmess[1, :], xsim[1, :], label=false, title="CH3OH")#, markersize=1)
    plot!([0, 1.1], [0, 1.1], label=false, c=:black, xlims=(minimum(xmess[1, :] * 0) * 0.8, maximum(xmess[1, :]) * 1.1), ylims=(minimum(xmess[1, :] * 0) * 0.8, maximum(xmess[1, :]) * 1.1))
    p2 = scatter(xmess[2, :], xsim[2, :], label=false, title="CO2")
    plot!([0, 1.3], [0, 1.3], label=false, c=:black, xlims=(minimum(xmess[2, :]) * 0.8, maximum(xmess[2, :]) * 1.1), ylims=(minimum(xmess[2, :]) * 0.8, maximum(xmess[2, :]) * 1.1))
    p3 = scatter(xmess[3, :], xsim[3, :], label=false, title="CO")
    plot!([0, 1.3], [0, 1.3], label=false, c=:black, xlims=(minimum(xmess[3, :]) * 0.8, maximum(xmess[3, :]) * 1.1), ylims=(minimum(xmess[3, :]) * 0.8, maximum(xmess[3, :]) * 1.1))
    p4 = scatter(xmess[4, :], xsim[4, :], label=false, title="H2")
    plot!([0.3, 1.7], [0.3, 1.7], label=false, c=:black, xlims=(minimum(xmess[4, :]) * 0.8, maximum(xmess[4, :]) * 1.1), ylims=(minimum(xmess[4, :]) * 0.8, maximum(xmess[4, :]) * 1.1))
    p5 = scatter(xmess[5, :], xsim[5, :], label=false, title="H2O")
    plot!([0.0, 0.115], [0.0, 0.115], c=:black, label=false, xlims=(minimum(xmess[5, :]) * 0.8, maximum(xmess[5, :]) * 1.1), ylims=(minimum(xmess[5, :]) * 0.8, maximum(xmess[5, :]) * 1.1))
    p6 = scatter(xmess[6, :], xsim[6, :], label=false, title="N2")
    plot!([0.1, 1.5], [0.1, 1.5], c=:black, label=false, xlims=(minimum(xmess[6, :]) * 0.8, maximum(xmess[6, :]) * 1.1), ylims=(minimum(xmess[6, :]) * 0.8, maximum(xmess[6, :]) * 1.1))
    return plot(p1, p2, p3, p4, p5, p6)
end
function plotReactorSteady!(p, cstr::Reactor)
    xsim, xmess = solveSS(cstr::Reactor)


    sum(((xsim[1:3, :] - xmess[1:3, :])) .^ 2)

    scatter!(p[1], xmess[1, :], xsim[1, :], c=:red, label=false, title="CH3OH")#, markersize=1)
    scatter!(p[2], xmess[2, :], xsim[2, :], c=:red, label=false, title="CO2")
    scatter!(p[3], xmess[3, :], xsim[3, :], c=:red, label=false, title="CO")
    scatter!(p[4], xmess[4, :], xsim[4, :], c=:red, label=false, title="H2")
    scatter!(p[5], xmess[5, :], xsim[5, :], c=:red, label=false, title="H2O")
    scatter!(p[6], xmess[6, :], xsim[6, :], c=:red, label=false, title="N2")

end

function plotDynamic(cstr::Reactor, file)

    ts, x_outs, sols, xins = simulateExperiments(cstr, file)
    idx = cstr.experiments.idx
    ## xout Plot
    p1 = plot(ts / 60, sols[1, :], label=false, title="CH3OH")
    plot!(ts / 60, x_outs["CH3OH"][idx:end], label=false, ylims=(minimum([minimum(sols[1, :]), minimum(x_outs["CH3OH"])]), maximum([maximum(sols[1, :]), maximum(x_outs["CH3OH"])])))

    p2 = plot(ts / 60, sols[2, :], label=false, title="CO2")
    plot!(ts / 60, x_outs["CO2"][idx:end], label=false, ylims=(minimum([minimum(sols[2, :]), minimum(x_outs["CO2"])]), maximum([maximum(sols[2, :]), maximum(x_outs["CO2"])])))

    p3 = plot(ts / 60, sols[3, :], label=false, title="CO")
    plot!(ts / 60, x_outs["CO"][idx:end], label=false, ylims=(minimum([minimum(sols[3, :]), minimum(x_outs["CO"])]), maximum([maximum(sols[3, :]), maximum(x_outs["CO"])])))

    p4 = plot(ts / 60, sols[4, :], label=false, title="H2")
    plot!(ts / 60, x_outs["H2"][idx:end], label=false, ylims=(minimum([minimum(sols[4, :]), minimum(x_outs["H2"])]), maximum([maximum(sols[4, :]), maximum(x_outs["H2"])])))

    p5 = plot(ts / 60, sols[5, :], label=false, title="H2O")
    plot!(ts / 60, x_outs["H2O"][idx:end], label=false, ylims=(minimum([minimum(sols[5, :]), minimum(x_outs["H2O"])]), maximum([maximum(sols[5, :]), maximum(x_outs["H2O"])])))

    p6 = plot(ts / 60, sols[6, :], label=false, title="N2")
    plot!(ts / 60, x_outs["N2"][idx:end], label=false, ylims=(minimum([minimum(sols[6, :]), minimum(x_outs["N2"])]), maximum([maximum(sols[6, :]), maximum(x_outs["N2"])])))
    return plot(p1, p2, p3, p4, p5, p6)
end

function plotDynamicPa(cstr::Reactor, files)

    tv, xoutv, solv, xinv, err = simulateAllExperiments(cstr, files)
    idx = cstr.experiments.idx
    p1 = scatter()
    for i in eachindex(files)
        scatter!(xoutv[i]["CH3OH"][idx:end], solv[i][1, :], label=false, title="CH3OH", markersize=3, color=:black, marker=:cross)
    end
    plot!([0, 1.1], [0, 1.1] * 1, label=false, xlims=(0, 0.05), ylims=(0, 0.05), color=:red)

    p2 = scatter()
    for i in eachindex(files)
        scatter!(xoutv[i]["CO2"][idx:end], solv[i][2, :], label=false, title="CO2", markersize=3, color=:black, marker=:cross)
    end
    plot!([0, 1.1], [0, 1.1] * 1, label=false, xlims=(0, 0.2), ylims=(0, 0.2), color=:red)

    p3 = scatter()
    for i in eachindex(files)
        scatter!(xoutv[i]["CO"][idx:end], solv[i][3, :], label=false, title="CO", markersize=3, markercolor=:black, marker=:cross)
    end
    plot!([0, 1.1], [0, 1.1] * 1, label=false, xlims=(0, 0.3), ylims=(0, 0.3), color=:red)

    p4 = scatter()
    for i in eachindex(files)
        scatter!(xoutv[i]["H2"][idx:end], solv[i][4, :], label=false, title="H2", markersize=3, markercolor=:black, marker=:cross)
    end
    plot!([0, 1.1], [0, 1.1] * 1, label=false, xlims=(0, 0.7), ylims=(0, 0.7), color=:red)

    p5 = scatter()
    for i in eachindex(files)
        scatter!(xoutv[i]["H2O"][idx:end], solv[i][5, :], label=false, title="H2O", markersize=3, markercolor=:black, marker=:cross)
    end
    plot!([0, 1.1], [0, 1.1] * 1, label=false, xlims=(0, 0.015), ylims=(0, 0.015), color=:red)

    p6 = scatter()
    for i in eachindex(files)
        scatter!(xoutv[i]["N2"][idx:end], solv[i][6, :], label=false, title="N2", markersize=3, markercolor=:black, marker=:cross)
    end
    plot!([0, 1.1], [0, 1.1] * 1, label=false, xlims=(0, 0.7), ylims=(0, 0.7), color=:red)

    plot(p1, p2, p3, p4, p5, p6)
end
function plotDynamicVollbrecht(cstr)
    ts, x_out, sol, xin = simulateVollbrechtDyn(cstr)

    p1 = plot(ts / 60, sol[1, :], label=false, ylabel=L"y_{CH_3OH}")
    scatter!(x_out[!, "dynamic_1"], x_out[!, "dynamic_2"], label=false, marker=:x, markersize=3, color=:black)
    p2 = plot(ts / 60, sol[2, :], label=false, ylabel=L"y_{CO_2}")
    scatter!(x_out[!, "dynamic_1"], x_out[!, "dynamic_3"], label=false, marker=:x, markersize=3, color=:black)
    p3 = plot(ts / 60, sol[3, :], ylabel=L"y_{CO}", label=false)
    scatter!(x_out[!, "dynamic_1"], x_out[!, "dynamic_4"], label=false, marker=:x, markersize=3, color=:black)
    p4 = plot(ts / 60, sol[4, :], label=false, ylabel=L"y_{H_2}")
    p5 = plot(ts / 60, sol[5, :], label=false, ylabel=L"y_{H_2O}", xlabel="t in min")
    scatter!(x_out[!, "dynamic_1"], x_out[!, "dynamic_6"], label=false, marker=:x, markersize=3, color=:black)
    p6 = plot(ts / 60, sol[6, :], label=false, ylabel=L"y_{N_2}", xlabel="t in min")
    plot(p1, p2, p3, p4, p5, p6, layout=(3, 2))
end


# cstr = Reactor(SetupVollbrecht(xin_mess=xin, xout_mess=xout, kinetic=Brilman(), files=dynamic))

# cstr.properties.qsat=0.0
# ts, x_out, sol, xin = simulateVollbrechtDyn(cstr);
# cstr.properties.qsat=0
# ts1, x_out1, sol1, xin1 = simulateVollbrechtDyn(cstr);


# l = @layout [a b ; b c;e f;b{0.1h}]

# p1 = plot(ts / 60, sol[1, :], label=false, ylabel=L"y_{CH_3OH}")
# #plot!(ts1 / 60, sol1[1, :], label=false,c=:red, ylabel=L"y_{CH_3OH}")
# scatter!(x_out[!, "dynamic_1"], x_out[!, "dynamic_2"], label=false, marker=:x, markersize=3, color=:black)
# p2=plot(ts / 60, sol[2, :], label=false, ylabel=L"y_{CO_2}")
# #plot!(ts1 / 60, sol1[2, :],c=:red, label=false, ylabel=L"y_{CO_2}")
# scatter!(x_out[!, "dynamic_1"], x_out[!, "dynamic_3"], label=false, marker=:x, markersize=3, color=:black)
# p3 = plot(ts / 60, sol[3, :], ylabel=L"y_{CO}", label=false)
# #plot!(ts1 / 60, sol1[3, :],c=:red, ylabel=L"y_{CO}", label=false)
# scatter!(x_out[!, "dynamic_1"], x_out[!, "dynamic_4"], label=false, marker=:x, markersize=3, color=:black)
# p4 = plot(ts / 60, sol[4, :], label=false, ylabel=L"y_{H_2}")
# #plot!(ts1 / 60, sol1[4, :],c=:red, label=false, ylabel=L"y_{H_2}")
# p5 = plot(ts / 60, sol[5, :], label=false, ylabel=L"y_{H_2O}", xlabel=L"t \ in \ min")
# #plot!(ts1 / 60, sol1[5, :],c=:red, label=false, ylabel=L"y_{H_2O}", xlabel=L"t \ in \ min")
# scatter!(x_out[!, "dynamic_1"], x_out[!, "dynamic_6"], label=false, marker=:x, markersize=3, color=:black)
# p6 = plot(ts / 60, sol[6, :], label=false, ylabel=L"y_{N_2}", xlabel=L"t \ in \ min")
# #plot!(ts1 / 60, sol1[6, :],c=:red, label=false, ylabel=L"y_{N_2}", xlabel=L"t \ in \ min")

# p7=plot(ts / 60, sol[5, :], label=L"Simulation:", ylabel=L"y_{H_2O}", xlabel=L"t \ in \ min",legend_columns=2)

# plot!([], [], grid=false, framestyle=:none, axis=false, background_color=:white,label=false,legend=:top)
# #plot!(ts1 / 60, sol1[5, :],c=:red,  label=L"Simulation: \ q_{sat}=0", ylabel=L"y_{H_2O}", xlabel=L"t \ in \ min")
# scatter!(x_out[!, "dynamic_1"], x_out[!, "dynamic_6"], label=L"Experiments", marker=:x, markersize=3, color=:black,xlims=(-10,0))

# plot(p1, p2, p3, p4, p5, p6,p7, layout=l,size = (600,400))



# # Deine Plots mit Labels versehen
# p1 = plot(ts / 60, sol[1, :], label="Simulation 1", ylabel=L"y_{CH_3OH}")
# plot!(ts1 / 60, sol1[1, :], label="Simulation 2", c=:red)
# scatter!(x_out[!, "dynamic_1"], x_out[!, "dynamic_2"], label="Messwerte", marker=:x, markersize=3, color=:black)

# p2 = plot(ts / 60, sol[2, :], label=false, ylabel=L"y_{CO_2}")
# plot!(ts1 / 60, sol1[2, :], c=:red, label=false)
# scatter!(x_out[!, "dynamic_1"], x_out[!, "dynamic_3"], label=false, marker=:x, markersize=3, color=:black)

# p3 = plot(ts / 60, sol[3, :], ylabel=L"y_{CO}", label=false)
# plot!(ts1 / 60, sol1[3, :], c=:red, label=false)
# scatter!(x_out[!, "dynamic_1"], x_out[!, "dynamic_4"], label=false, marker=:x, markersize=3, color=:black)

# p4 = plot(ts / 60, sol[4, :], label=false, ylabel=L"y_{H_2}")
# plot!(ts1 / 60, sol1[4, :], c=:red, label=false)

# p5 = plot(ts / 60, sol[5, :], ylabel=L"y_{H_2O}", xlabel=L"t \ in \ min", label=false)
# plot!(ts1 / 60, sol1[5, :], c=:red, label=false)
# scatter!(x_out[!, "dynamic_1"], x_out[!, "dynamic_6"], label=false, marker=:x, markersize=3, color=:black)

# p6 = plot(ts / 60, sol[6, :], ylabel=L"y_{N_2}", xlabel=L"t \ in \ min", label=false)
# plot!(ts1 / 60, sol1[6, :], c=:red, label=false)
# p7=plot(p5, p6, layout=(1,2), legend=:bottom)
# # Gesamte Plotanordnung mit externer Legende unten
# plot(p1, p2, p3, p4, p5, p6, layout=(4, 2), legend=:outerbottom)