
using DynaMeth, Plots, JLD2, LaTeXStrings

Tvec = 400:1:700
Tc0 = 492.05 #for position of heat removal curve (straight)                       
yin = [0.0, 0.25, 0.0, 0.74999, 0, 0.00001] #feed pure CO_2, minimal amount of nitrogen for equilibrium

#Kinetics
kinetics = [
    (:Seidel, Seidel()),
    (:Nestler, Nestler()),
    (:Bussche, Bussche()),
    (:Brilman, Brilman()),
    (:AutoCat, AutoCat())
]

#Heat production/removal and Continuation
results = Dict()
for (name, kinetic) in kinetics
    bif = SetupBifurcationAnalysis(kinetic=kinetic)
    cstr = Reactor(bif)
    cstr.properties.yin = yin
    cstr.properties.k = 0   #adiabatic
    if name == :AutoCat
        cstr.properties.H[3] = -41210.0
        param = vcat(load_object("data/AutoCat_SS_withoutPhi_GGWGraaf_e8.jld2")[1:end-4], 1, 1, 1, 1)
        set_kinetic_param!(cstr, param)
    end

    #Heat production/removal
    hp, hr, sol = solveHeatCascade(Tvec, Tc0, cstr)
    y = Yield(sol, cstr)

    #Continuation
    sol_cont = solveContinuationSingle(100, 550, cstr, ds=1)
    s, i, out = YieldBranches(sol_cont) #determine stable and instable branches

    results[name] = (hp=hp, hr=hr, sol=sol, y=y, cont=sol_cont, s=s, i=i, out=out)
end

### Equilibrium Yield
cstr_ref = Reactor(SetupBifurcationAnalysis(kinetic=Seidel()))
cstr_ref.properties.yin = yin
Y, T = equilibrium(cstr_ref.properties.yin, cstr_ref, Tvec=Tvec)

### Plots
begin
    left_margin = 5Plots.pt
    pal = get_color_palette(:auto, 1)

    ord = [(:Seidel, "Seidel", 1), (:AutoCat, "Kortuz", 2), (:Nestler, "Nestler", 3), (:Bussche, "Bussche", 4), (:Brilman, "van Schagen", 5)]

    p1 = plot(Tvec, results[:Seidel].y * 100, label="Seidel", lw=2, c=pal[1], size=(400, 400), left_margin=left_margin)
    for (k, lab, i) in ord[2:end]
        plot!(Tvec, results[k].y * 100, label=lab, lw=2, c=pal[i])
    end
    plot!(Tvec, Y * 100 .+ 0.0, label="Equilibrium", c=:black, lw=2, xlabel="T in K", ylabel="Yield in %", ls=:dash)
    plot!(Shape([(400, -100), (453.15, -100), (453.15, 100), (400, 100)]), c=:black, alpha=0.2, label=false, framestyle=:box,
        ylim=(0, 60), xlim=(400, 700), legendposition=:topright, legendfontsize=7)

    p2 = plot(Tvec, results[:Seidel].hp / cstr_ref.properties.mcat, label=L"Q_\mathrm{pr} \ \mathrm{Seidel}", lw=2, c=pal[1], size=(400, 400))
    for (k, lab, i) in ord[2:end]
        plot!(Tvec, results[k].hp / cstr_ref.properties.mcat, label=L"Q_\mathrm{pr} \ \mathrm{%$lab}", lw=2, c=pal[i])
    end
    plot!(Tvec, results[:Seidel].hr / cstr_ref.properties.mcat, label=L"Q_\mathrm{rm} \ \mathrm{Heat \ removal}", c=:red, lw=2, linestyle=:dash,
        xlabel="T in K", ylabel=L"\mathrm{Heat \ in \ } \frac{W}{m_\mathrm{cat}}", xlim=(400, 650), ylim=(-50, 100),
        legendposition=:bottomright, legendfontsize=7)
    plot!([Tc0, Tc0], [-50, 550], label=false, c=:red, lw=2, linestyle=:dash)
    plot!(Shape([(400, -50), (453.15, -50), (453.15, 600), (400, 600)]), c=:black, alpha=0.2, label=false, framestyle=:box)
    annotate!(559, 83, Plots.text("adiabatic", 8, "Computer Modern", :red, rotation=65))
    annotate!(485, 75, Plots.text("isothermal", 8, "Computer Modern", :red, rotation=90))

    p12 = plot(p1, p2, size=(700, 350), fontfamily="Computer Modern")

    p3 = plot(xlim=(450, 550), lw=2, xlabel=L"T_\mathrm{in} \  \mathrm{in \ K}", ylabel="Yield in %", ylim=(0, 20),
        framestyle=:box, left_margin=left_margin)
    for (k, _, i) in reverse(ord)
        plot!(results[k].cont, lw=2, c=pal[i], label=false)
    end

    p4 = plot(xlim=(450, 550), lw=2, xlabel=L"T_\mathrm{in} \  \mathrm{in \ K}", ylabel="Reactor temperature in K", ylims=(430, 550), framestyle=:box)
    for (k, _, i) in reverse(ord)
        plot!(results[k].cont, 3, lw=2, c=pal[i], label=false)
    end

    p34 = plot(p3, p4, size=(700, 350), fontfamily="Computer Modern")
    p_all = plot(p12, p34, layout=@layout([a; b]), size=(600, 600), fontfamily="Computer Modern")
end