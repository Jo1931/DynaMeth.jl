Base.@kwdef mutable struct Solution_3d
    Tm = []
    obj1m = []
    obj2m = []
    eigm = []
    ym = []
    ninm = []
    sol2d = []
end
function (p::Solution_3d)(sol::Solution_2d)
    push!(p.Tm, sol.Tv)
    push!(p.obj1m, sol.obj1)
    push!(p.obj2m, sol.obj2)
    push!(p.eigm, sol.eig)
    push!(p.ym, sol.y)
    push!(p.ninm, sol.nin)
    push!(p.sol2d, sol)
end
function area(p::Solution_3d)
    index_r = zeros(Int, length(p.Tm))
    index_l = zeros(Int, length(p.Tm))
    nin_r = zeros(length(p.Tm))
    nin_l = zeros(length(p.Tm))
    T_r = zeros(length(p.Tm))
    T_l = zeros(length(p.Tm))
    obj2_r = zeros(length(p.Tm))
    obj2_l = zeros(length(p.Tm))
    for j in eachindex(p.Tm)
        Tv = p.Tm[j]
        set = true
        for i in eachindex(Tv[1:end-1])
            if Tv[i] > Tv[i+1] && set == true
                set = false
                index_r[j] = i + 1
            end
            if Tv[i] < Tv[i+1] && set == false
                set = true
                index_l[j] = i + 1
            end
        end
        if index_r[j] > 0
            nin_r[j] = p.ninm[j][index_r[j]]
            nin_l[j] = p.ninm[j][index_l[j]]
            T_r[j] = p.Tm[j][index_r[j]]
            T_l[j] = p.Tm[j][index_l[j]]
        end
    end

    return T_r[T_r.>0], T_l[T_r.>0], nin_r[T_r.>0], nin_l[T_r.>0], obj2_r[T_r.>0], obj2_l[T_r.>0]
end
function Plots.plot(sol::Solution_3d, mcat; kwargs...)
    idx = 1
    p1 = plot(lw=2, xlabel=L"T_\mathrm{in} \  \mathrm{in \ K}", ylabel="Yield in %"; kwargs...)
    for idx = 1:length(sol.Tm)
        fac = 100
        st, inst, out = YieldBranches(sol.sol2d[idx])
        begin
            for i = 1:length(st)
                plot!(st[i][2, :], sol.ninm[idx][1:length(st[i][2, :])] / mcat, st[i][1, :] * fac, c=:black, label=false; kwargs...)
            end
            for i = 1:length(inst)
                plot!(inst[i][2, :], sol.ninm[idx][1:length(inst[i][2, :])] / mcat, inst[i][1, :] * fac, c=:black, ls=:dash, label=false, alpha=0.5; kwargs...)
            end
        end

    end
    p1
    T_r, T_l, nin_r, nin_l, obj2_r, obj2_l = area(sol)
    plot!(T_r, nin_r / mcat, obj2_r, color=:black; kwargs...)
    plot!(T_l, nin_l / mcat, obj2_l, color=:black; kwargs...)
    plot!(legend=false)
end
function plotEig8(sol::Solution_3d, idx1)
    eig8 = []
    for i in eachindex(sol.eigm[idx1])
        push!(eig8, sol.eigm[idx1][i][end])
    end
    print(length(eig8))
    t = sol.Tm[idx1]
    plot(t, real.(eig8) * 1, label=false)
    plot!(t, imag.(eig8), label=false)
end


function plot_area(sol::Solution_3d, mcat; kwargs...)

    T_r, T_l, nin_r, nin_l, obj2_r, obj2_l = area(sol)
    plot(T_r, nin_r / mcat, color=:black; kwargs...)
    plot!(T_l, nin_l / mcat, color=:black; kwargs...)
    plot!(legend=false)
end