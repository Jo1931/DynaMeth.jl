mutable struct SolutionFPO
    u
    x
    obj1
    obj2
    xin
    t
    z
    term
    model
end
function SolutionFPO(nall, mall, numP, nc)
    term = []
    model = []
    for i = 1:numP
        push!(term, TerminationStatusCode(0))
        push!(model, Dict())
    end
    SolutionFPO(zeros(15, numP),
        zeros(nc, nall, mall, numP),
        zeros(numP),
        zeros(numP),
        zeros(10, nall, numP),
        zeros(nall, numP),
        zeros(mall, numP),
        term,
        model)
end

function (sol::SolutionFPO)(model, i)
    x = value.(model[:x])
    u = value.(model[:u])
    ob1 = value(model[:obj1])
    ob2 = value(model[:obj2])
    xin1 = value.(model[:xin1])
    xin2 = value.(model[:xin2])
    xin3 = value.(model[:xin3])
    xin4 = value.(model[:xin4])
    xin5 = value.(model[:xin5])
    xin6 = value.(model[:xin6])
    xin7 = nothing
    try
        xin7 = value.(model[:xin7])
    catch
        xin7 = zeros(length(xin6))
    end
    xin8 = value.(model[:xin8])
    xin9 = value.(model[:xin9])
    xin10 = nothing
    try
        xin10 = value.(model[:xin10])
    catch
        xin10 = zero(xin9)
    end
    xin = hcat(xin1, xin2, xin3, xin4, xin5, xin6, xin7, xin8, xin9, xin10)
    t = value.(model[:tspan] * model[:tend])
    z = value.(model[:zspan] * model[:zend])
    term = termination_status(model)
    d = collect_results_JuMP(model)

    sol.u[:, i] = u
    sol.x[:, :, :, i] = x
    sol.obj1[i] = ob1
    sol.obj2[i] = ob2
    sol.xin[:, :, i] = xin'
    sol.t[:, i] = t
    if length(x[1, 1, :, 1]) == 1
    else
        sol.z[:, i] = z
    end
    sol.term[i] = term
    sol.model[i] = d
end

function Plots.plot(sol::T; kwargs...) where {T<:Union{SolutionFPO,NamedTuple,JLD2.ReconstructedMutable{:SolutionFPO,(:u, :x, :obj1, :obj2, :xin, :t, :z, :term),NTuple{8,Any}}}}
    plot(sol.obj2, sol.obj1, marker=:circle, linestyle=:auto; kwargs...)
end
function Plots.plot!(sol::T; kwargs...) where {T<:Union{SolutionFPO,NamedTuple,JLD2.ReconstructedMutable{:SolutionFPO,(:u, :x, :obj1, :obj2, :xin, :t, :z, :term),NTuple{8,Any}}}}
    plot!(sol.obj2, sol.obj1, marker=:circle, linestyle=:auto; kwargs...)
end

function plot_sol(sol; kwargs...)
    default(fontfamily="Computer Modern")
    p1 = plot(ylabel="methanol flowrate", xlabel="methanol yield", legend=:bottomleft, frame=:box, size=(400, 400))
    scatter!(p1, sol.obj2, sol.obj1; kwargs...)

    #scatter!(solSS.obj2,solSS.obj1, mc=:transparent,msc = :black, m =:x, label = "steady state")
    size_p = 500
    p2 = plot(ylabel=L"y_{CO_2,in}", xlabel="methanol yield", legend=:bottomleft, frame=:box, label=false)
    scatter!(sol.obj2, sol.u[1, :], label=false; kwargs...)
    p3 = plot(ylabel=L"y_{CO,in}", xlabel="methanol yield", legend=:bottomleft, frame=:box, label=false)
    scatter!(sol.obj2, sol.u[2, :], label=false; kwargs...)
    p4 = plot(ylabel=L"y_{H_2,in}", xlabel="methanol yield", legend=:bottomleft, frame=:box, label=false)
    scatter!(sol.obj2, sol.u[3, :], label=false; kwargs...)
    p6 = plot(ylabel=L"A_f", xlabel="methanol yield", legend=:bottomleft, frame=:box, label=false)
    scatter!(sol.obj2, sol.u[4, :], label=false; kwargs...)
    p7 = plot(ylabel=L"y_{A_{CO}}", xlabel="methanol yield", legend=:bottomleft, frame=:box, label=false)
    scatter!(sol.obj2, sol.u[5, :], label=false; kwargs...)
    p8 = plot(ylabel=L"A_{N2}", xlabel="methanol yield", legend=:bottomleft, frame=:box, label=false)
    scatter!(sol.obj2, sol.u[6, :], label=false; kwargs...)
    p9 = plot(ylabel=L"\tau", xlabel="methanol yield", legend=:bottomleft, frame=:box, label=false)
    scatter!(sol.obj2, sol.u[7, :], label=false; kwargs...)
    p10 = plot(ylabel=L"\Delta \phi", xlabel="methanol yield", legend=:bottomleft, frame=:box, label=false)
    scatter!(sol.obj2, sol.u[8, :], label=false; kwargs...)
    p5 = plot(ylabel=L"y_{N2_2,in}", xlabel="methanol yield", legend=:bottomleft, frame=:box, label=false)
    scatter!(sol.obj2, sol.u[9, :], label=false; kwargs...)
    return plot(p1, plot(p2, p3, p4, p5, p6, p7, p8, p9, p10), size=(1200, 600))
end
function plot_sol!(p, sol; kwargs...)
    scatter!(p[1], sol.obj2, sol.obj1, label=false; kwargs...)

    #scatter!(solSS.obj2,solSS.obj1, mc=:transparent,msc = :black, m =:x, label = "steady state")
    scatter!(p[2], sol.obj2, sol.u[1, :], label=false; kwargs...)
    scatter!(p[3], sol.obj2, sol.u[2, :], label=false; kwargs...)
    scatter!(p[4], sol.obj2, sol.u[3, :], label=false; kwargs...)
    scatter!(p[5], sol.obj2, sol.u[9, :], label=false; kwargs...)
    scatter!(p[6], sol.obj2, sol.u[4, :], label=false; kwargs...)
    scatter!(p[7], sol.obj2, sol.u[5, :], label=false; kwargs...)
    scatter!(p[8], sol.obj2, sol.u[6, :], label=false; kwargs...)
    scatter!(p[9], sol.obj2, sol.u[7, :], label=false; kwargs...)
    scatter!(p[10], sol.obj2, sol.u[8, :], label=false; kwargs...)
    return p
end

function Base.copy(sol::SolutionFPO)
    SolutionFPO(copy(sol.u), copy(sol.x), copy(sol.obj1), copy(sol.obj2), copy(sol.xin), copy(sol.t), copy(sol.z), copy(sol.term), copy(sol.model))
end
function convert_sol_to_tupel(sol::SolutionFPO)
    (
        u=copy(sol.u),
        x=copy(sol.x),
        obj1=copy(sol.obj1),
        obj2=copy(sol.obj2),
        xin=copy(sol.xin),
        t=copy(sol.t),
        z=copy(sol.z),
        term=copy(sol.term),
        model=copy(sol.model)
    )
end
function convert_tuple_to_sol(sol::NamedTuple)
    SolutionFPO(copy(sol.u), copy(sol.x), copy(sol.obj1), copy(sol.obj2), copy(sol.xin), copy(sol.t), copy(sol.z), copy(sol.term), copy(sol.model))
end



