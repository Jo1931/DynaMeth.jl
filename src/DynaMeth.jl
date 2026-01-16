module DynaMeth
#__precompile__(false)
using UnPack
using MAT
using CSV
using DataFrames
using Plots
using Dierckx
using NonlinearSolve
using DifferentialEquations
using LinearAlgebra
using ForwardDiff
using Sundials
using JuMP
using ComponentArrays
using LaTeXStrings
using JLD2
using Ipopt
using FileIO
using Statistics
using PGFPlotsX



abstract type AbstractReactor end

abstract type AbstractSetup <: AbstractReactor end
abstract type AbstractExperimental <: AbstractSetup end
abstract type AbstractFPO <: AbstractSetup end
abstract type AbstractBifurcation <: AbstractSetup end
abstract type AbstractMPC <: AbstractSetup end

abstract type AbstractKinetic <: AbstractReactor end
abstract type AbstractKineticComponent <: AbstractKinetic end
abstract type AbstractKineticSetup <: AbstractSetup end

abstract type AbstractParameterStruct <: AbstractKineticComponent end

function Dierckx.Spline1D(x::T) where {T<:Float64}
    t = [0.0, 4000000.0]
    y = [1.0, 1.0] * x
    return Spline1D(t, y, k=1)
end

function _run_scripts_in_subfolder(path; seen=Set{String}())
    for entry in sort(readdir(path))
        full = joinpath(path, entry)
        if isdir(full)
            _run_scripts_in_subfolder(full; seen=seen)
        elseif endswith(entry, ".jl")
            rp = realpath(full)
            rp in seen && continue
            push!(seen, rp)
            include(rp)
        end
    end
end

_protected_sqrt(x::T) where {T} =# real.(sqrt.(Complex.(x)))
    @inline myswish_0(x) = x / (1 + exp(-x))
@inline mysoftplus_0(x) = log(1 + exp(x))


dir = @__DIR__
path_n = dir * "/Reactor"

_run_scripts_in_subfolder(path_n)

let syms = Symbol[]
    for n in names(@__MODULE__; all=true, imported=false)
        s = String(n)
        startswith(s, "_") && continue
        startswith(s, "#") && continue      # <-- filtert #eval, ##meta, #..#..
        n in (:DynaMeth, :include, :eval) && continue  # extra safety

        isdefined(@__MODULE__, n) || continue
        v = getfield(@__MODULE__, n)
        (v isa Function || v isa DataType) || continue

        push!(syms, n)
    end
    eval(Expr(:export, syms...))
end


end