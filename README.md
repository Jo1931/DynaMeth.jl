# DynaMeth.jl

DynaMeth.jl is a Julia package developed as part of my PhD thesis on dynamic reactor operation for methanol synthesis.
It provides the core functions used to generate the results presented in the following publications:

- **Leipold, J., C. Seidel, D. Nikolic, A. Seidel-Morgenstern, and A. Kienle.** *Optimization of Methanol Synthesis under Forced Periodic Operation in Isothermal Fixed-Bed Reactors.* Computers & Chemical Engineering 175 (July 2023): 108285. https://doi.org/10.1016/j.compchemeng.2023.108285.

- **Leipold, J., D. Nikolic, A. Seidel-Morgenstern, and A. Kienle.** *Multi-Objective Optimization of Forced Periodic Operation of Methanol Synthesis in a Fixed-Bed Reactor.* In 34th European Symposium on Computer Aided Process Engineering (June 2024). https://doi.org/10.1016/B978-0-443-28824-1.50268-4.

- **Leipold, J., M. Jung, T. Keßler, and A. Kienle.** *Nonlinear Behavior of Methanol Synthesis Compared to CO2 Methanation.* Chemical Engineering & Technology 47, no. 3 (November 2024): 531–36. https://doi.org/10.1002/ceat.202300256.

- **Leipold, J., D. Nikolic, A. Seidel-Morgenstern, and A. Kienle.** *Optimization of Methanol Synthesis under Forced Periodic Operation in a Non-Isothermal Fixed-Bed Reactor.* Computers & Chemical Engineering 196 (May 2025): 109040. https://doi.org/10.1016/j.compchemeng.2025.109040.

- **Kortuz, W., J. Leipold, A. Kienle, and A. Seidel-Morgenstern.** *Kinetic Modeling of the Methanol-Assisted Autocatalytic Methanol Synthesis on Cu/ZnO/Al2O3.* Chemical Engineering Journal 518, (August 2025): 164505. https://doi.org/10.1016/j.cej.2025.164505.

- **Leipold, J., and A. Kienle.** *Nonlinear Behavior of Methanol Synthesis through CO2-Hydrogenation.* Chemical Engineering Journal 522 (October 2025): 167176. https://doi.org/10.1016/j.cej.2025.167176.

- **Kaps, L., W. Kortuz, J. Leipold, et al.** *Forced Periodic Reactor Operation Applied to Methanol Synthesis.* ChemCatChem n/a, no. n/a (n.d.): e01403. https://doi.org/10.1002/cctc.202501403.

## Installation

The package can be installed via:
```julia
using Pkg
Pkg.add(url="https://github.com/Jo1931/DynaMeth.jl")

```

## Basic Usage

The central object of this package is the Reactor struct.
A reactor is defined by its kinetic formulation (`kinetic`), model parameters (`properties`), and numerical methods (`methods`).

```julia
struct Reactor <: AbstractReactor
    properties
    kinetic
    methods
end
```
Depending on the chosen kinetic model and numerical methods, the package can be used to address different research questions.

### Parameter Estimation

A minimal example demonstrating how the reactor is generated and how model outputs are compared with experimental data is shown below.

```julia
using DynaMeth
using Plots

#Load Experimental Data
xin, xout = createSteadyStates("Volllbrecht")
dynamic = DynamicVollbrecht()

#Setup Reactor
vollbrecht = SetupVollbrecht(xin_mess=xin, xout_mess=xout, kinetic=Seidel(), files=dynamic)
cstr = Reactor(vollbrecht)

#Plot Simulation Results
plotReactorSteady(cstr)
plotDynamicVollbrecht(cstr)
```
A detailed example on parameter estimation is provided [here](https://github.com/Jo1931/DynaMeth.jl/blob/master/examples/Experiments.jl). Kortuz et al. ([August 2025](https://doi.org/10.1016/j.cej.2025.164505)).

### Bifurcation Analysis

```julia
#Setup Reactor
bifurcation = SetupBifurcationAnalysis(kinetic=Seidel())
cstr = Reactor(bifurcation)

#Temperature Range
Tvec = 250:700
Tc0 = 500   #cooling temperatur

#Caculate Heat Removal/Production Curve
hp, hr, sol = solveHeatCascade(Tvec, Tc0, cstr)

#Plot
plot_heat(hp, hr, Tvec, Tc0, cstr)
```
A detailed example on parameter estimation is provided [here](https://github.com/Jo1931/DynaMeth.jl/blob/master/examples/Bifurcation.jl).
Leipold et al. ([November 2024](https://doi.org/10.1002/ceat.202300256), [October 2025](https://doi.org/10.1016/j.cej.2025.167176))

### Forced periodic Operation

Leipold et al. ([July 2023](https://doi.org/10.1016/j.compchemeng.2023.108285), [June 2024](https://doi.org/10.1016/B978-0-443-28824-1.50268-4), [May 2025](https://doi.org/10.1016/j.cej.2025.167176))
Kaps et al. ([December 2025](https://doi.org/10.1002/cctc.202501403)) 

### Nonlinear Model Predictive Control

tbd

## Scope

This package contains the core models and methods for dynamic reactor simulation and optimization used in the associated research.
Scripts for reproducing specific figures and results will be provided in separate repositories corresponding to the individual journal publications.

## Status 

The package is under active development.

