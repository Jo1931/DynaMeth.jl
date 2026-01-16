# DynaMeth.jl

DynaMeth.jl is a Julia package developed as part of my PhD thesis on dynamic reactor operation for methanol synthesis.
It provides the core functions used to generate the results presented in my journal publications.

## Installation

Once registered, the package can be installed via:
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

tbd

### Bifurcation Analysis

tbd

### Forced periodic Operation

tbd

### Nonlinear Model Predictive Control

tbd

## Scope

This package contains the core models and methods for dynamic reactor simulation and optimization used in the associated research.
Scripts for reproducing specific figures and results are provided in separate repositories corresponding to the individual journal publications.

## Status 

The package is under active development.

