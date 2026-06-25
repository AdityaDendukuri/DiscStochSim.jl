# DiscStochSim.jl

[![Build Status](https://github.com/AdityaDendukuri/DiscStochSim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AdityaDendukuri/DiscStochSim.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/AdityaDendukuri/DiscStochSim.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/AdityaDendukuri/DiscStochSim.jl)

Flux-preserving adaptive Finite State Projection (FSP) for stochastic chemical reaction networks. Implementation for:

> A. Dendukuri, S. Chandrasekaran, and L. Petzold, *Flux-Preserving Adaptive Finite State Projection for Multiscale Stochastic Reaction Networks*, 2026.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/AdityaDendukuri/DiscStochSim.jl")
```

## Usage

```julia
using DiscStochSim, Catalyst

rn   = @reaction_network begin ... end
prob = FSPProblem(rn, u0, tspan, params; alpha=0.05, eps_flux=1e-6, eps_dt=1.0)
sol  = solve(prob, AdaptiveFSP())
```

Scripts reproducing all paper figures are in `examples/`.
