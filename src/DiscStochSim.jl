module DiscStochSim

using Catalyst
using SparseArrays
using LinearAlgebra
using Expokit
import CommonSolve

# Core types
include("reaction_system.jl")
include("boundary_conditions.jl")
include("state_space.jl")
include("generator.jl")

# Adaptive FSP solver
include("adaptive_fsp/problem.jl")
include("adaptive_fsp/solution.jl")
include("adaptive_fsp/solve.jl")

# Exports
export DiscreteStochasticSystem
export RectLatticeBoundaryCondition
export StateSpace, add_state!, remove_states!, expand!, compress!, renormalize!, get_global_ids
export build_generator, reconstruct_generator
export FSPProblem, FSPSolution, AdaptiveFSP
using CommonSolve: solve
export solve
export mean_trajectory, marginal

end
