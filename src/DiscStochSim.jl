module DiscStochSim

using Catalyst
using SparseArrays
using LinearAlgebra
using ExponentialUtilities
using Expokit
using UnicodePlots
using Random
using Distributions
import CommonSolve

# Core types
include("reaction_system.jl")
include("boundary_conditions.jl")
include("state_space.jl")
include("generator.jl")
include("else.jl")

# Adaptive FSP solver
include("adaptive_fsp/problem.jl")
include("adaptive_fsp/solution.jl")
include("adaptive_fsp/solve.jl")
include("adaptive_fsp/krylov_fsp.jl")
include("adaptive_fsp/terminal_progress.jl")

# Exports
export DiscreteStochasticSystem
export RectLatticeBoundaryCondition
export StateSpace, add_state!, remove_states!, expand!, expand_ssa!, expand_flux!, compress!, renormalize!, get_global_ids
export build_generator, reconstruct_generator
export ELSESubnetwork, else_law, else_step, else_trajectory
export else_population_step, else_population
export FSPProblem, FSPSolution, AdaptiveFSP, AdaptiveFSPDiagnostics, KrylovFSP, KrylovFSPDiagnostics, prune_threshold!
using CommonSolve: solve
export solve
export mean_trajectory, marginal, terminal_progress

end
