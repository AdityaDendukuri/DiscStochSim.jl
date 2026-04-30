module DiscStochSim

using Catalyst
using SparseArrays
using LinearAlgebra
using Statistics
using ExponentialUtilities
using Expokit
using UnicodePlots
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
include("adaptive_fsp/krylov_fsp.jl")
include("adaptive_fsp/terminal_progress.jl")

# RDME multigrid FSP
include("rdme/voxel_grid.jl")
include("rdme/rdme_model.jl")
include("rdme/operators.jl")
include("rdme/dynamic_pi.jl")
include("rdme/vcycle.jl")
include("rdme/solve.jl")

# Exports
export DiscreteStochasticSystem
export RectLatticeBoundaryCondition
export StateSpace, add_state!, remove_states!, expand!, expand_ssa!, expand_flux!, compress!, renormalize!, get_global_ids
export build_generator, reconstruct_generator
export FSPProblem, FSPSolution, AdaptiveFSP, AdaptiveFSPDiagnostics, KrylovFSP, KrylovFSPDiagnostics, prune_threshold!
using CommonSolve: solve
export solve
export mean_trajectory, marginal, terminal_progress

# RDME exports
export VoxelGrid, coarsen, build_hierarchy, diffusion_rate
export RDMEModel1D, build_rdme_system, build_intra_system, build_coarse_system, rdme_bc
export SchloglModel1D, build_schlogl_rdme_system, build_schlogl_coarse_system, build_schlogl_mixed_system
export build_schlogl_coarse_system_dynamic
export schlogl_rate_eq, schlogl_fixed_points
export compute_dynamic_pi, prolong_dynamic
export restrict, prolong, prolong_multiplicative, prolong_conditional, lumpability_ratio, per_molecule_flux, pair_admissibility, select_coarsening_mask, partial_restrict, partial_prolong
export select_coarsening_mask_mixed, adapt_mixed_state, mixed_coord_bounds
export _mixed_dim
export two_level_vcycle, two_level_vcycle_schlogl, two_level_vcycle_schlogl_injection
export two_level_vcycle_adaptive, two_level_step_mixed, MixedSolverCache
export RDMEMultigridFSP, RDMESolution, solve_rdme_multigrid
export mean_voxel_counts, marginal_voxel

end
