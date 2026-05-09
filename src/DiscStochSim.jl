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
export VoxelGrid, coarsen, coarsen_model, build_hierarchy, diffusion_rate
export Grid2D, coarsen2d, build_hierarchy2d, n_voxels
export CoarseningMap, CoarseningMap1D, CoarseningMap2D
export RDMEModel1D, build_rdme_system, build_intra_system, build_coarse_system, rdme_bc
export SchloglModel1D, build_schlogl_rdme_system, build_schlogl_coarse_system, build_schlogl_mixed_system
export build_schlogl_coarse_system_dynamic
export schlogl_rate_eq, schlogl_fixed_points
export compute_dynamic_pi, prolong_dynamic
export restrict, prolong, prolong_multiplicative, prolong_conditional
export restrict2s, prolong2s, prolong_multiplicative_2s
export lumpability_ratio, pair_admissibility
export multi_level_vcycle, multi_level_vcycle_2s, multi_level_vcycle_2d_2s
export BottleneckModel1D, build_bottleneck_system, build_intra_system_2s, bottleneck_bc
export build_bottleneck_system_2d, build_intra_system_2d
export RDMEMultigridFSP, RDMESolution, solve_rdme_multigrid
export mean_voxel_counts, marginal_voxel

end
