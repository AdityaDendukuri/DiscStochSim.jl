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

# Spatially adaptive FSP (new approach)
include("spatial_adaptive/bd_spatial_fsp.jl")
include("spatial_adaptive/flux_vcycle.jl")
# Schlögl as AbstractSpatialRDMEModel (must come after rdme_model.jl defines SchloglModel1D)

# RDME multigrid FSP
include("rdme/voxel_grid.jl")
include("rdme/catalyst_rdme.jl")
include("rdme/rdme_experiment.jl")
include("rdme/rdme_model.jl")
include("rdme/operators.jl")
include("rdme/dynamic_pi.jl")
include("rdme/patch_qsd.jl")
include("rdme/vcycle.jl")
include("rdme/mixed_vcycle.jl")
include("rdme/sp_transition.jl")
include("rdme/schlogl_mixed_1d.jl")
include("rdme/solve.jl")

# Graph-generalised RDME (depends on rdme_model, operators, solve)
include("rdme/graph_rdme.jl")

# Schlögl spatial FSP (must come after rdme_model.jl defines SchloglModel1D)
include("spatial_adaptive/schlogl_spatial_fsp.jl")

# Spatial CME: joint distribution over K voxels, multigrid coarsening
include("rdme/spatial_cme.jl")

# Tensor Train CME: MPS/TT representation of the joint RDME distribution
include("rdme/tensor_train_cme.jl")

# Row-TT solver: exact within-row joint + mean-field between rows
include("rdme/row_tt_solver.jl")

# Dynamic-K spatial CME: adaptive voxel activation/deactivation
include("rdme/dynamic_spatial_cme.jl")

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
export CoarseningMap, CoarseningMap1D, CoarseningMap2D, CoarseningMapFull, CoarseningMapRegions
export build_adaptive_cmap_1d, build_adaptive_cmap_2d
export fine_marginals_from_coarse
export RDMEModel1D, build_rdme_system, build_intra_system, build_coarse_system, rdme_bc
export SchloglModel1D, build_schlogl_rdme_system, build_schlogl_adaptive_system
export build_rdme_adaptive_2d
export remap_sp_transition
export build_schlogl_coarse_system, build_schlogl_mixed_system
export build_schlogl_coarse_system_dynamic
export schlogl_rate_eq, schlogl_fixed_points
export compute_dynamic_pi, prolong_dynamic
export restrict, prolong, prolong_product, prolong_multiplicative, prolong_conditional
export PatchPropensity, QSDProlongTable, build_qsd_table, prolong_qsd
export restrict2s, prolong2s, prolong_multiplicative_2s
export lumpability_ratio, pair_admissibility
export multi_level_vcycle, multi_level_vcycle_2s, multi_level_vcycle_2d_2s
export ActiveVoxelSet1D, n_active, left_inactive, right_inactive
export build_active_bottleneck_1d, mixed_vcycle_1d_2s
export update_inactive_means_1d!, adapt_active_1d, active_bc
export promote_left_1d, promote_right_1d, demote_left_1d, demote_right_1d
export build_active_schlogl_1d, update_inactive_means_schlogl_1d!
export promote_left_schlogl, promote_right_schlogl
export demote_left_schlogl, demote_right_schlogl
export adapt_schlogl_1d, schlogl_mixed_bc
export build_intra_schlogl_system, multi_level_vcycle_schlogl, evolve_active_schlogl!
export AnnihilationModel1D
export build_active_annihilation_1d, update_inactive_means_annihilation_1d!
export adapt_annihilation_1d
export PatchQSD, build_patch_qsd, prolong_patch_qsd
export patch_edges_1d, patch_edges_2d, patch_edges_from_cmap
export patch_qsd_correlation
# Graph RDME exports
export VoxelGraph, chain, grid_2d, grid_3d, degrees, adjacency_list
export build_schlogl_graph_system, solve_rdme_graph
export GraphRDMESolution, config_probs

# Spatial CME exports (joint distribution + multigrid on voxel axis)
export VoxelMesh, rect_voxel_mesh
export GraphRDMEModel, build_graph_rdme_system
export CoarseningPlan, make_coarsening_plan, restrict_state
export coarsen_mesh   # coarsen_model already exported above
export SpatialCMEState, set_voxel_ic!, spatial_vcycle!, mean_counts

# Row-TT Solver exports
export RowTTSolver, tt_step_row!, step!, mean_grid, max_bond_grid, row_joint_marginal, tt_norm_all

# Tensor Train CME exports
export CMETensorTrain, tt_norm, tt_marginals, tt_means, tt_compress!
export tt_evolve_local!, tt_evolve_all_local!, tt_evolve_pair!, tt_step_1d!, tt_step_2d!
export tt_bond_dims, tt_max_bond
export birth_death_generator, schlogl_generator

# Dynamic-K exports
export AbstractDynCME, DynCME
export VSTATE_ACTIVE, VSTATE_EQUIL, VSTATE_EMPTY
export marginal_k, marginal_mean_k, tv_vs_poisson
export marginalize_out, activate_voxel, deactivate_voxels
export evolve!, step!, activation_candidates
export fine_mean_grid, fine_status_grid, diagnose, place_molecules!

# Spatial adaptive FSP exports
export BirthDeathRDME, BranchingDeathRDME, SpatialRDME, SpatialFSP, set_ic!, step!, n_active
export RDMEExperiment, build_state, rect_grid
export ZeroIC, CenterIC, VoxelIC, RandomIC
export FluxGraph, build_flux_graph, update_flux!, evolve_flux_vcycle!
export n_equilibrated, n_empty, status_grid, mean_grid
export total_distribution, equilibration_fraction
export poisson_pmf_vec, jump_rate, ss_mean, voxel_mean

export BottleneckModel1D, build_bottleneck_system, build_intra_system_2s, bottleneck_bc
export build_bottleneck_system_2d, build_intra_system_2d
export RDMEMultigridFSP, RDMESolution, solve_rdme_multigrid
export mean_voxel_counts, marginal_voxel

end
