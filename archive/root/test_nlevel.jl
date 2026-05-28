using DiscStochSim
using LinearAlgebra

# Test simple birth-death on K=8 with 2 levels of coarsening (8 -> 4 -> 2)
K = 8
fine_grid = VoxelGrid(K, 1.0, 0)
u0 = CartesianIndex(ntuple(j -> j == 1 ? 10 : 0, Val(K))) # 10 molecules in first voxel

# High diffusion to make coarse models relevant
model = RDMEModel1D(10.0, 1.0, 0.1) 

alg = RDMEMultigridFSP(
    dt = 0.1,
    n_levels = 2,
    n_max = 30,
    save_every = 1,
    coarse_only = true,
    max_states = 100_000 
)

tspan = (0.0, 0.5)
rates = Float64[]

println("Starting 3-level RDME Multigrid test (K=8 -> 4 -> 2)...")
sol = solve_rdme_multigrid(model, fine_grid, u0, tspan, rates, alg)

println("Test completed successfully.")
println("Snapshots: ", length(sol.t))
println("State space sizes: ", sol.state_space_sizes)
