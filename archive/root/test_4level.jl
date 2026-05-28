using DiscStochSim
using LinearAlgebra

# Test simple birth-death on K=16 with 3 levels of coarsening (16 -> 8 -> 4 -> 2)
K = 16
fine_grid = VoxelGrid(K, 1.0, 0)
u0 = CartesianIndex(ntuple(j -> j == 1 ? 10 : 0, Val(K))) 

model = RDMEModel1D(10.0, 1.0, 0.1) 

alg = RDMEMultigridFSP(
    dt = 0.1,
    n_levels = 3,
    n_max = 20,
    save_every = 1,
    coarse_only = true,
    max_states = 10_000_000,
    weight_tol = 1e-12,
    binom_tol = 1e-4
)

tspan = (0.0, 0.3)
rates = Float64[]

println("Starting 4-level RDME Multigrid test (K=16 -> 8 -> 4 -> 2)...")
sol = solve_rdme_multigrid(model, fine_grid, u0, tspan, rates, alg)

println("Test completed successfully.")
println("Snapshots: ", length(sol.t))
println("State space sizes: ", sol.state_space_sizes)
