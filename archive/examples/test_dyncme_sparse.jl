"""
test_dyncme_sparse.jl

Demonstrates DynCME with bounded K':
  - Sparse regime: μ_ss = k_b/k_d = 0.2 (most voxels have 0 or 1 molecule)
  - n_max=2: joint state space stays tiny even at K'=8
  - Equil deactivation triggers quickly (T_eq = 1/k_d = 1 time unit)
  - Front propagates ≈ 1 voxel per time unit → K' bounded to front width ≈ 4-6

This shows the key idea: after restricting to the active front,
the joint CME is tractable for K > 100 fine voxels.
"""

using DiscStochSim, Printf

function run_sparse_test()
    K_fine = 30
    # Slow diffusion: D=0.5, sparse: k_b=0.2, k_d=1.0 → μ_ss=0.2
    mesh  = VoxelMesh(K_fine, [(i,i+1) for i in 1:K_fine-1], fill(0.5, K_fine-1))
    birth = fill(0.2, K_fine)
    death = 1.0    # μ_ss = 0.2 per voxel

    center = 15
    state = DynCME(mesh, birth, death;
                   init_voxels=[center], n_max=3, ε_equil=0.10, ε_expand=0.02)
    place_molecules!(state, center, 1)

    println("K_fine=$K_fine  μ_ss=0.2  D=0.5")
    println("─"^70)

    dt = 0.5
    for step in 1:30
        t1 = time()
        state = step!(state, dt; prune_tol=1e-10)
        elapsed = time() - t1
        t = step * dt
        n_eq  = count(==(VSTATE_EQUIL), state.vstate)
        n_emp = count(==(VSTATE_EMPTY),  state.vstate)
        K_act = length(state.active_vids)
        @printf("t=%5.1f  K'=%2d  eq=%2d  emp=%2d  states=%5d  (%.2fs)\n",
                t, K_act, n_eq, n_emp, length(state.sp), elapsed)
        length(state.sp) > 50000 && (println("State limit hit."); break)
    end

    println("\nFine grid means:")
    μ = fine_mean_grid(state)
    for i in 1:K_fine
        lab = state.vstate[i]==1 ? "ACT" : state.vstate[i]==2 ? "EQL" : "EMP"
        μ[i] > 0.001 || state.vstate[i] != VSTATE_EMPTY || continue
        @printf("  v%2d: μ=%.4f [%s]\n", i, μ[i], lab)
    end
    @printf("Total: %.4f  (target: %.1f)\n", sum(μ), K_fine * 0.2)
end

run_sparse_test()
