using DiscStochSim, Printf

function run_chain_test()
    K_fine = 12
    mesh  = VoxelMesh(K_fine, [(i,i+1) for i in 1:K_fine-1], fill(2.0, K_fine-1))
    birth = fill(1.0, K_fine); death = 0.5   # μ_ss = 2.0 per voxel
    center = 6

    state = DynCME(mesh, birth, death;
                   init_voxels=[center], n_max=6, ε_equil=0.12, ε_expand=0.15)
    place_molecules!(state, center, 1)

    dt = 0.5
    for step in 1:24
        t1 = time()
        state = step!(state, dt; prune_tol=1e-8)
        elapsed = time() - t1
        t = step * dt
        n_eq = count(==(VSTATE_EQUIL), state.vstate)
        @printf("t=%5.1f  K'=%2d  equil=%2d  states=%7d  (%.2fs)\n",
                t, length(state.active_vids), n_eq, length(state.sp), elapsed)
        length(state.sp) > 300000 && (println("State limit hit!"); break)
    end

    println("\nFine grid means:")
    μ = fine_mean_grid(state)
    for i in 1:K_fine
        lab = state.vstate[i]==1 ? "ACT" : state.vstate[i]==2 ? "EQL" : "EMP"
        @printf("  v%2d: μ=%.3f [%s]\n", i, μ[i], lab)
    end
    @printf("Total: %.3f  (target: %.1f)\n", sum(μ), K_fine * 2.0)
end

run_chain_test()
