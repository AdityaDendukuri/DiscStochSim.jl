using DiscStochSim, Printf

function test_birth()
    K_fine = 3
    mesh  = VoxelMesh(K_fine, [(1,2),(2,3)], fill(2.0, 2))
    birth = fill(1.0, K_fine); death = 0.5   # μ_ss = 2.0 per voxel

    # All 3 voxels active from t=0, 1 molecule at center (voxel 2)
    state = DynCME(mesh, birth, death; init_voxels=[1,2,3], n_max=8, ε_equil=0.15, ε_expand=0.1)
    place_molecules!(state, 2, 1)

    dt = 0.5
    for step in 1:10
        state = step!(state, dt; prune_tol=1e-10)
        μ = fine_mean_grid(state)
        total = sum(μ)
        ode_total = 3*2*(1-exp(-0.5*step*dt)) + exp(-0.5*step*dt)
        @printf("t=%.1f  sim=%.3f  ODE=%.3f  voxels=%s  K'=%d  states=%d\n",
                step*dt, total, ode_total, string(round.(μ,digits=2)),
                length(state.active_vids), length(state.sp))
    end
end

test_birth()
