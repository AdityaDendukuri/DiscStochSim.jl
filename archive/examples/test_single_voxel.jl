using DiscStochSim

function test_single_voxel()
    mesh  = VoxelMesh(1, Tuple{Int,Int}[], Float64[])
    model = GraphRDMEModel([1.0], [0.5])
    sys   = build_graph_rdme_system(model, mesh)
    println("Reactions: $(length(sys.stoichvecs))")
    println("Stoichs: $(sys.stoichvecs)")
    println("Prop at n=0 (birth): $(sys.propensities[1](CartesianIndex(0), Float64[], 0.0))")
    println("Prop at n=1 (death): $(sys.propensities[2](CartesianIndex(1), Float64[], 0.0))")

    sp = StateSpace{CartesianIndex{1}, Float64}()
    add_state!(sp, CartesianIndex(1), 1.0)
    bc = x -> Tuple(x)[1] <= 12
    dt = 0.5
    for step in 1:10
        expand!(sp, sys, bc; depth=1)
        A, = build_generator(sp, sys, Float64[], 0.0)
        sp.probs .= expv(dt, A, sp.probs; m=20)
        sp.probs .= max.(0.0, sp.probs); sp.probs ./= sum(sp.probs)
        prune_threshold!(sp, 1e-10)
        mu = sum(Tuple(ci)[1] * p for (ci,p) in zip(sp.states, sp.probs))
        println("t=$(step*dt)  mean=$(round(mu,digits=4))  states=$(length(sp))")
    end
    println("Analytical at t=5: $(round(2-exp(-2.5), digits=4))")
end

test_single_voxel()
