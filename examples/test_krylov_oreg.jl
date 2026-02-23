using DiscStochSim, Catalyst, Random

rn = @reaction_network begin
    @species X(t) Y(t) Z(t)
    @parameters k1 k2 k3 k4 k5
    k1, Y --> X
    k2, X + Y --> 0
    k3, X --> 2X + Z
    k4, 2X --> 0
    k5, Z --> Y
end

y1s, y2s, y3s = 500.0, 1000.0, 2000.0
mu1s, mu2s = 2000.0, 50000.0
rates = [mu1s/y2s, mu2s/(y1s*y2s), (mu1s+mu2s)/y1s, 2*mu1s/y1s^2, (mu1s+mu2s)/y3s]
prob = FSPProblem(rn, CartesianIndex(500, 1000, 2000), (0.0, 0.11), rates; bounds=(0, 50_000))

Random.default_rng(1)
t0 = time()
sol, diag = solve(prob, KrylovFSP(ε=1e-2, τ_init=1e-6, ℓ=5, r=1, save_interval=5000, max_iter=30000))
elapsed = time() - t0

println("Iters: $(diag.total_iters), rejections: $(diag.reject_count)")
println("Time reached: $(sol.t[end]), wall: $(round(elapsed; digits=1))s")
println("Accepted: $(diag.total_iters - diag.reject_count)")
if length(diag.size_log) > 0
    println("|S| range: [$(minimum(diag.size_log)), $(maximum(diag.size_log))]")
end
