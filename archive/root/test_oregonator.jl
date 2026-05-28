using DiscStochSim, Catalyst
import Random
Random.seed!(1)

function main()
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

    t_start = time()
    sol = solve(prob, AdaptiveFSP(ε_dt=1.0, prob_quantile=0.1, flux_tolerance=0.0, expansion_depth=1, save_interval=500))
    wall = time() - t_start
    @info "DONE" t=round(sol.t[end]; sigdigits=4) wall=round(wall; sigdigits=3) max_S=maximum(sol.state_space_sizes) snaps=length(sol.t)
end

main()
