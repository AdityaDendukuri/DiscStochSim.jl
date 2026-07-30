using DiscStochSim
using Test
using SparseArrays
using Catalyst
using Random

@testset "DiscStochSim.jl" begin

    @testset "StateSpace basics" begin
        sp = StateSpace{CartesianIndex{2}, Float64}()
        @test length(sp) == 0

        add_state!(sp, CartesianIndex(1, 0), 0.5)
        add_state!(sp, CartesianIndex(0, 1), 0.3)
        @test length(sp) == 2
        @test sp.probs == [0.5, 0.3]

        # Duplicate add is no-op
        add_state!(sp, CartesianIndex(1, 0), 0.9)
        @test length(sp) == 2
        @test sp.probs[1] == 0.5

        # Remove
        remove_states!(sp, BitVector([true, false]))
        @test length(sp) == 1
        @test sp.states == [CartesianIndex(0, 1)]
    end

    @testset "Renormalize" begin
        sp = StateSpace{CartesianIndex{2}, Float64}()
        add_state!(sp, CartesianIndex(0, 0), 0.3)
        add_state!(sp, CartesianIndex(1, 0), 0.2)
        s = renormalize!(sp)
        @test s ≈ 0.5
        @test sum(sp.probs) ≈ 1.0
    end

    @testset "BoundaryCondition" begin
        @test RectLatticeBoundaryCondition(CartesianIndex(5, 5), (0, 10)) == true
        @test RectLatticeBoundaryCondition(CartesianIndex(-1, 5), (0, 10)) == false
        @test RectLatticeBoundaryCondition(CartesianIndex(5, 11), (0, 10)) == false
        @test RectLatticeBoundaryCondition(CartesianIndex(3, 3), 10) == true
        @test RectLatticeBoundaryCondition(CartesianIndex(-1, 3), 10) == false
    end

    @testset "DiscreteStochasticSystem from Catalyst" begin
        rn = @reaction_network begin
            k1, A --> B
            k2, B --> A
        end
        model = DiscreteStochasticSystem(rn)
        @test length(model.stoichvecs) == 2
        @test length(model.propensities) == 2
    end

    @testset "Generator column sums zero" begin
        rn = @reaction_network begin
            k1, A --> B
            k2, B --> A
        end
        model = DiscreteStochasticSystem(rn)
        bc = x -> RectLatticeBoundaryCondition(x, (0, 10))
        sp = StateSpace{CartesianIndex{2}, Float64}()
        for i in 0:3, j in 0:3
            add_state!(sp, CartesianIndex(i, j), (i == 1 && j == 0) ? 1.0 : 0.0)
        end
        rates = [1.0, 0.5]
        A, _, _ = build_generator(sp, model, rates, 0.0)
        col_sums = vec(sum(A, dims=1))
        @test all(x -> abs(x) < 1e-12, col_sums)
    end

    @testset "General state-first ELSE" begin
        states = [:a, :b]
        R = sparse([-2.0 1.0; 1.0 -3.0])
        exits = Dict(
            :a => [:x => 1.0],
            :b => [:y => 2.0],
        )
        subnetwork = ELSESubnetwork(states, R, exits; name=:tiny)

        law = else_law(subnetwork, :a)
        @test law.occupation ≈ [0.6, 0.2]
        @test law.exit_probability ≈ [0.6, 0.4]
        @test law.conditional_time ≈ [2 / 3, 1.0]

        step = else_step(
            subnetwork,
            :a;
            rng=MersenneTwister(1),
        )
        @test step.preexit in states
        @test step.next_state in (:x, :y)

        population = else_population_step(
            subnetwork,
            Dict(:a => 10_000);
            rng=MersenneTwister(2),
        )
        @test sum(values(population.groups)) == 10_000
        @test size(population.M_counts) == (2, 1)
        @test population.occupation ≈ 10_000 .* law.occupation
        @test population.mean_waiting_time ≈ sum(law.occupation)

        left = ELSESubnetwork(
            [:left],
            sparse(reshape([-2.0], 1, 1)),
            Dict(:left => [:right => 2.0]);
            name=:left,
        )
        right = ELSESubnetwork(
            [:right],
            sparse(reshape([-1.0], 1, 1)),
            Dict(:right => [:left => 1.0]);
            name=:right,
        )
        subnetwork_for(state) = state == :left ? left : right

        together = else_population(
            subnetwork_for,
            :left,
            1_000,
            4;
            rng=MersenneTwister(3),
        )
        @test together.final_groups == Dict(:left => 1_000)
        @test Set(keys(together.occupations)) == Set((:left, :right))
        @test length(together.history) == 4
    end

    @testset "Reconstruction stays aligned after pruning" begin
        rn = @reaction_network begin
            kbirth, 0 --> A
            kdeath, A --> 0
        end
        model = DiscreteStochasticSystem(rn)
        rates = [2.0, 0.5]

        sp = StateSpace{CartesianIndex{1}, Float64}()
        for n in 0:5
            add_state!(sp, CartesianIndex(n), n == 2 ? 1.0 : 0.0)
        end

        A_old, _, _ = build_generator(sp, model, rates, 0.0)
        gids_old = get_global_ids(sp)

        # Emulate one step where states are pruned/reindexed, then expanded again.
        remove_states!(sp, BitVector([false, true, false, false, true, false]))
        add_state!(sp, CartesianIndex(6), 0.0)

        A_recon, _, _ = reconstruct_generator(sp, model, rates, 0.0, A_old, gids_old)
        A_full, _, _ = build_generator(sp, model, rates, 0.0)

        @test Matrix(A_recon) ≈ Matrix(A_full) atol=1e-12
    end

    @testset "Expand adds neighbors" begin
        rn = @reaction_network begin
            k1, A --> B
            k2, B --> A
        end
        model = DiscreteStochasticSystem(rn)
        bc = x -> RectLatticeBoundaryCondition(x, (0, 100))
        sp = StateSpace{CartesianIndex{2}, Float64}()
        add_state!(sp, CartesianIndex(5, 5), 1.0)
        expand!(sp, model, bc; depth=1)
        @test length(sp) >= 3  # original + at least 2 neighbors
    end

    @testset "FSPProblem + solve smoke test" begin
        rn = @reaction_network begin
            k1, A --> B
            k2, B --> A
        end
        rates = [1.0, 0.5]
        prob = FSPProblem(rn, CartesianIndex(5, 0), (0.0, 1.0), rates; bounds=(0, 20))
        sol, _ = solve(prob, AdaptiveFSP(ε_dt=0.5, prob_quantile=0.05,
                                        flux_tolerance=1e-6, save_interval=100))
        @test length(sol) >= 2
        traj = mean_trajectory(sol)
        @test size(traj, 2) == 2
        # Conservation: A + B should stay ≈ 5
        @test all(abs.(traj[:, 1] .+ traj[:, 2] .- 5.0) .< 0.5)
    end

    @testset "Marginal distribution" begin
        rn = @reaction_network begin
            k1, A --> B
            k2, B --> A
        end
        rates = [1.0, 0.5]
        prob = FSPProblem(rn, CartesianIndex(3, 0), (0.0, 0.5), rates; bounds=(0, 10))
        sol, _ = solve(prob, AdaptiveFSP(ε_dt=0.5, save_interval=50))
        m = marginal(sol, length(sol), 1)
        @test sum(m.probs) ≈ 1.0 atol=0.01
    end
end
