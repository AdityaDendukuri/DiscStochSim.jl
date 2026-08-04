@testset "FESS core" begin
    @testset "C++ reference formulas" begin
        network = FESSNetwork(:a)
        augment!(network, :b)
        subnetwork = FESSSubnetwork(network, [-2.0 1.0; 1.0 -3.0])
        entrance = [1.0, 0.0]

        @test subnetwork.Z ≈ [0.6 0.2; 0.2 0.4]
        @test escape_time(subnetwork, entrance) ≈ 0.8
        @test escape_probabilities(subnetwork, entrance) ≈ [0.6, 0.4]
        @test cut_time_losses(subnetwork, entrance) ≈ [0.8, 0.3]
        first = occupation(subnetwork, entrance)
        second = second_occupation(subnetwork, first)
        @test conditional_time(subnetwork, first, second, 1) ≈ 2 / 3
        @test conditional_time(subnetwork, first, second, 2) ≈ 1.0

        remaining, stats = shed!(
            subnetwork,
            entrance,
            [:a],
            FESSOptions(capacity=3, cut_threshold=0.4),
        )
        @test subnetwork.network.states == [:a]
        @test remaining == [1.0]
        @test stats.threshold_cuts == 1
        @test stats.actual_relative_loss ≈ 0.375
    end

    @testset "hard capacity" begin
        network = FESSNetwork(:a)
        augment!(network, :b)
        subnetwork = FESSSubnetwork(network, [-2.0 1.0; 1.0 -3.0])
        _, stats = shed!(
            subnetwork,
            [1.0, 0.0],
            [:a],
            FESSOptions(capacity=2, cut_threshold=0.0),
        )
        @test length(subnetwork.network) == 1
        @test stats.capacity_cuts == 1
    end

    @testset "single and multiple trajectories" begin
        rn = @reaction_network begin
            kbirth, 0 --> X
            kdeath, X --> 0
        end
        model = DiscreteStochasticSystem(rn)
        rates = [2.0, 0.5]
        initial = CartesianIndex(4)
        options = FESSOptions(
            capacity=12,
            cut_threshold=1.0e-3,
            expansion_depth=1,
        )

        redrawn = fess_step(
            initial,
            model,
            rates;
            options=options,
            rng=MersenneTwister(10),
        )
        @test !redrawn.terminated
        @test redrawn.network_size == 3
        @test redrawn.waiting_time > 0

        persistent = fess_fixed_steps(
            initial,
            4,
            model,
            rates;
            persistent=true,
            options=options,
            rng=MersenneTwister(11),
        )
        @test persistent.completed_steps == 4
        @test persistent.final_time > 0

        multiple = fess_ensemble(
            initial,
            3,
            0.0,
            model,
            rates;
            persistent=true,
            options=options,
            number_of_steps=4,
            seed=12,
        )
        @test multiple.shared_steps == 4
        @test multiple.step_counts == fill(4, 3)
        @test all(>(0), multiple.final_times)

        comparison = compare_fess_steps(
            initial,
            3,
            model,
            rates;
            n_trajectories=2,
            options=options,
            seed=13,
        )
        @test length(comparison.rows) == 4
        @test all(row -> row.steps == 3, comparison.rows)
        @test all(row -> row.average_final_time > 0, comparison.rows)
    end
end
