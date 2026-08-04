"""
    compare_fess_steps(initial_state, number_of_steps, model, rates; ...)

Run persistent and redrawn FESS for single and multiple trajectories. Each
method takes `number_of_steps`; the result reports the final physical time.
"""
function compare_fess_steps(
    initial_state,
    number_of_steps::Integer,
    model,
    rates;
    n_trajectories::Int=5,
    options::FESSOptions=FESSOptions(),
    seed::Integer=3100,
)
    persistent = fess_fixed_steps(
        initial_state, number_of_steps, model, rates;
        persistent=true, options=options, rng=MersenneTwister(seed),
    )
    redrawn = fess_fixed_steps(
        initial_state, number_of_steps, model, rates;
        persistent=false, options=options, rng=MersenneTwister(seed + 1),
    )
    multiple_persistent = fess_ensemble(
        initial_state, n_trajectories, 0.0, model, rates;
        persistent=true, options=options, seed=seed + 100,
        number_of_steps=number_of_steps,
    )
    multiple_redrawn = fess_ensemble(
        initial_state, n_trajectories, 0.0, model, rates;
        persistent=false, options=options, seed=seed + 200,
        number_of_steps=number_of_steps,
    )

    row(name, steps, times) = (
        method=name,
        steps=steps,
        trajectories=length(times),
        average_final_time=sum(times) / length(times),
        earliest_final_time=minimum(times),
        latest_final_time=maximum(times),
    )
    rows = [
        row("Persistent, single trajectory", persistent.completed_steps,
            [persistent.final_time]),
        row("Redrawn, single trajectory", redrawn.completed_steps,
            [redrawn.final_time]),
        row("Persistent, multiple trajectories",
            multiple_persistent.shared_steps, multiple_persistent.final_times),
        row("Redrawn, multiple trajectories",
            multiple_redrawn.shared_steps, multiple_redrawn.final_times),
    ]
    (
        rows=rows,
        single_persistent=persistent,
        single_redrawn=redrawn,
        multiple_persistent=multiple_persistent,
        multiple_redrawn=multiple_redrawn,
    )
end
