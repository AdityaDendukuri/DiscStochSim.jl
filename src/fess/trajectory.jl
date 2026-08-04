function _trajectory_step!(
    workspace,
    persistent,
    state,
    time,
    model,
    rates,
    options,
    rng,
)
    if persistent
        return fess_step!(
            workspace, state, model, rates;
            options=options, rng=rng, time=time, shed=true,
        )
    end
    fess_step(state, model, rates; options=options, rng=rng, time=time)
end

"""
    fess_trajectory(initial_state, final_time, model, rates; persistent, ...)

Simulate one FESS trajectory. Persistent trajectories reuse one `FESSWorkspace`;
redrawn trajectories construct a one-shell network at every step.
"""
function fess_trajectory(
    initial_state::S,
    final_time::Real,
    model,
    rates;
    persistent::Bool=true,
    options::FESSOptions=FESSOptions(),
    rng=Random.GLOBAL_RNG,
    max_steps::Int=100_000,
) where {S}
    final_time >= 0 || throw(ArgumentError("final time must be nonnegative"))
    workspace = FESSWorkspace(S)
    times = Float64[0.0]
    states = S[initial_state]
    step_times = Float64[]
    steps = FESSStep{S}[]
    state = initial_state
    time = 0.0

    for step_number in 1:max_steps
        time >= final_time && break
        push!(step_times, time)
        step = _trajectory_step!(
            workspace, persistent, state, time, model, rates, options, rng,
        )
        push!(steps, step)

        if step.terminated || time + step.waiting_time >= final_time
            times[end] < final_time && push!(times, Float64(final_time))
            length(states) < length(times) && push!(states, state)
            break
        end

        time += step.waiting_time
        state = step.next_state
        push!(times, time)
        push!(states, state)
        step_number == max_steps && error(
            "FESS trajectory exceeded $max_steps steps at time $time",
        )
    end
    FESSTrajectory(times, states, step_times, steps)
end

"""
    fess_fixed_steps(initial_state, number_of_steps, model, rates; persistent, ...)

Run an exact number of FESS macrosteps and report the physical time reached.
"""
function fess_fixed_steps(
    initial_state::S,
    number_of_steps::Integer,
    model,
    rates;
    persistent::Bool=true,
    options::FESSOptions=FESSOptions(),
    rng=Random.GLOBAL_RNG,
) where {S}
    number_of_steps >= 0 ||
        throw(ArgumentError("number of steps must be nonnegative"))
    workspace = FESSWorkspace(S)
    state = initial_state
    time = 0.0
    steps = FESSStep{S}[]

    for _ in 1:number_of_steps
        step = _trajectory_step!(
            workspace, persistent, state, time, model, rates, options, rng,
        )
        step.terminated && break
        push!(steps, step)
        state = step.next_state
        time += step.waiting_time
    end
    (final_time=time, final_state=state, steps=steps,
     completed_steps=length(steps))
end
