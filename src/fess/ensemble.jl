function _common_subnetwork!(
    workspace,
    entrances,
    entrance_counts,
    active_count,
    model,
    rates,
    options,
    persistent,
)
    augment!(
        workspace.network,
        entrances,
        model,
        rates;
        depth=options.expansion_depth,
        atol=options.atol,
    )
    subnetwork = FESSSubnetwork(workspace.network, model, rates)
    persistent || return subnetwork, _NO_SHEDDING

    distribution = zeros(length(workspace.network))
    for entrance in entrances
        distribution[workspace.network.index[entrance]] =
            entrance_counts[entrance] / active_count
    end
    _, shedding = shed!(subnetwork, distribution, entrances, options)
    subnetwork, shedding
end

function _common_occupations(subnetwork, entrances)
    right_hand_sides = zeros(length(subnetwork.network), length(entrances))
    columns = Dict{eltype(entrances),Int}()
    for (column, entrance) in enumerate(entrances)
        right_hand_sides[subnetwork.network.index[entrance], column] = 1.0
        columns[entrance] = column
    end
    first = max.(occupation(subnetwork, right_hand_sides), 0.0)
    second = second_occupation(subnetwork, first)
    columns, first, second
end

_state_values(state::CartesianIndex) = Float64[Tuple(state)...]
_state_values(state::Tuple) = Float64[state...]
_state_values(state::AbstractVector) = Float64[state...]

function _trajectory_moments(trajectory_times, trajectory_states, final_time, saveat)
    times = collect(0.0:saveat:final_time)
    (isempty(times) || times[end] < final_time) && push!(times, final_time)
    dimensions = length(_state_values(first(first(trajectory_states))))
    sums = [zeros(dimensions) for _ in times]
    sums_of_squares = [zeros(dimensions) for _ in times]

    for trajectory in eachindex(trajectory_times), time_index in eachindex(times)
        path_index = searchsortedlast(trajectory_times[trajectory], times[time_index])
        values = _state_values(trajectory_states[trajectory][path_index])
        sums[time_index] .+= values
        sums_of_squares[time_index] .+= values .^ 2
    end

    count = length(trajectory_times)
    means = [value ./ count for value in sums]
    variances = count == 1 ? [zeros(dimensions) for _ in times] : [
        max.((sums_of_squares[i] .- sums[i] .^ 2 ./ count) ./ (count - 1), 0.0)
        for i in eachindex(times)
    ]
    times, means, variances
end

"""
    fess_ensemble(initial_state, count, final_time, model, rates; ...)

Advance all unfinished trajectories once per step using one shared generator.
Each trajectory keeps its own state, clock, and RNG. The model is assumed to be
time homogeneous because trajectories at different physical times share `R`.
"""
function fess_ensemble(
    initial_state::S,
    count::Integer,
    final_time::Real,
    model,
    rates;
    persistent::Bool=true,
    options::FESSOptions=FESSOptions(),
    saveat::Real=0.1,
    seed::Integer=2026,
    number_of_steps::Union{Nothing,Integer}=nothing,
    max_steps::Int=100_000,
    max_network_size::Int=1_000,
) where {S}
    count > 0 || throw(ArgumentError("trajectory count must be positive"))
    final_time >= 0 || throw(ArgumentError("final time must be nonnegative"))
    saveat > 0 || throw(ArgumentError("save interval must be positive"))
    isnothing(number_of_steps) || number_of_steps >= 0 ||
        throw(ArgumentError("number of steps must be nonnegative"))
    _check(options)

    current_times = zeros(count)
    current_states = fill(initial_state, count)
    trajectory_times = [Float64[0.0] for _ in 1:count]
    trajectory_states = [S[initial_state] for _ in 1:count]
    rngs = [MersenneTwister(seed + i - 1) for i in 1:count]
    shared_workspace = FESSWorkspace(S)
    step_counts = zeros(Int, count)
    maximum_network_sizes = zeros(Int, count)
    network_sizes = Int[]
    shed_counts = Int[]
    active_counts = Int[]
    starting_state_counts = Int[]
    residuals = Float64[]

    wall_time = @elapsed begin
        step_limit = isnothing(number_of_steps) ? max_steps : number_of_steps
        for step_number in 1:step_limit
            candidates = isnothing(number_of_steps) ?
                findall(current_times .< final_time) : collect(1:count)
            isempty(candidates) && break
            active = Int[]
            for trajectory in candidates
                if _terminal(current_states[trajectory], model, rates, 0.0, options.atol)
                    if isnothing(number_of_steps)
                        current_times[trajectory] = final_time
                        push!(trajectory_times[trajectory], final_time)
                        push!(trajectory_states[trajectory], current_states[trajectory])
                    end
                else
                    push!(active, trajectory)
                end
            end
            isempty(active) && break

            counts = Dict{S,Int}()
            for trajectory in active
                state = current_states[trajectory]
                counts[state] = get(counts, state, 0) + 1
            end
            entrances = unique(current_states[active])
            workspace = persistent ? shared_workspace : FESSWorkspace(S)
            subnetwork, shedding = _common_subnetwork!(
                workspace, entrances, counts, length(active), model, rates,
                options, persistent,
            )
            length(subnetwork.network) <= max_network_size || error(
                "shared FESS network exceeded $max_network_size states",
            )
            columns, first, second = _common_occupations(subnetwork, entrances)
            network_size = length(subnetwork.network)
            push!(network_sizes, network_size)
            push!(shed_counts, shedding.removed)
            push!(active_counts, length(active))
            push!(starting_state_counts, length(entrances))
            push!(residuals, subnetwork.residual)

            for trajectory in active
                column = columns[current_states[trajectory]]
                next_state, _, waiting_time = _sample_exit(
                    subnetwork,
                    view(first, :, column),
                    view(second, :, column),
                    model,
                    rates,
                    rngs[trajectory];
                    atol=options.atol,
                )
                step_counts[trajectory] += 1
                maximum_network_sizes[trajectory] = max(
                    maximum_network_sizes[trajectory], network_size,
                )

                if isnothing(number_of_steps) &&
                   current_times[trajectory] + waiting_time >= final_time
                    current_times[trajectory] = final_time
                    push!(trajectory_times[trajectory], final_time)
                    push!(trajectory_states[trajectory], current_states[trajectory])
                else
                    current_times[trajectory] += waiting_time
                    current_states[trajectory] = next_state
                    push!(trajectory_times[trajectory], current_times[trajectory])
                    push!(trajectory_states[trajectory], next_state)
                end
            end
            isnothing(number_of_steps) && step_number == max_steps && error(
                "multiple-trajectory FESS exceeded $max_steps steps",
            )
        end
    end

    if isnothing(number_of_steps)
        times, means, variances = _trajectory_moments(
            trajectory_times, trajectory_states, final_time, saveat,
        )
    else
        times, means, variances = Float64[], Vector{Float64}[], Vector{Float64}[]
    end
    FESSEnsemble(
        times,
        means,
        variances,
        current_times,
        current_states,
        step_counts,
        maximum_network_sizes,
        network_sizes,
        shed_counts,
        active_counts,
        starting_state_counts,
        residuals,
        length(network_sizes),
        trajectory_times,
        trajectory_states,
        wall_time,
        persistent,
    )
end

function fess_summary(result::FESSEnsemble)
    shared_steps = result.shared_steps
    (
        trajectories=length(result.step_counts),
        total_trajectory_steps=sum(result.step_counts),
        average_steps_per_trajectory=
            sum(result.step_counts) / length(result.step_counts),
        maximum_network_size=maximum(result.maximum_network_sizes),
        shared_steps=shared_steps,
        trajectories_per_shared_step=shared_steps == 0 ? 0.0 :
            sum(result.step_counts) / shared_steps,
        total_states_removed=sum(result.step_removed_states),
        maximum_scaled_residual=isempty(result.residuals) ? 0.0 :
            maximum(result.residuals),
        wall_time=result.wall_time,
    )
end
