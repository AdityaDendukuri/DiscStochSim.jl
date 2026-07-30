"""
    ELSESubnetwork(states, R, exits; name=:subnetwork)

A finite ELSE subnetwork.

`R` is the killed generator on `states`. `exits[state]` is a collection of
`next_state => rate` pairs containing every transition that leaves the
subnetwork from `state`.
"""
struct ELSESubnetwork{S,M,F}
    name::Symbol
    states::Vector{S}
    index::Dict{S,Int}
    R::M
    boundary::Vector{Int}
    next_states::Vector{Vector{S}}
    exit_rates::Vector{Vector{Float64}}
    total_exit_rate::Vector{Float64}
    Z_solve::F
end

function ELSESubnetwork(
    states::AbstractVector{S},
    R::AbstractMatrix,
    exits::AbstractDict;
    name=:subnetwork,
) where {S}
    state_list = collect(states)
    state_index =
        Dict{S,Int}(state => i for (i, state) in enumerate(state_list))
    n = length(state_list)

    size(R) == (n, n) ||
        throw(ArgumentError("R must have one row and column per internal state"))

    boundary = Int[]
    next_states = Vector{Vector{S}}()
    exit_rates = Vector{Vector{Float64}}()
    total_exit_rate = zeros(n)

    for state in keys(exits)
        haskey(state_index, state) ||
            throw(ArgumentError("exit source $state is not an internal state"))
    end

    for (i, state) in enumerate(state_list)
        channels = get(exits, state, ())
        isempty(channels) && continue

        destinations = S[]
        rates = Float64[]

        for channel in channels
            destination = first(channel)
            rate = Float64(last(channel))
            rate >= 0.0 || throw(ArgumentError("exit rates must be nonnegative"))
            rate == 0.0 && continue
            push!(destinations, destination)
            push!(rates, rate)
        end

        isempty(rates) && continue
        push!(boundary, i)
        push!(next_states, destinations)
        push!(exit_rates, rates)
        total_exit_rate[i] = sum(rates)
    end

    column_loss = -vec(sum(R; dims=1))
    scale = max.(1.0, abs.(column_loss), abs.(total_exit_rate))
    balanced = abs.(column_loss .- total_exit_rate) .<= 1.0e-9 .* scale

    all(balanced) || throw(
        ArgumentError(
            "the supplied exit channels do not account for every column loss in R",
        ),
    )

    ELSESubnetwork(
        Symbol(name),
        state_list,
        state_index,
        R,
        boundary,
        next_states,
        exit_rates,
        total_exit_rate,
        lu(-R),
    )
end

"""
    else_law(subnetwork, start_state)

Compute the state-first ELSE law using two solves with the stored factorization:
`z = Z*e_a` and `z2 = Z*z`, where `Z = -inv(R)`.
"""
function else_law(subnetwork::ELSESubnetwork, start_state)
    start = get(subnetwork.index, start_state, 0)
    start > 0 || throw(ArgumentError("start state is outside the subnetwork"))

    p0 = zeros(length(subnetwork.states))
    p0[start] = 1.0

    occupation = subnetwork.Z_solve \ p0
    second_occupation = subnetwork.Z_solve \ occupation

    weights =
        subnetwork.total_exit_rate[subnetwork.boundary] .*
        occupation[subnetwork.boundary]
    weights = max.(weights, 0.0)
    mass = sum(weights)
    mass > 0.0 || throw(ArgumentError("the subnetwork has no reachable exit"))

    conditional_time = [
        occupation[index] > 0.0 ?
        second_occupation[index] / occupation[index] :
        0.0 for index in subnetwork.boundary
    ]

    (
        occupation=occupation,
        exit_probability=weights ./ mass,
        conditional_time=conditional_time,
    )
end

function _else_sample(probability, rng)
    target = rand(rng)
    total = 0.0

    for i in eachindex(probability)
        total += probability[i]
        total >= target && return i
    end

    lastindex(probability)
end

"""
    else_step(subnetwork, start_state; rng=Random.GLOBAL_RNG)

Sample the pre-exit state first, use its conditional mean exit time, and then
sample the physical transition leaving the subnetwork.
"""
function else_step(
    subnetwork::ELSESubnetwork,
    start_state;
    rng=Random.GLOBAL_RNG,
)
    law = else_law(subnetwork, start_state)
    boundary_position = _else_sample(law.exit_probability, rng)
    preexit_index = subnetwork.boundary[boundary_position]

    rates = subnetwork.exit_rates[boundary_position]
    channel_probability = rates ./ sum(rates)
    channel = _else_sample(channel_probability, rng)

    (
        subnetwork=subnetwork.name,
        waiting_time=law.conditional_time[boundary_position],
        preexit=subnetwork.states[preexit_index],
        next_state=subnetwork.next_states[boundary_position][channel],
        occupation=law.occupation,
    )
end

"""
    else_trajectory(subnetwork_for, initial_state, number_of_steps; rng)

Generate a state-first ELSE exit skeleton. `subnetwork_for(state)` selects the
subnetwork containing the current state.
"""
function else_trajectory(
    subnetwork_for,
    initial_state,
    number_of_steps;
    rng=Random.GLOBAL_RNG,
)
    state = initial_state
    time = 0.0
    times = Float64[time]
    states = [state]
    subnetworks = Symbol[]

    for _ in 1:number_of_steps
        subnetwork = subnetwork_for(state)
        step = else_step(subnetwork, state; rng=rng)

        time += step.waiting_time
        push!(times, time)
        push!(states, step.preexit)
        push!(times, time)
        push!(states, step.next_state)
        push!(subnetworks, step.subnetwork)

        state = step.next_state
    end

    (times=times, states=states, subnetworks=subnetworks)
end

"""
    else_population_step(subnetwork, groups; rng)

Advance grouped trajectory counts through one shared subnetwork. One
multinomial column is drawn for each occupied entry state; individual
trajectories are never looped over.
"""
function else_population_step(
    subnetwork::ELSESubnetwork{S},
    groups::AbstractDict;
    rng=Random.GLOBAL_RNG,
) where {S}
    sources = collect(groups)
    next_groups = Dict{S,Int}()
    occupation = zeros(length(subnetwork.states))
    M_counts = zeros(Int, length(subnetwork.boundary), length(sources))
    mean_waiting_time = 0.0
    total_population = sum(values(groups))

    for (column, (start_state, population)) in enumerate(sources)
        law = else_law(subnetwork, start_state)
        occupation .+= population .* law.occupation

        counts = rand(
            rng,
            Distributions.Multinomial(
                population,
                law.exit_probability,
            ),
        )
        M_counts[:, column] .= counts

        mean_waiting_time +=
            population *
            dot(law.exit_probability, law.conditional_time)

        for boundary_position in findall(>(0), counts)
            boundary_count = counts[boundary_position]
            rates = subnetwork.exit_rates[boundary_position]
            channel_probability = rates ./ sum(rates)
            channel_counts = rand(
                rng,
                Distributions.Multinomial(
                    boundary_count,
                    channel_probability,
                ),
            )

            for channel in findall(>(0), channel_counts)
                next_state =
                    subnetwork.next_states[boundary_position][channel]
                next_groups[next_state] =
                    get(next_groups, next_state, 0) +
                    channel_counts[channel]
            end
        end
    end

    (
        groups=next_groups,
        occupation=occupation,
        M_counts=M_counts,
        mean_waiting_time=mean_waiting_time / total_population,
    )
end

"""
    else_population(subnetwork_for, initial_state, population, number_of_steps; rng)

Advance a population of trajectories together. Trajectories are grouped first
by their current subnetwork and then by their entry state.
"""
function else_population(
    subnetwork_for,
    initial_state::S,
    population::Integer,
    number_of_steps;
    rng=Random.GLOBAL_RNG,
) where {S}
    groups = Dict{S,Int}(initial_state => population)
    occupations = Dict{Symbol,Vector{Float64}}()
    history = NamedTuple[]

    for step_number in 1:number_of_steps
        partitions = Dict{Symbol,Dict{S,Int}}()
        subnetworks = Dict{Symbol,Any}()

        for (state, count) in groups
            subnetwork = subnetwork_for(state)
            partition =
                get!(partitions, subnetwork.name, Dict{S,Int}())
            partition[state] = get(partition, state, 0) + count
            subnetworks[subnetwork.name] = subnetwork
        end

        next_groups = Dict{S,Int}()

        for name in sort!(collect(keys(partitions)))
            subnetwork = subnetworks[name]
            result = else_population_step(
                subnetwork,
                partitions[name];
                rng=rng,
            )

            total_occupation = get!(
                occupations,
                name,
                zeros(length(subnetwork.states)),
            )
            total_occupation .+= result.occupation

            for (state, count) in result.groups
                next_groups[state] =
                    get(next_groups, state, 0) + count
            end

            push!(
                history,
                (
                    step=step_number,
                    subnetwork=name,
                    sources=length(partitions[name]),
                    M_size=size(result.M_counts),
                    M_counts=result.M_counts,
                    mean_waiting_time=result.mean_waiting_time,
                ),
            )
        end

        groups = next_groups
    end

    (
        history=history,
        final_groups=groups,
        occupations=occupations,
    )
end
