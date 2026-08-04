function _probabilities(weights, atol)
    probability = max.(collect(Float64, weights), 0.0)
    total = sum(probability)
    total > atol || throw(ArgumentError("FESS probability has zero mass"))
    probability ./ total
end

function _draw_index(rng, probability)
    target = rand(rng)
    cumulative = 0.0
    for index in eachindex(probability)
        cumulative += probability[index]
        target <= cumulative && return index
    end
    lastindex(probability)
end

function _entrance_vector(network::FESSNetwork, state)
    entrance = zeros(length(network))
    entrance[network.index[state]] = 1.0
    entrance
end

function _prepare!(
    workspace::FESSWorkspace,
    entrances,
    model,
    rates,
    options::FESSOptions;
    time::Real=0.0,
    entrance_distribution=nothing,
    shed::Bool=true,
)
    _check(options)
    augment!(
        workspace.network,
        entrances,
        model,
        rates;
        depth=options.expansion_depth,
        time=time,
        atol=options.atol,
    )
    subnetwork = FESSSubnetwork(workspace.network, model, rates; time=time)
    shed || return subnetwork, _NO_SHEDDING

    distribution = isnothing(entrance_distribution) ?
        _entrance_vector(workspace.network, only(entrances)) :
        entrance_distribution
    _, stats = shed!(subnetwork, distribution, entrances, options)
    subnetwork, stats
end

function _sample_exit(
    subnetwork::FESSSubnetwork,
    occupation_vector,
    second_vector,
    model,
    rates,
    rng;
    time::Real=0.0,
    atol::Real=1.0e-12,
)
    probability = _probabilities(
        subnetwork.exit_rates .* occupation_vector,
        atol,
    )
    exit_index = _draw_index(rng, probability)
    occupation_vector[exit_index] > atol ||
        error("sampled FESS exit has zero occupation")
    waiting_time = conditional_time(
        subnetwork,
        occupation_vector,
        second_vector,
        exit_index,
    )
    isfinite(waiting_time) && waiting_time > atol ||
        error("FESS produced an invalid waiting time: $waiting_time")

    preexit = subnetwork.network.states[exit_index]
    propensities = _propensities(preexit, model, rates, time)
    reactions = Int[]
    outward_rates = Float64[]
    for reaction in eachindex(model.stoichvecs)
        destination = preexit + model.stoichvecs[reaction]
        if propensities[reaction] > atol &&
           !haskey(subnetwork.network.index, destination)
            push!(reactions, reaction)
            push!(outward_rates, propensities[reaction])
        end
    end
    isempty(reactions) &&
        error("sampled pre-exit state has no transition leaving the network")
    reaction = reactions[_draw_index(rng, _probabilities(outward_rates, atol))]
    preexit + model.stoichvecs[reaction], preexit, waiting_time
end

"""
    fess_step!(workspace, state, model, rates; options, rng, time, shed)

Take one state-first FESS step. Reusing `workspace` gives the persistent method;
passing a fresh workspace and `shed=false` gives the redrawn method.
"""
function fess_step!(
    workspace::FESSWorkspace{S},
    state::S,
    model,
    rates;
    options::FESSOptions=FESSOptions(),
    rng=Random.GLOBAL_RNG,
    time::Real=0.0,
    shed::Bool=true,
) where {S}
    _terminal(state, model, rates, time, options.atol) &&
        return FESSStep(
            state, state, 0.0, 0.0, true, length(workspace.network),
            _NO_SHEDDING, 0.0,
        )

    subnetwork, shedding = _prepare!(
        workspace,
        [state],
        model,
        rates,
        options;
        time=time,
        shed=shed,
    )
    entrance = _entrance_vector(subnetwork.network, state)
    first = max.(occupation(subnetwork, entrance), 0.0)
    second = second_occupation(subnetwork, first)
    next_state, preexit, waiting_time = _sample_exit(
        subnetwork,
        first,
        second,
        model,
        rates,
        rng;
        time=time,
        atol=options.atol,
    )
    FESSStep(
        next_state,
        preexit,
        waiting_time,
        sum(first),
        false,
        length(subnetwork.network),
        shedding,
        subnetwork.residual,
    )
end

function fess_step(
    state::S,
    model,
    rates;
    options::FESSOptions=FESSOptions(),
    rng=Random.GLOBAL_RNG,
    time::Real=0.0,
) where {S}
    workspace = FESSWorkspace(state)
    fess_step!(
        workspace,
        state,
        model,
        rates;
        options=options,
        rng=rng,
        time=time,
        shed=false,
    )
end
