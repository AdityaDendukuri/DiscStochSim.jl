mutable struct FESSSubnetwork{S}
    network::FESSNetwork{S}
    R::Matrix{Float64}
    Z::Matrix{Float64}
    exit_rates::Vector{Float64}
    residual::Float64
end

function FESSSubnetwork(network::FESSNetwork{S}, R::AbstractMatrix) where {S}
    n = length(network)
    size(R) == (n, n) ||
        throw(ArgumentError("generator size does not match the FESS network"))
    generator = Matrix{Float64}(R)
    factor = lu(-generator)
    Z = factor \ Matrix{Float64}(I, n, n)
    residual = norm((-generator) * Z - I, 1) /
        max(norm(generator, 1) * norm(Z, 1) + 1.0, eps(Float64))
    exit_rates = max.(-vec(sum(generator; dims=1)), 0.0)
    FESSSubnetwork(network, generator, Z, exit_rates, residual)
end

function FESSSubnetwork(
    network::FESSNetwork,
    model,
    rates;
    time::Real=0.0,
)
    n = length(network)
    R = zeros(n, n)
    for (column, state) in enumerate(network.states)
        propensities = _propensities(state, model, rates, time)
        for reaction in eachindex(model.stoichvecs)
            rate = propensities[reaction]
            rate == 0.0 && continue
            destination = state + model.stoichvecs[reaction]
            row = get(network.index, destination, 0)
            row == 0 || (R[row, column] += rate)
            R[column, column] -= rate
        end
    end
    FESSSubnetwork(network, R)
end

occupation(subnetwork::FESSSubnetwork, entrance) = subnetwork.Z * entrance
second_occupation(subnetwork::FESSSubnetwork, first) = subnetwork.Z * first
escape_time(subnetwork::FESSSubnetwork, entrance) =
    sum(occupation(subnetwork, entrance))

function escape_probabilities(subnetwork::FESSSubnetwork, entrance)
    subnetwork.exit_rates .* occupation(subnetwork, entrance)
end

function conditional_time(
    subnetwork::FESSSubnetwork,
    occupation_vector,
    second_occupation_vector,
    exit_index::Int,
)
    second_occupation_vector[exit_index] / occupation_vector[exit_index]
end

"""
    cut_time_losses(subnetwork, entrance)

C++ `cut_time_losses`, transposed for Julia's column-oriented generator.
Entry `j` is the loss in mean escape time caused by removing state `j`.
"""
function cut_time_losses(subnetwork::FESSSubnetwork, entrance)
    u = occupation(subnetwork, entrance)
    column_sums = vec(sum(subnetwork.Z; dims=1))
    u ./ diag(subnetwork.Z) .* column_sums
end

function _shed_state!(subnetwork::FESSSubnetwork, remove::Int)
    keep = trues(length(subnetwork.network))
    keep[remove] = false
    Z = subnetwork.Z
    subnetwork.Z = Z[keep, keep] -
        Z[keep, remove] * transpose(Z[remove, keep]) / Z[remove, remove]
    subnetwork.R = subnetwork.R[keep, keep]
    subnetwork.exit_rates = max.(-vec(sum(subnetwork.R; dims=1)), 0.0)
    _remove!(subnetwork.network, remove)
    subnetwork
end

"""
    shed!(subnetwork, entrance, protected, options)

Apply the C++ rule one state at a time. Capacity is hard; otherwise a state is
removed only when its individual relative loss is below `cut_threshold`.
"""
function shed!(
    subnetwork::FESSSubnetwork,
    entrance::AbstractVector,
    protected,
    options::FESSOptions,
)
    _check(options)
    length(protected) < options.capacity ||
        throw(ArgumentError("capacity must exceed the protected entrance count"))

    p = collect(Float64, entrance)
    baseline = escape_time(subnetwork, p)
    baseline > options.atol ||
        throw(ArgumentError("FESS network has no positive escape time"))
    estimated_loss = 0.0
    maximum_loss = 0.0
    threshold_cuts = 0
    capacity_cuts = 0

    while length(subnetwork.network) > length(protected)
        losses = cut_time_losses(subnetwork, p)
        minimum(losses) >= -100 * options.atol ||
            error("FESS cut-time loss has a significant negative entry")
        losses .= max.(losses, 0.0)
        for state in protected
            losses[subnetwork.network.index[state]] = Inf
        end

        remove = argmin(losses)
        isfinite(losses[remove]) || break
        current_time = escape_time(subnetwork, p)
        relative_loss = losses[remove] / current_time
        capacity_cut = length(subnetwork.network) >= options.capacity
        threshold_cut = options.cut_threshold > 0 &&
            relative_loss < options.cut_threshold
        (capacity_cut || threshold_cut) || break

        capacity_cut ? (capacity_cuts += 1) : (threshold_cuts += 1)
        maximum_loss = max(maximum_loss, relative_loss)
        estimated_loss += losses[remove]
        deleteat!(p, remove)
        _shed_state!(subnetwork, remove)
    end

    length(subnetwork.network) < options.capacity ||
        error("capacity cannot be met without removing a protected entrance")
    final_time = escape_time(subnetwork, p)
    stats = FESSSheddingStats(
        threshold_cuts + capacity_cuts,
        threshold_cuts,
        capacity_cuts,
        estimated_loss / baseline,
        abs(final_time - baseline) / baseline,
        maximum_loss,
    )
    p, stats
end

