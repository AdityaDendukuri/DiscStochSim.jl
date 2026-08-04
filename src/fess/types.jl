Base.@kwdef struct FESSOptions
    capacity::Int = 60
    cut_threshold::Float64 = 1.0e-3
    expansion_depth::Int = 1
    atol::Float64 = 1.0e-12
end

function _check(options::FESSOptions)
    options.capacity >= 2 ||
        throw(ArgumentError("FESS capacity must be at least two"))
    options.cut_threshold >= 0 ||
        throw(ArgumentError("FESS cut threshold must be nonnegative"))
    options.expansion_depth >= 0 ||
        throw(ArgumentError("FESS expansion depth must be nonnegative"))
    options.atol > 0 ||
        throw(ArgumentError("FESS tolerance must be positive"))
    options
end

struct FESSSheddingStats
    removed::Int
    threshold_cuts::Int
    capacity_cuts::Int
    estimated_relative_loss::Float64
    actual_relative_loss::Float64
    maximum_single_cut_loss::Float64
end

const _NO_SHEDDING = FESSSheddingStats(0, 0, 0, 0.0, 0.0, 0.0)

struct FESSStep{S}
    next_state::S
    preexit_state::S
    waiting_time::Float64
    mean_waiting_time::Float64
    terminated::Bool
    network_size::Int
    shedding::FESSSheddingStats
    residual::Float64
end

struct FESSTrajectory{S}
    times::Vector{Float64}
    states::Vector{S}
    step_times::Vector{Float64}
    steps::Vector{FESSStep{S}}
end

struct FESSEnsemble{S}
    times::Vector{Float64}
    mean_states::Vector{Vector{Float64}}
    variance_states::Vector{Vector{Float64}}
    final_times::Vector{Float64}
    final_states::Vector{S}
    step_counts::Vector{Int}
    maximum_network_sizes::Vector{Int}
    step_network_sizes::Vector{Int}
    step_removed_states::Vector{Int}
    active_trajectory_counts::Vector{Int}
    starting_state_counts::Vector{Int}
    residuals::Vector{Float64}
    shared_steps::Int
    trajectory_times::Vector{Vector{Float64}}
    trajectory_states::Vector{Vector{S}}
    wall_time::Float64
    persistent::Bool
end
