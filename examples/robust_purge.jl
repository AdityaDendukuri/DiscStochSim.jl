
function robust_purge!(
    X::Set{Element}, 
    p::Vector{T}, 
    model::Model,
    rates,
    t,
    prob_quantile::Number;
    flux_tolerance::Number = 1e-9 
) where {Element, T, Model}

    X_vec = collect(X)

	# 1. Candidate Selection: Find states with low probability mass
     idx = sortperm(p; rev=false)
    k = round(Int, prob_quantile * length(idx))            
    drop_idxs = k > 0 ? idx[1:k] : Int[]
	candidate_idxs = drop_idxs

    # Build new containers
    keep_mask = trues(length(p))
    keep_mask[drop_idxs] .= false

    flux_vector = zeros(T, length(p))
    local total_flux = 0
	if flux_tolerance > 0
    	for i in eachindex(X_vec)
        	current_state = X_vec[i]
        	weight = sum(prop(current_state, rates, t) 
						 for prop in model.propensities)
        	flux_vector[i] = p[i] * weight
    	end
    	total_flux = sum(flux_vector)
    	# Set an adaptive threshold based on the total system flux
   		flux_threshold = total_flux * flux_tolerance
		# Filter the candidates based on this adaptive threshold

    	final_idxs_to_prune = Set{Int}()
    	for idx in candidate_idxs
        	if flux_vector[idx] < flux_threshold
            	push!(final_idxs_to_prune, idx)
        	end
    	end
	else
		final_idxs_to_prune = candidate_idxs
	end
    # 3. Purge the final set of states
    new_p = [p[i] for i in eachindex(p) if i ∉ final_idxs_to_prune]
    states_to_remove = Set(X_vec[collect(final_idxs_to_prune)])
    new_X = setdiff(X, states_to_remove)
    return new_X, new_p, total_flux
end


