### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 453819dc-e207-46db-9c9b-9f10f7fd4daf
begin
	# numerical libraries
	using ExponentialUtilities, PROPACK, Arpack, SparseArrays, ExponentialUtilities, Expokit
	# output and plotting
	using ProgressLogging, JLD
	# modelling and statistics 
	using Catalyst, JumpProcesses, StatsBase, DifferentialEquations
	using Interpolations, LinearAlgebra
	# importing local fsp package
	using Revise, CairoMakie
	local_mod = include("../../src/DiscStochSim.jl")
	using .local_mod.DiscStochSim
end

# ╔═╡ cb9fa0aa-1588-4647-8b02-69fca1415974
rn = @reaction_network begin
    k₁, A --> B 
	k₂, B --> B + C
end

# ╔═╡ ae7ff050-c76d-4c0d-8c8a-23724f10485a
model = DiscreteStochasticSystem(rn);

# ╔═╡ 85c4ad39-d14a-4e73-bfa3-4a22fef4d537
def_params = begin
	# reaction rates 
	rates = [1e-6, 1e-1]
	# boundary 
	bounds = (0, 1e11) #(lower limit, upper limit)
	boundary_condition(x) = RectLatticeBoundaryCondition(x, bounds);
	# time interval and initial values
end

# ╔═╡ e3e0c398-9db0-420a-86b3-65220761028c
function robust_purge1!(
    X::Set{Element}, 
    p::Vector{T}, 
    model::Model,
    rates,
    t,
    prob_quantile::Number;
    flux_tolerance::Number = 1e-9 
) where {Element, T, Model}

    X_vec = collect(X)

    candidate_idxs = DiscStochSim.findLowestValuesPercent_naive(p, prob_quantile)

    flux_vector = zeros(T, length(p))

    for i in eachindex(X_vec)
        current_state = X_vec[i]
        weight = sum(prop(current_state, rates, t) for prop in model.propensities)
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

    # 3. Purge the final set of states
    new_p = [p[i] for i in eachindex(p) if i ∉ final_idxs_to_prune]
    states_to_remove = Set(X_vec[collect(final_idxs_to_prune)])
    new_X = setdiff(X, states_to_remove)

    return new_X, new_p, total_flux
end

# ╔═╡ 91bec3fc-65cc-4344-98d2-09f4944fc0e5
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
	candidate_idxs = DiscStochSim.findLowestValuesPercent_naive(p, prob_quantile)

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

# ╔═╡ ebd5f20c-88c2-46d1-982f-de864ba33727
md"""
### Initialize Finite State Space 
"""

# ╔═╡ 256f8d0d-7c2b-4415-969e-356ee271e908
function reconstruct_MasterEquation(
    S_new::Set{E},
    S_old_vec::Vector{E},         
    A_old::SparseMatrixCSC{T,Int}, 
    model, rates, bc::Function, t::Real;
    overlap_threshold::Real=0.3,
) where {E,T}

    # --- Setup ---
    S_new_vec = collect(S_new)
    n_new = length(S_new_vec)
    new_id = Dict{E,Int}(s => i for (i,s) in enumerate(S_new_vec))
    old_id  = Dict{E,Int}(s => i for (i,s) in enumerate(S_old_vec))

    # --- Overlap check (can be kept as a performance heuristic) ---
    n_retained = count(s -> haskey(old_id, s), S_new_vec)
    if n_new == 0 || n_retained / n_new < overlap_threshold
        return MasterEquation(S_new, model, rates, bc, t)
    end
    
    I = Int[]; J = Int[]; V = T[]
    @inbounds for (j_new, x) in enumerate(S_new_vec)
        col_sum = zero(T)
        
        # Build the entire column from scratch based on propensities
        for (k, ν) in enumerate(model.stoichvecs)
            y = x + ν
            # Check if the destination `y` is in our new state space
            i_new = get(new_id, y, 0)
            
            if i_new != 0
                α = model.propensities[k](x, rates, t)
                if α > 0
                    # Allow duplicates for sparse() to sum them up
                    push!(I, i_new); push!(J, j_new); push!(V, α)
                    col_sum += α
                end
            end
        end
        
        # Diagonal = minus the sum of all off-diagonals in this column
        if col_sum > 0
            push!(I, j_new); push!(J, j_new); push!(V, -col_sum)
        end
    end

    A = sparse(I, J, V, n_new, n_new)
    out_flow = spdiagm(0 => diag(A))
    in_flow  = A - out_flow
    
    return A, in_flow, out_flow
end

# ╔═╡ ca5f05c8-96ea-42a8-adb0-c095a223b4cf

fsp_sim_no_flux = begin
		# 4. Initialize p₀
	local U₀ = CartesianIndex(1,0,0)
	local 𝒮₀ = Set([U₀])
	local p₀ = zeros(length(𝒮₀))
	local p₀[FindElement(U₀,𝒮₀)] = 1.0
	# copy initial values 
	local pₜ = copy(p₀)
	local 𝒮ₜ = copy(𝒮₀)
	# time stepping loop
	local iter = 1
	local pₜ = p₀
	local δt = 1e2
	local T =0:δt:1e8
	
	# variables to store simulation observables
	size_𝒮ₜ1 = Int.(zeros(length(T))) # system sizes
	flux1 = zeros(length(T))
	flux_tol1 = zeros(length(T))
	local_ϵ1 = zeros(length(T))
	sol1 = []
	
	@progress for (iter, t) ∈ enumerate(T)
		
		# expand state space
	    𝒮ₜ, pₜ = expand!(𝒮ₜ, pₜ, model, boundary_condition, 1)
		global size_𝒮ₜ1[iter] = length(𝒮ₜ)
		local A, in_flow, out_flow = MasterEquation(𝒮ₜ, 
													 model,
													 rates,
													 boundary_condition, 
													 t)

		pₜ = expv(δt, A, pₜ)
    	push!(sol1, (𝒮ₜ, pₜ))
		
		# purge state space
		𝒮ₜ, pₜ, total_flux = robust_purge!(𝒮ₜ, 
									pₜ, 
									model,
									rates,
									δt,
									0.7; flux_tolerance=0.0)
	 	

		local_ϵ1[iter] = 1.0 - sum(pₜ)
    	pₜ ./= sum(pₜ)
		flux1[iter] = max(0.0, 1.0 - δt * total_flux)
	end
end


# ╔═╡ 6ad6db9b-f87c-4320-9a45-a18a170b5db7
fsp_sim_flux = begin
		# 4. Initialize p₀
	local U₀ = CartesianIndex(1,0,0)
	local 𝒮₀ = Set([U₀])
	local p₀ = zeros(length(𝒮₀))
	local p₀[FindElement(U₀,𝒮₀)] = 1.0
	# copy initial values 
	local pₜ = copy(p₀)
	local 𝒮ₜ = copy(𝒮₀)
	# time stepping loop
	local iter = 1
	local pₜ = p₀
	local δt = 10.0
	global savefreq=[]
	global T =0:δt:1e5
	
	# variables to store simulation observables
	size_𝒮ₜ = Int.(zeros(length(T))) # system sizes
	flux = zeros(length(T))
	flux_tol = zeros(length(T))
	local_ϵ = zeros(length(T))
	sol = []
	
	@progress for (iter, t) ∈ enumerate(T)
		
		# expand state space
	     𝒮ₜ, pₜ = expand!(𝒮ₜ, pₜ, model, boundary_condition, 1)
		local A, in_flow, out_flow = MasterEquation(𝒮ₜ, 
													 model,
													 rates,
													 boundary_condition, 
													 t)

		pₜ = expmv(δt, A, pₜ)

		if iter%savefreq == 0
			push!(size_𝒮ₜ[iter], length(𝒮ₜ))
    		push!(sol, (𝒮ₜ, pₜ))
			push!(times, )
		end
		
		# purge state space
		𝒮ₜ, pₜ, total_flux = robust_purge1!(𝒮ₜ, 
									pₜ, 
									model,
									rates,
									δt,
									0.9; flux_tolerance=0.000001)
	 	


		local_ϵ[iter] = 1.0 - sum(pₜ)
    	pₜ ./= sum(pₜ)
		#flux[iter] = max(0.0, 1.0 - δt * total_flux)
		pₜ ./= sum(pₜ)		
	end
	global afsp_mean = map(1:length(sol)) do i
		sum(collect.(Tuple.(sol[i][1])) .* sol[i][2])
	end 
	afsp_mean=hcat(afsp_mean...)'
end

# ╔═╡ 7ea7ecaf-7a46-4b9b-a0e0-56171ffd7c23
size_𝒮ₜ 

# ╔═╡ 18962357-9f6a-4668-82cb-7f133a130b5e
fsp_full = begin
	local U₀ = CartesianIndex(1,0,0)
	global 𝒮ₜ = Set{CartesianIndex{3}}()
	local iter = 1
	local δt = 1e4
	local T =0:δt:Integer(1e5)
	for i in 0:2
		for j in 0:2
			for k in 0:10000
				push!(𝒮ₜ, CartesianIndex(i, j, k))
			end
		end
	end
	idxmap = Dict([collect(𝒮ₜ)[i] => i for i in length(𝒮ₜ)])
	local pₜ = zeros(length(𝒮ₜ)) 
	local pₜ[FindElement(U₀,𝒮ₜ)] = 1.0
	local A, in_flow, out_flow = MasterEquation(𝒮ₜ, 
												model,
												rates,
												boundary_condition, 
												δt)
	global sol_exact = []
	@progress for t ∈ T
		p = expmv(t, A, pₜ)
		push!(sol_exact, p)
		iter+=1
	end
	global exact_mean = map(1:length(sol_exact)) do i
		    sum(collect.(Tuple.(𝒮ₜ)) .* sol_exact[i])
	end 
	exact_mean=hcat(exact_mean...)'
end

# ╔═╡ 7813ebfd-f674-44da-b2d5-dfed9e028457
begin
	# Distribution Comparison Functions - CairoMakie Version for Pluto
	# ============================================================================
	
	
	# ============================================================================
	# 1. MARGINAL DISTRIBUTIONS
	# ============================================================================
	
	"""
	Compute marginal distribution for a single species
	"""
	function compute_marginal(𝒮, p, species_idx)
	    # species_idx: 1=A, 2=B, 3=C
	    𝒮_vec = collect(𝒮)
	    
	    # Find range of values for this species
	    values = [state[species_idx] for state in 𝒮_vec]
	    min_val = minimum(values)
	    max_val = maximum(values)
	    
	    # Compute marginal probabilities
	    marginal = zeros(max_val - min_val + 1)
	    for (i, state) in enumerate(𝒮_vec)
	        idx = state[species_idx] - min_val + 1
	        marginal[idx] += p[i]
	    end
	    
	    return (values=collect(min_val:max_val), probs=marginal)
	end
	
	"""
	Compare marginal distributions between adaptive and exact
	"""
	function compare_marginals(time_idx_adaptive, time_idx_exact, species_idx; 
	                          species_name="Species")
	    
	    # Get solutions
	    𝒮_adaptive, p_adaptive = sol[time_idx_adaptive]
	    p_exact = sol_exact[time_idx_exact]
	    
	    # Compute marginals
	    marg_adaptive = compute_marginal(𝒮_adaptive, p_adaptive, species_idx)
	    marg_exact = compute_marginal(𝒮ₜ, p_exact, species_idx)
	    
	    # Compute error on overlapping support
	    min_val = max(marg_adaptive.values[1], marg_exact.values[1])
	    max_val = min(marg_adaptive.values[end], marg_exact.values[end])
	    
	    l1_error = 0.0
	    for val in min_val:max_val
	        idx_a = val - marg_adaptive.values[1] + 1
	        idx_e = val - marg_exact.values[1] + 1
	        
	        p_a = (1 <= idx_a <= length(marg_adaptive.probs)) ? marg_adaptive.probs[idx_a] : 0.0
	        p_e = (1 <= idx_e <= length(marg_exact.probs)) ? marg_exact.probs[idx_e] : 0.0
	        
	        l1_error += abs(p_a - p_e)
	    end
	    
	    # Account for probability mass outside overlapping region
	    for (i, val) in enumerate(marg_adaptive.values)
	        if val < marg_exact.values[1] || val > marg_exact.values[end]
	            l1_error += marg_adaptive.probs[i]
	        end
	    end
	    for (i, val) in enumerate(marg_exact.values)
	        if val < marg_adaptive.values[1] || val > marg_adaptive.values[end]
	            l1_error += marg_exact.probs[i]
	        end
	    end
	    
	    tv_distance = l1_error / 2.0
	    
	    println("\n$species_name Marginal Distribution:")
	    println("  L1 Error:     $l1_error")
	    println("  TV Distance:  $tv_distance")
	    println("  Support (Adaptive): [$(marg_adaptive.values[1]), $(marg_adaptive.values[end])]")
	    println("  Support (Exact):    [$(marg_exact.values[1]), $(marg_exact.values[end])]")
	    
	    return (adaptive=marg_adaptive, exact=marg_exact, l1=l1_error, tv=tv_distance)
	end
	
	"""
	Plot marginal distribution comparison using CairoMakie (black and white)
	"""
	function plot_marginal_comparison(time_idx_adaptive, time_idx_exact, species_idx;
	                                 species_name="Species")
	    
	    result = compare_marginals(time_idx_adaptive, time_idx_exact, species_idx; 
	                              species_name=species_name)
	    
	    fig = Figure(resolution=(600, 400))
	    ax = Axis(fig[1, 1], 
	             xlabel="$species_name Count", 
	             ylabel="Probability",
	             title="$species_name Marginal (TV=$(round(result.tv, digits=4)))")
	    
	    # Plot exact (solid line, circles)
	    lines!(ax, result.exact.values, result.exact.probs, 
	          label="Exact", linewidth=2, color=:black)
	    scatter!(ax, result.exact.values, result.exact.probs, 
	            color=:black, markersize=10, marker=:circle)
	    
	    # Plot adaptive (dashed line, squares)
	    lines!(ax, result.adaptive.values, result.adaptive.probs,
	          label="Adaptive", linewidth=2, color=:black, linestyle=:dash)
	    scatter!(ax, result.adaptive.values, result.adaptive.probs,
	            color=:black, markersize=10, marker=:rect)
	    
	    axislegend(ax, position=:rt)
	    
	    return fig
	end
	
	# ============================================================================
	# 2. JOINT DISTRIBUTIONS (2D)
	# ============================================================================
	
	"""
	Compute 2D joint distribution for two species
	"""
	function compute_joint_2d(𝒮, p, species_idx1, species_idx2)
	    𝒮_vec = collect(𝒮)
	    
	    # Find ranges
	    vals1 = [state[species_idx1] for state in 𝒮_vec]
	    vals2 = [state[species_idx2] for state in 𝒮_vec]
	    
	    min1, max1 = minimum(vals1), maximum(vals1)
	    min2, max2 = minimum(vals2), maximum(vals2)
	    
	    # Create 2D array
	    joint = zeros(max1 - min1 + 1, max2 - min2 + 1)
	    
	    for (i, state) in enumerate(𝒮_vec)
	        idx1 = state[species_idx1] - min1 + 1
	        idx2 = state[species_idx2] - min2 + 1
	        joint[idx1, idx2] += p[i]
	    end
	    
	    return (values1=collect(min1:max1), values2=collect(min2:max2), probs=joint)
	end
	
	"""
	Compare 2D joint distributions
	"""
	function compare_joint_2d(time_idx_adaptive, time_idx_exact, species_idx1, species_idx2;
	                         species_name1="Species 1", species_name2="Species 2")
	    
	    # Get solutions
	    𝒮_adaptive, p_adaptive = sol[time_idx_adaptive]
	    p_exact = sol_exact[time_idx_exact]
	    
	    # Compute joints
	    joint_adaptive = compute_joint_2d(𝒮_adaptive, p_adaptive, species_idx1, species_idx2)
	    joint_exact = compute_joint_2d(𝒮ₜ, p_exact, species_idx1, species_idx2)
	    
	    # Compute error on overlapping support
	    min1 = max(joint_adaptive.values1[1], joint_exact.values1[1])
	    max1 = min(joint_adaptive.values1[end], joint_exact.values1[end])
	    min2 = max(joint_adaptive.values2[1], joint_exact.values2[1])
	    max2 = min(joint_adaptive.values2[end], joint_exact.values2[end])
	    
	    l1_error = 0.0
	    
	    # Compare on overlapping region
	    for val1 in min1:max1
	        for val2 in min2:max2
	            idx1_a = val1 - joint_adaptive.values1[1] + 1
	            idx2_a = val2 - joint_adaptive.values2[1] + 1
	            idx1_e = val1 - joint_exact.values1[1] + 1
	            idx2_e = val2 - joint_exact.values2[1] + 1
	            
	            p_a = joint_adaptive.probs[idx1_a, idx2_a]
	            p_e = joint_exact.probs[idx1_e, idx2_e]
	            
	            l1_error += abs(p_a - p_e)
	        end
	    end
	    
	    # Add probability mass outside overlapping region (simplified)
	    total_mass_adaptive = sum(joint_adaptive.probs)
	    total_mass_exact = sum(joint_exact.probs)
	    
	    tv_distance = l1_error / 2.0
	    
	    println("\n$species_name1-$species_name2 Joint Distribution:")
	    println("  L1 Error:     $l1_error")
	    println("  TV Distance:  $tv_distance")
	    
	    return (adaptive=joint_adaptive, exact=joint_exact, l1=l1_error, tv=tv_distance)
	end
	
	"""
	Plot 2D joint distribution comparison using CairoMakie heatmaps (grayscale)
	"""
	function plot_joint_comparison(time_idx_adaptive, time_idx_exact, species_idx1, species_idx2;
	                              species_name1="Species 1", species_name2="Species 2")
	    
	    result = compare_joint_2d(time_idx_adaptive, time_idx_exact, species_idx1, species_idx2;
	                             species_name1=species_name1, species_name2=species_name2)
	    
	    fig = Figure(resolution=(1500, 400))
	    
	    # Plot exact
	    ax1 = Axis(fig[1, 1], 
	              xlabel=species_name2, ylabel=species_name1,
	              title="Exact Joint")
	    hm1 = heatmap!(ax1, result.exact.values2, result.exact.values1, result.exact.probs',
	                  colormap=:grays)
	    Colorbar(fig[1, 2], hm1)
	    
	    # Plot adaptive
	    ax2 = Axis(fig[1, 3],
	              xlabel=species_name2, ylabel=species_name1,
	              title="Adaptive Joint")
	    hm2 = heatmap!(ax2, result.adaptive.values2, result.adaptive.values1, result.adaptive.probs',
	                  colormap=:grays, colorrange=(0, maximum(result.exact.probs)))
	    Colorbar(fig[1, 4], hm2)
	    
	    # Plot difference on common grid
	    min1 = max(result.exact.values1[1], result.adaptive.values1[1])
	    max1 = min(result.exact.values1[end], result.adaptive.values1[end])
	    min2 = max(result.exact.values2[1], result.adaptive.values2[1])
	    max2 = min(result.exact.values2[end], result.adaptive.values2[end])
	    
	    diff = zeros(max1-min1+1, max2-min2+1)
	    for (i, val1) in enumerate(min1:max1)
	        for (j, val2) in enumerate(min2:max2)
	            idx1_e = val1 - result.exact.values1[1] + 1
	            idx2_e = val2 - result.exact.values2[1] + 1
	            idx1_a = val1 - result.adaptive.values1[1] + 1
	            idx2_a = val2 - result.adaptive.values2[1] + 1
	            
	            diff[i, j] = result.exact.probs[idx1_e, idx2_e] - result.adaptive.probs[idx1_a, idx2_a]
	        end
	    end
	    
	    ax3 = Axis(fig[1, 5],
	              xlabel=species_name2, ylabel=species_name1,
	              title="Difference (Exact - Adaptive)")
	    max_diff = maximum(abs.(diff))
	    hm3 = heatmap!(ax3, collect(min2:max2), collect(min1:max1), diff',
	                  colormap=:grays, colorrange=(-max_diff, max_diff))
	    Colorbar(fig[1, 6], hm3)
	    
	    return fig
	end
	
	# ============================================================================
	# 3. COMPREHENSIVE COMPARISON
	# ============================================================================
	
	"""
	Complete distribution comparison at a single time point
	"""
	function full_distribution_comparison(time_idx_adaptive, time_idx_exact)
	    
	    println("="^70)
	    println("DISTRIBUTION COMPARISON")
	    println("  Adaptive time index: $time_idx_adaptive")
	    println("  Exact time index:    $time_idx_exact")
	    println("="^70)
	    
	    # Marginals
	    println("\n--- MARGINAL DISTRIBUTIONS ---")
	    marg_A = compare_marginals(time_idx_adaptive, time_idx_exact, 1; species_name="A")
	    marg_B = compare_marginals(time_idx_adaptive, time_idx_exact, 2; species_name="B")
	    marg_C = compare_marginals(time_idx_adaptive, time_idx_exact, 3; species_name="C")
	    
	    # Joints
	    println("\n--- JOINT DISTRIBUTIONS ---")
	    joint_AB = compare_joint_2d(time_idx_adaptive, time_idx_exact, 1, 2; 
	                                species_name1="A", species_name2="B")
	    joint_BC = compare_joint_2d(time_idx_adaptive, time_idx_exact, 2, 3;
	                                species_name1="B", species_name2="C")
	    joint_AC = compare_joint_2d(time_idx_adaptive, time_idx_exact, 1, 3;
	                                species_name1="A", species_name2="C")
	    
	    # Means
	    println("\n--- MEAN VALUES ---")
	    mean_adaptive = afsp_mean[time_idx_adaptive, :]
	    mean_exact = exact_mean[time_idx_exact, :]
	    println("  Adaptive: A=$(mean_adaptive[1]), B=$(mean_adaptive[2]), C=$(mean_adaptive[3])")
	    println("  Exact:    A=$(mean_exact[1]), B=$(mean_exact[2]), C=$(mean_exact[3])")
	    println("  Errors:   A=$(mean_adaptive[1]-mean_exact[1]), B=$(mean_adaptive[2]-mean_exact[2]), C=$(mean_adaptive[3]-mean_exact[3])")
	    
	    println("\n" * "="^70)
	    
	    return (marginals=(A=marg_A, B=marg_B, C=marg_C),
	            joints=(AB=joint_AB, BC=joint_BC, AC=joint_AC))
	end
	
	"""
	Create comprehensive visualization with CairoMakie (black and white)
	"""
	function plot_full_comparison(time_idx_adaptive, time_idx_exact)
	    
	    # Get all marginal results
	    marg_A = compare_marginals(time_idx_adaptive, time_idx_exact, 1; species_name="A")
	    marg_B = compare_marginals(time_idx_adaptive, time_idx_exact, 2; species_name="B")
	    marg_C = compare_marginals(time_idx_adaptive, time_idx_exact, 3; species_name="C")
	    
	    # Create figure with all three marginals
	    fig = Figure(resolution=(1500, 400))
	    
	    # Species A
	    ax1 = Axis(fig[1, 1], xlabel="A Count", ylabel="Probability",
	              title="A Marginal (TV=$(round(marg_A.tv, digits=4)))")
	    lines!(ax1, marg_A.exact.values, marg_A.exact.probs, label="Exact", linewidth=2, color=:black)
	    scatter!(ax1, marg_A.exact.values, marg_A.exact.probs, color=:black, markersize=10, marker=:circle)
	    lines!(ax1, marg_A.adaptive.values, marg_A.adaptive.probs, label="Adaptive", 
	          linewidth=2, color=:black, linestyle=:dash)
	    scatter!(ax1, marg_A.adaptive.values, marg_A.adaptive.probs, 
	            color=:black, markersize=10, marker=:rect)
	    axislegend(ax1, position=:rt)
	    
	    # Species B
	    ax2 = Axis(fig[1, 2], xlabel="B Count", ylabel="Probability",
	              title="B Marginal (TV=$(round(marg_B.tv, digits=4)))")
	    lines!(ax2, marg_B.exact.values, marg_B.exact.probs, label="Exact", linewidth=2, color=:black)
	    scatter!(ax2, marg_B.exact.values, marg_B.exact.probs, color=:black, markersize=10, marker=:circle)
	    lines!(ax2, marg_B.adaptive.values, marg_B.adaptive.probs, label="Adaptive",
	          linewidth=2, color=:black, linestyle=:dash)
	    scatter!(ax2, marg_B.adaptive.values, marg_B.adaptive.probs,
	            color=:black, markersize=10, marker=:rect)
	    axislegend(ax2, position=:rt)
	    
	    # Species C
	    ax3 = Axis(fig[1, 3], xlabel="C Count", ylabel="Probability",
	              title="C Marginal (TV=$(round(marg_C.tv, digits=4)))")
	    lines!(ax3, marg_C.exact.values, marg_C.exact.probs, label="Exact", linewidth=2, color=:black)
	    scatter!(ax3, marg_C.exact.values, marg_C.exact.probs, color=:black, markersize=10, marker=:circle)
	    lines!(ax3, marg_C.adaptive.values, marg_C.adaptive.probs, label="Adaptive",
	          linewidth=2, color=:black, linestyle=:dash)
	    scatter!(ax3, marg_C.adaptive.values, marg_C.adaptive.probs,
	            color=:black, markersize=10, marker=:rect)
	    axislegend(ax3, position=:rt)
	    
	    return fig
	end
	
	# ============================================================================
	# USAGE IN PLUTO
	# ============================================================================
	
	#=
	# In your Pluto notebook, after running fsp_full and fsp_sim_flux:
	
	# Cell 1: Text comparison
	full_distribution_comparison(length(sol), length(sol_exact))
	
	# Cell 2: Plot all marginals (returns Figure directly)
	plot_full_comparison(length(sol), length(sol_exact))
	
	# Cell 3: Plot individual marginal
	plot_marginal_comparison(length(sol), length(sol_exact), 3; species_name="C")
	
	# Cell 4: Plot joint distribution
	plot_joint_comparison(length(sol), length(sol_exact), 2, 3; 
	                     species_name1="B", species_name2="C")
	
	# Cell 5: Time evolution of errors
	let
	    time_indices = 1:10:length(sol)
	    tv_A = Float64[]
	    tv_B = Float64[]
	    tv_C = Float64[]
	    
	    for i in time_indices
	        if i <= length(sol_exact)
	            result_A = compare_marginals(i, i, 1; species_name="A")
	            result_B = compare_marginals(i, i, 2; species_name="B")
	            result_C = compare_marginals(i, i, 3; species_name="C")
	            push!(tv_A, result_A.tv)
	            push!(tv_B, result_B.tv)
	            push!(tv_C, result_C.tv)
	        end
	    end
	    
	    fig = Figure(resolution=(800, 500))
	    ax = Axis(fig[1, 1], xlabel="Time Index", ylabel="TV Distance",
	             title="Marginal Distribution Errors", yscale=log10)
	    lines!(ax, time_indices[1:length(tv_A)], tv_A, label="A", linewidth=2)
	    lines!(ax, time_indices[1:length(tv_B)], tv_B, label="B", linewidth=2)
	    lines!(ax, time_indices[1:length(tv_C)], tv_C, label="C", linewidth=2)
	    axislegend(ax, position=:rt)
	    fig
	end
	=#
end

# ╔═╡ fd852a0c-24fe-4807-afb5-c9157d3c746a
begin
	local fig = Figure()
	local ax = Axis(fig[1, 1])
	#plot!(ax, exact_mean[:, 1])
	plot!(ax, exact_mean[:, 1], label="A")
	plot!(ax, exact_mean[:, 2], label="B")
	plot!(ax, exact_mean[:, 3] ./ maximum(exact_mean[:, 3]), label="C")

	axislegend(ax)
	fig
end

# ╔═╡ f541864a-38f7-4081-9a71-a1c6de95e5cf
global T1 = copy(T)

# ╔═╡ 6ba472e0-0c0f-4840-9640-fa4e848e2f8d
begin
    using LsqFit, Printf
    local δt = 1e4
    local T = 0:δt:Integer(1e5)
    
    # Match timepoints
    idxs = [findfirst(T1 .== t) for t in T]
    
    # SIAM theme
    siam_theme = Theme(
        fontsize = 10,
        font = "CMU Serif",
        Axis = (
            titlesize = 10,
            xlabelsize = 10,
            ylabelsize = 10,
            xticklabelsize = 9,
            yticklabelsize = 9,
            xgridvisible = false,
            ygridvisible = false,
            leftspinevisible = true,
            rightspinevisible = true,
            topspinevisible = true,
            bottomspinevisible = true,
            spinewidth = 0.9,
            xtickwidth = 0.9,
            ytickwidth = 0.9,
            xtickalign = 1,
            ytickalign = 1,
            xticksize = 3,
            yticksize = 3,
        ),
        Legend = (
            framevisible = false,
            patchsize = (20, 3),
            labelsize = 9,
            rowgap = 2,
            padding = (2, 2, 2, 2),
        ),
        Lines = (
            linewidth = 1.2,
        ),
        backgroundcolor = :white,
    )
    set_theme!(siam_theme)
    
    # Compute errors
    local errors_C = abs.(afsp_mean[idxs, 3] .- exact_mean[:, 3])
    local rel_errors_C = errors_C ./ exact_mean[:, 3]
    
    local fig = Figure(size=(504, 380))
    
    # Panel (a): Relative error
    local ax_a = Axis(fig[1, 1],
        xlabel = "Time (s)",
        ylabel = "Relative error (%)",
        title = "(a) Relative error in ⟨C⟩")
    
    lines!(ax_a, T, rel_errors_C .* 100,
        color = :black, linewidth = 1.2)
    
    # FIXED: Use lines! instead of hlines!
    lines!(ax_a, [T[1], T[end]], [0.01, 0.01], 
        color = (:gray, 0.5), linestyle = :dash, linewidth = 1)
    
    # Panel (b): Error rate
    local ax_b = Axis(fig[1, 2],
        xlabel = "Time (s)",
        ylabel = "Error rate (s⁻¹)",
        title = "(b) Instantaneous error rate")
    
    error_rate = diff(errors_C) ./ δt
    lines!(ax_b, T[1:end-1], error_rate,
        color = :black, linewidth = 1.2)
    
    # Panel (c): Mean trajectories
    local ax_c = Axis(fig[2, 1],
        xlabel = "Time (s)",
        ylabel = "Mean concentration",
        title = "(c) Mean trajectories")
    
    lines!(ax_c, T, exact_mean[:, 3],
        color = :black, linewidth = 1.5, linestyle = :solid, 
        label = "Full FSP")
    
    lines!(ax_c, T, afsp_mean[idxs, 3],
        color = (:black, 0.7), linewidth = 1.5, linestyle = :dash,
        label = "Adaptive FSP")
    
    lines!(ax_c, T, exact_mean[:, 1],
        color = (:gray50, 0.6), linewidth = 1.0, linestyle = :dot,
        label = "Species A")
    
    axislegend(ax_c, position = :lt, labelsize = 8)
    
    # Panel (d): Error scaling
    local ax_d = Axis(fig[2, 2],
        xlabel = "Mean population of C",
        ylabel = "Absolute error",
        title = "(d) Error scaling")
    
    scatter!(ax_d, exact_mean[:, 3], errors_C,
        color = (:black, 0.3), marker = :circle, markersize = 3)
    
    # Fit
    fit_idx = findfirst(exact_mean[:, 3] .> 50)
    lmodel(C, p) = p[1] .* C .+ p[2]
    fit = curve_fit(lmodel, 
                    exact_mean[fit_idx:end, 3], 
                    errors_C[fit_idx:end], 
                    [0.0001, 0.0])
    
    C_range = range(0, maximum(exact_mean[:, 3]), length=100)
    slope = coef(fit)[1]
    lines!(ax_d, C_range, lmodel(C_range, coef(fit)),
        color = :black, linewidth = 1.5, linestyle = :dash,
        label = @sprintf("ε = %.2e⟨C⟩ ", slope))
    
    text!(ax_d, 200, 0.055, 
          text = @sprintf("(%.2f%% relative)", slope*100),
          fontsize = 8, color = :gray30)
    
    axislegend(ax_d, position = :lt, labelsize = 8)
    
    # Layout
    colgap!(fig.layout, 10)
    rowgap!(fig.layout, 10)
    
    save("error_accumulation_analysis.pdf", fig, pt_per_unit=1)
    save("error_accumulation_analysis.png", fig, px_per_unit=2)
    
    # Diagnostics
    println("\n=== Error Analysis ===")
    println(@sprintf("Final absolute error: %.4f", errors_C[end]))
    println(@sprintf("Final relative error: %.3f%%", rel_errors_C[end] * 100))
    println(@sprintf("Error scaling: ε = %.2e × ⟨C⟩", slope))
    println(@sprintf("Mean error rate (early): %.2e s⁻¹", mean(error_rate[1:5])))
    println(@sprintf("Mean error rate (late):  %.2e s⁻¹", mean(error_rate[end-5:end])))
    println(@sprintf("State space size (final): %d", size_𝒮ₜ[idxs[end]]))
    
    fig
end

# ╔═╡ 3750408a-40e0-4ff4-b199-605886e3f0c2
begin
	local δt = 1e4
	local T =0:δt:Integer(1e5)
	

end

# ╔═╡ 83cf9310-51dc-458f-a82b-560d891388b4


# ╔═╡ 40a18845-9662-4028-ab7d-3480a9d5e9b1
idxs

# ╔═╡ 01a13951-2701-48bc-bb48-6c9594c89888
# ╠═╡ disabled = true
#=╠═╡
begin
	using LsqFit
	siam_theme = Theme(
	    fontsize = 10,
	    font = "CMU Serif",
	    
	    Axis = (
	        titlesize = 10,
	        xlabelsize = 10,
	        ylabelsize = 10,
	        xticklabelsize = 9,
	        yticklabelsize = 9,
	        
	        xgridvisible = false,
	        ygridvisible = false,
	        
	        leftspinevisible = true,
	        rightspinevisible = true,
	        topspinevisible = true,
	        bottomspinevisible = true,
	        
	        spinewidth = 0.75,
	        xtickwidth = 0.75,
	        ytickwidth = 0.75,
	        xtickalign = 1,
	        ytickalign = 1,
	        xticksize = 3,
	        yticksize = 3,
	    ),
	    
	    Legend = (
	        framevisible = false,
	        patchsize = (20, 3),
	        labelsize = 9,
	        rowgap = 2,
	        padding = (2, 2, 2, 2),
	    ),
	    
	    Lines = (
	        linewidth = 1.0,
	    ),
	    
	    backgroundcolor = :white,
	)
	
	set_theme!(siam_theme)
	
	# Create figure with 4 panels
	begin
	    local errors_C = abs.(afsp_mean[:, 3] .- exact_mean[:, 3])
	    Δt = 1.0
	    
	    local fig = Figure(size=(504, 380))  # SIAM double-column width
	    
	    # Panel (a): Error accumulation over time
	    local ax_a = Axis(fig[1, 1],
	        xlabel = "Time (s)",
	        ylabel = "Absolute error",
	        title = "(a) Mean error accumulation")
	    
	    lines!(ax_a, T, errors_C, 
	        color = :black, linewidth = 1.0)
	    
	    # Mark phase transition
	    vlines!(ax_a, [5000], 
	        color = :black, linestyle = :dash, linewidth = 0.75)
	    
	    text!(ax_a, 4200, 0.18, 
	        text = "Phase transition", 
	        fontsize = 8, rotation = π/2)
	    
	    # Panel (b): Error rate over time
	    local ax_b = Axis(fig[1, 2],
	        xlabel = "Time (s)",
	        ylabel = "Error rate (s⁻¹)",
	        title = "(b) Instantaneous error rate")
	    
	    error_rate = diff(errors_C) ./ Δt
	    # Smooth for visualization
	    window = 50
	    error_rate_smooth = [mean(error_rate[max(1,i-window):min(end,i+window)]) 
	                         for i in 1:length(error_rate)]
	    
	    lines!(ax_b, T[1:end-1], error_rate_smooth, 
	        color = :black, linewidth = 1.0)
	    
	    vlines!(ax_b, [5000], 
	        color = :black, linestyle = :dash, linewidth = 0.75)
	    
	    # Add phase rate annotations
	    text!(ax_b, 2500, 4.5e-5, 
	        text = "Phase 1\n3.8×10⁻⁶ s⁻¹", 
	        fontsize = 8, align = (:center, :top))
	    
	    text!(ax_b, 7500, 4.5e-5, 
	        text = "Phase 2\n4.2×10⁻⁵ s⁻¹", 
	        fontsize = 8, align = (:center, :top))
	    
	    # Panel (c): Mean trajectories comparison
	    local ax_c = Axis(fig[2, 1],
	        xlabel = "Time (s)",
	        ylabel = "Mean copy number",
	        title = "(c) Mean trajectories")
	    
	    # Species A
	    lines!(ax_c, T, exact_mean[:, 1], 
	        color = :black, linewidth = 1.0, label = "A (exact)")
	    lines!(ax_c, T, afsp_mean[:, 1], 
	        color = :black, linewidth = 1.0, linestyle = :dash, label = "A (adaptive)")
	    
	    # Species C
	    lines!(ax_c, T, exact_mean[:, 3], 
	        color = :gray40, linewidth = 1.0, label = "C (exact)")
	    lines!(ax_c, T, afsp_mean[:, 3], 
	        color = :gray40, linewidth = 1.0, linestyle = :dash, label = "C (adaptive)")
	    
	    vlines!(ax_c, [5000], 
	        color = :black, linestyle = :dot, linewidth = 0.75)
	    
	    axislegend(ax_c, position = :lt, nbanks = 2, labelsize = 8)
	    
	    # Panel (d): Error vs mean C (showing moment sensitivity)
	    local ax_d = Axis(fig[2, 2],
	        xlabel = "Mean C copy number",
	        ylabel = "Absolute error",
	        title = "(d) Error vs. mean C")
	    
	    # Color by time phase
	    phase1_mask = T .< 5000
	    phase2_mask = T .>= 5000
	    
	    scatter!(ax_d, exact_mean[phase1_mask, 3], errors_C[phase1_mask],
	        color = :black, marker = :circle, markersize = 2, 
	        label = "Phase 1 (t < 5000s)")
	    
	    scatter!(ax_d, exact_mean[phase2_mask, 3], errors_C[phase2_mask],
	        color = :gray40, marker = :xcross, markersize = 3,
	        label = "Phase 2 (t ≥ 5000s)")
	    
	    # Fit linear relationship
	    lmodel(C, p) = p[1] .* C .+ p[2]
	    fit = LsqFit.curve_fit(lmodel, exact_mean[1000:end, 3], errors_C[1000:end], [0.05, 0.0])
	    C_range = range(0, maximum(exact_mean[:, 3]), length=100)
	    
	    lines!(ax_d, C_range, lmodel(C_range, coef(fit)),
	        color = :black, linewidth = 1.0, linestyle = :dash,
	        label = "Linear fit")
	    
	    axislegend(ax_d, position = :lt, labelsize = 8)
	    
	    # Adjust layout
	    colgap!(fig.layout, 10)
	    rowgap!(fig.layout, 10)
	    
	    # Save figure
	    save("error_accumulation_analysis.eps", fig, pt_per_unit=1)
	    save("error_accumulation_analysis.png", fig, px_per_unit=2)
	    
	    fig
	end
end
  ╠═╡ =#

# ╔═╡ 035c4b08-7be6-4b89-b6cc-49374ee36418
# Critical insight: Check per-step error vs theoretical bound
begin
    local δt = 1e4
    local T = 0:δt:Integer(1e5)
    
    # Match timepoints
    local idxs = [findfirst(T1 .== t) for t in T]
    local errors_C = abs.(afsp_mean[idxs, 3] .- exact_mean[:, 3])
    local Δt = 1.0
    N = length(T) - 1
    
    # Average per-step error
    avg_per_step_obs = errors_C[end] / N
    
    # Theoretical bound: ε_n ≤ ε_flux for each step
    local ε_flux = 1e-4
    theoretical_bound = N * ε_flux
    
    println("\n=== GLOBAL ERROR COMPARISON ===")
    println("Number of steps:       $(N)")
    println("Time step size:        $(Δt) s")
    println("Flux tolerance:        $(ε_flux)")
    println("")
    println("Final observed error:  $(round(errors_C[end], sigdigits=5))")
    println("Theoretical bound:     $(round(theoretical_bound, sigdigits=5))")
    println("Observed/Bound ratio:  $(round(errors_C[end]/theoretical_bound, digits=4))")
    println("")
    println("Avg per-step observed: $(round(avg_per_step_obs, sigdigits=5))")
    println("Avg per-step theory:   $(ε_flux)")
    println("Ratio (obs/theory):    $(round(avg_per_step_obs/ε_flux, digits=4))")
    
    # This tells us if theory is conservative or tight
    if errors_C[end] < theoretical_bound
        println("\n✓ Theory is CONSERVATIVE (observed < bound)")
        println("  Bound is loose by factor: $(round(theoretical_bound/errors_C[end], digits=2))")
    else
        println("\n✗ Theory is VIOLATED (observed > bound)")
        println("  Theory underestimates by: $(round(errors_C[end]/theoretical_bound, digits=2))x")
    end
end

# ╔═╡ 32e6fa3b-dc9f-4cca-bfb1-380b4dddbf3c
function mean_sol(sol)
	sol_mean = map(1:length(sol)) do i
		    sum(collect.(Tuple.(sol[i][1])) .* sol[i][2])
	end 
	fsp_mean=hcat(sol_mean...)'
end


# ╔═╡ b9cf4dc2-9817-493b-940b-5e3677ceb92e
"""
Extracts the final marginal probability distribution for a single species.
This function creates a histogram by aggregating the probabilities of all states
that share the same population count for the specified species.
"""
function get_final_marginal_distribution(solution_data, species_index)
    states, probs = solution_data[end]
    dist_dict = Dict{Int, Float64}()
    
    for (state, prob) in zip(states, probs)
        count = state.I[species_index]
        dist_dict[count] = get(dist_dict, count, 0.0) + prob
    end
    
    counts = sort(collect(keys(dist_dict)))
    probabilities = [dist_dict[c] for c in counts]

    return counts, abs.(probabilities)
end

# ╔═╡ 690d5609-6bd0-4caa-b88e-a28e3b9da19c
sol[10]

# ╔═╡ e2e71c4e-4782-46bb-8560-9ffbae1145b3
begin
    # --- 1. Process Data ---
    num_species = 3 # A, B, C
    N_initial = 1   # initial population of A

    iters = 1:length(T)

    # Time vectors
    sol_t_correct = collect(T[iters])
    sol_t_failed  = collect(T[iters])  

    # Mean trajectories (ensure shape: time × species)
    μ_ok_raw   = mean_sol(sol[iters])
    μ_fail_raw = mean_sol(sol1[iters])

    μ_ok   = Matrix(μ_ok_raw)
    μ_fail = Matrix(μ_fail_raw)

    # If means are species × time, transpose to time × species
    if size(μ_ok, 1) == num_species && size(μ_ok, 2) == length(sol_t_correct)
        μ_ok = μ_ok'
    end
    if size(μ_fail, 1) == num_species && size(μ_fail, 2) == length(sol_t_failed)
        μ_fail = μ_fail'
    end
    @assert size(μ_ok,   1) == length(sol_t_correct)  "mean_sol(sol) must align with time"
    @assert size(μ_fail, 1) == length(sol_t_failed)   "mean_sol(sol1) must align with time"
    @assert size(μ_ok,   2) == num_species && size(μ_fail, 2) == num_species

    # State space sizes
    size_correct = size_𝒮ₜ[iters]
    size_failed  = [length(s[1]) for s in sol1[iters]]

    local_ϵ_ok   = local_ϵ[iters]
    local_ϵ_fail = local_ϵ1[iters]

    # Final marginals
    counts_A_correct, probs_A_correct = get_final_marginal_distribution(sol[iters], 1)
    counts_B_correct, probs_B_correct = get_final_marginal_distribution(sol[iters], 2)
    counts_C_correct, probs_C_correct = get_final_marginal_distribution(sol[iters], 3)

    counts_A_failed, probs_A_failed = get_final_marginal_distribution(sol1[iters], 1)
    counts_B_failed, probs_B_failed = get_final_marginal_distribution(sol1[iters], 2)
    counts_C_failed, probs_C_failed = get_final_marginal_distribution(sol1[iters], 3)

    # --- 2. Create Figure and Layout (SIAM style: smaller, grayscale) ---
    f = Figure(size = (1000, 650), fontsize = 11, fonts = (;regular = "Latin Modern Roman"))

    # Row titles (smaller, no bold for SIAM style)
    Label(f[1, 1:3], "Flux-Aware Pruning",
          fontsize = 12, tellwidth = false, padding = (0,0,15,0))
    Label(f[3, 1:3], "Probability-Only Pruning",
          fontsize = 12, tellwidth = false, padding = (0,0,15,0))

    # --- Row 1: Correct Simulation ---
    dist_grid_correct = f[2, 1] = GridLayout()
    ax_dist_A_correct = Axis(dist_grid_correct[1, 1], ylabel = "Probability", 
                             xticks = [0, 1], title = L"N_A")
    ax_dist_B_correct = Axis(dist_grid_correct[1, 2], xticks = [0, 1], title = L"N_B")
    
    # Fix for C: compute sensible xticks
    c_min = minimum(counts_C_correct)
    c_max = maximum(counts_C_correct)
    c_range = c_max - c_min
    # Choose 3-5 evenly spaced ticks
    if c_range <= 10
        c_ticks = [c_min, div(c_min + c_max, 2), c_max]
    else
        step = max(1, div(c_range, 4))
        c_ticks = collect(c_min:step:c_max)
    end
    ax_dist_C_correct = Axis(dist_grid_correct[1, 3], xticks = [c_max], title = L"N_C")

    ax_means_correct = Axis(f[2, 2], title = "Mean Trajectories", xscale = log10,
                            xgridvisible = false, ygridvisible = false)
    ax_size_correct  = Axis(f[2, 3], ylabel = "State Space Size", xscale = log10,
                            xgridvisible = false, ygridvisible = false)


    # --- Row 2: Failed Simulation ---
    dist_grid_failed = f[4, 1] = GridLayout()
    ax_dist_A_failed = Axis(dist_grid_failed[1, 1], ylabel = "Probability",
                            xlabel = L"N_A", xticks = [0, 1])
    ax_dist_B_failed = Axis(dist_grid_failed[1, 2], xlabel = L"N_B", xticks = [0, 1])
    ax_dist_C_failed = Axis(dist_grid_failed[1, 3], xlabel = L"N_C", xticks = [0, 1])

    ax_means_failed = Axis(f[4, 2], title = "Mean Trajectories",
                           xscale = log10, xlabel = "Time (s)",
                           xgridvisible = false, ygridvisible = false)
    ax_size_failed  = Axis(f[4, 3], title = "State Space Size",
                           xscale = log10, xlabel = "Time (s)",
                           xgridvisible = false, ygridvisible = false)
  

    # Link axes
    linkyaxes!(ax_means_correct, ax_means_failed)
    linkyaxes!(ax_size_correct,  ax_size_failed)
    linkyaxes!(ax_dist_A_correct, ax_dist_B_correct,
               ax_dist_A_failed, ax_dist_B_failed, ax_dist_C_failed)

    # --- 3. Plotting (SIAM style: grayscale, distinct line styles) ---
    bar_color = (:black)
    linestyles = [:solid, :dash, :dot]
    linewidths = [1.5, 1.5, 2.0]  # Slightly thicker dot for visibility
    species_labels = [L"A", L"B", L"C/\max(C)"]

    # Final distributions (correct)
    barplot!(ax_dist_A_correct, counts_A_correct, probs_A_correct, color = bar_color)
    barplot!(ax_dist_B_correct, counts_B_correct, probs_B_correct, color = bar_color)
    barplot!(ax_dist_C_correct, counts_C_correct, probs_C_correct, color = bar_color)

    # Final distributions (failed)
    barplot!(ax_dist_A_failed, counts_A_failed, probs_A_failed, color = bar_color)
    barplot!(ax_dist_B_failed, counts_B_failed, probs_B_failed, color = bar_color)
    barplot!(ax_dist_C_failed, counts_C_failed, probs_C_failed, color = bar_color)

    # Refined x-limits
    xlims!(ax_dist_A_correct, -0.5, 1.5)
    xlims!(ax_dist_B_correct, -0.5, 1.5)
    xlims!(ax_dist_A_failed,  -0.5, 1.5)
    xlims!(ax_dist_B_failed,  -0.5, 1.5)
    xlims!(ax_dist_C_failed,  -0.5, 1.5)
    xlims!(ax_dist_C_correct, c_min - 0.5, c_max + 0.5)

    # Mean trajectories (correct)
    valid_indices_correct = findall(t -> t > 0, sol_t_correct)
    # Plot species 1..(num_species-1) as-is
    for i in 1:(num_species - 1)
        lines!(ax_means_correct,
               sol_t_correct[valid_indices_correct],
               μ_ok[valid_indices_correct, i];
               color = :black, linestyle = linestyles[i], linewidth = linewidths[i],
               label = species_labels[i])
    end
    # Normalize last species for visibility
    denom = maximum(μ_ok[valid_indices_correct, end])
    last_series = denom > 0 ? μ_ok[valid_indices_correct, end] ./ denom :
                              zero.(μ_ok[valid_indices_correct, end])
    lines!(ax_means_correct,
           sol_t_correct[valid_indices_correct], last_series;
           color = :black, linestyle = linestyles[end], linewidth = linewidths[end],
           label = species_labels[end])
    axislegend(ax_means_correct, "Species", position = :lt, framevisible = false)

    # Mean trajectories (failed)
    valid_indices_failed = findall(t -> t > 0, sol_t_failed)
    for i in 1:num_species
        lines!(ax_means_failed,
               sol_t_failed[valid_indices_failed],
               μ_fail[valid_indices_failed, i];
               color = :black, linestyle = linestyles[i], linewidth = linewidths[i],
               label = species_labels[i])
    end
    axislegend(ax_means_failed, "Species", position = :rt, framevisible = false)

    # State space sizes + leakage (use gray for secondary axis)
    lines!(ax_size_correct,
           sol_t_correct[valid_indices_correct],
           size_correct[valid_indices_correct]; color = :black, linewidth = 1.5)
    

    lines!(ax_size_failed,
           sol_t_failed[valid_indices_failed],
           size_failed[valid_indices_failed]; color = :black, linewidth = 1.5)


    # --- 4. Final Layout Adjustments ---
    colgap!(f.layout, 20)
    rowgap!(f.layout, 8)
    rowgap!(f.layout, 2, 25) # extra space between the two main rows

    f
end

# ╔═╡ 838fb0ef-03b6-4501-b66c-07859355e5e2
begin
	local fig=Figure()
	a = Axis(fig[1, 1])
	lines!(a, μ_ok_raw[:, 1], label = "accurate value")
	lines!(a,Axis(fig[1, 1]), μ_ok_raw[:, 2], label = "pade approximation")
	lines!(a,Axis(fig[1, 1]), μ_ok_raw[:, 3], label = "taylor approximation")
	lines!(a, Axis(fig[1, 1]), exact_mean[:, 3], label = "empty val")
	axislegend(a)
	fig
end

# ╔═╡ d0ad5800-4b37-4722-b884-07af915a3d7f
save("bn.png", f)

# ╔═╡ c1a92cd0-ffed-43d7-81c5-021799578d0d
begin
	states1 = collect(sol[100][1])
	V1 = sol[100][2]
	P1 = zeros(600, 600)
	P1v = zeros(600, 600)
	P1[states1] .= 1;
	P1v[states1] .= V1;
	     
	
	states2 = sol[200][1]
	V2 = sol[100][2]
	P2 = zeros(600, 600)
	P2v = zeros(600, 600)
	P2[states2] .= 1;
	P2v[states2] .= V2;
	     
	
	statesm = states1 ∪ states2
	Pm = zeros(600, 600)
	Pmv = zeros(600, 600)
	Pm[statesm] .= 1;
	Pmv[states2] .= V2;
end

# ╔═╡ 5e389838-a3cc-46d2-986d-c407edc8e93a
begin
	
	# Initialize figure
	fig = Figure(size=(1000, 550))
	
	# Data arrays
	top = [P1, Pm, P2]
	bot = [P1v, Pmv, P2v]
	local T = LinRange(0, 30, length(sol))
	
	# Create 2x3 grid of axes
	axes = []
	for i in 1:2, j in 1:3
	    ax = Axis(fig[i, j])
	    push!(axes, ax)
	end
	
	# Plot contours
	for (i, idx) in enumerate(points)
	    # Top row - binary contours
	    ax_top = axes[i]
	    contour!(ax_top, top[i], color=:black, linewidth=2.5, levels=[1])
	    
	    # Bottom row - log contours  
	    ax_bot = axes[i + 3]
	    contour!(ax_bot, bot[i], color=:black, linewidth=2.5, levels=4)
	    
	    # Add colorbar only for rightmost column
	    if i == 3
	        Colorbar(fig[1, 4], colormap=:viridis, label=L"1 if active, 0 otherwise")
	        Colorbar(fig[2, 4], colormap=:viridis, label=L"\operatorname{log}_{10}(p_S(t))")
	    end
	end
	
	# Set axis properties
	for (i, ax) in enumerate(axes)
	    xlims!(ax, 400, 475)
	    ylims!(ax, 60, 135)
	    
	    # X labels for bottom row
	    if i > 3
	        ax.xlabel = L"|U|"
	    end
	    
	    # Y labels for left column
	    if i == 1 || i == 4
	        ax.ylabel = L"|V|"
	    end
	    
	    # Set ticks
	    ax.xticks = [150, 350, 550]
	end
	
	# Set font properties
	set_theme!(fontsize=12, fontfamily="Computer Modern")
	
	display(fig)
	save("toggle_switch_grid.pdf", fig)
end

# ╔═╡ 9867082b-f32b-4c92-adbf-866efeec8ced
lines(flux1)

# ╔═╡ 9e8a68b4-630f-4051-a764-dc44ccc8847e
# ╠═╡ disabled = true
#=╠═╡
begin
    # --- 1. Process Data ---
    num_species = 3 # A, B, C
    N_initial = 1 # The initial population of species A

    # Ensure time vectors are standard Arrays

	iters = 1:length(T)

    sol_t_correct =  collect(T[iters])
    sol_t_failed = collect(T[iters]) # Assuming a separate time vector for the failed sim

    # Calculate mean trajectories
    fsp_mean_correct = mean_sol(sol[iters])
    fsp_mean_failed = mean_sol(sol1[iters])

    # Calculate state space size evolution
    size_correct = size_𝒮ₜ[iters]
    size_failed = [length(s[1]) for s in sol1[iters]]
    
    # Extract final marginal distributions
    counts_A_correct, probs_A_correct = get_final_marginal_distribution(sol[iters], 1)
    counts_B_correct, probs_B_correct = get_final_marginal_distribution(sol[iters], 2)
    counts_C_correct, probs_C_correct = get_final_marginal_distribution(sol[iters], 3)

    counts_A_failed, probs_A_failed = get_final_marginal_distribution(sol1[iters], 1)
    counts_B_failed, probs_B_failed = get_final_marginal_distribution(sol1[iters], 2)
    counts_C_failed, probs_C_failed = get_final_marginal_distribution(sol1[iters], 3)

    # --- 2. Create Figure and Layout ---
    f = Figure(size = (1200, 800), fontsize=14)
    
    # Add row titles
    Label(f[1, 1:3], "Adaptive FSP With Flux Filter", fontsize=18, font=:bold, tellwidth=false, padding=(0,0,20,0))
    Label(f[3, 1:3], "No Flux Filter", fontsize=18, font=:bold, tellwidth=false, padding=(0,0,20,0))

    # --- Row 1: Correct Simulation ---
    dist_grid_correct = f[2, 1] = GridLayout()
    ax_dist_A_correct = Axis(dist_grid_correct[1, 1], ylabel="Probability", xticks=[0, 1])
    ax_dist_B_correct = Axis(dist_grid_correct[1, 2], xticks=[0, 1])
    ax_dist_C_correct = Axis(dist_grid_correct[1, 3], xticks=[minimum(counts_C_correct), maximum(counts_C_correct)]) 
    
    ax_means_correct = Axis(f[2, 2], title="Mean Trajectories", xscale=log10)
    ax_size_correct = Axis(f[2, 3], ylabel="State Space Size", xscale=log10)
	ax_flux_correct = Axis(f[2, 3], ylabel="Ergodcity after Pruning",  yaxisposition = :right, xscale=log10)

    # --- Row 2: Failed Simulation ---
    dist_grid_failed = f[4, 1] = GridLayout()
    ax_dist_A_failed = Axis(dist_grid_failed[1, 1], ylabel="Probability", xlabel=L"N_A", xticks=[0, 1])
    ax_dist_B_failed = Axis(dist_grid_failed[1, 2], xlabel=L"N_B", xticks=[0, 1])
    ax_dist_C_failed = Axis(dist_grid_failed[1, 3], xlabel=L"N_C", xticks=[0, 1])

    ax_means_failed = Axis(f[4, 2], title="Mean Trajectories", xscale=log10, xlabel="Time, t")
    ax_size_failed = Axis(f[4, 3], title="State Space Size", xscale=log10, xlabel="Time, t")
	ax_flux_failed = Axis(f[4, 3], ylabel="Ergodcity after Pruning",  yaxisposition = :right, xscale=log10)

    # Link axes for direct comparison
    linkyaxes!(ax_means_correct, ax_means_failed)
    linkyaxes!(ax_size_correct, ax_size_failed)
    linkyaxes!(ax_dist_A_correct, ax_dist_B_correct, ax_dist_A_failed, ax_dist_B_failed, ax_dist_C_failed)

    # --- 3. Plotting ---
    bar_color = (:black, 0.8)
    linestyles = [:solid, :dash, :dot]
    species_labels = ["A", "B", "C"]
    
    # Plot Correct Distributions
    barplot!(ax_dist_A_correct, counts_A_correct, probs_A_correct, color=bar_color)
    barplot!(ax_dist_B_correct, counts_B_correct, probs_B_correct, color=bar_color)
    barplot!(ax_dist_C_correct, counts_C_correct, probs_C_correct, color=bar_color)

    # Plot Failed Distributions
    barplot!(ax_dist_A_failed, counts_A_failed, probs_A_failed, color=bar_color)
    barplot!(ax_dist_B_failed, counts_B_failed, probs_B_failed, color=bar_color)
    barplot!(ax_dist_C_failed, counts_C_failed, probs_C_failed, color=bar_color)

    # Set refined x-limits to focus on the single tick
    xlims!(ax_dist_A_correct, 0, 1)
    xlims!(ax_dist_B_correct, 0, 1)
    #xlims!(ax_dist_C_correct, N_initial - 50, N_initial + 50)
    xlims!(ax_dist_A_failed, 0, 1)
    xlims!(ax_dist_B_failed, 0, 1)
    xlims!(ax_dist_C_failed, 0, 1)

    # Plot Correct Mean Trajectories
    valid_indices_correct = findall(t -> t > 0, sol_t_correct)
    for i in 1:num_species-1
        lines!(ax_means_correct, sol_t_correct[valid_indices_correct], fsp_mean_correct[valid_indices_correct, i], color=:black, linestyle=linestyles[i], label=species_labels[i])
    end
	lines!(ax_means_correct, sol_t_correct[valid_indices_correct], fsp_mean_correct[valid_indices_correct, end] ./ maximum(fsp_mean_correct[valid_indices_correct, end]), color=:black, linestyle=linestyles[end], label=species_labels[end])
    axislegend(ax_means_correct, "Species", position=:cc)

    # Plot Failed Mean Trajectories
    valid_indices_failed = findall(t -> t > 0, sol_t_failed)
    for i in 1:num_species
        lines!(ax_means_failed, sol_t_failed[valid_indices_failed], fsp_mean_failed[valid_indices_failed, i], color=:black, linestyle=linestyles[i], label=species_labels[i])
    end
    axislegend(ax_means_failed, "Species", position=:rc)

    # Plot State Space Sizes
    lines!(ax_size_correct, sol_t_correct[valid_indices_correct], size_correct[valid_indices_correct], color=:black)
	lines!(ax_flux_correct, sol_t_correct[valid_indices_correct], local_ϵ[valid_indices_correct],  linestyle=:dash,color=:black)
    lines!(ax_size_failed, sol_t_failed[valid_indices_failed], size_failed[valid_indices_failed], color=:black)
	lines!(ax_flux_failed, sol_t_failed[valid_indices_failed], local_ϵ1[valid_indices_failed], color=:black, linestyle=:dash)
    
    # --- 4. Final Layout Adjustments ---
    colgap!(f.layout, 25)
    rowgap!(f.layout, 10)
    rowgap!(f.layout, 2, 30) # Add extra space between the two main rows

    f
end

  ╠═╡ =#

# ╔═╡ 8d9d0907-dc0d-4ad2-9faf-fbd015eb7234
full_distribution_comparison(length(sol), length(sol_exact))

# ╔═╡ 6f4ca781-97e6-4c70-a200-ea230eebd449
plot_full_comparison(length(sol), length(sol_exact))

# ╔═╡ afe2ea2f-b3e8-49f9-a42a-bbfc60d53c04
plot_joint_comparison(length(sol), length(sol_exact), 2, 3; 
                     species_name1="B", species_name2="C")

# ╔═╡ e7b739ac-a39e-494f-8ba1-4e554f9ac06c


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Arpack = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
Catalyst = "479239e8-5488-4da2-87a7-35f2df7eef83"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
Expokit = "a1e7a1ef-7a5d-5822-a38c-be74e1bb89f4"
ExponentialUtilities = "d4d017d3-3776-5f7e-afef-a10c40355c18"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
JLD = "4138dd39-2aa7-5051-a626-17a0bb65d9c8"
JumpProcesses = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
PROPACK = "b169e327-5944-5131-97a6-5d3d3f0a476a"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
Arpack = "~0.5.4"
CairoMakie = "~0.13.10"
Catalyst = "~15.0.8"
DifferentialEquations = "~7.16.1"
Expokit = "~0.2.0"
ExponentialUtilities = "~1.27.0"
Interpolations = "~0.15.1"
JLD = "~0.13.5"
JumpProcesses = "~9.16.1"
LsqFit = "~0.15.1"
PROPACK = "~0.5.0"
ProgressLogging = "~0.1.5"
Revise = "~3.8.1"
StatsBase = "~0.34.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.0"
manifest_format = "2.0"
project_hash = "fd48092b9ab7d6f9dfbe9b356a8cb3f417041fa8"

[[deps.ADTypes]]
git-tree-sha1 = "27cecae79e5cc9935255f90c53bb831cc3c870d7"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.18.0"
weakdeps = ["ChainRulesCore", "ConstructionBase", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.AlmostBlockDiagonals]]
deps = ["ConcreteStructs"]
git-tree-sha1 = "743abe5e5fe8cff96dad4123f263c0d8eee281c0"
uuid = "a95523ee-d6da-40b5-98cc-27bc505739d5"
version = "0.1.10"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "d81ae5489e13bc03567d4fbbb06c546a5e53c857"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.22.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "355ab2d61069927d4247cd69ad0e1f140b31e30d"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.12.0"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "SIMD", "TranscodingStreams"]
git-tree-sha1 = "a8f503e8e1a5f583fbef15a8440c8c7e32185df2"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "4826c9fe6023a87029e54870ad1a9800c7ea6623"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.10.1"

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"
    CliqueTreesExt = "CliqueTrees"

    [deps.BandedMatrices.weakdeps]
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "bca794632b8a9bbe159d56bf9e31c422671b35e0"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.3.2"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.BlockArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "f3132231bfb207e2793aa4f5f933f39c38d2cd60"
uuid = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
version = "1.9.0"
weakdeps = ["Adapt", "BandedMatrices"]

    [deps.BlockArrays.extensions]
    BlockArraysAdaptExt = "Adapt"
    BlockArraysBandedMatricesExt = "BandedMatrices"

[[deps.Blosc]]
deps = ["Blosc_jll"]
git-tree-sha1 = "310b77648d38c223d947ff3f50f511d08690b8d5"
uuid = "a74b3585-a348-5f62-a45c-50e91977d574"
version = "0.7.3"

[[deps.Blosc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Lz4_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "535c80f1c0847a4c967ea945fca21becc9de1522"
uuid = "0b7ba130-8d10-5ba8-a3d6-c5182647fed9"
version = "1.21.7+0"

[[deps.BoundaryValueDiffEq]]
deps = ["ADTypes", "BoundaryValueDiffEqAscher", "BoundaryValueDiffEqCore", "BoundaryValueDiffEqFIRK", "BoundaryValueDiffEqMIRK", "BoundaryValueDiffEqMIRKN", "BoundaryValueDiffEqShooting", "DiffEqBase", "FastClosures", "ForwardDiff", "LinearAlgebra", "Reexport", "SciMLBase"]
git-tree-sha1 = "d6ec33e4516b2e790a64128afdb54f3b536667a7"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "5.18.0"

    [deps.BoundaryValueDiffEq.extensions]
    BoundaryValueDiffEqODEInterfaceExt = "ODEInterface"

    [deps.BoundaryValueDiffEq.weakdeps]
    ODEInterface = "54ca160b-1b9f-5127-a996-1867f4bc2a2c"

[[deps.BoundaryValueDiffEqAscher]]
deps = ["ADTypes", "AlmostBlockDiagonals", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield"]
git-tree-sha1 = "47c833c459738a3f27c5b458ecf7832a4731ef4d"
uuid = "7227322d-7511-4e07-9247-ad6ff830280e"
version = "1.8.0"

[[deps.BoundaryValueDiffEqCore]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "ForwardDiff", "LineSearch", "LinearAlgebra", "Logging", "NonlinearSolveFirstOrder", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseConnectivityTracer", "SparseMatrixColorings"]
git-tree-sha1 = "b7b4d8cc80f116eab2eb6124dba58ea7aef31b85"
uuid = "56b672f2-a5fe-4263-ab2d-da677488eb3a"
version = "1.11.1"

[[deps.BoundaryValueDiffEqFIRK]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "325e6981a414cfa5181218936c23f0e16dee8f08"
uuid = "85d9eb09-370e-4000-bb32-543851f73618"
version = "1.9.0"

[[deps.BoundaryValueDiffEqMIRK]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "da6ae5e564ad06ced4d7504929c58130558007dd"
uuid = "1a22d4ce-7765-49ea-b6f2-13c8438986a6"
version = "1.9.0"

[[deps.BoundaryValueDiffEqMIRKN]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "609c2d03ea024df0d475fee483b93cf0e87c29d6"
uuid = "9255f1d6-53bf-473e-b6bd-23f1ff009da4"
version = "1.8.0"

[[deps.BoundaryValueDiffEqShooting]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "ba9bd1f31b58bfd5e48a56da0a426bcbd3462546"
uuid = "ed55bfe0-3725-4db6-871e-a1dc9f42a757"
version = "1.9.0"

[[deps.BracketingNonlinearSolve]]
deps = ["CommonSolve", "ConcreteStructs", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "a9014924595b7a2c1dd14aac516e38fa10ada656"
uuid = "70df07ce-3d50-431d-a3e7-ca6ddb60ac1e"
version = "1.3.0"
weakdeps = ["ChainRulesCore", "ForwardDiff"]

    [deps.BracketingNonlinearSolve.extensions]
    BracketingNonlinearSolveChainRulesCoreExt = ["ChainRulesCore", "ForwardDiff"]
    BracketingNonlinearSolveForwardDiffExt = "ForwardDiff"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Preferences", "Static"]
git-tree-sha1 = "f3a21d7fc84ba618a779d1ed2fcca2e682865bab"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.7"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "9bd45574379e50579a78774334f4a1f1238c0af5"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.13.10"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.Catalyst]]
deps = ["Combinatorics", "DataStructures", "DiffEqBase", "DocStringExtensions", "DynamicPolynomials", "DynamicQuantities", "EnumX", "Graphs", "JumpProcesses", "LaTeXStrings", "Latexify", "LinearAlgebra", "MacroTools", "ModelingToolkit", "Parameters", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SymbolicUtils", "Symbolics", "Unitful"]
git-tree-sha1 = "83a85091233ac1d45ae7cbec0c01eb1968b5e5f7"
uuid = "479239e8-5488-4da2-87a7-35f2df7eef83"
version = "15.0.8"

    [deps.Catalyst.extensions]
    CatalystBifurcationKitExtension = "BifurcationKit"
    CatalystCairoMakieExtension = "CairoMakie"
    CatalystGraphMakieExtension = ["GraphMakie", "NetworkLayout", "Makie"]
    CatalystHomotopyContinuationExtension = "HomotopyContinuation"
    CatalystStructuralIdentifiabilityExtension = "StructuralIdentifiability"

    [deps.Catalyst.weakdeps]
    BifurcationKit = "0f109fa4-8a5d-4b75-95aa-f515264e7665"
    CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
    GraphMakie = "1ecd5474-83a3-4783-bb4f-06765db800d2"
    HomotopyContinuation = "f213a82b-91d6-5c5d-acf7-10f1c761b327"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    NetworkLayout = "46757867-2c16-5918-afeb-47bfcb05e46a"
    StructuralIdentifiability = "220ca800-aa68-49bb-acd8-6037fa93a544"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "980f01d6d3283b3dbdfd7ed89405f96b7256ad57"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "2.0.1"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonMark]]
deps = ["PrecompileTools"]
git-tree-sha1 = "351d6f4eaf273b753001b2de4dffb8279b100769"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.9.1"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.Compiler]]
git-tree-sha1 = "382d79bfe72a406294faca39ef0c3cef6e6ce1f1"
uuid = "807dbc54-b67e-4c79-8afb-eafe4df6f2e1"
version = "0.1.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "5620ff4ee0084a6ab7097a27ba0c19290200b037"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.4"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "FastBroadcast", "ForwardDiff", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqRosenbrock", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack", "SymbolicIndexingInterface"]
git-tree-sha1 = "38a6a29851deae145184e2010eb683d56988fb93"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.59.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "FastPower", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "1b1e070e57681d1176d99a5ec455717e24686612"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.183.2"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseForwardDiffExt = ["ForwardDiff"]
    DiffEqBaseGTPSAExt = "GTPSA"
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseMooncakeExt = "Mooncake"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseSparseArraysExt = "SparseArrays"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffEqCallbacks]]
deps = ["ConcreteStructs", "DataStructures", "DiffEqBase", "DifferentiationInterface", "LinearAlgebra", "Markdown", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "448d118c4a7763b041fc41cf91c85d23773f1dc2"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "4.10.1"

    [deps.DiffEqCallbacks.extensions]
    DiffEqCallbacksFunctorsExt = "Functors"

    [deps.DiffEqCallbacks.weakdeps]
    Functors = "d9f16b24-f501-4c13-a1f2-28368ffc5196"

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "516d553f5deee7c55b2945b5edf05b6542837887"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.24.1"

    [deps.DiffEqNoiseProcess.extensions]
    DiffEqNoiseProcessReverseDiffExt = "ReverseDiff"

    [deps.DiffEqNoiseProcess.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "JumpProcesses", "LinearAlgebra", "LinearSolve", "NonlinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "afdc7dfee475828b4f0286d63ffe66b97d7a3fa7"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.16.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "c8d85ecfcbaef899308706bebdd8b00107f3fb43"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.6.54"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = "GPUArraysCore"
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.DispatchDoctor]]
deps = ["MacroTools", "Preferences"]
git-tree-sha1 = "fc34127e78323c49984e1a146d577d0f890dd2b4"
uuid = "8d63f2c5-f18a-4cf2-ba9d-b3f60fc568c8"
version = "0.4.26"
weakdeps = ["ChainRulesCore", "EnzymeCore"]

    [deps.DispatchDoctor.extensions]
    DispatchDoctorChainRulesCoreExt = "ChainRulesCore"
    DispatchDoctorEnzymeCoreExt = "EnzymeCore"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3bc002af51045ca3b47d2e1787d6ce02e68b943a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.122"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "c249d86e97a7e8398ce2068dce4c078a1c3464de"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.16"
weakdeps = ["Makie", "Random"]

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"
    DomainSetsRandomExt = "Random"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "ca693f8707a77a0e365d49fe4622203b72b6cf1d"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.3"

[[deps.DynamicQuantities]]
deps = ["DispatchDoctor", "PrecompileTools", "TestItems", "Tricks"]
git-tree-sha1 = "57c48a46e27d67208ad51c564b078a90bbd0dc2c"
uuid = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
version = "1.10.0"

    [deps.DynamicQuantities.extensions]
    DynamicQuantitiesLinearAlgebraExt = "LinearAlgebra"
    DynamicQuantitiesMeasurementsExt = "Measurements"
    DynamicQuantitiesScientificTypesExt = "ScientificTypes"
    DynamicQuantitiesUnitfulExt = "Unitful"

    [deps.DynamicQuantities.weakdeps]
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    ScientificTypes = "321657f4-b219-11e9-178b-2701a2544e81"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.EnzymeCore]]
git-tree-sha1 = "f91e7cb4c17dae77c490b75328f22a226708557c"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.15"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.ExplicitImports]]
deps = ["Markdown", "Pkg", "PrecompileTools", "TOML"]
git-tree-sha1 = "fde76669757deacce495be6018d17ffe9d70f214"
uuid = "7d51a73a-1435-4ff3-83d9-f097790105c7"
version = "1.13.2"

[[deps.Expokit]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "b0313f5f1825aabb1adb186c837e0d339d02f351"
uuid = "a1e7a1ef-7a5d-5822-a38c-be74e1bb89f4"
version = "0.2.0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "cae251c76f353e32d32d76fae2fea655eab652af"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.27.0"
weakdeps = ["StaticArrays"]

    [deps.ExponentialUtilities.extensions]
    ExponentialUtilitiesStaticArraysExt = "StaticArrays"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "eaa040768ea663ca695d442be1bc97edfe6824f2"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "6.1.3+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FastAlmostBandedMatrices]]
deps = ["ArrayInterface", "ArrayLayouts", "BandedMatrices", "ConcreteStructs", "LazyArrays", "LinearAlgebra", "MatrixFactorizations", "PrecompileTools", "Reexport"]
git-tree-sha1 = "9482a2b4face8ade73792c23a54796c79ed1bcbf"
uuid = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
version = "0.1.5"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "ab1b34570bcdf272899062e1a56285a53ecaae08"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.3.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "0044e9f5e49a57e88205e8f30ab73928b05fe5b6"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "1.1.0"

[[deps.FastPower]]
git-tree-sha1 = "e47c70bf430175e077d1955d7f04923504acc74c"
uuid = "a4df4552-cc26-4903-aec0-212e50a0e84b"
version = "1.2.0"

    [deps.FastPower.extensions]
    FastPowerEnzymeExt = "Enzyme"
    FastPowerForwardDiffExt = "ForwardDiff"
    FastPowerMeasurementsExt = "Measurements"
    FastPowerMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    FastPowerMooncakeExt = "Mooncake"
    FastPowerReverseDiffExt = "ReverseDiff"
    FastPowerTrackerExt = "Tracker"

    [deps.FastPower.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "d60eb76f37d7e5a40cc2e7c36974d864b82dc802"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.1"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "173e4d8f14230a7523ae11b9a3fa9edb3e0efd78"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.14.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FindFirstFunctions]]
git-tree-sha1 = "544bdc2902fa966900d354510be775f8668a57bd"
uuid = "64ca27bc-2ba2-4a57-88aa-44e436879224"
version = "1.4.2"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "9340ca07ca27093ff68418b7558ca37b05f8aeb1"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.29.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "ba6ce081425d0afb2bedd00d9884464f764a9225"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.2.2"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "a694e2a57394e409f7a11ee0977362a9fafcb8c7"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.6"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "1f5a80f4ed9f5a4aada88fc2db456e637676414b"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.10"

    [deps.GeometryBasics.extensions]
    GeometryBasicsGeoInterfaceExt = "GeoInterface"

    [deps.GeometryBasics.weakdeps]
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "50c11ffab2a3d50192a228c313f05b5b5dc5acb2"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.0+0"

[[deps.Glob]]
git-tree-sha1 = "97285bbd5230dd766e9ef6749b80fc617126d496"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.1"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "7a98c6502f4632dbe9fb1973a4244eaa3324e84d"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.1"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "93d5c27c8de51687a2c70ec0716e6e76f298416f"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.2"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.H5Zblosc]]
deps = ["Blosc", "HDF5"]
git-tree-sha1 = "d778420e524bcf56066e8c63c7aa315ae7269da2"
uuid = "c8ec2601-a99c-407f-b158-e79c03c2f5f7"
version = "0.1.2"

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "MPIPreferences", "Mmap", "Preferences", "Printf", "Random", "Requires", "UUIDs"]
git-tree-sha1 = "e856eef26cf5bf2b0f95f8f4fc37553c72c8641c"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.17.2"

    [deps.HDF5.extensions]
    MPIExt = "MPI"

    [deps.HDF5.weakdeps]
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"

[[deps.HDF5_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "LibCURL_jll", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "OpenSSL_jll", "TOML", "Zlib_jll", "libaec_jll"]
git-tree-sha1 = "e94f84da9af7ce9c6be049e9067e511e17ff89ec"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.14.6+0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XML2_jll", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "3d468106a05408f9f7b6f161d9e7715159af247b"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.12.2+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "bf0210c01fb7d67c31fed97d7c1d1716b98ea689"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.1"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "5fbb102dcb8b1a858111ae81d56682376130517d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.11"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD]]
deps = ["Compat", "FileIO", "H5Zblosc", "HDF5", "Printf"]
git-tree-sha1 = "e42f32690d41f758e126a48ee43459ef91179d1f"
uuid = "4138dd39-2aa7-5051-a626-17a0bb65d9c8"
version = "0.13.5"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "06ea418d0c95878c8f3031023951edcf25b9e0ef"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.2.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4255f0032eafd6451d707a51d5f0248b8a165e4d"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.3+0"

[[deps.JuliaFormatter]]
deps = ["CommonMark", "Glob", "JuliaSyntax", "PrecompileTools", "TOML"]
git-tree-sha1 = "f512fefd5fdc7dd1ca05778f08f91e9e4c9fdc37"
uuid = "98e50ef6-434e-11e9-1051-2b60c6c9e899"
version = "2.1.6"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "277779adfedf4a30d66b64edc75dc6bb6d52a16e"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.10.6"

[[deps.JuliaSyntax]]
git-tree-sha1 = "937da4713526b96ac9a178e2035019d3b78ead4a"
uuid = "70703baa-626e-46a2-a12c-08ffd08c73b4"
version = "0.4.10"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "f8da88993c914357031daf0023f18748ff473924"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.16.1"
weakdeps = ["FastBroadcast"]

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "ba51324b894edaf1df3ab16e2cc6bc3280a2f1a7"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.10"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "d1fc961038207e43982851e57ee257adc37be5e8"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.10.2"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "79ee64f6ba0a5a49930f51c86f60d7526b5e12c8"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.8.0"

    [deps.LazyArrays.extensions]
    LazyArraysBandedMatricesExt = "BandedMatrices"
    LazyArraysBlockArraysExt = "BlockArrays"
    LazyArraysBlockBandedMatricesExt = "BlockBandedMatrices"
    LazyArraysStaticArraysExt = "StaticArrays"

    [deps.LazyArrays.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "56513a09b8e0ae6485f34401ea9e2f31357958ec"
uuid = "2d8b4e74-eb68-11e8-0fb9-d5eb67b50637"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.11.1+1"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3acf07f130a76f87c041cfb2ff7d7284ca67b072"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.2+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2a7a12fc0a4e7fb773450d17975322aa77142106"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.2+0"

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "97d502765cc5cf3a722120f50da03c2474efce04"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.4"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "4adee99b7262ad2a1a4bbbc59d993d24e55ea96f"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.4.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LinearOperators]]
deps = ["FastClosures", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "TimerOutputs"]
git-tree-sha1 = "db137007d2c4ed948aa5f2518a2b451851ea8bda"
uuid = "5c8ed15e-5a4c-59e4-a42b-c7e8811fb125"
version = "2.11.0"

    [deps.LinearOperators.extensions]
    LinearOperatorsAMDGPUExt = "AMDGPU"
    LinearOperatorsCUDAExt = "CUDA"
    LinearOperatorsChainRulesCoreExt = "ChainRulesCore"
    LinearOperatorsJLArraysExt = "JLArrays"
    LinearOperatorsLDLFactorizationsExt = "LDLFactorizations"
    LinearOperatorsMetalExt = "Metal"

    [deps.LinearOperators.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    JLArrays = "27aeb0d3-9eb9-45fb-866b-73c2ecf80fcb"
    LDLFactorizations = "40e66cde-538c-5869-a4ad-c39174c6795b"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "GPUArraysCore", "InteractiveUtils", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "OpenBLAS_jll", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "9112d1a01dad29464dbec00e4ac10b499dfe5bcb"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "3.45.0"

    [deps.LinearSolve.extensions]
    LinearSolveAMDGPUExt = "AMDGPU"
    LinearSolveBLISExt = ["blis_jll", "LAPACK_jll"]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveCUSOLVERRFExt = ["CUSOLVERRF", "SparseArrays"]
    LinearSolveCliqueTreesExt = ["CliqueTrees", "SparseArrays"]
    LinearSolveEnzymeExt = "EnzymeCore"
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveFastLapackInterfaceExt = "FastLapackInterface"
    LinearSolveForwardDiffExt = "ForwardDiff"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolveMooncakeExt = "Mooncake"
    LinearSolvePardisoExt = ["Pardiso", "SparseArrays"]
    LinearSolveRecursiveFactorizationExt = "RecursiveFactorization"
    LinearSolveSparseArraysExt = "SparseArrays"
    LinearSolveSparspakExt = ["SparseArrays", "Sparspak"]

    [deps.LinearSolve.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    CUSOLVERRF = "a8cc9031-bad2-4722-94f5-40deabb4245c"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    FastLapackInterface = "29a986be-02c6-4525-aec4-84b980013641"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    LAPACK_jll = "51474c39-65e3-53ba-86ba-03b1b862ec14"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Sparspak = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
    blis_jll = "6136c539-28a5-5bf0-87cc-b183200dce32"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoweredCodeUtils]]
deps = ["CodeTracking", "Compiler", "JuliaInterpreter"]
git-tree-sha1 = "e24491cb83551e44a69b9106c50666dea9d953ab"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.4.4"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "Printf", "StatsAPI"]
git-tree-sha1 = "f386224fa41af0c27f45e2f9a8f323e538143b43"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.15.1"

[[deps.Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "191686b1ac1ea9c89fc52e996ad15d1d241d1e33"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.10.1+0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "9341048b9f723f2ae2a72a5269ac2f15f80534dc"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.3.2+0"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "c105fe467859e7f6e9a852cb15cb4301126fac07"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.11"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "e214f2a20bdd64c04cd3e4ff62d3c9be7e969a59"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.5.4+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "1d7d16f0e02ec063becd7a140f619b2ffe5f2b11"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.22.10"

[[deps.MakieCore]]
deps = ["ColorTypes", "GeometryBasics", "IntervalSets", "Observables"]
git-tree-sha1 = "c3159eb1e3aa3e409edbb71f4035ed8b1fc16e23"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.9.5"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "a370fef694c109e1950836176ed0d5eabbb65479"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.6"

[[deps.MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "3bb3cf4685f1c90f22883f4c4bb6d203fa882b79"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "3.1.3"
weakdeps = ["BandedMatrices"]

    [deps.MatrixFactorizations.extensions]
    MatrixFactorizationsBandedMatricesExt = "BandedMatrices"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "54e2fdc38130c05b42be423e90da3bade29b74bd"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.4"
weakdeps = ["SparseArrays"]

    [deps.MaybeInplace.extensions]
    MaybeInplaceSparseArraysExt = "SparseArrays"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bc95bf4149bf535c09602e3acdf950d9b4376227"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.ModelingToolkit]]
deps = ["ADTypes", "AbstractTrees", "ArrayInterface", "BlockArrays", "ChainRulesCore", "Combinatorics", "CommonSolve", "Compat", "ConstructionBase", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "DiffRules", "DifferentiationInterface", "Distributed", "Distributions", "DocStringExtensions", "DomainSets", "DynamicQuantities", "EnumX", "ExprTools", "FindFirstFunctions", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "Graphs", "InteractiveUtils", "JuliaFormatter", "JumpProcesses", "Latexify", "Libdl", "LinearAlgebra", "MLStyle", "Moshi", "NaNMath", "NonlinearSolve", "OffsetArrays", "OrderedCollections", "PrecompileTools", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SCCNonlinearSolve", "SciMLBase", "SciMLStructures", "Serialization", "Setfield", "SimpleNonlinearSolve", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "SymbolicUtils", "Symbolics", "URIs", "UnPack", "Unitful"]
git-tree-sha1 = "c323b230fe37eb029771a0670a07a2cf805c2ac2"
uuid = "961ee093-0014-501f-94e3-6117800e7a78"
version = "9.83.0"

    [deps.ModelingToolkit.extensions]
    MTKBifurcationKitExt = "BifurcationKit"
    MTKCasADiDynamicOptExt = "CasADi"
    MTKDeepDiffsExt = "DeepDiffs"
    MTKFMIExt = "FMI"
    MTKInfiniteOptExt = "InfiniteOpt"
    MTKLabelledArraysExt = "LabelledArrays"

    [deps.ModelingToolkit.weakdeps]
    BifurcationKit = "0f109fa4-8a5d-4b75-95aa-f515264e7665"
    CasADi = "c49709b8-5c63-11e9-2fb2-69db5844192f"
    DeepDiffs = "ab62b9b5-e342-54a8-a765-a90f495de1a6"
    FMI = "14a09403-18e3-468f-ad8a-74f8dda2d9ac"
    InfiniteOpt = "20393b10-9daf-11e9-18c9-8db751c92c57"
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "53f817d3e84537d84545e0ad749e483412dd6b2a"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.5.20"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "fade91fe9bee7b142d332fc6ab3f0deea29f637b"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.9"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "22df8573f8e7c593ac205455ca088989d0a2c7a0"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.7"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "25a6638571a902ecfb1ae2a18fc1575f86b1d4df"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.10.0"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "NonlinearSolveBase", "NonlinearSolveFirstOrder", "NonlinearSolveQuasiNewton", "NonlinearSolveSpectralMethods", "PrecompileTools", "Preferences", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SparseArrays", "SparseMatrixColorings", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "d2ec18c1e4eccbb70b64be2435fc3b06fbcdc0a1"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "4.10.0"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = ["NLsolve", "LineSearches"]
    NonlinearSolvePETScExt = ["PETSc", "MPI"]
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"
    NonlinearSolveSundialsExt = "Sundials"

    [deps.NonlinearSolve.weakdeps]
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    LineSearches = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    PETSc = "ace2c81b-2b5f-4b1e-a30d-d662738edfe0"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"
    Sundials = "c3572dad-4567-51f8-b174-8c6c989267f4"

[[deps.NonlinearSolveBase]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "Compat", "ConcreteStructs", "DifferentiationInterface", "EnzymeCore", "FastClosures", "LinearAlgebra", "Markdown", "MaybeInplace", "Preferences", "Printf", "RecursiveArrayTools", "SciMLBase", "SciMLJacobianOperators", "SciMLOperators", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "1d42a315ba627ca0027d49d0efb44e3d88db24aa"
uuid = "be0214bd-f91f-a760-ac4e-3421ce2b2da0"
version = "1.14.0"
weakdeps = ["BandedMatrices", "DiffEqBase", "ForwardDiff", "LineSearch", "LinearSolve", "SparseArrays", "SparseMatrixColorings"]

    [deps.NonlinearSolveBase.extensions]
    NonlinearSolveBaseBandedMatricesExt = "BandedMatrices"
    NonlinearSolveBaseDiffEqBaseExt = "DiffEqBase"
    NonlinearSolveBaseForwardDiffExt = "ForwardDiff"
    NonlinearSolveBaseLineSearchExt = "LineSearch"
    NonlinearSolveBaseLinearSolveExt = "LinearSolve"
    NonlinearSolveBaseSparseArraysExt = "SparseArrays"
    NonlinearSolveBaseSparseMatrixColoringsExt = "SparseMatrixColorings"

[[deps.NonlinearSolveFirstOrder]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConcreteStructs", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLJacobianOperators", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "3f1198ae5cbf21e84b8251a9e62fa1f888f3e4cb"
uuid = "5959db7a-ea39-4486-b5fe-2dd0bf03d60d"
version = "1.7.0"

[[deps.NonlinearSolveQuasiNewton]]
deps = ["ArrayInterface", "CommonSolve", "ConcreteStructs", "DiffEqBase", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLOperators", "StaticArraysCore"]
git-tree-sha1 = "40dfaf1bf74f1f700f81d0002d4dd90999598eb2"
uuid = "9a2c21bd-3a47-402d-9113-8faf9a0ee114"
version = "1.8.0"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveQuasiNewton.extensions]
    NonlinearSolveQuasiNewtonForwardDiffExt = "ForwardDiff"

[[deps.NonlinearSolveSpectralMethods]]
deps = ["CommonSolve", "ConcreteStructs", "DiffEqBase", "LineSearch", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "84de5a469e119eb2c22ae07c543dc4e7f7001ee7"
uuid = "26075421-4e9a-44e1-8bd1-420ed7ad02b2"
version = "1.3.0"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveSpectralMethods.extensions]
    NonlinearSolveSpectralMethodsForwardDiffExt = "ForwardDiff"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "567515ca155d0020a45b05175449b499c63e7015"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.29+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML", "Zlib_jll"]
git-tree-sha1 = "ec764453819f802fc1e144bfe750c454181bd66d"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "5.0.8+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["Compat", "EnumX", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "61942645c38dd2b5b78e2082c9b51ab315315d10"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.13.2"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c392fc5dd032381919e3b22dd32d6443760ce7ea"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.5.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "ExplicitImports", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqAdamsBashforthMoulton", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqExplicitRK", "OrdinaryDiffEqExponentialRK", "OrdinaryDiffEqExtrapolation", "OrdinaryDiffEqFIRK", "OrdinaryDiffEqFeagin", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqHighOrderRK", "OrdinaryDiffEqIMEXMultistep", "OrdinaryDiffEqLinear", "OrdinaryDiffEqLowOrderRK", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqNordsieck", "OrdinaryDiffEqPDIRK", "OrdinaryDiffEqPRK", "OrdinaryDiffEqQPRK", "OrdinaryDiffEqRKN", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqSSPRK", "OrdinaryDiffEqStabilizedIRK", "OrdinaryDiffEqStabilizedRK", "OrdinaryDiffEqSymplecticRK", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "Static", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "92cbff9044dc0a035b859de3778a9d0bfe73bdea"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.102.0"

[[deps.OrdinaryDiffEqAdamsBashforthMoulton]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqLowOrderRK", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "09aae1486c767caa6bce9de892455cbdf5a6fbc8"
uuid = "89bda076-bce5-4f1c-845f-551c83cdda9a"
version = "1.5.0"

[[deps.OrdinaryDiffEqBDF]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqSDIRK", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "ce8db53fd1e4e41c020fd53961e7314f75e4c21c"
uuid = "6ad6398a-0878-4a85-9266-38940aa047c8"
version = "1.10.1"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "FastBroadcast", "FastClosures", "FastPower", "FillArrays", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleUnPack", "Static", "StaticArrayInterface", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "e579c9a4f9102e82da3d97c349a74d6bc11cf8dc"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "1.30.0"

    [deps.OrdinaryDiffEqCore.extensions]
    OrdinaryDiffEqCoreEnzymeCoreExt = "EnzymeCore"
    OrdinaryDiffEqCoreMooncakeExt = "Mooncake"

    [deps.OrdinaryDiffEqCore.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"

[[deps.OrdinaryDiffEqDefault]]
deps = ["ADTypes", "DiffEqBase", "EnumX", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "PrecompileTools", "Preferences", "Reexport", "SciMLBase"]
git-tree-sha1 = "7d5ddeee97e1bdcc848f1397cbc3d03bd57f33e7"
uuid = "50262376-6c5a-4cf5-baba-aaf4f84d72d7"
version = "1.8.0"

[[deps.OrdinaryDiffEqDifferentiation]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqCore", "SciMLBase", "SciMLOperators", "SparseMatrixColorings", "StaticArrayInterface", "StaticArrays"]
git-tree-sha1 = "320b5f3e4e61ca0ad863c63c803f69973ba6efce"
uuid = "4302a76b-040a-498a-8c04-15b101fed76b"
version = "1.16.1"
weakdeps = ["SparseArrays"]

    [deps.OrdinaryDiffEqDifferentiation.extensions]
    OrdinaryDiffEqDifferentiationSparseArraysExt = "SparseArrays"

[[deps.OrdinaryDiffEqExplicitRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "4c0633f587395d7aaec0679dc649eb03fcc74e73"
uuid = "9286f039-9fbf-40e8-bf65-aa933bdc4db0"
version = "1.4.0"

[[deps.OrdinaryDiffEqExponentialRK]]
deps = ["ADTypes", "DiffEqBase", "ExponentialUtilities", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "3b81416ff11e55ea0ae7b449efc818256d9d450b"
uuid = "e0540318-69ee-4070-8777-9e2de6de23de"
version = "1.8.0"

[[deps.OrdinaryDiffEqExtrapolation]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastPower", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "9e1b11cf448a2c1bca640103c1c848a20aa2f967"
uuid = "becaefa8-8ca2-5cf9-886d-c06f3d2bd2c4"
version = "1.9.0"

[[deps.OrdinaryDiffEqFIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastGaussQuadrature", "FastPower", "LinearAlgebra", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "b968d66de3de5ffcf18544bc202ca792bad20710"
uuid = "5960d6e9-dd7a-4743-88e7-cf307b64f125"
version = "1.16.0"

[[deps.OrdinaryDiffEqFeagin]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "815b54211201ec42b8829e0275ab3c9632d16cbe"
uuid = "101fe9f7-ebb6-4678-b671-3a81e7194747"
version = "1.4.0"

[[deps.OrdinaryDiffEqFunctionMap]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "fe750e4b8c1b1b9e1c1319ff2e052e83ad57b3ac"
uuid = "d3585ca7-f5d3-4ba6-8057-292ed1abd90f"
version = "1.5.0"

[[deps.OrdinaryDiffEqHighOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "42096f72136078fa02804515f1748ddeb1f0d47d"
uuid = "d28bc4f8-55e1-4f49-af69-84c1a99f0f58"
version = "1.5.0"

[[deps.OrdinaryDiffEqIMEXMultistep]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "a5dcd75959dada0005b1707a5ca9359faa1734ba"
uuid = "9f002381-b378-40b7-97a6-27a27c83f129"
version = "1.7.0"

[[deps.OrdinaryDiffEqLinear]]
deps = ["DiffEqBase", "ExponentialUtilities", "LinearAlgebra", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "925fc0136e8128fd19abf126e9358ec1f997390f"
uuid = "521117fe-8c41-49f8-b3b6-30780b3f0fb5"
version = "1.6.0"

[[deps.OrdinaryDiffEqLowOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "3cc4987c8e4725276b55a52e08b56ded4862917e"
uuid = "1344f307-1e59-4825-a18e-ace9aa3fa4c6"
version = "1.6.0"

[[deps.OrdinaryDiffEqLowStorageRK]]
deps = ["Adapt", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "StaticArrays"]
git-tree-sha1 = "e6bd0a7fb6643a57b06a90415608a81aaf7bd772"
uuid = "b0944070-b475-4768-8dec-fb6eb410534d"
version = "1.7.0"

[[deps.OrdinaryDiffEqNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "PreallocationTools", "RecursiveArrayTools", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "StaticArrays"]
git-tree-sha1 = "b05226afc8fa6b8fc6f2258a89987b4f5bd0db4e"
uuid = "127b3ac7-2247-4354-8eb6-78cf4e7c58e8"
version = "1.14.1"

[[deps.OrdinaryDiffEqNordsieck]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "c90aa7fa0d725472c4098096adf6a08266c2f682"
uuid = "c9986a66-5c92-4813-8696-a7ec84c806c8"
version = "1.4.0"

[[deps.OrdinaryDiffEqPDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "9d599d2eafdf74ab26ea6bf3feb28183a2ade143"
uuid = "5dd0a6cf-3d4b-4314-aa06-06d4e299bc89"
version = "1.6.0"

[[deps.OrdinaryDiffEqPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "Reexport", "SciMLBase"]
git-tree-sha1 = "8e35132689133255be6d63df4190b5fc97b6cf2b"
uuid = "5b33eab2-c0f1-4480-b2c3-94bc1e80bda1"
version = "1.4.0"

[[deps.OrdinaryDiffEqQPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "63fb643a956b27cd0e33a3c6d910c3c118082e0f"
uuid = "04162be5-8125-4266-98ed-640baecc6514"
version = "1.4.0"

[[deps.OrdinaryDiffEqRKN]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "a31c41f9dbea7c7179c6e544c25c7e144d63868c"
uuid = "af6ede74-add8-4cfd-b1df-9a4dbb109d7a"
version = "1.5.0"

[[deps.OrdinaryDiffEqRosenbrock]]
deps = ["ADTypes", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "f34bc2f58656843596d09a4c4de8c20724ebc2f1"
uuid = "43230ef6-c299-4910-a778-202eb28ce4ce"
version = "1.18.1"

[[deps.OrdinaryDiffEqSDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "20caa72c004414435fb5769fadb711e96ed5bcd4"
uuid = "2d112036-d095-4a1e-ab9a-08536f3ecdbf"
version = "1.7.0"

[[deps.OrdinaryDiffEqSSPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "StaticArrays"]
git-tree-sha1 = "3bce87977264916bd92455754ab336faec68bf8a"
uuid = "669c94d9-1f4b-4b64-b377-1aa079aa2388"
version = "1.7.0"

[[deps.OrdinaryDiffEqStabilizedIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqStabilizedRK", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "75abe7462f4b0b2a2463bb512c8a5458bbd39185"
uuid = "e3e12d00-db14-5390-b879-ac3dd2ef6296"
version = "1.6.0"

[[deps.OrdinaryDiffEqStabilizedRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "7e94d3d1b3528b4bcf9e0248198ee0a2fd65a697"
uuid = "358294b1-0aab-51c3-aafe-ad5ab194a2ad"
version = "1.4.0"

[[deps.OrdinaryDiffEqSymplecticRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "e8dd5ab225287947016dc144a5ded1fb83885638"
uuid = "fa646aed-7ef9-47eb-84c4-9443fc8cbfa8"
version = "1.7.0"

[[deps.OrdinaryDiffEqTsit5]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "778c7d379265f17f40dbe9aaa6f6a2a08bc7fa3e"
uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a"
version = "1.5.0"

[[deps.OrdinaryDiffEqVerner]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "185578fa7c38119d4318326f9375f1cba0f0ce53"
uuid = "79d7bb75-1356-48c1-b8c0-6832512096c2"
version = "1.6.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "d922b4d80d1e12c658da7785e754f4796cc1d60d"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.36"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.PROPACK]]
deps = ["LinearAlgebra", "LinearOperators", "PROPACK_jll"]
git-tree-sha1 = "f500ecef3edcf1c8889fc861a0debcdb6c714c57"
uuid = "b169e327-5944-5131-97a6-5d3d3f0a476a"
version = "0.5.0"

[[deps.PROPACK_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg", "libblastrampoline_jll"]
git-tree-sha1 = "b6141cf260569331be16207e0886e71f71ccd5b1"
uuid = "2c99449f-f194-5e92-8ff0-13814ff41dcf"
version = "0.2.1+0"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1f7f9bbd5f7a2e5a9f7d96e51c9754454ea7f60b"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.4+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.PoissonRandom]]
deps = ["LogExpFunctions", "Random"]
git-tree-sha1 = "67afbcbe9e184d6729a92a022147ed4cf972ca7b"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.7"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "6f7cd22a802094d239824c57d94c8e2d0f7cfc7d"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.18"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "PrecompileTools"]
git-tree-sha1 = "07626c22057b3562d6f5f92b3402dbd5bfe82de8"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.30"

    [deps.PreallocationTools.extensions]
    PreallocationToolsForwardDiffExt = "ForwardDiff"
    PreallocationToolsReverseDiffExt = "ReverseDiff"
    PreallocationToolsSparseConnectivityTracerExt = "SparseConnectivityTracer"

    [deps.PreallocationTools.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "d95ed0324b0799843ac6f7a6a85e65fe4e5173f0"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.5"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "8b3fc30bc0390abdce15f8822c889f669baed73d"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "dbe5fd0b334694e905cb9fda73cd8554333c46e2"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.1"

[[deps.RandomNumbers]]
deps = ["Random"]
git-tree-sha1 = "c6ec94d2aaba1ab2ff983052cf6a606ca5985902"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.6.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "51bdb23afaaa551f923a0e990f7c44a4451a26f1"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.39.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsKernelAbstractionsExt = "KernelAbstractions"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTablesExt = ["Tables"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "20ccb7e2501e9da93fe8450d01aeabf16a5f0c82"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.8.1"
weakdeps = ["Distributed"]

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "86a8a8b783481e1ea6b9c91dd949cb32191f8ab4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.15"

[[deps.SCCNonlinearSolve]]
deps = ["CommonSolve", "PrecompileTools", "Reexport", "SciMLBase", "SymbolicIndexingInterface"]
git-tree-sha1 = "5595105cef621942aceb1aa546b883c79ccbfa8f"
uuid = "9dfe8606-65a1-4bb3-9748-cb89d1561431"
version = "1.4.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "1a1b2ed5178814d32374d98f808163c24d47f1d8"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.110.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMLStyleExt = "MLStyle"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    MLStyle = "d8e11817-5142-5d16-987a-aa16d5891078"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLJacobianOperators]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DifferentiationInterface", "FastClosures", "LinearAlgebra", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "a273b291c90909ba6fe08402dd68e09aae423008"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.11"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "c1053ba68ede9e4005fc925dd4e8723fcd96eef8"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.9.0"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLPublic]]
git-tree-sha1 = "ed647f161e8b3f2973f24979ec074e8d084f1bee"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.0"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "566c4ed301ccb2a44cbd5a27da5f885e0ed1d5df"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "818554664a2e01fc3784becb2eb3a82326a604b6"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "09d986e27a606f172c5b6cffbd8b8b2f10bf1c75"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "2.7.0"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolveDiffEqBaseExt = "DiffEqBase"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveTrackerExt = "Tracker"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffEqBase = "2b5f629d-d688-5b77-993f-72d75c75574e"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SparseConnectivityTracer]]
deps = ["ADTypes", "DocStringExtensions", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "182990067a09adf950274f97f38f68c76f81d2d0"
uuid = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
version = "0.6.21"

    [deps.SparseConnectivityTracer.extensions]
    SparseConnectivityTracerDataInterpolationsExt = "DataInterpolations"
    SparseConnectivityTracerLogExpFunctionsExt = "LogExpFunctions"
    SparseConnectivityTracerNNlibExt = "NNlib"
    SparseConnectivityTracerNaNMathExt = "NaNMath"
    SparseConnectivityTracerSpecialFunctionsExt = "SpecialFunctions"

    [deps.SparseConnectivityTracer.weakdeps]
    DataInterpolations = "82cc6244-b520-54b8-b5a6-8a565e85f1d0"
    LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "Random", "SparseArrays"]
git-tree-sha1 = "d3f3b7bb8a561b5ff20ee7cf9574841ee4e4e137"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.22"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsCUDAExt = "CUDA"
    SparseMatrixColoringsCliqueTreesExt = "CliqueTrees"
    SparseMatrixColoringsColorsExt = "Colors"

    [deps.SparseMatrixColorings.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "49440414711eddc7227724ae6e570c7d5559a086"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.3.1"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "b8693004b385c842357406e3af647701fe783f98"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.15"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "a136f98cefaf3e2924a66bd75173d1c891ab7453"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.7"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.SteadyStateDiffEq]]
deps = ["ConcreteStructs", "DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NonlinearSolveBase", "Reexport", "SciMLBase"]
git-tree-sha1 = "1919936d42482832a4b5d0b678513ad6756c6f00"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "2.6.0"

[[deps.StochasticDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FastPower", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SparseArrays", "StaticArrays", "UnPack"]
git-tree-sha1 = "c3a55a2a1e180e249a0550d30a58c700487aa7ef"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.81.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "83151ba8065a73f53ca2ae98bc7274d817aa30f2"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.8"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "a2c37d815bf00575332b7bd0389f771cb7987214"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.2"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "cd47aa083c9c7bdeb7b92de26deb46d6a33163c9"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.5.1"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "7c7a7ee705724b3c80d5451ac49779db36c6f758"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.28.0"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "SuiteSparse_jll", "libblastrampoline_jll"]
git-tree-sha1 = "91db7ed92c66f81435fe880947171f1212936b14"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.3+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "94c58884e013efff548002e8dc2fdd1cb74dfce5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.46"

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

    [deps.SymbolicIndexingInterface.weakdeps]
    PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "f75c7deb7e11eea72d2c1ea31b24070b713ba061"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.3"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "ExproniconLite", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "a85b4262a55dbd1af39bb6facf621d79ca6a322d"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "3.32.0"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "LaTeXStrings", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "OffsetArrays", "PrecompileTools", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "8206e177903a41519145f577cb7f3793f3b7c960"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.57.0"

    [deps.Symbolics.extensions]
    SymbolicsD3TreesExt = "D3Trees"
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"
    SymbolicsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Symbolics.weakdeps]
    D3Trees = "e3df1716-f71e-5df9-9e2d-98e193103c45"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TestItems]]
git-tree-sha1 = "42fd9023fef18b9b78c8343a4e2f3813ffbcefcb"
uuid = "1c621080-faea-4a02-84b6-bbd5e436b8fe"
version = "1.0.0"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "d969183d3d244b6c33796b5ed01ab97328f2db85"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.5"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "98b9352a24cb6a2066f9ababcc6802de9aed8ad8"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.6"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "83360bda12f61c250835830cc40b64f487cc2230"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.25.1"
weakdeps = ["ConstructionBase", "ForwardDiff", "InverseFunctions", "LaTeXStrings", "Latexify", "Printf"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    PrintfExt = "Printf"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "80d3930c6347cfce7ccf96bd3bafdf079d9c0390"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.9+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaec_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1aa23f01927b2dac46db77a56b31088feee0a491"
uuid = "477f73a3-ac25-53e9-8cc3-50b2fa2566f0"
version = "1.1.4+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.13.1+1"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "07b6a107d926093898e82b3b1db657ebe33134ec"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.50+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.5.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"
"""

# ╔═╡ Cell order:
# ╠═453819dc-e207-46db-9c9b-9f10f7fd4daf
# ╠═7813ebfd-f674-44da-b2d5-dfed9e028457
# ╠═cb9fa0aa-1588-4647-8b02-69fca1415974
# ╠═ae7ff050-c76d-4c0d-8c8a-23724f10485a
# ╠═85c4ad39-d14a-4e73-bfa3-4a22fef4d537
# ╠═e3e0c398-9db0-420a-86b3-65220761028c
# ╠═91bec3fc-65cc-4344-98d2-09f4944fc0e5
# ╟─ebd5f20c-88c2-46d1-982f-de864ba33727
# ╠═256f8d0d-7c2b-4415-969e-356ee271e908
# ╠═ca5f05c8-96ea-42a8-adb0-c095a223b4cf
# ╠═6ad6db9b-f87c-4320-9a45-a18a170b5db7
# ╠═7ea7ecaf-7a46-4b9b-a0e0-56171ffd7c23
# ╠═18962357-9f6a-4668-82cb-7f133a130b5e
# ╠═fd852a0c-24fe-4807-afb5-c9157d3c746a
# ╠═f541864a-38f7-4081-9a71-a1c6de95e5cf
# ╠═3750408a-40e0-4ff4-b199-605886e3f0c2
# ╠═83cf9310-51dc-458f-a82b-560d891388b4
# ╠═40a18845-9662-4028-ab7d-3480a9d5e9b1
# ╠═6ba472e0-0c0f-4840-9640-fa4e848e2f8d
# ╠═01a13951-2701-48bc-bb48-6c9594c89888
# ╠═035c4b08-7be6-4b89-b6cc-49374ee36418
# ╠═838fb0ef-03b6-4501-b66c-07859355e5e2
# ╠═32e6fa3b-dc9f-4cca-bfb1-380b4dddbf3c
# ╠═b9cf4dc2-9817-493b-940b-5e3677ceb92e
# ╠═690d5609-6bd0-4caa-b88e-a28e3b9da19c
# ╠═e2e71c4e-4782-46bb-8560-9ffbae1145b3
# ╠═d0ad5800-4b37-4722-b884-07af915a3d7f
# ╠═c1a92cd0-ffed-43d7-81c5-021799578d0d
# ╠═5e389838-a3cc-46d2-986d-c407edc8e93a
# ╠═9867082b-f32b-4c92-adbf-866efeec8ced
# ╠═9e8a68b4-630f-4051-a764-dc44ccc8847e
# ╠═8d9d0907-dc0d-4ad2-9faf-fbd015eb7234
# ╠═6f4ca781-97e6-4c70-a200-ea230eebd449
# ╠═afe2ea2f-b3e8-49f9-a42a-bbfc60d53c04
# ╠═e7b739ac-a39e-494f-8ba1-4e554f9ac06c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
