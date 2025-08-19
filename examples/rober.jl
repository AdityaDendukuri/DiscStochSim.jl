### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ d71daa7e-5fa2-11f0-0a79-ab2680039a47
begin
	# numerical libraries
	using ExponentialUtilities, PROPACK, Arpack, SparseArrays
	# output and plotting
	using ProgressLogging, JLD, CairoMakie, Expokit
	# modelling and statistics 
	using Catalyst, JumpProcesses, StatsBase, DifferentialEquations
	using Interpolations, LinearAlgebra, NaNMath
	# importing local fsp package
	using Revise
	local_mod = include("../src/DiscStochSim.jl")
	using .local_mod.DiscStochSim
end

# ╔═╡ eef14334-d43c-49d2-b40b-fe16a3747fb2
rn = @reaction_network begin
    k1, A --> B
    k2, 2B --> C + B
    k3, C + B --> A + C
end

# ╔═╡ 068fa5b2-629c-44a7-aaa1-f1e2ed68c840
model = DiscreteStochasticSystem(rn);

# ╔═╡ 60d1b852-a518-42a3-be5b-f8466fd8d9c5
rober_params = begin
	Ω = 1e5
	k1 = 0.04                     # X → Y
	k2 = 1.0e4 / Ω                # Y + Z → X + Z
	k3 = 3.0e7 / Ω  
    rates = [k1, k2, k3]
    bounds = (0, Ω+5) 
    global boundary_condition(x) = RectLatticeBoundaryCondition(x, bounds)
end

# ╔═╡ 08430a94-9c7f-4616-9141-7846bd7e2c7a

function dobrushin(A::SparseMatrixCSC{Float64,Int}, τ::Float64=1.0)
    n = size(A, 1)
    γ_min = 1.0

    for i in 1:n-1
        eᵢ = zeros(n); eᵢ[i] = 1.0
        Pi = expv(τ, A, eᵢ)

        for j in i+1:n
            eⱼ = zeros(n); eⱼ[j] = 1.0
            Pj = expv(τ, A, eⱼ)

            γ_ij = sum(min(Pi[k], Pj[k]) for k in 1:n)
            γ_min = min(γ_min, γ_ij)
        end
    end

    return 1 - γ_min
end


# ╔═╡ 7daa7457-7516-4566-8d6f-82b886c6eb05
function jump_matrix(A::AbstractMatrix{T}) where T<:Real
    Q̄ = copy(A)
    for i in 1:size(Q̄, 1)
        Q̄[i, i] = 0.0
    end
    return Q̄
end

# ╔═╡ 655b46ff-bec5-4dfe-9cfc-45e3418d8be5
fsp_sim_rober = begin
    # --- Initialization ---
    tf = 1e4
    ϵ_dt = 10.0# Tolerance for flux-based adaptive dt

    # Initial condition for Robertson problem
    U₀ = CartesianIndex(Int(Ω), 0, 0)
    global 𝒮ₜ = Set([U₀])
    pₜ = zeros(length(𝒮ₜ))
    pₜ[FindElement(U₀, 𝒮ₜ)] = 1.0

    # Time and diagnostics
    local t = 0.0
    sol_t = [t]
    sol_S_size = [length(𝒮ₜ)]
    sol = [(copy(𝒮ₜ), copy(pₜ))]
    flx = Float64[]
    tresh = Float64[]

	global 𝒮ₜ, pₜ = expand!(𝒮ₜ, pₜ, model, boundary_condition, 10)
	global A = MasterEquation(𝒮ₜ, model, rates, boundary_condition, t)
	global in_flux = Vector(-diag(A) .* pₜ)
	δt = 1e-5

	iter = 0
    while t < tf
		
        # Expand state space 
		global 𝒮ₜ, pₜ = expand!(𝒮ₜ, pₜ, in_flux, model, rates, t, boundary_condition, 1; flux_tol=1e-5)
        
		global in_flow, out_flow = ComputeFlow(𝒮ₜ, 
											   model, 
											   rates, 
											   boundary_condition, 
											   t) 
	
        A = MasterEquation(in_flow, out_flow)
		global in_flux, out_flux = Vector(in_flow * pₜ) , Vector(-out_flow * pₜ)
		total_flux = sum(in_flux)

		# based on out-flux, compute δt
		global δt = (total_flux > 0.0) ? (ϵ_dt / total_flux) : (tf - t)
        δt = min(δt, tf - t, 1e6) 

		# evolve probability distribution
		pₜ = expv(δt, A, pₜ)

		
        𝒮ₜ, pₜ = purge!(𝒮ₜ, pₜ, in_flux, model, rates, t, 0.1; flux_tolerance=1e-9)
		
        if sum(pₜ) <= 0
            @warn "All probability mass lost in pruning – skipping"
        end

        # 7. Normalize
        pₜ ./= sum(pₜ)

        # 8. Advance time
        t += δt

        # note .. im storing values in logarithemic intervals 
		# (too many snapshots to save)
		if iter % 1e3 == 0
        	push!(sol, (copy(𝒮ₜ), copy(pₜ)))
        	push!(sol_t, t)
        	push!(sol_S_size, length(𝒮ₜ))
        	push!(tresh, total_flux)
		end
		global iter = iter + 1
    end
end


# ╔═╡ ac8c249b-0284-4355-b410-8d47efcc829b
begin
	sol_mean = map(1:length(sol)) do i
	    sum(collect.(Tuple.(sol[i][1])) .* sol[i][2])
	end 
	fsp_mean=hcat(sol_mean...)'

end

# ╔═╡ 1a3f749e-6dce-42e8-9b2e-bd096e24e665
begin
	
	# --- Define a Camera-Ready Theme for SIAM Publications ---
	siam_theme = Theme(
	    # Set a professional, clean font (Computer Modern is standard for TeX/LaTeX)
	    font = "Computer Modern",
	    fontsize = 18, # Base font size for the figure
	    
	    Axis = (
	        # Axis labels and title settings
	        titlefont = :bold,
	        titlesize = 22,
	        xlabelsize = 20,
	        ylabelsize = 20,
	        xticklabelsize = 16,
	        yticklabelsize = 16,
	        
	        # Grid settings for a clean look
	        xgridstyle = :dash,
	        ygridstyle = :dash,
	        xgridwidth = 0.7,
	        ygridwidth = 0.7,
	        xgridcolor = :gray80,
	        ygridcolor = :gray80,
			xscale=log10,
	
	        # Spine (box around the plot) settings
	        leftspinevisible = true,
	        rightspinevisible = false,
	        topspinevisible = false,
	        bottomspinevisible = true,
	    ),
	
	    Legend = (
	        framevisible = false,
	        patchsize = (30, 10),
	        labelsize = 16,
	    ),
	
	    Lines = (
	        color = :black,
	        linewidth = 2.0,
	    )
	)
	
	# --- Generate the 4-Panel Figure ---
	begin
	    # Apply the custom theme for all subsequent plots
	    set_theme!(siam_theme)
	
	    # Prepare the data from your simulation results
	    # Ensure all vectors are the same length for plotting
	    local num_steps = length(sol_t)
	    local time_points = sol_t[1:num_steps]
	    global dt_values = vcat([0.0], sol_t[2:end] .- sol_t[1:end-1])
	    local state_space_sizes = sol_S_size[1:num_steps]
	    local total_fluxes = tresh[1:num_steps-1] # Assuming flx has one less element
	    local flux_thresholds = tresh[1:num_steps-2] # Assuming tresh has one less element
	    
	    # Ensure fsp_mean has dimensions (timesteps, species)
	    mean_A = fsp_mean[1:num_steps, 1]
	    mean_B = fsp_mean[1:num_steps, 2]
	    mean_C = fsp_mean[1:num_steps, 3]
	
	    # Create the figure with a 2x2 grid layout
	    local fig = Figure(size = (1200, 800))
	
	    # --- Top Row ---
	    
	    # Panel 1: Mean Trajectories
	    ax11 = Axis(fig[1, 1], title = "Mean Trajectories", xlabel = "Simulation Time, t", ylabel = "Mean Population")
	    lines!(ax11, time_points, mean_A, label = "A (y₁)")
	    lines!(ax11, time_points, mean_B .* 4.e3, label = "B (y₂) * 2.5e5", linestyle = :dash)
	    lines!(ax11, time_points, mean_C, label = "C (y₃)", linestyle = :dot)
	    axislegend(ax11, position = :rc)
	
	    # Panel 2: State Space Size
	    ax12 = Axis(fig[1, 2], title = "State Space Size", xlabel = "Simulation Time, t", ylabel = "|S|")
	    lines!(ax12, time_points, state_space_sizes)
	    
	    # --- Bottom Row ---
	
	    # Panel 3: Adaptive Time Step (dt)
	    ax21 = Axis(fig[2, 1], title = "Adaptive Time Step", xlabel = "Simulation Time, t", ylabel = "δt")
	    lines!(ax21, time_points, dt_values)
	
	    # Panel 4: Flux Diagnostics
	    ax22 = Axis(fig[2, 2], title = "Flux Diagnostics", xlabel = "Simulation Time, t", ylabel = "Flux Value")
	    #lines!(ax22, time_points[2:end], total_fluxes, label = "Total Flux")
	    lines!(ax22, time_points[3:end], flux_thresholds, label = "Flux Threshold", linestyle = :dash)
	    axislegend(ax22, position = :rc)
	    
	    # Link all x-axes for synchronized zooming/panning
	    linkxaxes!(ax11, ax12, ax21, ax22)
	
	    # Add spacing between rows for clarity
	    rowgap!(fig.layout, 1, 60)
	
	    # Display the final figure
	    fig
	end
end

# ╔═╡ Cell order:
# ╠═d71daa7e-5fa2-11f0-0a79-ab2680039a47
# ╠═eef14334-d43c-49d2-b40b-fe16a3747fb2
# ╠═068fa5b2-629c-44a7-aaa1-f1e2ed68c840
# ╠═60d1b852-a518-42a3-be5b-f8466fd8d9c5
# ╠═08430a94-9c7f-4616-9141-7846bd7e2c7a
# ╠═7daa7457-7516-4566-8d6f-82b886c6eb05
# ╠═655b46ff-bec5-4dfe-9cfc-45e3418d8be5
# ╠═ac8c249b-0284-4355-b410-8d47efcc829b
# ╠═1a3f749e-6dce-42e8-9b2e-bd096e24e665
