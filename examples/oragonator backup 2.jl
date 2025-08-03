### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ f61fe968-f0ad-4d2c-9af5-c27f9224a1cc
begin
	# numerical libraries
	using ExponentialUtilities, Expokit, PROPACK, Arpack, SparseArrays
	# output and plotting
	using ProgressLogging, JLD
	# modelling and statistics 
	using Catalyst, JumpProcesses, StatsBase, DifferentialEquations
	using Interpolations
	# importing local fsp package
	using Revise
	local_mod = include("../src/DiscStochSim.jl")
	using .local_mod.DiscStochSim
end

# ╔═╡ 9dd29549-5238-4314-a3dd-96f4a9e385d9
oregonator_rn = @reaction_network begin
    @species X(t) Y(t) Z(t)
    @parameters k1_eff k2 k3_eff k4 k5

    # The 5 reactions of the standard Oregonator
    k1_eff, Y --> X
    k2, X + Y --> 0
    k3_eff, X --> 2X + Z
    k4, 2X --> 0
    k5, Z --> Y
end

# ╔═╡ 4e84db3b-5377-4f42-ae38-c7f867fd5827
model = DiscreteStochasticSystem(oregonator_rn)

# ╔═╡ 64061d79-cbfe-4219-a7a8-cdc57ffb9344
oregonator_gillespie_params = begin
    # Rate constants calculated directly from Gillespie (1977), equations 61a-e and 62
    k1_eff = 2.0       # Corresponds to c₁*X₁ in the paper
    k2 = 0.1         # Corresponds to c₂
    k3_eff = 104.0     # Corresponds to c₃*X₂
    k4 = 0.016       # Corresponds to c₄
    k5 = 26.0        # Corresponds to c₅*X₃
    
    # The final rates vector for the solver
    rates = [k1_eff, k2, k3_eff, k4, k5]

    # Bounds must be large enough to contain the oscillations shown in the paper's Fig. 24
    bounds = (0, 10000) 
    boundary_condition(x) = RectLatticeBoundaryCondition(x, bounds)
end

# ╔═╡ 9f592c0e-313e-48dd-bcbe-22e4b41ad55a
fsp_sim_oregonator = begin

    tf = 50.0   # Simulate for enough time to see several oscillations
    ϵ_dt = 10.0 # Tolerance for adaptive dt calculation

    # Initial state (X, Y) near the unstable steady state
    U₀ = CartesianIndex(500, 1000, 2000)
    𝒮ₜ = Set([U₀])
    pₜ = zeros(length(𝒮ₜ))
    pₜ[FindElement(U₀, 𝒮ₜ)] = 1.0

    # Initialize current time and dynamic solution arrays
    local t = 0.0
    sol_t = [t]
    sol_S_size = [length(𝒮ₜ)]
    sol = [(copy(𝒮ₜ), copy(pₜ))]

    # --- Main Adaptive Loop ---
    
    while t < tf

        # 1. Calculate adaptive δt based on total expected system activity
        X_vec = collect(𝒮ₜ)
        total_flux = 0.0
        for i in eachindex(X_vec)
            state = X_vec[i]
            weight = maximum([prop(state, rates, t) for prop in model.propensities])
            total_flux += pₜ[i] * weight
        end

        δt = (total_flux > 0.0) ? (ϵ_dt / total_flux) : (tf - t)
        δt = min(δt, tf - t)

        # 2. Expand state space
        global 𝒮ₜ, pₜ = expand!(𝒮ₜ, pₜ, model, rates, t, boundary_condition, 2)

        # 3. Build Master Equation and Evolve over the calculated δt
        A = MasterEquation(𝒮ₜ, model, rates, boundary_condition, t)
        pₜ = expv(δt, A, pₜ)
        
        # 4. Purge state space using the robust flux-based method
        𝒮ₜ, pₜ = purge!(𝒮ₜ, pₜ, model, rates, t, 0.6, 1e-6)
        
        # 5. Renormalize probability
        pₜ ./= sum(pₜ)
        
        # 6. Update time and store results
        t += δt
        
        push!(sol, (copy(𝒮ₜ), copy(pₜ)))
        push!(sol_t, t)
        push!(sol_S_size, length(𝒮ₜ))
    end
    
end

# ╔═╡ 500e6bd4-2fc6-4b6b-8923-5faf534e15cf
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	# Initial state (X, Y) near the unstable steady state
	local U₀ = CartesianIndex(10, 100, 10)
	local 𝒮ₜ = Set([U₀])
	local 𝒮ₜ= expand!(𝒮ₜ, model, boundary_condition,4)
	@show length(𝒮ₜ)
	𝒮ₜ
end
  ╠═╡ =#

# ╔═╡ a4822e24-f8d9-40ac-974c-364dc93c4eac
begin
	sol_mean = map(1:size(sol)[1]) do i
	    sum(collect.(Tuple.(sol[i][1])) .* sol[i][2])
	end 
	fsp_mean=hcat(sol_mean...)'
end

# ╔═╡ dcf8bd9b-d8ae-4a30-bae3-e6d4b08925b2
begin
	using CairoMakie
	
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
	
	        # Spine (box around the plot) settings
	        leftspinevisible = true,
	        rightspinevisible = false,
	        topspinevisible = false,
	        bottomspinevisible = true,
	    ),
	    
	    Axis3 = (
	        protrusions = (40, 20, 20, 20), # Adjust space for labels
	        xlabelsize = 18,
	        ylabelsize = 18,
	        zlabelsize = 18,
	        xticklabelsize = 14,
	        yticklabelsize = 14,
	        zticklabelsize = 14,
	        titlefont = :bold,
	        titlesize = 22,
	        xgridcolor = :gray80,
	        ygridcolor = :gray80,
	        zgridcolor = :gray80,
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
	
	# --- Generate the 6-Panel Figure ---
	begin
	    # Apply the custom theme for all subsequent plots
	    set_theme!(siam_theme)
	
	    # Prepare the data from your simulation results
	    # Ensure fsp_mean has dimensions (timesteps, species)
	    time_points = sol_t
	    dt_values = vcat([0.0], sol_t[2:end] .- sol_t[1:end-1])
	    state_space_sizes = sol_S_size
	    mean_y1 = fsp_mean[:, 1]
	    mean_y2 = fsp_mean[:, 2]
	    mean_y3 = fsp_mean[:, 3]
	
	    # Create the figure with a 2x3 grid layout
	    fig = Figure(size = (1200, 800))
	
	    # --- Top Row: Time Series Plots ---
	    
	    # Panel 1: Adaptive Time Step (dt)
	    ax11 = Axis(fig[1, 1], title = "Adaptive Time Step", xlabel = "Simulation Time, t", ylabel = "δt")
	    lines!(ax11, time_points, dt_values)
	
	    # Panel 2: State Space Size
	    ax12 = Axis(fig[1, 2], title = "State Space Size", xlabel = "Simulation Time, t", ylabel = "|S|")
	    lines!(ax12, time_points, state_space_sizes)
	    
	    # Panel 3: Mean Trajectories
	    ax13 = Axis(fig[1, 3], title = "Mean Trajectories", xlabel = "Simulation Time, t", ylabel = "Mean Population")
	    lines!(ax13, time_points, fsp_mean[:, 1], label = "X (y₁)")
	    lines!(ax13, time_points, fsp_mean[:, 2], label = "Y (y₂)", linestyle = :dash)
	    lines!(ax13, time_points, fsp_mean[:, 3], label = "Z (y₃)", linestyle = :dot)
	    axislegend(ax13, position = :rc) # Place legend in the right-center
	
	    # Link the x-axes of the top row for synchronized zooming/panning
	    linkxaxes!(ax11, ax12, ax13)
	
	    # --- Bottom Row: 3D Limit Cycle Perspectives ---
	
	    # Panel 4: Perspective 1 (Isometric View)
	    ax21 = Axis3(fig[2, 1], title = "Limit Cycle (Isometric View)", xlabel = "X", ylabel = "Y", zlabel = "Z", 
	                   azimuth = 0.25pi, elevation = 0.15pi)
	    lines!(ax21, mean_y1, mean_y2, mean_y3)
	
	    # Panel 5: Perspective 2 (Top-Down View, X-Y Plane)
	    ax22 = Axis3(fig[2, 2], title = "Limit Cycle (X-Y Plane)", xlabel = "X", ylabel = "Y", zlabel = "Z",
	                   azimuth = 0.5pi, elevation = 0.5pi) # Looking straight down the Z-axis
	    lines!(ax22, mean_y1, mean_y2, mean_y3)
	
	    # Panel 6: Perspective 3 (Side View, Y-Z Plane)
	    ax23 = Axis3(fig[2, 3], title = "Limit Cycle (Y-Z Plane)", xlabel = "X", ylabel = "Y", zlabel = "Z",
	                   azimuth = 0.0pi, elevation = 0.0pi) # Looking straight along the X-axis
	    lines!(ax23, mean_y1, mean_y2, mean_y3)
	
	    # Add spacing between rows for clarity
	    rowgap!(fig.layout, 1, 60)
	
	    # Display the final figure
	    fig
	end
	
end

# ╔═╡ 5c736465-5bfc-45e6-8445-266e3eaffd47
plot(sol_S_size, title="state space size")

# ╔═╡ f0de5260-35fb-4b19-ada4-13d7764151f8
# ╠═╡ skip_as_script = true
#=╠═╡
plot(log.(collect(0:0.1:1)))
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═f61fe968-f0ad-4d2c-9af5-c27f9224a1cc
# ╠═9dd29549-5238-4314-a3dd-96f4a9e385d9
# ╠═4e84db3b-5377-4f42-ae38-c7f867fd5827
# ╠═64061d79-cbfe-4219-a7a8-cdc57ffb9344
# ╠═9f592c0e-313e-48dd-bcbe-22e4b41ad55a
# ╠═500e6bd4-2fc6-4b6b-8923-5faf534e15cf
# ╠═a4822e24-f8d9-40ac-974c-364dc93c4eac
# ╠═dcf8bd9b-d8ae-4a30-bae3-e6d4b08925b2
# ╠═5c736465-5bfc-45e6-8445-266e3eaffd47
# ╠═f0de5260-35fb-4b19-ada4-13d7764151f8
