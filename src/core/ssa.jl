using Random
using Interpolations
using Statistics

# Define the reaction rates for the Lotka-Volterra model
function reaction_rates(X, Y, params)
    r1, r2, r3 = params
    # Reaction 1: Prey reproduction (A → 2A)
    rate1 = r1 * X
    # Reaction 2: Prey consumption by predators (A + B → B + B)
    rate2 = r2 * X * Y
    # Reaction 3: Predator death (B → ∅)
    rate3 = r3 * Y
    return [rate1, rate2, rate3]
end

# Define the stoichiometry for the reactions
function reaction_stoichiometries()
    # Reaction 1: Prey birth (A → 2A)
    r1 = [-1, 0]
    # Reaction 2: Prey consumed by predator (A + B → B + B)
    r2 = [-1, 1]
    # Reaction 3: Predator death (B → ∅)
    r3 = [0, -1]
    return [r1, r2, r3]
end

# Gillespie algorithm to simulate the system
function gillespie(X0, Y0, params, T_max)
    # Initial conditions
    X, Y = X0, Y0
    t = 0.0  # Initial time
    traj_t = [t]  # List to store time values
    traj_u = [[X, Y]]  # List to store the states (populations), using vectors
    
    # Get the reaction stoichiometries
    stoichs = reaction_stoichiometries()
    
    while t < T_max
        rates = reaction_rates(X, Y, params)  # Calculate the reaction rates
        total_rate = sum(rates)  # Total rate of reactions
        
        # If no reactions left, stop the simulation
        if total_rate == 0
            break
        end
        
        # Time to next reaction
        tau = -log(rand()) / total_rate
        t += tau
        
        # Randomly choose which reaction occurs
        r = rand() * total_rate
        cumulative_rate = 0.0
        reaction_idx = 0
        for i in 1:3
            cumulative_rate += rates[i]
            if r < cumulative_rate
                reaction_idx = i
                break
            end
        end
        
        # Update the populations based on the selected reaction
        if reaction_idx == 1
            X += 1  # Prey birth
        elseif reaction_idx == 2
            X -= 1  # Prey consumed
            Y += 1  # Predator population increases
        elseif reaction_idx == 3
            Y -= 1  # Predator dies
        end
        
        # Store the new state as a vector
        push!(traj_t, t)
        push!(traj_u, [X, Y])  # Store the state as a vector [X, Y]
    end

    # Return as a structure with time (t) and states (u)
    return (t=traj_t, u=traj_u)
end


# Function to interpolate and compute the mean over a uniform grid
function mean_trajectory(trajectories, T_min, T_max, num_points)
    # Create the uniform grid
    uniform_time = LinRange(T_min, T_max, num_points)
    n_species = length(trajectories[1].u[1])
    
    # Interpolate each trajectory onto the uniform grid
    traj_mean = zeros(num_points, n_species)
    for traj in trajectories        
        interp = LinearInterpolation(traj.t, traj.u, extrapolation_bc=Linear())
        traj_mean += reduce(hcat,[interp(t) for t in uniform_time])' ./ length(trajectories)
    end
    
    # Return as vectors (mean prey and predator populations at each time)
    return uniform_time, traj_mean
end

export gillespie, mean_trajectory

