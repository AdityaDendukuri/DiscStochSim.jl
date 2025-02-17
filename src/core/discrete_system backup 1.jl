using Catalyst

struct DiscreteStochasticSystem{ElementType}
  stoichvecs::Array{ElementType}
  propensities::Array{Function}
end

"""
    DiscreteStochasticSystem(reaction_system::ReactionSystem)

TBW
"""
function DiscreteStochasticSystem(reaction_system::ReactionSystem)
  X = Catalyst.get_species(reaction_system)
  P = Catalyst.get_ps(reaction_system)

  # for each reaction, define a forward reaction
  stoichvecs = map(reaction_system |> netstoichmat |> eachcol) do stoichvec
    CartesianIndex(stoichvec...)
  end

  # for each reaction, define a propensity function
  propensities = map(reaction_system |> Catalyst.get_rxs) do reaction
    eqn = jumpratelaw(reaction, combinatoric_ratelaw=true)
    build_function(eqn, X, P, expression=Val{false})
  end

  DiscreteStochasticSystem{CartesianIndex}(stoichvecs, propensities)
end

# getters 
stoichvecs(system::DiscreteStochasticSystem) = system.stoichvecs
propensities(system::DiscreteStochasticSystem) = system.propensities

export DiscreteStochasticSystem
