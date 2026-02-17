"""
    FSPProblem(rn, u0, tspan, rates; bounds=(0, 50_000))

Defines a Chemical Master Equation problem for the Finite State Projection method.

# Arguments
- `rn::ReactionSystem`: Catalyst reaction network.
- `u0::CartesianIndex`: Initial state.
- `tspan::Tuple{Real,Real}`: Time interval `(t0, tf)`.
- `rates::Vector`: Reaction rate parameters.
- `bounds`: Domain bounds passed to `RectLatticeBoundaryCondition`.
"""
struct FSPProblem{E, T<:Real}
    model::DiscreteStochasticSystem{E}
    u0::E
    tspan::Tuple{T,T}
    rates::Vector{T}
    bc::Function
end

function FSPProblem(rn::ReactionSystem, u0::CartesianIndex, tspan::Tuple,
                    rates::Vector; bounds=(0, 50_000))
    model = DiscreteStochasticSystem(rn)
    T = promote_type(eltype(rates), typeof(tspan[1]), typeof(tspan[2]))
    tspan_T = (T(tspan[1]), T(tspan[2]))
    rates_T = T.(rates)
    bc = x -> RectLatticeBoundaryCondition(x, bounds)
    FSPProblem{CartesianIndex, T}(model, u0, tspan_T, rates_T, bc)
end
