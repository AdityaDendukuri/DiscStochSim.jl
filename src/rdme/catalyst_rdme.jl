"""
    SpatialRDME

Wraps a Catalyst `ReactionSystem` with diffusion coefficients to create
an `AbstractSpatialRDMEModel` compatible with `SpatialFSP`.

Currently supports single-species networks with linear (birth + death) reactions.
Reactions are defined with Catalyst's `@reaction_network` macro; rates are
extracted automatically and diffusion is specified separately.

# Example
```julia
rn = @reaction_network begin
    k_b, ∅ → A
    k_d, A → ∅
end

model = SpatialRDME(rn;
    diffusion = (A = 1.0,),
    params    = (k_b = 1.0, k_d = 0.5),
    dx        = 1.0)

state = SpatialFSP(model, 100, 100; n_max=16, ε_expand=0.15, ε_equil=0.07)
set_ic!(state, CartesianIndex(51, 51), 1)
```
"""
struct SpatialRDME <: AbstractSpatialRDMEModel
    rn      :: ReactionSystem   # Catalyst network (kept for introspection)
    k_b     :: Float64          # total zeroth-order birth rate per voxel
    k_d     :: Float64          # total first-order death rate per molecule
    D       :: Float64          # diffusion coefficient for the (single) species
    dx      :: Float64          # voxel side length
end

"""
    SpatialRDME(rn; diffusion, params, dx=1.0)

Construct a `SpatialRDME` from a Catalyst `ReactionSystem`.

Arguments
---------
- `rn`        : Catalyst reaction network (single species, linear reactions)
- `diffusion` : NamedTuple or Dict mapping species symbol → diffusion coefficient D
- `params`    : NamedTuple mapping parameter names → numerical values
- `dx`        : voxel side length (default 1.0)

The constructor automatically identifies birth reactions (net stoich +1, zeroth-order)
and death reactions (net stoich −1, first-order) and extracts their rates.
"""
function SpatialRDME(rn::ReactionSystem;
                     diffusion,
                     params::NamedTuple,
                     dx::Float64 = 1.0)

    specs = Catalyst.get_species(rn)
    length(specs) == 1 ||
        error("SpatialRDME currently supports single-species networks only. " *
              "Got $(length(specs)) species: $(nameof.(specs))")

    # Resolve the diffusion coefficient for the single species
    species_sym = Symbol(Catalyst.getname(first(specs)))
    D = if diffusion isa NamedTuple
        haskey(diffusion, species_sym) ?
            Float64(diffusion[species_sym]) :
            error("Diffusion coefficient for :$species_sym not found in diffusion=$diffusion")
    else
        Float64(diffusion[species_sym])
    end

    # Parameter vector in the order Catalyst expects
    ps_syms   = Symbol.(Catalyst.getname.(Catalyst.get_ps(rn)))
    p_vals    = [Float64(params[s]) for s in ps_syms]

    # Build propensity functions via Catalyst
    sys = DiscreteStochasticSystem(rn)

    # Classify reactions by net stoichiometry and integrate rates
    k_b = 0.0
    k_d = 0.0
    for (prop, stoich) in zip(sys.propensities, sys.stoichvecs)
        net = first(Tuple(stoich))          # net change in the single species
        if net > 0
            # Birth-type: evaluate propensity at n=0 to get zeroth-order rate
            k_b += prop(CartesianIndex(0), p_vals, 0.0)
        elseif net < 0
            # Death-type: evaluate propensity at n=1; divide by n to get per-molecule rate
            k_d += prop(CartesianIndex(1), p_vals, 0.0)   # = k_d × 1 for mass-action
        end
    end

    k_b > 0 || k_d > 0 ||
        error("No birth or death reactions found in the network")

    SpatialRDME(rn, k_b, k_d, D, dx)
end

# ── AbstractSpatialRDMEModel interface ────────────────────────────────────────
jump_rate(m::SpatialRDME)              = m.D / m.dx^2
local_birth_rate(m::SpatialRDME, n)   = m.k_b
local_total_birth_rate(m::SpatialRDME, n, n_sites::Int=1) = m.k_b * n_sites
has_finite_local_ss(::SpatialRDME)    = true
ss_mean(m::SpatialRDME)               = m.k_d > 0 ? m.k_b / m.k_d : Inf
