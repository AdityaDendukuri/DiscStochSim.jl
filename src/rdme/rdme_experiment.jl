"""
RDME Experiment API
===================

Declarative interface for setting up spatial stochastic simulations.

# IC specification
    ZeroIC()                          # all voxels empty
    CenterIC(:A => 1)                 # one A at grid centre
    CenterIC(:A => 1, :B => 3)        # multi-species at centre
    VoxelIC((51,51) => (:A=>1,))      # specific voxel(s)
    VoxelIC((10,10)=>(:A=>5,), (90,90)=>(:B=>5,))

# Grids
    rect_grid(K)           # K×K rectangle  →  (K,K)
    rect_grid(Kx, Ky)      # Kx×Ky rectangle → (Kx,Ky)
    (any VoxelGraph from graph_rdme.jl)

# Experiment builder
    exp = RDMEExperiment(rn;
        grid      = rect_grid(100),
        diffusion = (A=1.0,),
        params    = (k_b=1.0, k_d=0.5),
        ic        = CenterIC(:A=>1))

# Running
    state = build_state(exp; n_max=16, ε_expand=0.15, ε_equil=0.07)
    step!(state, 0.3)
"""

# ── IC Specification ──────────────────────────────────────────────────────────

abstract type AbstractICSpec end

"""All voxels start empty (zero molecules of every species)."""
struct ZeroIC <: AbstractICSpec end

"""
Place molecules of one or more species at the geometric centre of the grid.

    CenterIC(:A => 1)
    CenterIC(:A => 1, :B => 3)
"""
struct CenterIC <: AbstractICSpec
    counts :: Dict{Symbol,Int}      # species → count
    CenterIC(pairs...) = new(Dict(Symbol(k) => v for (k,v) in pairs))
end

"""
Place molecules at specific (row,col) voxels.

    VoxelIC((51,51) => (:A=>1,))
    VoxelIC((10,10)=>(:A=>5,:B=>0), (90,90)=>(:B=>5,:A=>0))

Each key is a (row,col) tuple; each value is a NamedTuple or Tuple of Pairs
giving species=>count.
"""
struct VoxelIC <: AbstractICSpec
    assignments :: Vector{Pair{Tuple{Int,Int}, Dict{Symbol,Int}}}
    function VoxelIC(pairs...)
        assignments = map(pairs) do (pos, sc)
            Tuple(pos) => Dict(Symbol(k) => v for (k,v) in pairs(sc))
        end
        new(collect(assignments))
    end
end

# Convenience: single-voxel shorthand
VoxelIC(pos::Tuple{Int,Int}, sc) = VoxelIC(pos => sc)

"""
Place N molecules of one species randomly across all voxels.

    RandomIC(:A, 10; seed=42)
"""
struct RandomIC <: AbstractICSpec
    species :: Symbol
    count   :: Int
    seed    :: Int
    RandomIC(species, count; seed=0) = new(Symbol(species), count, seed)
end

# ── Grid helpers ──────────────────────────────────────────────────────────────

"""
    rect_grid(K)        →  (K, K)
    rect_grid(Kx, Ky)   →  (Kx, Ky)

Shorthand for a rectangular grid specification.  Pass the result as the
`grid` argument of `RDMEExperiment`.
"""
rect_grid(K::Int) = (K, K)
rect_grid(Kx::Int, Ky::Int) = (Kx, Ky)

# ── RDMEExperiment ────────────────────────────────────────────────────────────

"""
    RDMEExperiment

Declarative specification for a spatial stochastic simulation on a rectangular
or graph-based grid.

Fields
------
- `network`   : Catalyst `ReactionSystem`
- `grid`      : `(Kx, Ky)` for rectangular, or a `VoxelGraph` for arbitrary topology
- `diffusion` : `NamedTuple` mapping species symbol → diffusion coefficient D
- `params`    : `NamedTuple` mapping parameter name → numerical value
- `ic`        : Initial condition specification (see `ZeroIC`, `CenterIC`, `VoxelIC`)

For spatially heterogeneous birth rates, supply `birth_region` as a
`Dict{Symbol,Function}` where the function `f(i,j)::Bool` returns true for
voxels that produce the given species.

# Example
```julia
rn = @reaction_network begin
    k_b, ∅ → A
    k_d, A → ∅
end

exp = RDMEExperiment(rn;
    grid      = rect_grid(100),
    diffusion = (A = 1.0,),
    params    = (k_b = 1.0, k_d = 0.5),
    ic        = CenterIC(:A => 1))

state = build_state(exp; n_max=16, ε_expand=0.15, ε_equil=0.07)
```
"""
struct RDMEExperiment{G}
    network       :: ReactionSystem
    grid          :: G                          # (Kx,Ky) or VoxelGraph
    diffusion     :: NamedTuple
    params        :: NamedTuple
    ic            :: AbstractICSpec
    birth_region  :: Union{Nothing, Dict{Symbol,Function}}  # spatially heterogeneous birth
end

function RDMEExperiment(rn::ReactionSystem;
                        grid,
                        diffusion,
                        params::NamedTuple,
                        ic::AbstractICSpec = ZeroIC(),
                        birth_region = nothing)
    diffusion_nt = diffusion isa NamedTuple ? diffusion :
                   NamedTuple{Tuple(keys(diffusion))}(Tuple(values(diffusion)))
    RDMEExperiment{typeof(grid)}(rn, grid, diffusion_nt, params, ic, birth_region)
end

# ── build_state: rectangular grid, single species ─────────────────────────────

"""
    build_state(exp::RDMEExperiment; kwargs...) → SpatialFSP

Construct and initialise a `SpatialFSP` state from an `RDMEExperiment`.
All keyword arguments are forwarded to `SpatialFSP(model, Kx, Ky; ...)`.

Handles:
- Homogeneous birth/death (SpatialRDME → SpatialFSP)
- Spatially heterogeneous birth (birth_region mask)
- Intuitive IC types (ZeroIC, CenterIC, VoxelIC, RandomIC)
"""
function build_state(exp::RDMEExperiment{Tuple{Int,Int}};
                     n_max::Int         = 40,
                     ε_expand::Float64  = 0.3,
                     ε_equil::Float64   = 0.06,
                     use_block_vcycle::Bool = true,
                     kwargs...)

    Kx, Ky = exp.grid

    # Build model: heterogeneous rates require a masked birth model
    model = if isnothing(exp.birth_region)
        SpatialRDME(exp.network;
                    diffusion = exp.diffusion,
                    params    = exp.params)
    else
        _masked_rdme(exp)
    end

    state = SpatialFSP(model, Kx, Ky;
                       n_max, ε_expand, ε_equil,
                       use_block_vcycle, kwargs...)

    _apply_ic!(state, exp.ic, Kx, Ky)
    state
end

# Apply IC ────────────────────────────────────────────────────────────────────

function _apply_ic!(state::SpatialFSP, ic::ZeroIC, Kx, Ky)
    # No IC — default already zero
end

function _apply_ic!(state::SpatialFSP, ic::CenterIC, Kx, Ky)
    # CenterIC only makes sense for single-species SpatialFSP (uses first species)
    isempty(ic.counts) && return
    # SpatialFSP.set_ic! takes a voxel index and count
    cx, cy = Kx÷2+1, Ky÷2+1
    # For single-species SpatialFSP, just use the first (and only expected) species count
    count = first(values(ic.counts))
    set_ic!(state, CartesianIndex(cx, cy), count)
end

function _apply_ic!(state::SpatialFSP, ic::VoxelIC, Kx, Ky)
    for ((r, c), species_counts) in ic.assignments
        count = first(values(species_counts))   # single-species: take first
        set_ic!(state, CartesianIndex(r, c), count)
    end
end

function _apply_ic!(state::SpatialFSP, ic::RandomIC, Kx, Ky)
    rng = ic.seed > 0 ? Random.MersenneTwister(ic.seed) : Random.default_rng()
    for _ in 1:ic.count
        r = rand(rng, 1:Kx)
        c = rand(rng, 1:Ky)
        # Add one molecule; SpatialFSP may not support cumulative set_ic
        set_ic!(state, CartesianIndex(r, c), 1)
    end
end

# Masked RDME for heterogeneous birth ────────────────────────────────────────

struct MaskedRDME <: AbstractSpatialRDMEModel
    base   :: SpatialRDME
    mask   :: Function      # (i,j) → Bool: true if birth allowed here
end

jump_rate(m::MaskedRDME)                     = jump_rate(m.base)
local_birth_rate(m::MaskedRDME, n)            = local_birth_rate(m.base, n)
local_total_birth_rate(m::MaskedRDME, n, ns) = local_total_birth_rate(m.base, n, ns)
has_finite_local_ss(m::MaskedRDME)            = has_finite_local_ss(m.base)
ss_mean(m::MaskedRDME)                        = ss_mean(m.base)

# k_d forwarded for SpatialFSP internals
Base.getproperty(m::MaskedRDME, s::Symbol) =
    s === :k_d ? m.base.k_d : getfield(m, s)

function _masked_rdme(exp::RDMEExperiment)
    base = SpatialRDME(exp.network; diffusion=exp.diffusion, params=exp.params)
    # Use the first species' birth region (single-species for now)
    specs = Catalyst.get_species(exp.network)
    sym   = Symbol(Catalyst.getname(first(specs)))
    mask_fn = get(exp.birth_region, sym, (i,j)->true)
    MaskedRDME(base, mask_fn)
end

# ── Pretty printing ───────────────────────────────────────────────────────────

function Base.show(io::IO, exp::RDMEExperiment)
    specs = join(string.(Symbol.(Catalyst.getname.(Catalyst.get_species(exp.network)))), ", ")
    rxs   = length(Catalyst.get_rxs(exp.network))
    grid_str = exp.grid isa Tuple ? "$(exp.grid[1])×$(exp.grid[2]) rectangle" :
                                    "VoxelGraph($(exp.grid.n_voxels) voxels)"
    print(io, "RDMEExperiment{$specs}  grid=$grid_str  rxs=$rxs  ic=$(typeof(exp.ic))")
end

function Base.show(io::IO, m::SpatialRDME)
    rxs = length(Catalyst.get_rxs(m.rn))
    print(io, "SpatialRDME(rxs=$rxs, k_b=$(m.k_b), k_d=$(m.k_d), D=$(m.D))")
end
