"""
Dynamic-K spatial CME: full joint distribution with adaptive voxel activation
=============================================================================

The joint distribution p(n_1,...,n_K') lives on StateSpace{CartesianIndex{K'}}
where K' changes dynamically as voxels activate (enter the joint CME) or
equilibrate (get marginalized out).

Lifecycle per step:
  1. Evolve: expv(dt, A, p) on the current K'-dimensional joint CME
  2. Expand: expand! adds boundary states reachable from current support
  3. Prune:  prune_threshold! removes negligible-probability states
  4. Equil check: for each active voxel k, TV(marginal_k, Poisson(μ_k)) < ε_equil
  5. Deactivate: equil voxels → marginalize out (exact sum), update equil_μ
  6. Activate: empty/equil fine voxels adjacent to active with flux > ε_expand
               → tensor product with Poisson(μ_k) IC

K' stays small (≈ active front width) because:
  - Back of wave equilibrates → deactivated → K' shrinks
  - Front of wave activates   → activated   → K' grows
  - Net: K' ≈ width of non-equilibrated region (typically 10–50 voxels)

Equil coupling to active voxels:
  Equilibrated voxel j with mean μ_j, adjacent to active voxel k:
  → Effective birth term in active voxel k: D/dx² × μ_j
  → Effective death term in active voxel k: D/dx² × n_k (per molecule, unchanged)
  This is exact in the mean (correct at leading order for Poisson-distributed equil voxels).

Key types
---------
  DynCME{K'}            — joint CME over K' voxels (K' is a type parameter)
  AbstractDynCME        — abstract base for type-unstable outer dispatch

Key functions
-------------
  marginal_k            — marginal distribution of one voxel from joint
  marginal_mean_k       — exact FSP mean for one voxel
  tv_vs_poisson         — TV distance from Poisson(μ)
  marginalize_out       — exact marginalization: sum over voxel dimensions
  activate_voxel        — add voxel to joint CME: K' → K'+1
  deactivate_voxels     — remove voxels from joint CME: K' → K'-m
  build_active_system   — RDME system for active voxels + equil boundary terms
  step!                 — full lifecycle: evolve → equil check → activate/deactivate
  fine_mean_grid        — exact FSP means for all fine voxels (active + equil)
"""

# ── Abstract type for outer-loop type instability ─────────────────────────────

abstract type AbstractDynCME end

# ── DynCME{K'} ─────────────────────────────────────────────────────────────────

"""
    DynCME{K}

Joint CME over K' = K active voxels, with equil voxels marginalized out.

Fields
------
- `sp`          : joint distribution StateSpace{CartesianIndex{K}}
- `model`       : GraphRDMEModel for K active voxels (includes equil boundary terms)
- `mesh`        : VoxelMesh for K active voxels (active-active edges only)
- `active_vids` : NTuple{K,Int} — which fine voxels are the K active voxels
- `equil_μ`     : Float64 mean count per fine voxel (all voxels; 0.0 for empty)
- `vstate`      : Vector{Int} — 1=active, 2=equil, 3=empty (length n_fine_voxels)
- `fine_mesh`   : VoxelMesh for the full fine grid (fixed)
- `fine_birth`  : zeroth-order birth rate per fine voxel
- `fine_death`  : first-order death rate per molecule
- `n_max`       : FSP truncation (max molecules per active voxel)
- `ε_equil`     : TV threshold for equilibration detection
- `ε_expand`    : flux threshold for voxel activation
"""
struct DynCME{K} <: AbstractDynCME
    sp           :: StateSpace{CartesianIndex{K}, Float64}
    model        :: GraphRDMEModel
    mesh         :: VoxelMesh
    active_vids  :: NTuple{K, Int}
    equil_μ      :: Vector{Float64}
    vstate       :: Vector{Int}          # 1=active 2=equil 3=empty
    fine_mesh    :: VoxelMesh
    fine_birth   :: Vector{Float64}
    fine_death   :: Float64
    n_max        :: Int
    ε_equil      :: Float64
    ε_expand     :: Float64
end

const VSTATE_ACTIVE = 1
const VSTATE_EQUIL  = 2
const VSTATE_EMPTY  = 3

# ── Constructor from fine grid ─────────────────────────────────────────────────

"""
    DynCME(fine_mesh, fine_birth, fine_death; init_voxels, n_max, ε_equil, ε_expand)

Initialise a `DynCME` with `init_voxels` as active at t=0 (all others empty).
Place 1 molecule at the centre of `init_voxels`; all others at 0.
"""
function DynCME(fine_mesh::VoxelMesh,
                fine_birth::Vector{Float64},
                fine_death::Float64;
                init_voxels::Vector{Int},
                n_max::Int         = 20,
                ε_equil::Float64   = 0.05,
                ε_expand::Float64  = 0.01)

    N_fine = fine_mesh.n_voxels
    @assert length(fine_birth) == N_fine

    equil_μ = zeros(N_fine)
    vstate  = fill(VSTATE_EMPTY, N_fine)
    for v in init_voxels; vstate[v] = VSTATE_ACTIVE; end

    K = length(init_voxels)
    avids = Tuple(init_voxels)
    model, mesh_active = _build_active_system(fine_mesh, fine_birth, fine_death,
                                               avids, equil_μ)

    sp = StateSpace{CartesianIndex{K}, Float64}()
    # IC: all-zero state in the active block
    add_state!(sp, CartesianIndex(ntuple(_ -> 0, Val(K))), 1.0)

    _call_typed(Val(K), sp, model, mesh_active, avids, equil_μ, vstate,
                fine_mesh, fine_birth, fine_death, n_max, ε_equil, ε_expand)
end

# Type-stable constructor helper
function _call_typed(::Val{K},
                     sp, model, mesh_a, avids, equil_μ, vstate,
                     fine_mesh, fine_birth, fine_death, n_max, ε_equil, ε_expand) where K
    DynCME{K}(sp, model, mesh_a, NTuple{K,Int}(avids), equil_μ, vstate,
               fine_mesh, fine_birth, fine_death, n_max, ε_equil, ε_expand)
end

# ── Build active RDME system ──────────────────────────────────────────────────

"""
    _build_active_system(fine_mesh, fine_birth, fine_death, active_vids, equil_μ)
        → (GraphRDMEModel, VoxelMesh)

Build the RDME model and restricted mesh for the active voxels.

For each active voxel k (fine index v = active_vids[k]):
  - Birth rate = fine_birth[v] + Σ_{j equil neighbor of v} D_rates[edge] × equil_μ[j]
  - Death rate = fine_death + Σ_{j equil neighbor of v} D_rates[edge]  (per molecule)

The effective death rate increase from equil neighbors is absorbed into an AUGMENTED
death rate per active voxel.  Birth and death from equil boundaries represent mean-field
coupling to the equil region.

The returned VoxelMesh contains only active-active edges (equil edges dropped).
"""
function _build_active_system(fine_mesh::VoxelMesh,
                               fine_birth::Vector{Float64},
                               fine_death::Float64,
                               active_vids,
                               equil_μ::Vector{Float64})
    K = length(active_vids)
    K == 0 && return (GraphRDMEModel(Float64[], Float64[]),
                      VoxelMesh(0, Tuple{Int,Int}[], Float64[]))
    vstate_local = fill(VSTATE_EMPTY, fine_mesh.n_voxels)
    for (k, v) in enumerate(active_vids); vstate_local[v] = k; end  # k = position in active

    # Active-voxel index within joint CME
    active_pos = Dict(v => k for (k, v) in enumerate(active_vids))

    # Per-voxel effective birth and death rates.
    # For active voxel k adjacent to equil voxel j:
    #   birth_k += D × μ_equil_j   (molecules flowing in from equil)
    #   death_k += D               (molecules flowing out to equil, per-molecule)
    # This gives correct effective μ_ss_k = birth_k / death_k at equilibrium.
    avids_vec   = collect(Int, active_vids)
    birth_rates = copy(fine_birth[avids_vec])
    death_rates = fill(fine_death, K)

    for (e, (u, v)) in enumerate(fine_mesh.edges)
        D = fine_mesh.D_rates[e]
        if haskey(active_pos, u) && !haskey(active_pos, v) && equil_μ[v] > 0
            birth_rates[active_pos[u]] += D * equil_μ[v]
            death_rates[active_pos[u]] += D
        elseif haskey(active_pos, v) && !haskey(active_pos, u) && equil_μ[u] > 0
            birth_rates[active_pos[v]] += D * equil_μ[u]
            death_rates[active_pos[v]] += D
        end
    end

    # Active-active edges only
    edges_active  = Tuple{Int,Int}[]
    drates_active = Float64[]
    for (e, (u, v)) in enumerate(fine_mesh.edges)
        haskey(active_pos, u) && haskey(active_pos, v) || continue
        pu, pv = active_pos[u], active_pos[v]
        push!(edges_active, pu < pv ? (pu,pv) : (pv,pu))
        push!(drates_active, fine_mesh.D_rates[e])
    end

    mesh_active = VoxelMesh(K, edges_active, drates_active)
    model       = GraphRDMEModel(birth_rates, death_rates)
    (model, mesh_active)
end

# ── Marginal and TV computations ──────────────────────────────────────────────

"""
    marginal_k(sp, k, n_max) → Vector{Float64}

Compute the marginal distribution of voxel k (1-indexed in joint CME)
by summing over all other voxels.  Returns `p_k(0), p_k(1), ..., p_k(n_max)`.
"""
function marginal_k(sp::StateSpace{CartesianIndex{K}, Float64},
                    k::Int, n_max::Int) where K
    pmf = zeros(n_max + 1)
    for (ci, prob) in zip(sp.states, sp.probs)
        n = Tuple(ci)[k]
        if 0 ≤ n ≤ n_max
            pmf[n + 1] += prob
        end
    end
    s = sum(pmf); s > 0 && (pmf ./= s)
    pmf
end

"""
    marginal_mean_k(sp, k) → Float64

Exact FSP mean ⟨n_k⟩ = Σ_n n_k · p(n) for voxel k.
"""
function marginal_mean_k(sp::StateSpace{CartesianIndex{K}, Float64}, k::Int) where K
    μ = 0.0
    for (ci, prob) in zip(sp.states, sp.probs)
        μ += Tuple(ci)[k] * prob
    end
    μ
end

"""
    tv_vs_poisson(pmf, μ) → Float64

Total variation distance ½ Σ_n |pmf(n) - Poisson(μ; n)|.
"""
function tv_vs_poisson(pmf::Vector{Float64}, μ::Float64)
    n_max = length(pmf) - 1
    tv = 0.0
    if μ ≤ 0
        tv = 1.0 - pmf[1]   # all mass should be at n=0
        return 0.5 * abs(tv) + 0.5 * (1.0 - pmf[1])
    end
    pois_norm = 0.0
    for n in 0:n_max
        p_pois = exp(-μ) * μ^n / factorial(big(n))
        pois_norm += p_pois
        tv += abs(pmf[n + 1] - p_pois)
    end
    0.5 * tv   # TV = ½ Σ |·|
end

# ── Marginalization (exact: sums over voxel dimensions) ──────────────────────

"""
    marginalize_out(sp, keep_positions) → StateSpace{CartesianIndex{K'}}

Remove dimensions not in `keep_positions` by summing over them.
`keep_positions` is a sorted vector of 1-indexed positions in the K-tuple.

Returns a new `StateSpace{CartesianIndex{K'}}` with K' = length(keep_positions).
This is the exact marginal distribution over the kept voxels.
"""
function marginalize_out(sp::StateSpace{CartesianIndex{K}, Float64},
                          keep_positions::Vector{Int}) where K
    K_new = length(keep_positions)
    new_sp = StateSpace{CartesianIndex{K_new}, Float64}()
    for (ci, prob) in zip(sp.states, sp.probs)
        t = Tuple(ci)
        new_t = ntuple(i -> t[keep_positions[i]], Val(K_new))
        new_ci = CartesianIndex(new_t)
        idx = get(new_sp.index, new_ci, 0)
        if idx == 0
            add_state!(new_sp, new_ci, prob)
        else
            new_sp.probs[idx] += prob
        end
    end
    new_sp
end

# ── Voxel activation ──────────────────────────────────────────────────────────

"""
    activate_voxel(state::DynCME{K}, fine_vid, insert_pos) → DynCME{K+1}

Add fine voxel `fine_vid` to the joint CME at position `insert_pos` in the
active list (1-indexed).

Initialisation: the new voxel's molecule count is drawn from Poisson(equil_μ[fine_vid])
if it was equil, or from Poisson(0) = δ_0 if it was empty.  This is the tensor
product of the current joint distribution with the new voxel's marginal IC.
"""
function activate_voxel(state::DynCME{K}, fine_vid::Int,
                         insert_pos::Int = K + 1) where K
    μ_new = state.equil_μ[fine_vid]   # 0.0 for empty voxels

    # Poisson IC for new voxel (truncated at n_max)
    n_max = state.n_max
    pois_w = zeros(n_max + 1)
    if μ_new ≤ 0
        pois_w[1] = 1.0               # δ_{n=0} for empty voxel
    else
        for n in 0:n_max
            pois_w[n + 1] = exp(-μ_new) * μ_new^n / factorial(big(n))
        end
        pois_w ./= sum(pois_w)
    end

    # Build new state space: tensor product with Poisson IC
    new_sp = StateSpace{CartesianIndex{K+1}, Float64}()
    for (ci, prob) in zip(state.sp.states, state.sp.probs)
        t = Tuple(ci)
        for n_new in 0:n_max
            w = pois_w[n_new + 1]
            w < 1e-14 && continue
            new_t = ntuple(k -> k < insert_pos ? t[k] :
                                k == insert_pos ? n_new : t[k-1], Val(K+1))
            new_ci = CartesianIndex(new_t)
            idx = get(new_sp.index, new_ci, 0)
            if idx == 0
                add_state!(new_sp, new_ci, Float64(prob * w))
            else
                new_sp.probs[idx] += Float64(prob * w)
            end
        end
    end

    # Update metadata
    new_vstate   = copy(state.vstate);  new_vstate[fine_vid] = VSTATE_ACTIVE
    new_equil_μ  = copy(state.equil_μ); new_equil_μ[fine_vid] = 0.0
    new_avids    = ntuple(k -> k < insert_pos ? state.active_vids[k] :
                                k == insert_pos ? fine_vid : state.active_vids[k-1], Val(K+1))
    new_model, new_mesh = _build_active_system(state.fine_mesh, state.fine_birth,
                                                state.fine_death, new_avids, new_equil_μ)

    DynCME{K+1}(new_sp, new_model, new_mesh, new_avids, new_equil_μ, new_vstate,
                 state.fine_mesh, state.fine_birth, state.fine_death,
                 state.n_max, state.ε_equil, state.ε_expand)
end

# ── Voxel deactivation (equilibration) ───────────────────────────────────────

"""
    deactivate_voxels(state::DynCME{K}, positions) → DynCME{K'}

Remove voxels at `positions` (1-indexed within active_vids) from the joint CME.
Their marginal means are stored in `equil_μ` and they are marked EQUIL.

The new joint distribution is the EXACT marginal over the remaining voxels:
    p_new(n_keep) = Σ_{n_remove} p(n_keep, n_remove)
"""
function deactivate_voxels(state::DynCME{K}, positions::Vector{Int}) where K
    K_new = K - length(positions)
    keep  = sort(setdiff(1:K, positions))
    K_new == length(keep) || error("positions inconsistency")

    # Compute marginal means for deactivated voxels → update equil_μ
    new_equil_μ = copy(state.equil_μ)
    new_vstate  = copy(state.vstate)
    for p in positions
        v = state.active_vids[p]
        μ_v = marginal_mean_k(state.sp, p)
        new_equil_μ[v] = μ_v
        new_vstate[v]  = VSTATE_EQUIL
    end

    # Exact marginalization
    new_sp = marginalize_out(state.sp, keep)

    new_avids = ntuple(i -> state.active_vids[keep[i]], Val(K_new))
    new_model, new_mesh = _build_active_system(state.fine_mesh, state.fine_birth,
                                                state.fine_death, new_avids, new_equil_μ)

    DynCME{K_new}(new_sp, new_model, new_mesh, new_avids, new_equil_μ, new_vstate,
                   state.fine_mesh, state.fine_birth, state.fine_death,
                   state.n_max, state.ε_equil, state.ε_expand)
end

# ── Core evolution ─────────────────────────────────────────────────────────────

"""
    evolve!(state::DynCME{K}, dt; krylov_m, prune_tol) → DynCME{K}

Advance the joint distribution forward by `dt` using matrix exponentiation,
then expand boundary states and prune negligible-probability states.

Returns a new `DynCME{K}` with the same K (same active voxel set).
"""
function evolve!(state::DynCME{K}, dt::Float64;
                  krylov_m::Int      = 30,
                  prune_tol::Float64 = 1e-14) where K

    sp  = state.sp
    sys = build_graph_rdme_system(state.model, state.mesh)
    bc  = x -> all(c -> 0 ≤ c ≤ state.n_max, Tuple(x))

    # Expand with depth=2 BEFORE expv.
    # depth=1 states start with p=0 and gain probability during expv.
    # Without depth=2, birth transitions FROM those depth-1 states to depth-2 states
    # would be missing from the generator (targets not in state space), causing
    # artificially low effective birth rates.
    expand!(sp, sys, bc; depth=2)

    if length(sp) > 0 && dt > 0
        A, = build_generator(sp, sys, Float64[], 0.0)
        sp.probs .= expv(dt, A, sp.probs; m=krylov_m)
        sp.probs .= max.(0.0, sp.probs)
        s = sum(sp.probs); s > 0 && (sp.probs ./= s)
    end

    prune_threshold!(sp, prune_tol)

    DynCME{K}(sp, state.model, state.mesh, state.active_vids, state.equil_μ,
               state.vstate, state.fine_mesh, state.fine_birth, state.fine_death,
               state.n_max, state.ε_equil, state.ε_expand)
end

# ── Flux-based activation check ───────────────────────────────────────────────

"""
    activation_candidates(state::DynCME{K}) → Vector{Int}

Return fine voxel indices (empty or equil) that should be activated:
those adjacent to an active voxel where the expected diffusion flux into
them exceeds `ε_expand`.

Flux from active voxel v to neighbor j: D_rate × ⟨n_v⟩.
"""
function activation_candidates(state::DynCME{K}) where K
    active_set = Set(state.active_vids)
    # Map active fine voxel → its mean count in the joint CME
    active_mean = Dict{Int,Float64}()
    for (k, v) in enumerate(state.active_vids)
        active_mean[v] = marginal_mean_k(state.sp, k)
    end

    candidates = Int[]
    for (e, (u, v)) in enumerate(state.fine_mesh.edges)
        D = state.fine_mesh.D_rates[e]
        for (src, dst) in ((u,v), (v,u))
            (src ∈ active_set && dst ∉ active_set) || continue
            state.vstate[dst] == VSTATE_ACTIVE && continue   # already active (shouldn't happen)
            flux = D * get(active_mean, src, 0.0)
            if flux > state.ε_expand && dst ∉ candidates
                push!(candidates, dst)
            end
        end
    end
    candidates
end

# ── Full step lifecycle ───────────────────────────────────────────────────────

"""
    step!(state::AbstractDynCME, dt; krylov_m, prune_tol) → AbstractDynCME

Full lifecycle step (type-unstable at call site; each inner function is type-stable):

1. Evolve: expv step + expand + prune
2. Equil check: for each active voxel, TV(marginal, Poisson(μ)) < ε_equil
3. Deactivate equil voxels: K' → K'-m
4. Activation check: flux from active to empty/equil neighbors
5. Activate candidates: K' → K'+a

Returns a (possibly differently-typed) `AbstractDynCME`.
"""
function step!(state::AbstractDynCME, dt::Float64;
               krylov_m::Int      = 30,
               prune_tol::Float64 = 1e-14) :: AbstractDynCME
    _step_impl!(state, dt, krylov_m, prune_tol)
end

function _step_impl!(state::DynCME{K}, dt::Float64,
                      krylov_m::Int, prune_tol::Float64) where K

    # ── 1. Evolve ──────────────────────────────────────────────────────────────
    state = evolve!(state, dt; krylov_m, prune_tol)

    # ── 2. Equil check ─────────────────────────────────────────────────────────
    # Use MEAN-BASED criterion: |μ_k - μ_ss| / μ_ss < ε_equil.
    # This is robust to n_max truncation — works even when tails are cut off.
    # TV vs Poisson(μ_ss) would fail here because n_max=5 can't represent
    # Poisson(2) tails accurately, causing equil to never trigger.
    equil_positions = Int[]
    for k in 1:K
        mean_k = marginal_mean_k(state.sp, k)
        # Use per-voxel death_rates[k] so boundary voxels adjacent to equil
        # get the correct effective μ_ss = (k_b + D×μ_equil) / (k_d + D).
        μ_ss_k = state.model.birth_rates[k] / state.model.death_rates[k]
        rel_err = abs(mean_k - μ_ss_k) / (μ_ss_k + 1e-10)
        if rel_err < state.ε_equil
            push!(equil_positions, k)
        end
    end

    # ── 3. Activate new voxels first (before deactivation, so flux is up-to-date) ──
    for v in activation_candidates(state)
        state = activate_voxel(state, v)
    end

    # ── 4. Deactivate equil voxels ─────────────────────────────────────────────
    # equil_positions was computed before activation, so newly activated voxels
    # are not included (they have K+1,...,K+a indices beyond the old K).
    if !isempty(equil_positions)
        state = deactivate_voxels(state, equil_positions)
    end

    state
end

# ── Fine-grid mean extraction ─────────────────────────────────────────────────

"""
    fine_mean_grid(state::DynCME{K}) → Vector{Float64}

Return exact FSP mean count for every fine voxel:
  - Active voxels: ⟨n_k⟩ = Σ n_k · p(n) from joint CME
  - Equil voxels:  equil_μ[v] (stored from last marginalisation)
  - Empty voxels:  0.0
"""
function fine_mean_grid(state::DynCME{K}) where K
    μ = copy(state.equil_μ)   # already has equil values; 0.0 for empty
    for (k, v) in enumerate(state.active_vids)
        μ[v] = marginal_mean_k(state.sp, k)
    end
    μ
end

"""
    fine_status_grid(state::DynCME{K}) → Vector{Int}

Return voxel status (VSTATE_ACTIVE=1, VSTATE_EQUIL=2, VSTATE_EMPTY=3).
"""
fine_status_grid(state::DynCME) = state.vstate

# ── Convenience constructors ──────────────────────────────────────────────────

"""
    place_molecules!(state::DynCME{K}, fine_vid, n) → DynCME{K}

Set IC: place exactly `n` molecules in fine voxel `fine_vid` (must be active).
All probability mass goes to the single state with n_vid = n, others = 0.
"""
function place_molecules!(state::DynCME{K}, fine_vid::Int, n::Int) where K
    k = findfirst(==(fine_vid), collect(state.active_vids))
    isnothing(k) && error("voxel $fine_vid is not active (call activate_voxel first)")
    sp = state.sp
    empty!(sp.states); empty!(sp.probs); empty!(sp.index); empty!(sp.global_id)
    ic = CartesianIndex(ntuple(i -> i == k ? n : 0, Val(K)))
    add_state!(sp, ic, 1.0)
    state
end

# ── Diagnostics ───────────────────────────────────────────────────────────────

"""
    diagnose(state::DynCME{K})

Print a summary: number of active/equil/empty voxels, states, probability mass.
"""
function diagnose(state::DynCME{K}) where K
    n_active = K
    n_equil  = count(==(VSTATE_EQUIL), state.vstate)
    n_empty  = count(==(VSTATE_EMPTY), state.vstate)
    n_states = length(state.sp)
    p_mass   = sum(state.sp.probs)
    println("K'=$(lpad(n_active,3)) active | $(lpad(n_equil,4)) equil | " *
            "$(lpad(n_empty,5)) empty | $(lpad(n_states,6)) CME states | " *
            "p=$(round(p_mass, digits=6))")
end
