"""
    remap_sp_transition(sp_old, old_cmap, new_cmap, ::Val{new_K}; weight_tol) -> StateSpace

Remap an FSP state from one spatial CoarseningMap to another.

For each old coarse region r, its n_r molecules are redistributed to whichever
new regions absorb its fine voxels:

  * Full overlap (all of r → one new region): exact, no approximation.
  * Partial overlap (r splits among several new regions): sequential Binomial
    split with probabilities = fine-voxel overlap fractions. Exact under the
    diffusion-equilibrium assumption (molecules uniformly distributed within r).

Typical transitions for the 2D BD wavefront:
  - Front voxels settle → count folds into settled background (full overlap → exact).
  - Empty background gains new front voxels → background count splits; the
    background fraction going to the new front region is proportional to
    the number of voxels crossing the threshold (usually O(1) out of O(K²)).
"""
function remap_sp_transition(
    sp_old::StateSpace{CartesianIndex{old_K}, Float64},
    old_cmap::CoarseningMap,
    new_cmap::CoarseningMap,
    ::Val{new_K};
    weight_tol::Float64 = 1e-14
) where {old_K, new_K}
    old_cmap.n_coarse == old_K ||
        error("old_cmap.n_coarse=$(old_cmap.n_coarse) ≠ old_K=$old_K")
    new_cmap.n_coarse == new_K ||
        error("new_cmap.n_coarse=$(new_cmap.n_coarse) ≠ new_K=$new_K")
    old_cmap.n_fine == new_cmap.n_fine ||
        error("Fine grid mismatch: $(old_cmap.n_fine) ≠ $(new_cmap.n_fine)")

    # Build fine-voxel overlap counts: overlap[r_old][r_new] = # fine voxels shared
    K_fine = old_cmap.n_fine
    overlap = [Dict{Int,Int}() for _ in 1:old_K]
    for k in 1:K_fine
        r_old = old_cmap.fine_to_coarse[k]
        r_new = new_cmap.fine_to_coarse[k]
        overlap[r_old][r_new] = get(overlap[r_old], r_new, 0) + 1
    end

    # splits[r] = sorted [(new_r, fraction)] for old region r
    splits = Vector{Vector{Tuple{Int,Float64}}}(undef, old_K)
    for r in 1:old_K
        total = isempty(overlap[r]) ? 1 : sum(values(overlap[r]))
        splits[r] = sort!([(k, v / total) for (k, v) in overlap[r]]; by = first)
    end

    # Precompute log-factorial table up to the largest molecule count in sp_old
    n_max_total = isempty(sp_old.states) ? 0 :
                  maximum(sum(Tuple(s)) for s in sp_old.states)
    logfact = _sp_logfact(n_max_total)

    sp_new = StateSpace{CartesianIndex{new_K}, Float64}()
    n_buf  = zeros(Int, new_K)

    for i in eachindex(sp_old.states)
        p = sp_old.probs[i]
        p > weight_tol || continue
        n_old = Tuple(sp_old.states[i])
        fill!(n_buf, 0)
        _remap_add!(sp_new, splits, n_old, n_buf, 1, p, logfact, weight_tol)
    end

    renormalize!(sp_new)
    sp_new
end

# ── helpers ────────────────────────────────────────────────────────────────────

function _sp_logfact(n_max::Int)
    t = zeros(max(n_max, 0) + 1)
    for i in 1:length(t)-1
        t[i+1] = t[i] + log(Float64(i))
    end
    t
end

# Recurse over old regions; for each, distribute its n_r molecules to new regions.
function _remap_add!(
    sp_new::StateSpace{CartesianIndex{new_K}, Float64},
    splits::Vector{Vector{Tuple{Int,Float64}}},
    n_old::NTuple{old_K, Int},
    n_buf::Vector{Int},
    old_r::Int,
    weight::Float64,
    logfact::Vector{Float64},
    weight_tol::Float64
) where {old_K, new_K}
    if old_r > old_K
        s = CartesianIndex(NTuple{new_K, Int}(n_buf))
        idx = get(sp_new.index, s, 0)
        if idx == 0
            add_state!(sp_new, s, weight)
        else
            sp_new.probs[idx] += weight
        end
        return
    end

    n_r   = n_old[old_r]
    split = splits[old_r]

    if length(split) == 1
        # Entire region → one new region: exact transfer, no enumeration.
        r_new = split[1][1]
        n_buf[r_new] += n_r
        _remap_add!(sp_new, splits, n_old, n_buf, old_r + 1, weight, logfact, weight_tol)
        n_buf[r_new] -= n_r
    else
        # Sequential Binomial split: enumerate Multinomial outcomes.
        _remap_split!(sp_new, splits, n_old, n_buf, old_r, n_r, split, 1,
                      weight, logfact, weight_tol)
    end
end

# Sequential Binomial decomposition of Multinomial(n_rem, split[idx:end]).
# At each level, Binom(n_rem, k, p_this / p_remaining) gives the conditional
# probability of placing k molecules in split[idx], given the remaining pool.
function _remap_split!(
    sp_new::StateSpace{CartesianIndex{new_K}, Float64},
    splits::Vector{Vector{Tuple{Int,Float64}}},
    n_old::NTuple{old_K, Int},
    n_buf::Vector{Int},
    old_r::Int,
    n_rem::Int,
    split::Vector{Tuple{Int,Float64}},
    idx::Int,
    weight::Float64,
    logfact::Vector{Float64},
    weight_tol::Float64
) where {old_K, new_K}
    if idx == length(split)
        # Last destination: all remaining go here (no enumeration needed).
        r_new = split[idx][1]
        n_buf[r_new] += n_rem
        _remap_add!(sp_new, splits, n_old, n_buf, old_r + 1, weight, logfact, weight_tol)
        n_buf[r_new] -= n_rem
        return
    end

    r_new, f_this = split[idx]
    f_rest = sum(f for (_, f) in split[idx+1:end])   # prob mass for remaining destinations
    f_tot  = f_this + f_rest
    p_this = f_tot > 0.0 ? f_this / f_tot : 0.0      # conditional Binomial probability

    lf = logfact   # shorthand
    for k in 0:n_rem
        # Binom(n_rem, k, p_this)
        log_bw = lf[n_rem+1] - lf[k+1] - lf[n_rem-k+1]
        bw = if p_this <= 0.0
            k == 0 ? 1.0 : 0.0
        elseif p_this >= 1.0
            k == n_rem ? 1.0 : 0.0
        else
            exp(log_bw + k * log(p_this) + (n_rem - k) * log(1.0 - p_this))
        end
        w = weight * bw
        w > weight_tol || continue

        n_buf[r_new] += k
        _remap_split!(sp_new, splits, n_old, n_buf, old_r, n_rem - k, split, idx + 1,
                      w, logfact, weight_tol)
        n_buf[r_new] -= k
    end
end
