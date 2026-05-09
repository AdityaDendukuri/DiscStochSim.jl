---
name: project_multigrid_method
description: Working multigrid approach, fixes applied, known limitations, and next steps
type: project
---

Working n-level multigrid V-cycle for Schlögl RDME.

## What's working

- `multi_level_vcycle` runs for Schlögl without state-space explosion
- `prolong_conditional` (not Binomial prolong) for Schlögl prolongation
- `_filter_positive` before restrict eliminates zero-coverage coarse states
- `prolong_conditional` step 3 now uses `haskey` check (not positional `n_pre+1`)
- State space grows at ~1000 states/step at the bistable front — manageable
- τ_pre > 0 (intra-pair diffusion pre-smooth) essential for within-pair equilibration
- Splitting demo (t=12, K=4, D=5): P(low) grows to 37%, P(high) to 10%

## Key fixes applied (multigrid-fsp branch)

### Fix 1: prolong_conditional step 3 ordering fix (operators.jl)
Old: `for i in (n_pre+1):length(sp_coarse_post)` — assumed expanded states at tail
New: `for i in 1:length ... haskey(sp_coarse_pre.index, x) && continue`
Why: for n_levels≥2, sp_coarse_post from recursive call has different ordering

### Fix 2: _filter_positive before restrict (vcycle.jl)
Prune zero-prob boundary fine states before restricting to get sp_coarse_pre.
Why: expand! adds p=0 boundary states → they restrict to zero-coverage coarse states
→ those coarse states receive prob after solve but prolong_conditional can't assign it back
(they're "covered" by sp_coarse_pre but have no positive-prob fine counterpart)
Result: eliminated n_zero_pre issue; mult mass now correctly equals covered mass

### Fix 3: prolong_conditional not prolong_multiplicative for Schlögl (vcycle.jl)
prolong_multiplicative drops ~30% mass (probability that moved to expanded coarse states)
prolong_conditional recovers 95-99% via nearest-neighbor conditional extension
State space grows bounded: ~1K per step at bistable front

### Fix 4: size guards (vcycle.jl)
Added max_states check after coarse expand! (not just at V-cycle entry)
Binomial prolong for RDME gets max_states guard passed through

## Structural limitations

**Voxel 1/3 asymmetry**: prolong_conditional keeps n1_j (first voxel of pair) fixed at
pre-smoothed value, shifts n2_j = nc_new - n1_j. When nc drops to low stable state (62),
n2 = 62 - 72 = -10 → skipped. Voxels 1 and 3 remain in unstable basin in simulations.
**Fix**: need dynamic-π (pi_table) prolongation to distribute nc symmetrically.

**expand_depth=0 required**: with expand_depth=1, fine state space grows 3× per step
at the bistable front → hits max_states quickly. Prolong_conditional alone drives ~1000
states/step which is acceptable.

**Concentrated IC**: V-cycle poorly handles IC with all molecules in one voxel.
Prolong_conditional can only extend by ±1 nc per pair per step → takes O(40-90) steps
to reach stable states. Use coarse_only mode for early time with concentrated IC.

## Scaling results

K=4, n_levels=1, bistable front IC (148,148,31,31):
- |S| stabilizes at ~1105 with expand_depth=0, prune_tol=1e-10
- Means stable: high pair stays at 148, low pair slowly increases

K=4, n_levels=1, unstable IC (72,72,72,72), D=5:
- |S| grows linearly: 15 → 115K at t=12
- P(low basin) voxels 2,4: 0% → 37% at t=12
- P(high basin) voxels 2,4: 0% → 10% at t=12
- Max states 500K hit at ~t=20 (needs prune_tol tuning for longer runs)

## State-space explosion root causes (all now fixed)

1. Binomial prolong of corrections: nc≈300 at bistable front → 26^4=457K states per coarse state
2. Zero-coverage coarse states: mass dropped silently (looked covered but wasn't)
3. Binomial prolong of corrections also had mass inflation (P(δ_neg) subtracted from non-existent fine states → zeroed → mass inflated to 1.97)

## Next development priorities

1. Dynamic-π prolongation wired into multi_level_vcycle (needs pi_table parameter)
   → fixes voxel 1/3 asymmetry, enables correct bistable splitting in all voxels
2. K=8 n_levels=2 test with bistable IC — state space unknown, needs testing
3. Comparison with direct fine FSP for accuracy validation (K=4 is feasible)
4. Coarse_only path (_solve_rdme_coarse_only) also needs prolong_conditional for Schlögl
   (currently uses Binomial _multi_prolong → explodes for large nc)
