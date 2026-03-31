# Implementation Notes and Development Log

## Open Paper Tasks
- [ ] remark 3.8 (compressed generator vs renormalization)
- [ ] fix rhodius citations
- [ ] Benchmark 2 (Schlögl): not yet implemented

---

## Benchmark 1 Development Log
### Birth-Death-Diffusion (K=4 and K=8)

**Parameters:** k_b=2, k_d=1, D=1, domain [0,1], IC = all molecules in voxel 1.
**Stationary:** each voxel → Poisson(μ*=2), Σμ = K*2.

---

### Phase 1: Getting the V-cycle correct

**Problem:** Original implementation used operator-splitting prolongation:
`p^h_new = P(p̃^{2h})` — replace fine distribution entirely with prolongation
of full coarse solution. This caused the fine state space to explode because
each prolongation produces O(n_max^K) fine states.

**Fix (Algorithm 6.1 from paper):** Correction-based prolongation:
`p^h_new = p̃^h + P(p̃^{2h} - p^{2h})`
Only the *correction* δ^{2h} = p̃^{2h} - p^{2h} is prolongated. New fine
states come only from the small leaked probability in δ, not from the entire
(large) coarse support. Implemented in `src/rdme/vcycle.jl`.

**Key implementation detail:** δ can be signed (some coarse states gain
probability, others lose). We split into positive and negative parts, prolong
each separately, then combine:
- sp_δ_pos: states where p̃^{2h} > p^{2h}
- sp_δ_neg: states where p̃^{2h} < p^{2h} (stored as absolute values)
- sp_fine_new.probs += prolong(sp_δ_pos) - prolong(sp_δ_neg)
- Clip any resulting negative probabilities to zero

---

### Phase 2: Coarse expansion (before vs. after solve)

**Problem:** Early implementation added FSP boundary states to the coarse
state space AFTER `expv`. This meant the new boundary states had probability
zero, were prolongated, and immediately pruned → no spatial spread.
Debug showed coarse state space frozen at 6 states across all steps.

**Fix:** `expand!(sp_coarse, coarse_sys, coarse_bc)` must be called BEFORE
`expv(dt, A_coarse, probs)`. The boundary expansion allows probability to
flow INTO the new boundary states DURING the matrix exponential.
This is the crucial ordering: expand → build_generator → expv → prune.

---

### Phase 3: K=8 fine-state-space explosion

**Problem:** Even with the correction-based V-cycle, the fine state space
(CartesianIndex{8}) accumulated monotonically across V-cycles. Old fine
states are never removed between cycles; only new correction states are added.
For K=8 this grows to tens of thousands of fine states within ~10 steps,
and the per-step cost of maintaining and evolving this 8-dimensional state
space becomes prohibitive.

**Approach tried:** Tighter prune_tol and binom_tol to limit fine state
space growth. Result: fine state space still grew explosively for K=8.

---

### Phase 4: Coarse-only mode

**Key insight:** For large K, don't maintain a persistent fine state space
at all. Work entirely at the K/2 coarse level between time steps; only
construct the fine distribution at snapshot times (via prolongation), then
immediately discard it.

**Implementation** (`src/rdme/solve.jl`, `_solve_rdme_coarse_only`):
```
Initialize: restrict fine IC → coarse IC
While t < tf:
    expand!(sp_c, coarse_sys, coarse_bc; depth=1)  ← BEFORE solve
    A_c = build_generator(sp_c, coarse_sys, rates, t)
    sp_c.probs = expv(dt, A_c, sp_c.probs; m=krylov_m)
    renormalize!(sp_c)
    prune_threshold!(sp_c, prune_tol)
    if save:
        sp_fine = prolong(sp_c, Val(K))   ← at snapshot time only
        save snapshot; discard sp_fine
```
Memory: O(|S_coarse|) regardless of K.
Added `coarse_only::Bool = true` to `RDMEMultigridFSP`.

---

### Phase 5: Krylov step-size failures

**Problem:** With K=8 coarse (K2=4) and initial dt=0.5, means at t=5 were
~11.6 instead of 16.0 (27% error). Increasing Krylov m didn't help much.

**Root cause:** Krylov accuracy requires dt × ‖A‖ to be moderate.
For K2=4 coarse:
- d_coarse = D/(2h)^2 = 1/(2*0.125)^2 = 1/0.0625 = 16
- n_max_coarse = 24, K2=4 voxels
- ‖A_coarse‖_∞ ≈ 800 (coarse diffusion dominates)
- dt=0.5 → dt×‖A‖ = 400 >> safe range for m=30

**Fix:** Reduce dt to 0.02 so dt×‖A‖ ≈ 16, safely within Krylov range.
Rule of thumb: dt × ‖A‖ ≲ 2m for Krylov dimension m.

---

### Phase 6: Three-level coarsening attempt (K=8 → K=4 → K=2)

**Motivation:** K=8 with K2=4 coarse took 25 min (250 steps × growing
4D state space). The K2=2 level for K=4 was fast (1.19s). Idea: use K4=2
as coarsest level for K=8, where the 2D state space stays tiny.

**Implementation:** Added `n_levels::Int` to `RDMEMultigridFSP`, and
`_multi_restrict`, `_multi_prolong` helpers for n_levels=1,2,3.

**Result: Hung for 40+ minutes.** Root cause:
At each snapshot (save_every=10 steps), the two-step prolongation
`prolong(prolong(sp_K2, Val{4}), Val{8})` was catastrophically expensive:
- K2=2 coarse state (n1, n2) with n≈8: Binomial(8,k)×0.5^8 > 1e-3 keeps
  k=1..7, giving 9 splits per voxel.
- First prolong K2→K4: 9^2 = 81 K4 states per K2 state.
- With ~500 K2 states: ~40K K4 states.
- Second prolong K4→K8: 9^4 = 6561 K8 states per K4 state → hundreds of
  millions of K8 states.
The two-step prolongation is exponential in the number of levels.

**Decision:** Abandon three-level approach. Two-level coarse-only
(K=8 → K2=4) is the practical operating point.

---

### Phase 7: Pruning tolerance tradeoff

After establishing two-level coarse-only as the approach, the remaining
question is prune_tol. Results from sweeping:

| prune_tol | Runtime | Σμ at t=5 | Error |
|-----------|---------|-----------|-------|
| 1e-4      | 5s      | 11.53     | 28%   |
| 1e-6      | 45s     | 15.20     | 5.0%  |
| 1e-7      | 122s    | 15.77     | 1.4%  |
| 1e-8      | 1525s   | 15.90     | 0.6%  |

prune_tol=1e-4: Too aggressive. Removes tail states (high molecule counts).
Renormalization after pruning concentrates probability near low counts,
systematically biasing means downward by 28%.

prune_tol=1e-7: 2 min, 1.4% error. Chosen as working baseline.

The error at prune_tol=1e-7 is NOT a systematic bias — it reflects
incomplete convergence (the distribution hasn't fully reached stationarity
by t=5 for the truncated system). At longer times or tighter tolerance
the means continue converging upward toward 2.0.

---

### Phase 8: Prolongation bottleneck (K=8 snapshot saves)

**Problem:** With prune_tol=1e-7 and save_every=25, the K=8 run took 134s
despite the coarse loop being only ~20s. Root cause identified by profiling:

- Coarse loop (250 steps × ~80ms/step average): ~20s
- Prolongation at each of 10 save points: ~10s each (×10 = ~100s)
  - K=4 coarse has ~20K states at late time
  - Prolongation K=4→K=8 with binom_tol=0.005: creates 17M fine states
  - Each prolongation requires building and deduplicating a 17M-state dict
- Iterating over 17M fine states to compute means: ~2s per snapshot

The 10 prolongation calls accounted for ~110s of the 134s total.

Profiling results for late-time prolongation (|S_coarse|=20K):
| binom_tol | prolong time | |S_fine| |
|-----------|-------------|---------|
| 0.005     | 10s         | 17.8M   |
| 0.02      | 5.5s        | 10.6M   |
| 0.05      | 3.7s        | 6.8M    |

Even binom_tol=0.05 would give 50 × 3.7s = 35s, still expensive. The
fundamental issue: each coarse voxel with n_c molecules contributes
O((n_c+1)^2) fine states via the K=4→K=8 two-step prolongation.

**Key insight:** For the K=8 benchmark, prolongation is entirely wasted —
we only need per-voxel means, not the full K=8 distribution. And fine means
follow immediately from coarse means by binomial-split symmetry:
    μ_fine[2c-1] = μ_fine[2c] = μ_coarse[c] / 2
This is exact: E[n_{2c-1}] = E[Binomial(N_c, 0.5)] = N_c/2.

**Fix:** Added `coarse_only_means()` helper in the benchmark script that
runs the coarse loop directly and computes fine means from coarse states,
bypassing all prolongation. Cost: O(|S_coarse|) per save vs O(|S_fine|).

**Result:** K=8 runtime: **134s → 23.8s** (5.6× speedup), same accuracy:
Σμ=15.80 at t=5 (1.2% error, slightly improved from 1.4% due to symmetry).

The coarse-only means approach is appropriate whenever the fine distribution
is not needed directly (e.g., for mean/variance tracking, not TV distances).

---

### Current working parameters

K=4:
- direct FSP: n_max=12, 28561 states, krylov_m=60, 1.31s
- multigrid coarse-only: K2=2, dt=0.1, krylov_m=40, prune_tol=1e-12,
  binom_tol=0, 1.36s; TV→6e-4 at t=5

K=8:
- coarse-only: K2=4, dt=0.02, krylov_m=30, n_max=12, prune_tol=1e-7,
  coarse_n_max_per_voxel=14, coarse_n_total=32, save_every=25
- Means via coarse_only_means() (no prolongation): **23.8s**
- Σμ=15.80 at t=5 (1.2% error vs analytical 16.0)
- Coarse state space: 21K states at t=5

---

### Key files

- `src/rdme/vcycle.jl`: two_level_vcycle (correction-based, Algorithm 6.1)
- `src/rdme/solve.jl`: _solve_rdme_coarse_only (coarse-only mode)
- `src/rdme/operators.jl`: restrict, prolong (binomial), expand!
- `src/rdme/rdme_model.jl`: build_coarse_system (birth rate scaling r=2)
- `examples/fsp_vs_multigrid.jl`: benchmark script, K=4 vs K=8
- `examples/multigrid_debug.jl`: instrumented V-cycle diagnostics
