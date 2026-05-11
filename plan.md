# DiscStochSim.jl — Scientific Plan

## The Problem

The Reaction-Diffusion Master Equation (RDME) tracks the joint probability
distribution p(n₁, ..., nₖ, t) over molecule counts across K spatial voxels.
The state space is n_max^(S·K) — exponential in K. Even K=8, S=1 is intractable
with naive FSP.

Two sources of computational hardness:
1. **State space size** — exponential in K·S
2. **Stiffness** — fast diffusion forces tiny timesteps in explicit methods

---

## The Scientific Niche

The multigrid FSP approach targets a specific, well-defined regime:

| | Fast diffusion | Slow / arbitrary diffusion |
|---|---|---|
| **Low copy numbers (n ≈ 1–10)** | **Multigrid FSP** ← this work | MPS / tensor networks |
| **High copy numbers** | Chemical Langevin / SDE | Mean-field PDE |

**Why this quadrant matters:** SSA is stiff (must simulate every diffusion hop),
mean-field PDEs are wrong (copy numbers too low for Gaussian approximation),
standard FSP is intractable (state space too large). Multigrid FSP is the only
tractable exact method here.

**Why SSA fails specifically:** For D_A = 1.0 on an 8×8 grid with dx = 1/8,
the diffusion rate is D/dx² = 64 hops/molecule/time. With ~2 molecules/voxel
across 64 voxels, SSA spends ~8000 diffusion events per unit time to simulate
~10 reactions. Efficiency: ~0.1%. Multigrid FSP uses dt = 1.0 — a 1024× larger
timestep — because fast diffusion is implicit in the Binomial coarsening.

**The compelling biological case:** Signaling molecules in small cells — 1–10
copies of a transcription factor diffusing rapidly in a nucleus. The probability
distribution over spatial locations determines activation probability of a
distant target gene. This is a quantity SSA requires many trajectories to
estimate and mean-field gets qualitatively wrong at low copy numbers.

---

## The Core Algorithm

**Key insight:** When diffusion between two adjacent voxels dominates reactions,
their joint distribution equilibrates to Bin(n̄, ½) — exactly. Replacing the
pair (n₁, n₂) by their total n̄ reduces the per-pair state space from n_max²
to n_max. Applied to all K voxels: n_max^K → n_max^(K/2). Square root of the
original state space.

**Coarsening is valid when** α_j = |μ_{2j-1} − μ_{2j}| / (μ_{2j-1} + μ_{2j}) ≈ 0
(voxels balanced). Fails at bistable fronts where α_j ≈ 1.

**The adaptive mixed-resolution algorithm (talk, slides 9–10):**
- Maintain mixed state p^M: coarse groups (n̄_j) where α_j < α_lo, fine groups
  (n_{2j-1}, n_{2j}) where α_j > α_hi
- Evolve ṗ^M = A^M p^M directly — no compress→reconstruct cycle
- Reassign (refine/coarsen) only when front crosses a group boundary
- The mixed generator A^M couples coarse and fine blocks at group boundaries

**Why no reconstruction:** Reconstruction after every step forces Binomial splits
on all groups including the front, where the true distribution is bimodal. Evolving
A^M directly avoids this — coarse groups carry the Binomial assumption implicitly,
fine groups track the true within-group distribution exactly.

---

## Scalability

For truly large spatial systems, the multigrid hierarchy reduces to a coarse-only
solve with fine-level mean recovery:

```
Fine grid: 64×64 (CartesianIndex{8192}) — intractable, never stored
                    ↕ two levels of 4× coarsening (analytical)
Coarse solve: 2×2 (CartesianIndex{8}) — ~20k states, adaptive FSP
                    ↕ mean broadcasting at save times
Output: E[nA(x,y)], E[nB(x,y)] at full fine resolution
```

**The O(1) scaling result:** The coarse solve cost is independent of fine grid
size. A 16×16 grid and a 64×64 grid have the same coarsest-level computation.
The multigrid hierarchy justifies the approximation and provides an error
certificate: run at 2×2 and 4×4, compare means — agreement validates coarsening.

**Limitation:** Coarse-only gives first moments (spatial mean field), not the
full joint distribution. For spatial correlations, noise spectra, or rare event
probabilities at the fine level, the full distribution is needed and is
fundamentally intractable for large K. This is an honest boundary of the method.

---

## What the Method Does NOT Address

- Arbitrary diffusion rates (coarsening requires d >> reaction rates)
- High copy numbers (mean-field PDE is adequate there)
- Full fine-level joint distribution for large systems (exponential wall is real)
- Multi-species with many species per voxel (state space is n_max^(S·K))

For arbitrary diffusion + low copy numbers: MPS/tensor network methods are the
right tool (bond dimension χ captures spatial correlations; multigrid is the χ=1
MPS limit under fast diffusion).

---

## Current State of the Code

### Working
- `examples/birthdeath_stiff.jl` — 1D birth-death V-cycle, 1024× CFL speedup
- `examples/rdme_birthdeath_1d.jl` — multi-level benchmark, Poisson validation
- `examples/bottleneck_rdme.jl` — 1D two-species A→B, 8× speedup, TV=0.009

### In progress
- `examples/fig_2d_bottleneck.jl` — 2D coarse-only FSP with mean visualization
  - Runs correctly (36.8s, 20k states at t=60) but PNG save not yet verified
  - Fix needed: remove `display(fig)` (headless), tighten bc to `2*N_MAX`

### Missing (talk slides 11, 12, 15, 16 are empty)
- **Slide 11:** TV distance over time — birth-death and Schlögl front,
  naive grouping vs adaptive mixed-resolution
- **Slide 12:** Per-voxel means accuracy + speedup vs K scaling curve
- **Slide 15:** Self-consistency check — mixed-resolution ODE vs Monte Carlo
  on same generator (means agree to <1%, TV ~ N^{-1/2})
- **Slide 16:** Interface orientation figure K=8 Schlögl — coarse P(n̄₂)
  loses orientation; fine P(n₃, n₄) retains bimodal structure

---

## Priority Order

1. **Fix and verify `fig_2d_bottleneck.jl`** — remove display(), tighten bc,
   confirm PNG saves correctly. This is the 2D coarse-only result.

2. **Generate slide 11 figures** — TV distance comparison for birth-death and
   Schlögl. These validate the core adaptive criterion claim.

3. **Generate slide 12 figures** — speedup vs K scaling. This is the
   quantitative performance claim.

4. **Generate slides 15–16** — self-consistency and interface orientation.
   These are appendix/supporting material.

5. **Paper write-up** — framing as above: fast diffusion + low copy numbers
   is the target regime; multigrid is the only tractable exact method there.

---

## Open Questions

- Rigorous error bound replacing the heuristic imbalance criterion α_j
- Extension to multi-species (S > 2): state space is n_max^(S·K), coarsening
  scales by patch_size^S — still a square root reduction per level
- Deeper hierarchies (m > 2 levels): each level halves the exponent; for K=64
  you need 5 levels to reach CartesianIndex{8}
- Connection to MPS: multigrid coarsening = bond dimension χ=1 MPS under fast
  diffusion. Could a variable-χ MPS extend the approach to arbitrary diffusion?
