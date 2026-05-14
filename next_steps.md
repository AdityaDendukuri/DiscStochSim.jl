# Next Steps — Adaptive FSP for 2D RDME (Start Over)

## What we want to build

An adaptive FSP that solves the 2D birth-death RDME on a 20×20 grid efficiently. The key insight: the RDME is a set of per-voxel CMEs (one per voxel) coupled by off-diagonal diffusion terms. This is spiritually identical to the user's existing adaptive CME paper (DendukuriPetzoldAdaptiveFSP_2025), just extended to voxels.

---

## The right mental model

**RDME generator = block-diagonal CME generators + diffusion coupling**

```
A = diag(A_1, A_2, ..., A_K) + D_diffusion
```

where `A_j` is the per-voxel birth-death CME and `D_diffusion` is the off-diagonal hopping.

For K_eff tracked voxel groups, the generator is:
- Birth: λ·M_j per group j (M_j = # fine voxels in group)
- Death: μ·n_j per group j
- Diffusion hop j→k: D·n_boundary(j,k)/M_j · n_j  [exact for M_j=1; valid for M_j>1 when group is in diffusion equilibrium]

For individual voxels (M=1): rate = D·n_j (exact, same as fine RDME).

**Adaptive FSP (from the paper, extended to voxels):**
1. Start sparse: IC = point mass, |S|=1
2. Expand: add states reachable via reactions + diffusion hops
3. Prune: remove low-probability states (quantile-based or threshold)
4. Merge voxels that are in diffusion equilibrium (exact restriction)
5. Extend to new voxels when diffusion flux exits tracked region (n_new=0, exact)

---

## Key design decisions learned this session

### 1. Mean-field for spatial decisions, FSP for distribution

The mean-field ODE gives E[n_ij(t)] exactly for linear systems. Use it to:
- Determine which voxels are at the wavefront (coarsening decisions)
- Compute generator boundary fractions (D·mf[boundary_voxel]/mf_total_region)

The FSP gives P(n_1,...,n_K, t) — variance, bimodality, full distribution. The ODE cannot do this.

**Never** use FSP-derived uniform spread (n_region/M) for coarsening decisions when regions are not in equilibrium — this creates a self-defeating feedback loop.

### 2. Individual front voxels + global merged backgrounds

Design B that works:
- **Settled background** (:hi, n̄ ≥ θ_hi): one global merged region — equilibrated, Multinomial valid
- **Empty background** (:lo, n̄ ≤ θ_lo): one global merged region — n≈0, trivially valid
- **Active front**: each voxel its own FSP dimension — exact, no Binomial needed
- **Excess front** (if > max_front): one bulk region — approximate, but small and near-equilibrium

**Why this is correct:**
- Merged regions: the generator rate D·n_adj/M is exact when the region IS in diffusion equilibrium. Only merge when that condition holds (verified via mf-state gradient).
- Individual front voxels: rate = D·n_j, always exact.
- Extension from empty: n_new=0, exact since empty background truly has n≈0.

**Generator boundary fractions** (critical fix in `build_rdme_adaptive_2d`):
Use D · (Σ mf[boundary_k] / mf_total_region) instead of D · n_adj/M when mf_state shows non-uniform distribution within a merged region. This gives the correct expected diffusion rate without approximation error for linear systems.

### 3. State-space transition (remap_sp_transition)

When K_eff changes, map each old FSP state to the new K_eff-dimensional state:
- **Fine voxel stays fine**: copy count directly (1:1 mapping)
- **Fine voxel settles → background**: fold count into background (exact marginalization, probs accumulate for same new state)
- **New fine voxel from background**: initialize n=0 (exact since background is near-empty)
- **Background → background**: use **majority-destination** rule (most voxels of old region go to which new region)

The majority-destination rule (not first-representative voxel) is critical when regions partially change labels between rebuilds.

### 4. No V-cycle needed

The V-cycle (restrict → coarser expv → prolong) was our original approach. It was trying to solve:
- Large K_eff (many merged super-voxels) → expensive Level-1 expv
- Use coarser Level-2 to accelerate

But this caused the Binomial prolongation explosion (Part B: large-n Binomial splits for settled background pairs).

**With individual front voxels, the V-cycle is unnecessary:**
- K_eff ≤ max_front + 3 ≈ 8-10
- State space ≤ ~5k-50k states
- Direct expv is fast at this size
- The Binomial explosion is avoided entirely

### 5. State space size

With individual_front=True and global merged backgrounds:
- K_eff = 2 backgrounds + F individual front + 0 or 1 bulk front
- F ≤ max_front (typically 4-8)
- State space: n_settled × n_empty × 6^F (approximately)
- With n_settled ≈ 18 values, n_empty ≈ 22 values, F=6: 18×22×6^6 ≈ 20M states — too large!
- The two large background dimensions × fine voxels is the explosion

**Key constraint:** The product of background dimensions × fine voxel dimensions is the bottleneck. For 2 backgrounds with 18×22 values each and 6 fine voxels with 6^6: the state space is ~20M, too large.

**Practical working range:** max_front=4-6 with one or two backgrounds gives ~5k-50k states at peak.

---

## What to implement (clean start)

### File structure
- Keep `src/rdme/rdme_model.jl` with `build_rdme_adaptive_2d` (CLEAN version: exact `n_adj/M` rate, no mf_state parameter)
- Keep `src/rdme/voxel_grid.jl` with `build_adaptive_cmap_2d(individual_front=true)` (global :hi/:lo merging)
- Keep `src/rdme/operators.jl` with mode-expansion `_prolong_vw_split!` (for any future use)
- Rewrite `examples/fig_2d_birthdeath_adaptive.jl` cleanly

### Algorithm (clean)

```
t=0:
  - Build cmap from mf_state (initial IC) with individual_front=true
  - Build sys from cmap (exact generator)
  - Set IC as point mass
  - expand!(sp, sys, bc; depth=2)

each step:
  - expand!(sp, sys, bc; depth=1)    # grow state space with reactions+diffusion
  - expv(dt, A, sp.probs)            # direct FSP step
  - prune_threshold!(sp, tol)        # remove low-probability states
  - advance mf_state ODE             # for coarsening decisions

every rebuild_every steps:
  - new_cmap = build_adaptive_cmap_2d(mf_state, ..., individual_front=true, max_front=N)
  - if topology changed:
      new_sys = build_rdme_adaptive_2d(model, grid, new_cmap, Val(new_Ke))
      new_sp  = remap_sp_transition(sp, old_cmap, new_cmap, Val(new_Ke))
      expand!(new_sp, new_sys, new_bc; depth=1)
  - else:
      rebuild sys (generator rates may change as mf_state evolves)
```

### Parameters to calibrate
- `max_front`: 4-6 (larger = more accurate but larger state space)
- `rebuild_every`: 10-20 steps (more frequent = better tracking but more rebuilds)
- `θ_lo_2d`: 0.5*n_ss = 1.0 (must be above global birth equilibrium to keep empty region stable)
- `θ_hi_2d`: 0.85*n_ss = 1.7 (settled threshold)
- `prune_tol`: 1e-9

### Correctness check
For linear birth-death: FSP E[n_region_j] = Σ mf_state[voxel in region j]. The error measures discrepancy between FSP dynamics and the true mean-field. Should be <5 molecules max error for well-chosen coarsening.

---

## What NOT to do (lessons learned)

1. **No V-cycle** — the Binomial Part B explodes for large-n merged regions. Unnecessary with individual front voxels.

2. **No Multinomial prolongation** at K_eff transitions — use `remap_sp_transition` (exact for individual voxels and settled backgrounds, n=0 for new front voxels).

3. **No FSP-derived voxel means for coarsening decisions** — creates a self-defeating loop where wrong FSP dynamics → wrong coarsening → worse FSP dynamics. Use mf_state instead.

4. **No gradient-based coarsening** — the gradient criterion identifies currently-active voxels but doesn't guarantee regions are in diffusion equilibrium. Threshold-based (:hi/:lo/:front) correctly identifies equilibrated regions by construction.

5. **Don't merge disconnected settled regions** — the two corners should stay as ONE global settled dimension (not two separate connected components) to keep the state space tractable.

6. **Don't make the empty background too small** — θ_lo must be above the global birth equilibrium (≈0.18 at t=1) or the empty region disappears and all voxels become :front.

---

## Current state of the branch

The `multigrid-fsp` branch has many iterations of changes to:
- `src/rdme/operators.jl`: mode-expansion Binomial (√n instead of n per split) — keep
- `src/rdme/voxel_grid.jl`: individual_front coarsening, majority-destination in `build_coarse_pairing` — partly keep
- `src/rdme/vcycle.jl`: K2>=K1 fallback in `adaptive_vcycle_step!` — can ignore (V-cycle not used)
- `examples/fig_2d_birthdeath_adaptive.jl`: many iterations — rewrite cleanly

The clean rewrite should start from scratch on the example, using only the library functions that are needed (no V-cycle calls).
