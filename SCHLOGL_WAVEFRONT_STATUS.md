# Schlogl Wavefront Fix — Current Status

## Goal
Make the `SpatialFSP` correctly simulate bistable Schlogl wavefront propagation on a 2D grid.
The wavefront should expand from an initial equil-hi droplet into an equil-lo background.

## What Was Working Before This Session
- BranchingDeath model wavefront: works with multigrid (validated, wave speeds L=1,2,3 measured)
- Schlogl 2D droplet figure: existed in `examples/schlogl_2d_droplet.jl` but produced **frozen** droplets

## Root Causes Found (Three Bugs)

### Bug 1: Bidirectional cross-basin flux (FIXED)
**File:** `src/spatial_adaptive/bd_spatial_fsp.jl`, step 4 (~line 2731)  
**Problem:** The old code had an `elseif` case that let equil-lo voxels generate flux BACK to equil-hi, causing the hi region to collapse.  
**Fix:** Removed the elseif; flux is now unidirectional (equil-hi → equil-lo only).  
**Status:** DONE

### Bug 2: Multigrid classifies all Schlogl active voxels as "interior" (FIXED)
**File:** `src/spatial_adaptive/bd_spatial_fsp.jl`, `_step_multigrid!` (~line 2083)  
**Problem:** For Schlogl (no empty voxels), ALL active voxels were classified as "interior" → multigrid total-count restriction loses bistable structure.  
**Fix:** Added Schlogl-specific edge classification: active voxels adjacent to equil of the opposite basin → "edge" (per-voxel expv).  
**Status:** DONE

### Bug 3: `_expv_voxel` support too narrow + krylov_m too small (PARTIALLY FIXED, STUCK)
**File:** `src/spatial_adaptive/bd_spatial_fsp.jl`, per-voxel fallback loop in `step!` (~line 2612)  
**Problem:** Two sub-problems:
1. `n_include = n_un = 63` → support [30, 71] — too narrow for strong inflow (wavefront voxel with 1 equil-hi neighbor has in_rate=556; probability should reach n_hi=164 but escapes the truncated window)
2. `krylov_m = 30` → insufficient for λ·dt ≈ 278 — `expv` returns sum=0.034 instead of sum≈1.0; normalised remnant has mean≈57 < n_un → snap fires but snaps to equil-LO not equil-hi

**Fixes applied:**
- Changed `n_include = n_hi` (extended support to n_hi+8 ≈ 172)
- Added time-splitting: `n_sub = ceil(max_rate * dt / krylov_m)` sub-steps each with λ·dt_sub ≤ krylov_m

**Verification:** Direct test with m=30 vs m=100 on Q matrix shows:
- m=30: sum=0.034 (broken)
- m=100: sum=1.0 (correct, mean≈135 after dt=0.5 from delta(n_lo))

**Status:** Code written, but compilation/testing TIMED OUT in the last run. Needs verification.

### Bug 4 (consequence of Bug 3): Deactivation criterion (FIXED, depends on Bug 3)
**File:** `src/spatial_adaptive/bd_spatial_fsp.jl`, step 3a (~line 2695)  
**Problem:** `p[n_un+1] < ε_equil` never fires for wavefront voxels stuck at mean-field fixed point n≈68.  
**Fix:** Added mean-based snap: `abs(mean_p - n_un) > 4.0` → snap to nearest basin (delta(n_hi) or delta(n_lo)).  
**Status:** DONE — but only effective once Bug 3 is fixed (otherwise expv gives wrong mean)

## Current State of `bd_spatial_fsp.jl`

All four fixes are in the source file. The time-splitting change (Bug 3) was the LAST edit made before the session got stuck.

Key lines to verify:
- **~2612-2640**: per-voxel loop with time-splitting + `n_include=n_hi`
- **~2695-2729**: step 3a with mean-based snap criterion (`snap_margin=4.0`)
- **~2083-2099**: multigrid edge classification (Schlogl front voxels → edge)
- **~2731**: cross-basin flux (unidirectional, fixed)

## What Needs Testing

1. Force recompile: `rm ~/.julia/compiled/v1.12/DiscStochSim/*.{ji,so}` then `using DiscStochSim`
2. Run 15×15 test: start with R=2 equil-hi circle, rest equil-lo, step 20 times with dt=0.5
3. Expected: equil-hi count grows from 13 → 25 → 41 → ... (wave expands ring by ring)
4. If still frozen: check `n_sub` calculation — `max_rate = in_rate + out_coeff * n_hi = 556 + 8*164 = 1868`, `n_sub = ceil(1868*0.5/30) = ceil(31.1) = 32` sub-steps. Each sub-step dt_sub = 0.5/32 ≈ 0.016, λ·dt_sub ≈ 1868*0.016 ≈ 30 ✓

## Test Script
```julia
using DiscStochSim, Printf
D = 2.0; K = 15
model = SchloglModel1D(D, 0.028, 3.2e-4, 21.0, 1.0)
n_lo, n_un, n_hi = schlogl_fixed_points(model)
s = SpatialFSP(model, K, K; n_max=220, n_max_eq=220,
               ε_expand=0.15, ε_equil=0.12, use_block_vcycle=false)
ci, cj = K÷2+1, K÷2+1
for i in 1:K, j in 1:K
    r = sqrt((i-ci)^2 + (j-cj)^2)
    p = zeros(221); p[r <= 2.0 ? n_hi+1 : n_lo+1] = 1.0
    s.equil_dists[CartesianIndex(i,j)] = p
end
n_hi_count(s) = count(voxel_mean(p) > n_un for p in values(s.equil_dists))
@printf("t=0: equil-hi=%d\n", n_hi_count(s))
for step in 1:20
    step!(s, 0.5)
    @printf("step %2d: active=%3d equil-hi=%3d\n", step, length(s.dists), n_hi_count(s))
end
```

## Key Files
- `src/spatial_adaptive/bd_spatial_fsp.jl` — all four fixes applied
- `src/spatial_adaptive/schlogl_spatial_fsp.jl` — Schlogl generator (unchanged, correct)
- `examples/schlogl_2d_droplet.jl` — the figure script (needs working wavefront first)
- `paper/talk_multigrid_v5.tex` — slides (rewritten slide 3, 6, 7 to match current method)
