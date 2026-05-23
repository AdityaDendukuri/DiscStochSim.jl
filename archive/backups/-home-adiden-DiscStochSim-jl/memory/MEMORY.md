# DiscStochSim.jl — Session Memory

## Project overview
Julia package implementing adaptive FSP for stochastic chemical kinetics.
Key methods: `AdaptiveFSP` (FP-FSP), `KrylovFSP` (baseline).
See `src/adaptive_fsp/solve.jl`, `src/adaptive_fsp/krylov_fsp.jl`.

## Separate comparison scripts (3-script pattern per system)
Each system has: `examples/<system>_fpfsp.jl`, `examples/<system>_krylov.jl`, `examples/<system>_plot.jl`
Results serialized to `examples/results/<system>_{fpfsp,krylov}.jls`.
Output PDFs: `examples/output/comparison/<system>_comparison.{png,pdf}` → copied to `paper/`.

### Oregonator (verified, paper-matching results)
- FP-FSP: `oregonator_fpfsp.jl` → AdaptiveFSP(ε_dt=1.0, prob_quantile=0.1, flux_tolerance=1.0, expansion_depth=1, save_interval=1000, flux_method=:total, **expand_method=:ssa**)
  - Results: 36s, 26149 iters, max|S|=1057, dt∈[1.6e-6, 1.9e-5] (paper updated to match)
- KrylovFSP: `oregonator_krylov.jl` → INLINED old algorithm (non-absorbing gen, expmv, no τ growth, ℓ=5 lookback)
  - Parameters: ε=1e-2, τ_init=1e-6, ℓ=5, r=1, max_iter=3000
  - Results: 73.6s, 3000 iters, 0 rejections, τ=1e-6 (constant), t=0.003 (2.73%), max|S|=13294

**Key fix**: old solve.jl used `expand_ssa!` (1 SSA neighbor per state); current uses `expand!` (all stoichiometric).
Added `expand_method::Symbol = :stoich` param to `AdaptiveFSP`; Oregonator uses `:ssa` to replicate old behavior.
**Stoich for Oregonator is structurally impossible**: growth factor = (1-α)×6. Need α≥0.83 to prevent explosion, but that collapses |S| to 1. SSA (growth factor (1-α)×2) stabilizes at ~436 mean states with α=0.1. Do not attempt to switch Oregonator to stoich.
Old KrylovFSP (non-absorbing, expmv): τ-lock via NaN detection, NOT FSP mass conservation. τ never grows.
New KrylovFSP (absorbing, Algorithm 2): used for Robertson.

### Bottleneck results (current, verified by running unified_comparison_plot.jl)
- FP-FSP: 1.0s, ⟨C⟩=4.978, |S|=1001 ✓
- KrylovFSP (absorbing gen): ⟨C⟩≈0.001, ⟨B⟩=0.010 (correct!), 7 iters, τ: 100→3700
  - droptol=1e-10 retains B=1 correctly; failure is WAVEFRONT TRUNCATION
  - τ doubles each step; after 7 steps |S|=14 (only C≤12), but true dist. spans C=[0,1000]
  - With absorbing generator: probability to C≥13 is absorbed/lost → ⟨C⟩≈0 in projection
  - Previous "verified" ⟨C⟩=0.12 was from old non-absorbing KrylovFSP; current absorbing gives ≈0

### Robertson (verified results from examples/results/*.jls)
- **Correct params**: `flux_method=:maximum, expand_method=:stoich, ε_dt=0.01, prob_quantile=0.4, flux_tolerance=1e-9`
- FP-FSP (tf=1e4): 1144s, 662,472 iters, dt∈[1.1e-9, 0.026], burst-and-relax ✓
  - Final: A=506.6, B≈0, C=9493.4, |S|_final=442
- KrylovFSP (absorbing gen, max_iter=500): 500 total iters, 336 rejections (67%), 1.92s, t=227.8 (2.3%)
  - Accepted τ range: [2e-11, 51.2] — NOT "locked" at narrow range; τ fluctuates widely
  - 67% rejection rate from FSP failures when B≥2 states hit projection boundary
  - Paper narrative: "high rejection rate prevents efficient progress, reaches 2.3% in 500 iters"
  - τ_log has 827 entries (every attempt including rejected); time_log also 827 entries
  - To show accepted-only Krylov τ: filter by `diff(time_log) .> 0`
- Data: examples/results/robertson_{fpfsp,krylov}.jls
- Scripts: robertson.jl, robertson_krylov.jl

## KrylovFSP implementations (TWO distinct versions used in paper)
1. **Old (non-absorbing, `examples/oregonator_krylov.jl` inline)**: `build_generator` (conservative), `expmv` from Expokit, no τ growth, lookback window ℓ=5, FSP: sum(p)≥1-ε. τ-lock = NaN in expmv for large propensities.
2. **New (absorbing, `src/adaptive_fsp/krylov_fsp.jl`)**: Algorithm 2 Vo & Sidje 2016, `build_generator(...;absorbing=true)`, `expv` from ExponentialUtilities, τ doubles on accept, time-scaled FSP. Used for Robertson.

## AdaptiveFSP parameters (src/adaptive_fsp/solve.jl)
- `flux_method`: `:total` (default) uses Φ=Σp·w; `:maximum` uses Φ=max(p·w). Robertson REQUIRES `:maximum`: in slow post-burst phase, Σp·w≈300 (inflated by tiny-p/huge-w B states) but max(p·w)≈0.03, giving 10000× larger dt with `:maximum`. Using `:total` with ε_dt=0.01 gives ~400M steps; `:maximum` gives ~662K.
- `expand_method`: `:stoich` for Robertson/bottleneck; `:ssa` for Oregonator (avoids state space inflation)
- `expansion_depth`: 1 for all systems (increasing d for SSA explodes state space and gives wrong results)
- `dt_max`: set to `1/λ_downstream` for slow-start systems (Bottleneck: `dt_max=1/k2=10.0`); `Inf` otherwise
- Time step formula: `dt = min(ε_dt/Φ, expansion_depth * dt_max, tf-t)` where Φ depends on flux_method

## Paper tex (paper/ex_article.tex)
### Comparison section (simplified this session)
- Sec 5.6 only (Sec 5.7 "Corrected Efficiency" removed); single unified 3-panel figure
- Figure: unified_comparison.pdf → loaded from pre-computed .jls files + Bottleneck inline
- Script: examples/unified_comparison_plot.jl (fast, ~30s, NO inline Oregonator/Robertson)
- One summary table (tab:comparison_summary), one figure (fig:unified_comparison)
- Oregonator (panel a): |S| vs time — FP-FSP compact, Krylov explodes
- Bottleneck (panel b): ⟨C⟩ vs time — FP-FSP accurate, Krylov ≈0 (wavefront truncated)
- Robertson (panel c): accepted dt/τ vs time — FP-FSP full horizon, Krylov stops at t=228

## Source files modified (this and previous sessions)
- `src/adaptive_fsp/solve.jl` — added flux_method, expand_method params; max_propensity dt cap
- `src/adaptive_fsp/krylov_fsp.jl` — complete rewrite (faithful Algorithm 2)
- `src/DiscStochSim.jl` — added KrylovFSP exports + krylov_fsp.jl include
- `src/state_space.jl` — added `prune_threshold!` function

## Plot style
B&W SIAM. Analytical = thin dotted line + star5 markers (markersize=9).
Axis labels: Unicode ⟨C⟩ (U+27E8/U+27E9) — NOT LaTeX \langle (GR backend unsupported).

## Seeding
Use `Random.seed!(1)` (not `Random.default_rng(1)` which is invalid syntax).
