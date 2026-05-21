# Adaptive Spatial FSP — Detailed Methods Notes

*DiscStochSim.jl, 2026-05-20*

---

## 1. Problem Setting

We simulate the **2D spatial RDME** on a $K_x \times K_y$ lattice of voxels.  Each voxel $k$ holds an integer molecule count $n_k \geq 0$ of a single species.  The dynamics combine:

- **Local reactions** with propensities that depend only on $n_k$
- **Diffusion** between nearest-neighbour voxels: molecule hops at rate $D \cdot n_k$ (from voxel $k$ to each of its 4 neighbours)

The exact probability distribution is $P(\mathbf{n}, t)$ over all joint configurations $\mathbf{n} = (n_1, \ldots, n_{K_x K_y})$, which obeys a high-dimensional master equation.

**Core approximation:** outside the active wave front, each voxel's distribution is factorised (product approximation) and either frozen at a background state or updated algebraically.  Only voxels at the wave front are tracked with the full per-voxel CME.

---

## 2. Chemical Models

### 2.1 BirthDeathRDME

Linear reactions:
$$\emptyset \xrightarrow{k_b} X, \qquad X \xrightarrow{k_d} \emptyset$$

Unique local steady state: $n^* = k_b / k_d$.  Background state: Poisson$(n^*)$.

### 2.2 BranchingDeathRDME

$$X \xrightarrow{k_b} 2X, \qquad X \xrightarrow{k_d} \emptyset$$

Branching birth introduces nonlinear growth; background treated as saturated (mean near $n_\text{max}$).

### 2.3 SchloglModel1D (Bistable)

Four reactions in a single voxel:
$$
\emptyset \xrightarrow{c_3} X, \quad
X \xrightarrow{c_4} \emptyset, \quad
2X \xrightarrow{c_1 \binom{n}{2}} 3X, \quad
3X \xrightarrow{c_2 \binom{n}{3}} 2X
$$

Propensities: $c_3$, $c_4 n$, $c_1 n(n-1)/2$, $c_2 n(n-1)(n-2)/6$.

With standard parameters (c1=0.028, c2=3.2e-4, c3=19.5, c4=1.0) the isolated voxel has **three fixed points**:
- $n_\text{lo} \approx 31$ — stable (low state)  
- $n_\text{un} \approx 72$ — unstable (barrier)  
- $n_\text{hi} \approx 162$ — stable (high state)

The thermodynamic bias toward each state is captured by the **quasi-potential sum**:
$$\Delta\Phi = \sum_{k=n_\text{lo}}^{n_\text{hi}-1} \ln\frac{b(k)}{d(k+1)}$$
where $b(n)$, $d(n)$ are birth and death propensities.  $\Delta\Phi > 0$ means high state preferred; $\Delta\Phi < 0$ means low state preferred; $\Delta\Phi \approx 0$ is Maxwell coexistence.

---

## 3. SpatialFSP Data Structure

```julia
mutable struct SpatialFSP
    model        :: AbstractSpatialRDMEModel
    K_x, K_y     :: Int                  # grid dimensions

    dists        :: Dict{CartesianIndex{2}, Vector{Float64}}   # ACTIVE voxels: full p_k
    equil_dists  :: Dict{CartesianIndex{2}, Vector{Float64}}   # BACKGROUND voxels: δ or Poisson

    block_states :: Dict{CartesianIndex{2}, StateSpace}        # V-cycle: cached joint 2×2 states
    pending_flux :: Dict{CartesianIndex{2}, Float64}           # accumulated cross-basin flux

    t            :: Float64
    n_max        :: Int        # truncation for active voxel distributions
    n_max_eq     :: Int        # truncation for background (= n_max for Schlögl; auto for BD)
    ε_expand     :: Float64    # flux threshold for background → active transition
    ε_equil      :: Float64    # convergence threshold for active → background transition
    ε_prune      :: Float64    # zero-out entries below this in p_k
    krylov_m_eq  :: Int        # Krylov dimension for background equil updates
    ...
end
```

**Default tolerances:**
| Parameter | Default | Role |
|---|---|---|
| `n_max` | 40 | Truncation of active distribution |
| `ε_expand` | 0.3 | Min flux to promote background → active |
| `ε_equil` | 0.06 | Convergence criterion for demotion |
| `ε_prune` | 1e-14 | Zero-threshold for probability entries |
| `krylov_m_eq` | 8 | Krylov dim for background expv |
| `block_max_states` | 50,000 | Max states in V-cycle joint space |

---

## 4. Voxel States

**Active** (`dists`): full CME distribution tracked.  The voxel is at or near the transition zone (wave front).

**Background** (`equil_dists`): compressed background distribution.
- BD models: Poisson($\mu_\text{loc}$) updated at each step.
- Schlögl: $\delta_{n_\text{lo}}$ or $\delta_{n_\text{hi}}$ — **frozen** (not updated each step).

There is no "empty" class in the Schlögl droplet setup — every voxel is initialised in one of the two basins.

---

## 5. Standard Step Function

`step!(state, Δt)` — one time step without the V-cycle (per-voxel product approximation):

### Step 1 — Compute means
$$\mu_k = \langle n_k \rangle = \sum_n n \, p_k(n) \quad \forall k \in \texttt{dists} \cup \texttt{equil\_dists}$$
Stored in a `Dict` called `means`.

### Step 2 — Evolve active voxels
For each active voxel $k$, compute mean-field rates from stored neighbour means:
$$\lambda_k^\text{in} = D \sum_{j \sim k} \mu_j, \qquad \gamma_k = D |\mathcal{N}(k)|$$
Build per-voxel generator $Q_k$ (see §7) and solve:
$$p_k(t + \Delta t) = e^{Q_k \Delta t}\, p_k(t)$$
via `expv(Δt, Q_k, p_k; m=krylov_m)` from ExponentialUtilities.jl.

Support of $p_k$ is adaptively clipped to where $p_k > \varepsilon_\text{prune}$, extended by `margin=8` on each side.  For Schlögl, $n_\text{un}$ is forced into the support (via `n_include=n_un`) so that cross-basin transitions are not missed numerically.

### Step 3 — Update background
- **BD models:** boundary background voxels (adjacent to an active voxel) update via `expv`; interior background voxels update algebraically to their local Poisson SS.
- **Schlögl:** background voxels are **frozen** (no update).

### Step 3a — Early deactivation [Schlögl only, fires BEFORE cross-basin flux]
For each active voxel, check the outflux through the unstable state:
$$p_k(n_\text{un}) < \varepsilon_\text{equil} \implies \text{move } k \text{ to background}$$
This fires **before** Steps 4–5 so that a newly-deactivated equil-hi voxel immediately generates cross-basin flux in the same step, enabling single-step wave-front propagation.

*Why p(n_un) as criterion:* This is the probability mass at the bistable barrier — the "outflux" through the saddle point.  When it drops below threshold the distribution is fully committed to one basin.  It is the exact bistable analog of the TV-distance criterion used for BD.

### Step 4 — Flux accumulation
**BD models:** for each active/background voxel $k$ with $\mu_k > 0$, accumulate flux toward empty neighbours:
$$\phi_j \mathrel{+}= D \cdot \mu_k \cdot \Delta t \quad \forall j \in \mathcal{N}(k) \text{ that are empty}$$

**Schlögl models:** for each **background** voxel $k$ with $\mu_k > n_\text{un}$ (equil-hi), accumulate cross-basin flux toward adjacent background voxels $j$ with $\mu_j < n_\text{un}$ (equil-lo):
$$\phi_j \mathrel{+}= D \cdot (\mu_k - n_\text{un}) \cdot \Delta t$$
Only **settled background** voxels generate flux.  Active (transitioning) voxels do not.

This is the **unified flux mechanism** — the bistable analog of the BD rule:

| Model | Source | Target | Flux quantity |
|---|---|---|---|
| BD | any active/equil | empty | $D \cdot \mu_k \cdot \Delta t$ |
| Schlögl | equil-hi | equil-lo | $D \cdot (\mu_k - n_\text{un}) \cdot \Delta t$ |

Both use the same activation threshold `ε_expand`.

### Step 5 — Activate target voxels
Any target voxel with $\phi_j \geq \varepsilon_\text{expand}$ is promoted from background to active, starting from its current stored distribution:
- BD empty: starts from $\delta_{n=0}$
- Schlögl equil-lo/hi: starts from stored $\delta_{n_\text{lo}}$ or $\delta_{n_\text{hi}}$

### Step 6 — Deactivate active voxels [BD only]
For BD, check each active voxel for convergence to its local Poisson SS:
$$\text{TV}(p_k, \text{Pois}(\mu_\text{loc})) < 0.08 \;\text{ AND }\; \frac{|\mu_k - \mu_\text{loc}|}{\mu_\text{loc}} < \varepsilon_\text{equil}$$
Schlögl deactivation was handled in Step 3a.

---

## 6. V-cycle Step (for BD Wave Front)

`step!(state, Δt, fg)` — ADI pair-based V-cycle.  Used for the BD 2D case where adjacent active voxels are strongly correlated and the product approximation significantly underestimates the wave-front speed.

The V-cycle processes pairs of adjacent active voxels jointly, alternating row-adjacent and column-adjacent pairs each step (ADI pattern).

### V-cycle stages for a pair $(v_1, v_2)$:

**1. Restrict — build joint product state:**
$$p_{v_1 v_2}(n_1, n_2) = p_{v_1}(n_1) \cdot p_{v_2}(n_2)$$
stored as a sparse `StateSpace{CartesianIndex{2}}`.  Expand by one reaction ("rstep") to capture states reachable in one hop.

**2. Pre-smooth ($\tau_\text{pre} = 0.35 \Delta t$):**
Evolve the joint state under the **intra-pair diffusion only** (molecules hopping between $v_1$ and $v_2$):
$$\dot{p}_{12} = Q_\text{intra}\, p_{12}$$
where $Q_\text{intra}$ has transitions $(n_1, n_2) \to (n_1-1, n_2+1)$ at rate $D n_1$ and reverse.

**3. Restrict — total count:**
Project to the **total molecule count** distribution:
$$p_N = \sum_{n_1+n_2=N} p_{12}(n_1, n_2)$$

**4. Coarse solve ($\Delta t$):**
Evolve $p_N$ under a single-voxel CME for a "2-site" voxel:
- Effective in-rate: $\lambda^\text{in} = D \sum_{j \notin \{v_1,v_2\}} \mu_j$ (mean-field from external neighbours)
- Effective out-rate: $\gamma = D \cdot (\text{out-degree of } v_1 + v_2) / 2$
- Reactions: same as single-voxel but `n_sites=2` (e.g. birth rate $= 2 c_b$)

$$p_N(t + \Delta t) = e^{Q_{2\text{-site}} \Delta t}\, p_N(t)$$

**5. Prolong — product-convolution:**

This is the step where different approaches diverge.  The goal is to distribute probability across fine configurations $(n_1, n_2)$ consistent with the updated coarse totals $p_N^\text{new}$.

There are four prolongation approaches of increasing sophistication.  Which one is used depends on the problem:

**Product-convolution prolongation** — unified approach for all model types:
$$p(n_1, n_2 \mid N) \;\propto\; m_1(n_1) \cdot m_2(N - n_1)$$
where $m_1, m_2$ are the per-voxel marginals extracted from the **post-pre-smooth fine joint state** `sp`.

**Why this is the right prolongation for all systems:**

In PDE multigrid, prolongation exploits *smoothness* — polynomials of degree $p$ reconstruct a solution exactly at $p+1$ points. In the CME, the analog is *quasi-stationarity*: after pre-smoothing for $\tau_\text{pre}$, the within-pair joint distribution relaxes toward its quasi-stationary distribution (QSD), at which point the joint approximately factorises as $p(n_1,n_2) \approx m_1(n_1) \cdot m_2(n_2)$. Product-convolution is therefore exact in the limit $\tau_\text{pre} \to \infty$, for any marginal shape.

**Contains previous approaches as special cases:**

For Poisson marginals $m_k \approx \text{Pois}(\mu_k)$:
$$p(n_1 \mid N) \propto \frac{\mu_1^{n_1}}{n_1!} \cdot \frac{\mu_2^{N-n_1}}{(N-n_1)!} = \binom{N}{n_1}\!\left(\frac{\mu_1}{\mu_1{+}\mu_2}\right)^{n_1}\!\left(\frac{\mu_2}{\mu_1{+}\mu_2}\right)^{N-n_1}$$
which is exactly $\text{Binomial}(N, \mu_1/(\mu_1+\mu_2))$. So for BD wave fronts with unimodal near-Poisson marginals, product-convolution gives the same result as the old mean-weighted Binomial — no regression.

**Comparison across all approaches:**

| Approach | Bistable interface pair | BD wave front | Handles new N? |
|---|---|---|---|
| Symmetric Binomial $(N, 0.5)$ | ✗ Unimodal at saddle $N/2$ | ✗ Wrong weights | ✓ |
| Mean-weighted Binomial $(N, \mu_1/\mu)$ | ✗ Unimodal at wrong point | ✓ Correct | ✓ |
| `prolong_multiplicative` | ✓ Preserves existing asymmetry | ✓ Correct | ✗ Needs seed for new N |
| **`prolong_product`** (current) | ✓ Correct bimodal | ✓ Reduces to Binomial | ✓ Uniform treatment |

**Implementation:** `_evolve_pair_vcycle` now uses product-convolution. Marginals $m_1, m_2$ are extracted from `sp` (post-pre-smooth fine state), the convolution $Z_N = \sum_{n_1} m_1(n_1) m_2(N-n_1)$ is computed for each total $N$, and all fine states are populated uniformly.

The older approaches (`prolong_multiplicative`, `prolong_multiplicative_2s`) remain in `operators.jl` for use by the standalone RDME multigrid solver (`vcycle.jl`) which has a different calling convention (operates on full joint `StateSpace` objects, not per-voxel marginals).

**6. Post-smooth ($\tau_\text{post} = 0.15 \Delta t$):**
Final intra-pair diffusion smoothing of the prolongated state.

**7. Marginalise:**
$$p_{v_1}^{\rm new}(n_1) = \sum_{n_2} p_{12}^\text{new}(n_1, n_2), \qquad p_{v_2}^{\rm new}(n_2) = \sum_{n_1} p_{12}^\text{new}(n_1, n_2)$$

If the joint state exceeds `block_max_states` at any point, the pair is skipped and each voxel falls back to the independent per-voxel `_expv_voxel`.

### Why V-cycle for BD but not Schlögl?

For BD, the wave front propagates because a high-concentration voxel transfers molecules to an empty neighbour — a **joint** event that the product approximation misses.  The V-cycle captures the exact probability of this correlated hop.

For Schlögl, the front advances via per-voxel barrier crossing driven by the asymmetric quasi-potential.  The equil-only cross-basin flux mechanism (Steps 3a–4–5) already implements this correctly at the per-voxel level.  The V-cycle would add cost without proportionate physical gain.

---

## 7. Per-Voxel Generator

The generator $Q_k$ for voxel $k$ with in-rate $\lambda^\text{in}$ and per-molecule out-rate $\gamma$:

$$Q_k[n, n+1] = \lambda^\text{in} + b(n) \quad \text{(upward: diffusion in + local birth)}$$
$$Q_k[n, n-1] = d(n) + \gamma n \quad \text{(downward: local death + diffusion out)}$$
$$Q_k[n, n] = -(Q_k[n,n+1] + Q_k[n,n-1])$$

where $b(n)$, $d(n)$ are the local reaction propensities.

**For BD:** $b(n) = k_b$, $d(n) = k_d n$ → dense or tridiagonal, generator `_single_voxel_generator`.

**For Schlögl:** $b(n) = c_3 + c_1 n(n-1)/2$, $d(n) = c_4 n + c_2 n(n-1)(n-2)/6$ → **tridiagonal** (birth-death chain), returned as `LinearAlgebra.Tridiagonal(dl, d, du)`.

Using `Tridiagonal` reduces per matrix-vector product cost from $O(n^2)$ to $O(n)$ in the Krylov `expv`, which is critical given the large spectral radius ($\|A\| \cdot \Delta t \approx 270$ for Schlögl with $D=2$).

The support of the generator is adaptively clipped to the active support of $p_k$ plus a margin, via:
```julia
lo_mol = max(0, findfirst(>(ε_prune), p) - 1 - margin)
hi_mol = min(n_max, findlast(>(ε_prune), p) - 1 + margin)
# For Schlögl: also force n_un into support
lo_mol = min(lo_mol, n_un - margin)
hi_mol = max(hi_mol, n_un + margin)
```

---

## 8. Wave Front Mechanism Summary

### BD:
```
[empty voxel] ← flux from active high-concentration neighbour (Step 4)
    ↓ [activated when flux ≥ ε_expand] (Step 5)
[active voxel, starts from δ₀]
    ↓ [CME evolves: molecules accumulate, spread toward Poisson(μ_loc)] (Step 2)
    ↓ [TV distance converges] (Step 6)
[equil voxel, Poisson(μ_loc)]
    ↓ [generates flux toward next empty layer] (Step 4)
```

### Schlögl:
```
[equil-hi voxel, δ_{n_hi}, P(hi)=1]
    ↓ [generates cross-basin flux to adjacent equil-lo] (Step 4)
[equil-lo activated, starts from δ_{n_lo}] (Step 5)
    ↓ [CME evolves: driven past n_un by diffusion from equil-hi neighbour] (Step 2)
    ↓ [p(n_un) drops below ε_equil] (Step 3a)
[new equil-hi voxel, δ_{n_hi}, P(hi)=1]
    ↓ [generates cross-basin flux to next equil-lo layer] (Step 4)
```

The two mechanisms are structurally identical: a settled voxel in basin A radiates flux that activates a neighbouring voxel in basin B, which transitions and becomes settled in basin A, then radiates to the next layer.

---

## 9. Visualization: P(hi) Colormap

For bistable models, colour each voxel by the probability of being in the **high state**:

$$P_\text{hi}(k) = \sum_{n > n_\text{un}} p_k(n)$$

- Equil-lo ($\delta_{n_\text{lo}}$): $P_\text{hi} = 0$ → **deep red**
- Equil-hi ($\delta_{n_\text{hi}}$): $P_\text{hi} = 1$ → **deep blue**
- Active/transitioning: $P_\text{hi} \in (0,1)$ → **white**

**Why not mean-based?** Colouring by $\langle n_k \rangle$ makes active voxels with mean $\in (n_\text{un}, n_\text{hi})$ appear "almost blue" — indistinguishable from equil-hi.  This creates the visual impression of scattered isolated high-state nucleation events, which is physically misleading.  These are actually the stochastic transition zone (voxels mid-crossing of the bistable barrier), not spurious nucleation.  The P(hi) map correctly represents them as the white uncertainty ring between the definite red/blue regions.

---

## 10. Key Files

| File | Contents |
|---|---|
| `src/spatial_adaptive/bd_spatial_fsp.jl` | `SpatialFSP` struct; `step!` (standard and V-cycle variants); `_expv_voxel`; `_evolve_pair_vcycle`; `_step_adi_pair!`; all helper functions |
| `src/spatial_adaptive/schlogl_spatial_fsp.jl` | `SchloglModel1D` dispatch: `_single_voxel_generator_range` (returns `Tridiagonal`), `_voxel_down_rate`, `local_birth_rate`, `schlogl_fixed_points` |
| `src/rdme/rdme_model.jl` | Model structs: `BirthDeathRDME`, `BranchingDeathRDME`, `SchloglModel1D`; `schlogl_fixed_points` |
| `examples/anim_bdmultigrid_20x20.jl` | BD 2D wave front animation (V-cycle) |
| `examples/schlogl_2d_droplet.jl` | Schlögl bistable critical droplet demo (P(hi) colormap) |

---

## 11. Parameter Guide for Schlögl Critical Droplet

| Parameter | Effect | Recommended |
|---|---|---|
| `c3` | Thermodynamic bias. `c3=19.9` ≈ Maxwell coexistence; `c3=21.0` strongly biases high state | 20.4–21.0 |
| `D` | Diffusion. Larger → faster wave, wider stochastic front | 1.5–2.0 |
| `ε_expand` | Flux threshold for activation. Too small → cascade; too large → front stalls | 0.15 |
| `ε_equil` | Deactivation threshold (p(n_un) criterion). Too large → spurious deactivation | 0.12 |
| `interface_width` | Initial IC: ring half-width around droplet boundary | 2–3 |
| `T_END` | Larger allows fuller saturation; simulation time ∝ T_END | 80–200 |

Quasi-potential as a function of c3 (D-independent):
| c3 | ΔΦ | Interpretation |
|---|---|---|
| 19.5 | −1.01 | Low state strongly preferred — all droplets collapse |
| 19.9 | ≈ 0 | Maxwell coexistence — Rc → ∞ |
| 20.4 | +0.18 | High state weakly preferred — slow, clean front |
| 20.6 | +0.43 | Moderate bias — intermediate wave speed |
| 21.0 | +0.93 | Strong bias — fast but wide stochastic front |
