# ELSE on the toggle switch — results

ELSE (Exact Linear System Embedding, Fornace) applied to bistable chemical master
equations, explained from the two methods you already use — SSA and FSP — and
validated end to end, from a toy toggle up to the Tian–Burrage toggle at its real
parameters.

---

## 1. What ELSE is, next to SSA and FSP

A reaction system with two basins. The trajectory sits in one basin for a very
long time, then rarely hops to the other. The question: **when does it hop, and
where does it cross?**

| method | what it does | cost |
|---|---|---|
| **SSA** | simulate every reaction until the hop | one event at a time — millions of in-basin reactions before the rare escape |
| **FSP** | truncate the state space, solve `dp/dt = A·p` for the whole distribution in time | one matrix step per `Δt`, out to some horizon `tf` |
| **ELSE** | restrict the generator to **one basin** `J` and solve **once** | a single linear solve on `|J|` states — cost set by the basin size, **not** by how rare the hop is |

The waste SSA suffers is the millions of fast in-basin reactions before the one
that escapes. ELSE deletes exactly that: it takes the basin as a block and, with
one linear solve, returns the exact escape **time** and exit **state**, then jumps
the trajectory there. The deeper the basin (the rarer the hop), the more SSA
suffers and the more ELSE saves — its cost doesn't change.

FSP is the bridge: ELSE reuses the **same generator** FSP builds, but instead of
evolving it forward in time it integrates it to escape in closed form.

---

## 2. The method — one matrix

Restrict the generator to the basin `J`; every reaction leaving `J` is absorbed
(the "transient subnetwork," Fornace Fig. 1). Call it `A_JJ` (column convention,
`dp/dt = A_JJ p`), with per-state escape rate `Δ_j = Σ_{y∉J} rate(j→y)`.

The one object ELSE is built on is the **fundamental matrix**

```math
Z = -\mathbf{A}_{JJ}^{-1} = \int_0^\infty e^{t\mathbf{A}_{JJ}}\,dt \ \ge 0
\qquad\text{(Fornace Eq. 8).}
```

`Z_ij` = expected time spent in state `j` before escaping, started at `i`. We only
ever need `Z p0 = -A_JJ⁻¹ p0` — one solve — and everything reads off it:

| quantity | formula | Fornace |
|---|---|---|
| mean escape time | `1ᵀ Z p0` | Eq. 15 |
| moments / variance | `E[τ^m] = m! · 1ᵀ Z^m p0` | Eq. 14 |
| exit-state law | `π_exit ∝ Δ ∘ (Z p0)` | Eq. 10 |
| full escape law | survival `Γ_J = 1ᵀ e^{tA_JJ} p0`, density `f = Δᵀ e^{tA_JJ} p0` | (†) |

The escape density `f` is the boundary outflux `Φ_out` — *what FSP treats as
truncation error to suppress is, read off the transient subnetwork, exactly ELSE's
answer.*

---

## 3. Validation on a toy toggle

Symmetric toggle (`b=20, K=10, n=3, d=1`), basin `J = {x ≥ y}` (`|J| = 1081`),
escape = crossing the separatrix `x = y`. Reference: plain Gillespie.

**ELSE reproduces the whole switch-time distribution, not just the mean.** One
eigendecomposition of `A_JJ` is shared across the batch; each switch is then one
draw (invert the survival + pick the exit state). Gillespie needs ~**1381 reaction
events per switch**; ELSE needs one shared solve.

- mean switch time: ELSE **31.11**, Gillespie **31.25 ± 0.4** — statistically exact.
- exit-state distribution: total variation **0.009** vs 20k Gillespie crossings.

![switch-time CDF](figs/toggle_sampler_cdf.png)
![exit state](figs/toggle_sampler_exit.png)

---

## 4. Scaling the sampler — the resolvent (Laplace transform)

The eigendecomposition above is **dense** and won't scale. The escape law lives in
the **resolvent** `(sI − A_JJ)⁻¹`, which is the Laplace transform of the semigroup:

```math
(sI - \mathbf{A}_{JJ})^{-1}\,p_0 = \int_0^\infty e^{-st}\,e^{t\mathbf{A}_{JJ}}p_0\,dt .
```

`s = 0` gives `Z` (mean, exit). Sweeping `s` gives the whole distribution: the
escape-density transform is `f̂(s) = Δᵀ (sI − A_JJ)⁻¹ p0`. We evaluate it at a fixed
set of Bromwich nodes `s_k = a + i·kπ/T` (**one sparse solve each, shared across all
times and all samples**), invert by the Fourier-series method to get `f(t)`,
integrate to the CDF, and sample. No eigendecomposition — only sparse solves.

*(Detail that matters: invert the density `f`, not the survival. Started deep in
the well `Δᵀp0 = 0` and many derivatives vanish, so `f` is smooth with `f(0)=0` and
`f̂(s)` decays fast — a few hundred fixed nodes suffice. Keep `a·T` mild so the
`e^{at}` inversion factor doesn't amplify ripple into a spurious tail.)*

On the toy (501 shared resolvent solves) the resolvent density lands exactly on the
eigendecomposition, and the sampled CDF on Gillespie (mean **31.09**).

![escape density via resolvent](figs/toggle_laplace_density.png)
![sampled CDF](figs/toggle_laplace_cdf.png)

---

## 5. The real toggle (Tian & Burrage)

The actual model, at its actual parameters — `U`/`V` mutual repression, Hill
coefficient 3, production `≈420`, degradation `≈1`, `u0 = (85,5)`. Basin `J = U-high
well = {U ≥ V}` on a grid large enough to contain the well (mode `≈ (378,27)`), so
the only escape is the `U=V` crossing.

| | value |
|---|---|
| basin size `|J|` | **232,221 states** |
| generator build | 0.56 s |
| **one sparse solve** `-(A_JJ \ p0)` | **1.7 s** |
| escape flux across the separatrix `U=V` | **100.0 %** (no spurious grid-edge escape) |
| mean crossing coordinate `U (= V)` | **≈ 132** (the saddle) |
| **mean hop time `U`-high → `V`-high** | **3.4 × 10¹⁴** |

That last number is the point:

| method | can it reach the hop? |
|---|---|
| **your FSP** (`tf = 30`) | no — the hop is ~**10¹³×** past the horizon; the `V`-high basin stays empty |
| **SSA** | no — ~**10¹⁴ reaction events per trajectory**, times many trajectories |
| **ELSE** | **yes — exact mean in one 1.7 s solve on 232k states** |

ELSE delivers a quantity that is flatly inaccessible to SSA and FSP, on the real
model, in under two seconds.

---

## 6. Scope and next steps

- **Mean + exit, toy and real:** done, exact, one sparse solve each. ✓
- **Full distribution, toy:** done via the resolvent (sparse solves, no
  eigendecomposition), validated against both references. ✓
- **Full distribution, real:** the same resolvent route, but ~hundreds of *complex*
  232k-state solves — minutes, not seconds. Feasible; not yet run.
- **Compactness (the open lever):** `|J| = 232k` is the entire `U≥V` half-grid, but
  the occupation is a tight blob at the mode plus a thin channel to the saddle. A
  much smaller `J` — a box around the well reaching just past the separatrix —
  should give the same `3.4×10¹⁴` far cheaper. This is the experiment that decides
  how cheap ELSE really is here.

**Code:** `examples/else_toggle.jl` (walk-through notebook), `examples/toggle_laplace.jl`
(resolvent sampler), `examples/toggle_real_else.jl` (real toggle).
