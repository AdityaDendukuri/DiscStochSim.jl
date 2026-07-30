# ELSE on a stiff CME (prototype for §4.3)

A chemical master equation example for ELSE: a fast reversible module gating a
rare slow escape — the canonical CME source of the timescale separation ELSE
targets. Code: `examples/fast_slow_else.jl`.

## System

$$A \underset{k_r}{\overset{k_f}{\rightleftharpoons}} B \ (\text{fast}), \qquad B \xrightarrow{k_s} C \ (\text{slow}), \qquad A+B+C = N.$$
Start with $N$ copies of $A$. The $(A,B)$ manifold ($C=0$, a 1D chain in the
$A$-count) equilibrates over many fast flips, then rarely drains one molecule to
$C$. The observable is the **first-escape MFPT** to that drain.

## ELSE

Before the first escape the dynamics live entirely on the $C=0$ manifold, with
$B\to C$ leaving it. Take $R = A_{JJ}$ over the $N{+}1$ manifold states with the
escape absorbed. The fundamental matrix $Z = -R^{-1}$ (ELSE Eq. 8) gives, in one
linear solve,
$$\text{MFPT} = \mathbf 1^\top Z\,e_{\text{init}}, \qquad
\text{exit distribution} \propto (k_s\,b)\odot (Z\,e_{\text{init}}),$$
i.e. the exact expected time to first escape and the exact state you escape from.
Cost is one $(N{+}1)$-state solve — **independent of the stiffness $k_f/k_s$**.

## Results ($N=20$, $k_s=0.05$)

**Exact, certified by FSP.** At $k_f/k_s = 100$: ELSE MFPT $= 2.102495$, flux-FSP
(survival integral $\int P(C{=}0,t)\,dt$ on the full 231-state CME) $= 2.102495$,
relative difference $8\times10^{-9}$. ELSE's subnetwork embedding is exact on the
CME.

**Unbiased where QSSA is not.** The fast-equilibrium (QSSA) MFPT $= 1/(k_s N/2) = 2.0$
is a *constant*; the true MFPT is not. QSSA is badly wrong at low/moderate
stiffness (at $k_f/k_s = 2$, exact $= 4.98$ vs QSSA $2.0$) and only becomes valid
as $k_f/k_s \to \infty$. ELSE is exact at every stiffness. SSA agrees with ELSE
within Monte-Carlo error throughout.

![mfpt](figs/fastslow_mfpt.png)

**Cost is flat in the stiffness SSA chokes on.** SSA must simulate every fast
$A\rightleftharpoons B$ flip: mean reaction events to one escape grow from 11 to
409 across the sweep, while ELSE stays at one $21$-state solve. Crossover at
$k_f/k_s \approx 6$; beyond it ELSE is cheaper, and the gap widens with stiffness.

![cost](figs/fastslow_cost.png)

| $k_f/k_s$ | ELSE MFPT (exact) | SSA MFPT | QSSA | SSA events | ELSE work |
|---|---|---|---|---|---|
| 2   | 4.98 | 5.01 ±0.06 | 2.0 | 11 | 21 |
| 10  | 2.87 | 2.93 ±0.04 | 2.0 | 30 | 21 |
| 50  | 2.20 | 2.23 ±0.04 | 2.0 | 112 | 21 |
| 200 | 2.05 | 2.04 ±0.04 | 2.0 | 409 | 21 |

## Scope (honest)

- The comparison to SSA/QSSA is the real content: ELSE is **exact** (unlike QSSA)
  at a cost **independent of stiffness** (unlike SSA). That is ELSE's pitch, shown
  on a CME.
- At *this* scale ELSE and FSP are both small linear solves, so FSP here is the
  exactness **certificate**, not a competitor — the two do not race. ELSE's unique
  advantage over a direct solve appears only in the regime FSP cannot reach:
  large or many-basin CMEs where the global law is intractable but each
  subnetwork solve is cheap, stitched along a trajectory. That scaling claim
  needs a larger system; this example establishes exactness + the stiffness
  argument that motivates it.
- Contribution scope: this is the worked chemical example plus its FSP-certified
  reference that §4.3 is missing — not a new method.
