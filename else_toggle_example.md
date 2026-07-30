# ELSE on a multistable CME: the genetic toggle switch (scaling side of §4.3)

A second CME example for ELSE, complementary to the 1D fast/slow chain
(`else_cme_example.md`). Here the escape is a genuine **2D basin escape** across
a **separatrix**, so ELSE's identity is exercised on *both* of its outputs — the
first-passage time (WHEN) and the exit distribution (WHERE). Code:
`examples/toggle_else.jl`.

## System

Symmetric genetic toggle switch: two mutually repressing genes/proteins,
$$\varnothing \xrightarrow{\ h(Y)\ } X,\quad
  \varnothing \xrightarrow{\ h(X)\ } Y,\quad
  X \xrightarrow{d} \varnothing,\quad
  Y \xrightarrow{d} \varnothing,\qquad
  h(z) = \frac{b\,K^{n}}{K^{n}+z^{n}}.$$
For $b=20,\,K=10,\,d=1$ and cooperativity $n\gtrsim 2$ the system is bistable:
an X-high/Y-low basin and an X-low/Y-high basin, separated by the diagonal
$x=y$. Within a basin the 2D fluctuations equilibrate fast; a rare noise
excursion crosses the separatrix and commits to the other basin.

## ELSE

Take the X-high basin $J=\{x\ge y\}$ on the grid $0\!:\!G$, with the separatrix
crossing $x<y$ absorbed: $R=A_{JJ}$. The fundamental matrix $Z=-R^{-1}$ (ELSE
Eq. 8) gives, in one linear solve on the $|J|$ basin states,
$$\text{MFPT} = \mathbf 1^\top Z\,e_{\text{init}}, \qquad
  \text{exit distribution} \propto (\text{cross-rate})\odot(Z\,e_{\text{init}}),$$
i.e. the exact expected time to first cross the separatrix **and** the exact
state you cross *from*. Crossings happen only from the diagonal $x=y$ (a Y-birth
$y\!\to\!y{+}1$ or an X-death $x\!\to\!x{-}1$), so the exit distribution is a
genuine distribution over the saddle region. Cost is one $|J|$-state solve —
**independent of the barrier height**.

## Transient vs recurrent (why the reference is SSA, not FSP)

In the fast/slow chain, $C$ is a *true absorbing sink* ($B\to C$ irreversible),
so the FSP survival integral $\int P(C{=}0,t)\,dt$ on the closed CME **is** the
exact MFPT. The toggle is different: it is **recurrent**. On the closed CME
$P(x\ge y)$ relaxes to its stationary plateau ($\approx\tfrac12$ by symmetry),
not to $0$, because trajectories re-cross — the survival integral diverges. The
separatrix crossing is a first-passage on a recurrent chain, whose exact
reference is **SSA first-passage** (equivalently, the FSP of the basin with the
crossing made absorbing — but that operator *is* the $R$ ELSE inverts, so it
cannot be an independent check). This is exactly the transient/recurrent split
of the resolvent theory in `else_theory.md`: $s=0$ regular (absorbing, fundamental
matrix) vs $s=0$ pole (recurrent, deviation matrix).

## Results ($b=20,\,K=10,\,d=1$, grid $0\!:\!45$, $|J|=1081$)

**Exact — MFPT.** At $n=3$ (init at the X-high mode $(20,2)$): ELSE MFPT
$=31.107$; SSA first-passage over $2\times10^4$ trajectories $=31.14\pm0.21$,
a $0.14\sigma$ discrepancy — statistically indistinguishable. ELSE's basin
embedding is exact on the CME.

**Exact — exit distribution.** ELSE's exit distribution $(\text{cross-rate})\odot
(Z e_{\text{init}})$ over the separatrix matches the empirical crossing histogram
from $2\times10^4$ SSA crossings to **total variation $0.009$** (sampling-limited).
The distribution peaks at $x\approx 9\text{–}10$, right at the deterministic saddle
$x=y\approx 9.7$ — ELSE says *where* you leave, exactly.

![exit](figs/toggle_exit.png)

**Cost is flat in the barrier SSA chokes on.** Deepening the barrier via the Hill
cooperativity $n$ — with the grid and basin held **fixed** — leaves ELSE at one
$1081$-state solve while the SSA event count to a single crossing grows from 36
to 3533. MFPT agreement holds within Monte-Carlo error throughout. Crossover at
$n\approx 2.8$; beyond it ELSE is cheaper and the gap widens with the barrier.

![cost](figs/toggle_cost.png)

| $n$ | ELSE MFPT (exact) | SSA MFPT | SSA events | ELSE work |
|---|---|---|---|---|
| 2.0 | 0.872 | 0.849 ±0.08 | 36   | 1081 |
| 2.5 | 16.38 | 16.32 ±0.38 | 721  | 1081 |
| 3.0 | 31.11 | 31.01 ±0.73 | 1369 | 1081 |
| 3.5 | 52.62 | 51.58 ±1.3  | 2256 | 1081 |
| 4.0 | 79.62 | 81.76 ±2.1  | 3533 | 1081 |

## Scope (honest)

- This is the 2D companion to the 1D fast/slow example: it upgrades the escape
  from a single chain to a genuine basin with a separatrix, and adds the **exit
  distribution** — ELSE's second output — as a certified artifact. The hard axis
  is **metastability** (barrier height) rather than reversible stiffness, and the
  exact reference is **SSA first-passage** because the chain is recurrent.
- As with the 1D example, at *this* scale ELSE and a direct basin solve are both
  small; ELSE's unique advantage over a global solve appears only when the joint
  law is intractable but each basin stays local — e.g. a network with many
  weakly-coupled bistable modules whose $2^M$ attractors make the global FSP grid
  exponential while each module's escape is one small ELSE solve. That
  multi-basin scaling is the future direction in `else_multiscale_fsp.md`; this
  example establishes exactness of *both* ELSE outputs on a real multistable CME
  plus the metastability argument that motivates it.
- Contribution scope: worked chemical example + independent SSA-certified
  reference for §4.3 — not a new method.
