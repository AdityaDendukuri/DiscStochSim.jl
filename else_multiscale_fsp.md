# Multiscale FSP for stiff CMEs — observation, dead-ends, future work

Notes from probing whether ELSE's ideas can accelerate a *distribution-resolving*
FSP on stiff systems. Two concrete ideas were tested and refuted; the reason
points to the direction worth pursuing. Recorded as future work.

## Observation: transient FSP is stiffness-bound, like SSA

Fast/slow CME $A \rightleftharpoons B$ (fast, $k_f=k_r=5$), $B \to C$ (slow,
$k_s=0.05$), $N=20$. First-escape MFPT (`examples/fast_slow_else.jl`):

| method | state space | steps / events | wall time |
|---|---|---|---|
| ELSE (fundamental-matrix solve) | 21 | 1 solve | 20 µs |
| flux-FSP (exact, un-pruned) | 231 (full closed simplex) | 15 677 time-steps | 24 s |

flux-FSP is exact but ~$10^6\times$ slower — **not** because of state count (231
is tiny), but because the adaptive step $dt = \varepsilon_{dt}/\Phi$ uses the
*internal* gross flux $\Phi = \mathbf 1^\top\!\operatorname{diag}$-rate. That flux
is dominated by the fast $A\rightleftharpoons B$ churn, which never quiets — even
at equilibrium the gross flipping rate stays $\approx k_f\langle A\rangle + k_r\langle B\rangle \approx 100$.
Measured over the run: $\Phi \in [23,100]$, $dt \in [0.002, 0.008]$. So the step
controller keys on the fast timescale, not the slow observable — the definition of
stiffness. SSA pays this tax in events; transient FSP pays it in time-steps.

## Dead-end 1: boundary-flux step control (refuted)

*Idea.* Drive $dt$ by the **boundary** outflux $\Phi_{\text{out}} = \mathbf 1^\top A_{J'J}p_J$
(the actual truncation-error driver, Cor. 4.7) rather than internal $\Phi$, so
fast internal recirculation stops forcing tiny steps. `build_generator` already
returns `exit_rate_boundary`; the docstring even names this rule.

*Refuted.* The matrix exponential is **not** free over a big stiff step. Krylov
`expv` error grows with $dt\cdot\lVert A\rVert$. On an open stiff system
($A\rightleftharpoons B$, $B\to B+C$), computing $\langle C\rangle(40)$ (true value
$19.95$):

| method | $\langle C\rangle(40)$ |
|---|---|
| dense $\exp(40A)\,p$ / analytic | 19.95 |
| `expv`, Krylov dim $m=30$ (default) | 4.31 |
| `expv`, $m=80$ | 11.48 |
| `expv`, $m=150$ | 19.35 |

Accuracy needs $m \propto dt\cdot\lVert A\rVert$, i.e. **work $\propto T\lVert A\rVert$
regardless of how the interval is chunked.** Big steps just trade the flux-heuristic
wall for the `expv`-accuracy wall $dt\lVert A\rVert \lesssim m$, which is also set
by the fast rate. You cannot step around the fast modes. (An earlier closed-system
"100×" from larger steps was per-step *overhead* amortization — fewer generator
rebuilds / snapshots — not less `expv` work; it held accuracy only while
$dt\lVert A\rVert$ stayed modest.)

## Dead-end 2: exact reinjection boundary closure (refuted)

Closing the FSP boundary with a Schur complement over an exterior shell,
$A_{JJ} + A_{JJ'}(-A_{J'J'})^{-1}A_{J'J}$, routed via ELSE's fundamental matrix.
This is classical *stochastic complementation* (Meyer 1989), not ELSE and not new.
A partial shell buys nothing over simply enlarging the interior and marginalizing;
exactness only appears with the full exterior — i.e. the full solve. For computing
a distribution there is no boundary shortcut. (Details + code:
`examples/boundary_closure.jl`.)

## Future work: multiscale FSP by fast-subnetwork elimination

Both dead-ends share a lesson: **you cannot step or close your way around fast
modes — you must remove them from the propagated operator.** That is exactly what
ELSE does for trajectories (analytic escape via $Z=-A_{JJ}^{-1}$), and its honest
FSP analog is a *multiscale FSP*:

- Partition into fast and slow reaction blocks. Project onto the fast block's
  quasi-stationary manifold and propagate only the slow dynamics.
- Use the fast block's fundamental matrix $Z_{\text{fast}} = -A_{\text{fast}}^{-1}$
  to carry the **exact** fast-fluctuation corrections — i.e. QSSA *plus* the terms
  QSSA drops. (QSSA alone was 5% biased on the first-escape example: exact $4.98$
  vs QSSA $2.0$ at low stiffness.)
- The reduced slow generator is non-stiff, so the FSP step is no longer bounded by
  the fast rate.

Open questions before this is worth building: the construction cost of the reduced
operator, its accuracy vs full FSP across the timescale gap, and the regime where
the gap is large enough to pay for the reduction. This is the genuinely
ELSE-inspired direction; the two tested shortcuts are not.
