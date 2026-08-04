# FESS next steps

## Implementation status

The C++-matched baseline now lives in `src/fess/`; `examples/notebook.jl` is a
thin experiment and plotting client.

- `FESSNetwork` and `FESSWorkspace` own persistent states and their indices.
- `FESSSubnetwork` owns `R`, `Z`, and the boundary exit rates, matching the
  structure of the C++ reference subnetwork.
- `augment!`, `shed!`, `escape_time`, `escape_probabilities`,
  `cut_time_losses`, and `conditional_time` are small independent operations.
- The full fundamental matrix is obtained by one factorized solve against the
  identity. Shedding then uses the same rank-one inverse downdate as the C++
  implementation.
- Shedding removes one state at a time and applies the C++ threshold separately
  to each removal.
- `fess_step!`, `fess_trajectory`, and `fess_ensemble` share one state-first
  exit and conditional-time implementation.
- `fess_ensemble` builds one shared generator per step, solves every unique
  active entrance together, and retains an independent clock and RNG for every
  trajectory.
- `compare_fess_steps` runs the four FESS variants for the same number of steps
  and compares the final physical time they reach.

As in the C++ implementation, capacity is a hard constraint and can force a
cut whose relative loss exceeds the threshold. The forward--adjoint and
cumulative-budget method below remains a proposed alternative to benchmark
against the C++-matched baseline.

## Immediate baseline

1. Keep the redrawn and persistent single-trajectory implementations as the primary comparison.
2. Use identical network-expansion and exit-sampling rules so persistence and shedding are the only differences.
3. Compare ensembles rather than individual paths using cycle period, peak amplitude, mean, variance, FESS steps, wall time, and peak memory.
4. Retain the asynchronous group implementation as a diagnostic, but do not use it for the Lotka--Volterra performance benchmark. The measured solve reuse is approximately `1.0003`, so grouping provides essentially no benefit in this example.

## 1. Remove explicit matrix inverses

Factorize the killed generator once,

```math
-A = LU,
```

and obtain the occupation vectors through linear solves,

```math
u=-A^{-1}p,
\qquad
u^{(2)}=-A^{-1}u.
```

This is more stable and usually faster than constructing `inv(A)`. Use multiple right-hand sides when several entrance states share a network.

## 2. Add numerical accuracy checks

For every solve, monitor the scaled residual

```math
\eta=
\frac{\lVert Au+p\rVert}
{\lVert A\rVert\lVert u\rVert+\lVert p\rVert}.
```

Apply iterative refinement when the residual is too large, and rebuild the factorization if refinement fails. Also record pivot growth and an inexpensive condition estimate.

## 3. Replace inverse-based shedding with forward--adjoint sensitivity

Solve

```math
Au=-p,
\qquad
A^\top v=-\mathbf 1.
```

Use `u` and `v` to estimate how much removing a state changes the residence time. A simple first indicator is

```math
s_j=|u_jv_j|.
```

This requires only one forward and one adjoint solve and avoids forming the dense fundamental matrix.

## 4. Shed against a cumulative error budget

Sort removable states by their estimated effect and remove the largest set satisfying

```math
\sum_{j\in R}\delta T_j
\leq
\varepsilon_{\mathrm{shed}}T.
```

After shedding, recompute the residence time and verify the actual relative change. This prevents many individually small removals from accumulating into a large error.

## 5. Reuse factorizations after network changes

When a small block of states is added, write

```math
A_{\mathrm{new}}=
\begin{bmatrix}
A & B\\
C & D
\end{bmatrix}
```

and reuse the factorization of `A` through the Schur complement

```math
S=D-CA^{-1}B.
```

Factor only the small block `S`. Periodically perform a complete sparse refactorization to control numerical drift and fill-in.

## 6. Improve sparse ordering

Maintain separate logical state indices and numerical factorization indices. Test approximate minimum-degree and nested-dissection orderings, placing recently added or likely-to-be-shed boundary states last. Measure `nnz(L + U)` rather than only the number of states.

## 7. Eliminate historical states instead of deleting them

Partition the persistent space into an active core `C` and historical halo `H`, then form

```math
A_{\mathrm{eff}}
=A_{CC}-A_{CH}A_{HH}^{-1}A_{HC}.
```

This compresses the halo while retaining its algebraic effect. Direct deletion instead turns every transition into a removed state into an immediate network exit, which can change the FESS dynamics.

## 8. Adapt expansion using boundary indicators

Rank boundary states using their outward occupation flux,

```math
\phi_j=u_j\Delta_j,
```

and expand where this indicator is largest. Stop when the unresolved boundary contribution is below an expansion tolerance rather than always using a fixed graph depth.

## 9. Add rebuild criteria

Rebuild the generator and its factorization when any of the following occurs:

- the scaled residual exceeds its tolerance;
- the update count exceeds a fixed limit;
- factorization fill becomes excessive;
- the condition estimate deteriorates;
- an incremental update costs more than a fresh factorization.

## 10. Benchmark each change independently

For every implementation stage, report:

- error in mean residence time;
- cycle-period and peak-amplitude error relative to SSA;
- number and size of FESS steps;
- generator-build, factorization, solve, shedding, and sampling time;
- `nnz(A)` and `nnz(L + U)`;
- residual and condition estimates;
- peak memory and persistent-state count.

After validating the C++-matched baseline, benchmark forward--adjoint shedding
and a cumulative shedding budget as a separate alternative. This keeps changes
to the intended FESS construction distinguishable from numerical optimizations.
