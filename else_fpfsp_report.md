# State-first ELSE and Flux-Preserving FSP

## Scope

This report develops the exact linear system embedding (ELSE) method as a
replacement for reaction-by-reaction stochastic simulation. It then places
ELSE in the notation of finite state projection (FSP) and flux-preserving FSP
(FP-FSP).

The implementation discussed here is the state-first method:

1. choose a finite subnetwork;
2. factorize its killed generator;
3. use the fundamental matrix to sample the state from which escape occurs;
4. use a conditional mean for the elapsed time;
5. sample the physical reaction that leaves the subnetwork;
6. restart from the resulting state.

The method does not compute a matrix exponential. It also does not simulate
hundreds or millions of trajectories independently. For an ensemble, it
constructs multinomial count columns while reusing the same linear
factorization.

The implementation is divided into:

- the model-independent code in `src/else.jl`;
- the Burrage toggle example in
  `examples/burrage_else_state_first_notebook.jl`.

## 1. Chemical master equation

Consider a continuous-time Markov chain with states

$$
x\in\mathcal X\subseteq\mathbb Z_{\geq 0}^d.
$$

Reaction channel $k$ has stoichiometric change $\nu_k$ and propensity
$\alpha_k(x)$. A firing produces

$$
x\longmapsto x+\nu_k.
$$

With a column probability vector $p(t)$, the chemical master equation is

$$
\frac{d}{dt}p(t)=A p(t).
$$

The generator entries follow the column convention

$$
A_{y x}
=
\begin{cases}
\alpha_k(x), & y=x+\nu_k,\\[1mm]
-\displaystyle\sum_k\alpha_k(x), & y=x,\\
0, & \text{otherwise}.
\end{cases}
$$

Thus $A_{yx}>0$ is the rate of the transition from source state $x$ to
destination state $y$, and

$$
\mathbf 1^\top A=0.
$$

At state $x$, the stochastic simulation algorithm (SSA) uses the total
activity

$$
w(x)=\sum_k\alpha_k(x)=-A_{xx}.
$$

SSA samples one holding time and one reaction at a time. This becomes expensive
when a trajectory makes many internal reactions before a rare transition of
interest.

## 2. A finite ELSE subnetwork

Choose a finite set

$$
J=\{x_1,\ldots,x_n\}\subset\mathcal X
$$

containing the current state. The set may be:

- a metastable basin;
- a fast reaction class;
- a local neighborhood;
- a union of disconnected components;
- any finite state collection for which all exits can be enumerated.

ELSE does not require the set to be a basin. The set affects computational
efficiency, not the algebra of the first-exit calculation.

### 2.1 Killed generator

Let $R$ be the restriction of the full generator to $J$, retaining the full
physical diagonal:

$$
R=A_{JJ}.
$$

If a reaction leaves $J$, its rate remains part of the diagonal loss in $R$,
but its destination does not appear as an off-diagonal entry of $R$.
Consequently,

$$
\mathbf 1^\top R\leq 0.
$$

This is the killed generator. It must not be replaced by a reflected or
mass-conserving generator.

### 2.2 Exit-channel table

For every internal state $b\in J$ and external destination $y\notin J$, define

$$
B_{yb}
$$

as the total physical rate of transitions from $b$ to $y$. In a reaction
network,

$$
B_{yb}
=
\sum_{k:b+\nu_k=y}\alpha_k(b).
$$

The total escape rate from internal state $b$ is

$$
\Delta_b=\sum_{y\notin J}B_{yb}.
$$

The exit-channel table and killed generator must balance:

$$
-\mathbf 1^\top R
=
\Delta^\top.
$$

Equivalently, for every column $b$,

$$
-\sum_{j\in J}R_{jb}
=
\sum_{y\notin J}B_{yb}.
$$

This condition is important. It ensures that every loss of probability from
$R$ corresponds to a known physical exit. The general Julia constructor checks
this equality.

## 3. The fundamental matrix

Assume the process eventually leaves $J$. Define

$$
Z=-R^{-1}.
$$

The implementation does not form the inverse. It factorizes $-R$ once and
performs linear solves.

For an entry state $a\in J$, let $e_a$ be the corresponding coordinate vector.
The occupation vector is obtained from

$$
(-R)z=e_a,
$$

or equivalently

$$
z=Ze_a.
$$

Its component $z_b$ is the expected time spent in state $b$ before leaving the
subnetwork:

$$
z_b
=
\mathbb E_a
\left[
\int_0^\tau
\mathbf 1_{\{X(t)=b\}}\,dt
\right],
$$

where $\tau$ is the first-exit time from $J$.

The mean exit time is therefore

$$
\mathbb E_a[\tau]
=
\mathbf 1^\top z.
$$

### 3.1 Exit-state probabilities

While the process occupies state $b$, it escapes at rate $\Delta_b$. Integrated
over the complete visit to the subnetwork,

$$
\Pr(B_{\mathrm{pre}}=b\mid a)
=
\Delta_b z_b.
$$

The probabilities sum to one:

$$
\sum_{b\in J}\Delta_bz_b
=
\Delta^\top Z e_a
=
1.
$$

Only states with $\Delta_b>0$ need to be included in the categorical or
multinomial distribution. These are the boundary states of the subnetwork.

### 3.2 Physical exit destination

After sampling the pre-exit state $b$, sample its physical exit destination
$y$ using

$$
\Pr(Y_{\mathrm{exit}}=y\mid B_{\mathrm{pre}}=b)
=
\frac{B_{yb}}{\Delta_b}.
$$

The complete joint exit law is consequently

$$
\Pr(B_{\mathrm{pre}}=b,Y_{\mathrm{exit}}=y\mid a)
=
B_{yb}z_b.
$$

This formulation is general. It does not refer to a toggle switch, a diagonal
separatrix, or exactly two basins.

## 4. State-first elapsed time

The state-first approximation samples the pre-exit state before assigning an
elapsed time. Compute a second occupation vector from another solve:

$$
(-R)z^{(2)}=z,
$$

so that

$$
z^{(2)}=Zz=Z^2e_a.
$$

For a sampled pre-exit state $b$,

$$
\mathbb E_a[\tau\mid B_{\mathrm{pre}}=b]
=
\frac{z_b^{(2)}}{z_b}.
$$

The cancellation of $\Delta_b$ can be seen from

$$
\mathbb E_a[
\tau\mathbf 1_{\{B_{\mathrm{pre}}=b\}}
]
=
\Delta_b z_b^{(2)}
$$

and

$$
\Pr(B_{\mathrm{pre}}=b\mid a)
=
\Delta_bz_b.
$$

The simulated ELSE leg uses

$$
\widehat\tau_b
=
\frac{z_b^{(2)}}{z_b}.
$$

### 4.1 What is exact

The following quantities are exact up to numerical linear-solve error:

- expected occupation $z$;
- mean first-exit time $\mathbf 1^\top z$;
- pre-exit-state probabilities $\Delta_bz_b$;
- physical exit-channel probabilities;
- conditional mean time $z_b^{(2)}/z_b$.

### 4.2 What is approximated

The algorithm does not sample the complete conditional distribution

$$
\tau\mid B_{\mathrm{pre}}=b.
$$

It replaces that random time by its conditional mean. Therefore:

- the exit-state sequence has the correct transition law;
- elapsed times are correct in conditional expectation;
- ensemble mean occupations and mean passage times are preserved;
- exit-time variances and other fine temporal statistics are not preserved.

This is the intended bias-variance trade:

- temporal variance is reduced by replacing each conditional distribution with
  one value;
- the relevant conditional expectations remain correct.

No claim should be made that the resulting time skeleton has the full SSA
pathwise timing law.

## 5. General single-trajectory algorithm

### Inputs

A general ELSE step needs:

1. the internal state list $J$;
2. the killed generator $R$;
3. all external destinations and rates $B_{yb}$;
4. an entry state $a\in J$;
5. a random-number generator.

### Precomputation

For a reused subnetwork:

1. verify the balance between $R$ and the exit table;
2. factorize $-R$;
3. record the boundary indices and their outgoing channels.

The factorization is reused for every entry state in that subnetwork.

### One state-first ELSE step

Starting from $a$:

1. solve

   $$
   (-R)z=e_a;
   $$

2. solve

   $$
   (-R)z^{(2)}=z;
   $$

3. form the boundary-state distribution

   $$
   \pi_b=\Delta_bz_b;
   $$

4. sample the pre-exit state

   $$
   b\sim\operatorname{Categorical}(\pi);
   $$

5. assign the conditional mean duration

   $$
   \widehat\tau=\frac{z_b^{(2)}}{z_b};
   $$

6. sample a physical destination with probabilities

   $$
   \rho_{y\mid b}=\frac{B_{yb}}{\Delta_b};
   $$

7. update

   $$
   t\leftarrow t+\widehat\tau,
   \qquad
   x\leftarrow y;
   $$

8. select a subnetwork containing $y$ and repeat.

### Pseudocode

```text
ELSETrajectory(subnetwork_for, initial_state, number_of_steps)
    state <- initial_state
    time  <- 0

    repeat number_of_steps times
        E <- subnetwork_for(state)

        z  <- solve(-E.R, point_mass(state))
        z2 <- solve(-E.R, z)

        boundary_probability <- E.exit_rate .* z
        b <- categorical(boundary_probability)

        waiting_time <- z2[b] / z[b]
        y <- categorical(E.exit_channels[b])

        time  <- time + waiting_time
        state <- y
    end
```

There is no assumption that the next subnetwork differs from the previous one.
The destination state itself determines what happens next.

## 6. Simulating a population together

Suppose $n_a$ trajectories currently enter the same subnetwork at state $a$.
Computing and sampling every trajectory independently would repeat identical
linear algebra.

Instead, construct one probability column for entry state $a$:

$$
\pi_{ba}
=
\Delta_bZ_{ba}.
$$

Then draw all pre-exit-state counts together:

$$
M_{\cdot a}
\sim
\operatorname{Multinomial}
\left(
n_a,\{\pi_{ba}\}_{b\in\partial J}
\right).
$$

Here:

- column $a$ corresponds to one occupied entry state;
- row $b$ corresponds to one pre-exit boundary state;
- $M_{ba}$ is the number of trajectories assigned to $b$.

If several entry states are occupied, construct one column per occupied entry
state:

$$
M
=
\begin{bmatrix}
M_{\cdot a_1} &
M_{\cdot a_2} &
\cdots
\end{bmatrix}.
$$

The number of random draws is controlled by the number of occupied entry
states and boundary states, not directly by the number of trajectories.

### 6.1 Sampling the physical exit channels

For the $M_{ba}$ trajectories assigned to pre-exit state $b$, draw another
small multinomial across its physical exit channels:

$$
C_{\cdot b}
\sim
\operatorname{Multinomial}
\left(
M_{ba},
\left\{\frac{B_{yb}}{\Delta_b}\right\}_y
\right).
$$

The destination counts are accumulated into a new dictionary:

$$
n_y^{\mathrm{next}}
=
\sum_{a,b}C_{yba}.
$$

Trajectories that reach the same destination are immediately regrouped.

### 6.2 Reusing the linear algebra

For each subnetwork, the expensive operation is the factorization of $-R$.
After that:

- each occupied entry state requires two triangular solves;
- each population column requires one multinomial draw;
- each occupied boundary state requires one small exit-channel multinomial.

For a population of one million trajectories, the code does not execute one
million ELSE loops. It manipulates integer counts.

### 6.3 Populations spread over several subnetworks

At the beginning of a macro-step, group states by

```julia
subnetwork_for(state)
```

Each group uses the factorization belonging to that subnetwork. The resulting
destination counts are merged afterward. This supports:

- more than two basins;
- disconnected subnetworks;
- different entry states in the same subnetwork;
- trajectories that remain in the same named subnetwork;
- model-specific or adaptive subnetwork-selection rules.

No Boolean variable is used to force an alternation.

### 6.4 Ensemble occupation

For an entry group of size $n_a$, its expected accumulated occupation is

$$
n_aZe_a.
$$

Summing this vector over all entry groups and macro-steps gives the expected
time spent in every internal state by the complete population.

This occupation can be visualized directly. If different basins have very
different mean residence times, a globally normalized heatmap may hide the
shorter-lived basin. Normalizing each basin by its own peak displays basin
shape, but no longer displays their relative time weights. The figure must say
which normalization is being used.

## 7. The reusable Julia interface

The general implementation is centered on:

```julia
ELSESubnetwork(states, R, exits; name=:subnetwork)
```

The arguments are:

- `states`: ordered internal state list;
- `R`: killed generator in the same ordering;
- `exits[state]`: `next_state => rate` pairs;
- `name`: identifier used when grouping populations and occupations.

The exported operations are:

```julia
else_law(subnetwork, start_state)
else_step(subnetwork, start_state; rng)
else_trajectory(subnetwork_for, initial_state, number_of_steps; rng)
else_population_step(subnetwork, groups; rng)
else_population(subnetwork_for, initial_state, population, number_of_steps; rng)
```

`else_law` performs the two solves and returns:

- `occupation`;
- `exit_probability`;
- `conditional_time`.

`else_step` samples one state-first macro-event.

`else_population_step` creates the multinomial matrix $M$ for grouped entry
states sharing one subnetwork.

`else_population` first partitions a mixed population by subnetwork, advances
each partition, and merges the destination counts.

## 8. FP-FSP mathematics

FSP and ELSE begin with the same CME generator, but they ask different
questions.

### 8.1 Standard killed FSP

For an active set $J$, standard FSP evolves

$$
\frac{d}{dt}p_J(t)=A_{JJ}p_J(t).
$$

Because transitions leaving $J$ are killed,

$$
\mathbf 1^\top p_J(t)\leq 1.
$$

The missing mass is the probability that has escaped the finite projection by
time $t$:

$$
\varepsilon_J(t)
=
1-\mathbf 1^\top p_J(t).
$$

FSP expands $J$ until this leakage is sufficiently small.

### 8.2 Mass-conserving finite generator

FP-FSP propagates a normalized approximation on an active set. If

$$
\Delta_b=-\sum_{j\in J}(A_{JJ})_{jb},
$$

then the mass-conserving compressed generator is

$$
\widetilde A_J
=
A_{JJ}+\operatorname{diag}(\Delta).
$$

It satisfies

$$
\mathbf 1^\top\widetilde A_J=0.
$$

This compressed operator is suitable for the FP-FSP distribution update. It is
not suitable for ELSE because it suppresses the escape that ELSE needs to
sample.

### 8.3 State activity and total flux

For probability $p_b(t)$ at state $b$, define its reaction activity

$$
\Phi_b(t)=p_b(t)w_b.
$$

The total flux is

$$
\Phi_{\mathrm{total}}(t)
=
\sum_{b\in J}p_b(t)w_b.
$$

It is the expected SSA event rate under the current distribution.

Boundary outflux is different:

$$
\Phi_{\mathrm{out}}(t)
=
\sum_{b\in J}p_b(t)\Delta_b.
$$

Total flux counts all reaction activity. Boundary outflux counts only
transitions leaving the active set.

### 8.4 Adaptive FP-FSP

The FP-FSP procedure repeatedly:

1. expands the active set along reaction-reachable directions;
2. constructs the finite generator;
3. selects a timestep from the total flux;
4. propagates the finite distribution;
5. prunes low-probability states while protecting important flux carriers.

A representative flux-adaptive timestep is

$$
\Delta t
=
\frac{\varepsilon_{\Delta t}}
{\Phi_{\mathrm{total}}}.
$$

For pruning, probability alone is insufficient. A low-probability state may
have large $p_bw_b$ and carry important current between regions. FP-FSP
therefore protects states with significant flux.

If a set of states with total probability mass $m$ is removed and the retained
distribution is renormalized, the restriction and renormalization contribution
to the $L^1$ error is

$$
2m.
$$

Boundary interaction, pruning, and numerical propagation error together
determine the FP-FSP approximation error.

## 9. ELSE in FP-FSP language

The direct correspondence is:

| FP-FSP concept | State-first ELSE interpretation |
|---|---|
| Active set $J$ | One reusable subnetwork |
| Killed restriction $A_{JJ}$ | ELSE generator $R$ |
| Boundary rate $\Delta_b$ | Rate of a macro-exit from state $b$ |
| Integrated finite-state mass | Expected occupation $z$ |
| Boundary leakage | Probability assigned to an ELSE exit |
| Expansion of $J$ | A potentially larger macro-step |
| Pruning of $J$ | More transitions treated as macro-exits |
| State flux | Guide for retaining useful internal connector states |

The crucial difference is interpretive:

- FSP regards escape from $J$ as truncation leakage to be reduced;
- ELSE regards escape from $J$ as the next simulated macro-event.

### 9.1 Integrated total flux

For entry state $a$, the expected occupation is

$$
z=Ze_a.
$$

The integrated reaction activity is

$$
\mathcal N_J
=
\sum_{b\in J}w_bz_b.
$$

This is the expected number of SSA reaction firings represented by one ELSE
leg. Meanwhile,

$$
\sum_b\Delta_bz_b=1,
$$

because there is one final exit.

Thus $\mathcal N_J$ is a useful acceleration diagnostic:

- if $\mathcal N_J$ is close to one, ELSE skips little SSA work;
- if $\mathcal N_J$ is very large, one ELSE leg replaces many SSA events;
- the benefit must still exceed the factorization and solve cost.

### 9.2 FP-FSP-guided subnetwork selection

The integrated analogue of instantaneous FP-FSP state flux is

$$
\Psi_b=w_bz_b.
$$

It is the expected number of reactions originating from state $b$ during one
ELSE leg. This suggests a general adaptive construction:

1. initialize $J$ around the current state;
2. build $R$ and the exit table;
3. solve for $z$;
4. expand toward states with large occupation or boundary interaction;
5. protect connector states with large $\Psi_b$;
6. stop when additional skipped SSA work is not worth the larger factorization.

Changing $J$ changes where macro-events are placed. If every exit is accounted
for, it does not change the pre-exit-state law for that chosen first-exit
problem.

## 10. Burrage toggle as an application

The toggle-specific notebook defines

$$
\emptyset
\xrightarrow{
\eta\left(
\alpha_1+\frac{\beta_1K_1^3}{K_1^3+V^3}
\right)}
U,
$$

$$
U\xrightarrow{d_1+d_U}\emptyset,
$$

$$
\emptyset
\xrightarrow{
\eta\left(
\alpha_2+\frac{\beta_2K_2^3}{K_2^3+U^3}
\right)}
V,
$$

and

$$
V\xrightarrow{d_2}\emptyset.
$$

The example chooses

$$
J_U=\{(u,v):u\geq v\}
$$

and

$$
J_V=\{(u,v):v\geq u\},
$$

within a large computational grid.

These choices are not part of ELSE. They are only the model-specific
subnetwork definition.

The notebook supplies:

1. `toggle_transitions(state)`, which enumerates physical destinations and
   rates;
2. `make_basin`, which creates $J$, builds $R$, and records every exit;
3. `basin_for(state, basins)`, which chooses the next subnetwork;
4. plotting code that projects a state to its $(U,V)$ coordinates.

Everything involving the $Z$ solves, categorical sampling, multinomial matrix
$M$, regrouping, and occupation accumulation is in the general implementation.

The two-basin heatmap normalizes the two occupations separately:

$$
\widetilde z_U
=
\frac{z_U}{\max z_U},
\qquad
\widetilde z_V
=
\frac{z_V}{\max z_V}.
$$

This displays both basin shapes. It is not a globally normalized stationary
distribution and should not be used to compare the relative residence times of
the two basins.

## 11. Correctness and limitations

### 11.1 Required conditions

For the state-first construction:

1. $R$ must retain the complete physical diagonal;
2. all transitions leaving $J$ must be included in the exit-channel table;
3. the process must eventually leave $J$;
4. the entry state must belong to $J$;
5. the next step must begin at the sampled physical destination;
6. the factorization and solves must be numerically accurate.

### 11.2 Computational truncation

If a model uses a finite grid to represent an otherwise unbounded state space,
reactions crossing the outer grid are also exits. They must be included in the
exit table.

If such an outer-grid exit is sampled, the next subnetwork must contain that
destination. Otherwise the grid must be enlarged. Silently renormalizing away
outer-grid loss would change the exit probabilities.

### 11.3 Timing limitation

The state-first method preserves conditional mean times, not the complete
first-exit-time distribution. It is suited to:

- mean trajectories;
- ensemble occupations;
- mean switching times;
- discrete basin or exit sequences;
- large population transport between subnetworks.

It is not sufficient by itself for:

- exit-time variances;
- temporal quantiles;
- high-frequency path observations;
- exact SSA event histories.

Recovering those quantities requires additional temporal machinery, such as a
transform-based convolution or conditional time sampling. That machinery is
separate from the simple state-first algorithm.

## 12. Summary

The general ELSE object is

$$
\mathcal E=(J,R,B).
$$

Its reusable factorization supplies

$$
z=Ze_a,
\qquad
z^{(2)}=Z^2e_a,
\qquad
Z=-R^{-1}.
$$

The state-first law is

$$
\Pr(B_{\mathrm{pre}}=b\mid a)
=
\Delta_bz_b,
$$

the conditional mean time is

$$
\mathbb E[\tau\mid B_{\mathrm{pre}}=b,a]
=
\frac{z_b^{(2)}}{z_b},
$$

and the physical exit law is

$$
\Pr(Y_{\mathrm{exit}}=y\mid B_{\mathrm{pre}}=b)
=
\frac{B_{yb}}{\Delta_b}.
$$

For a grouped population, one multinomial column replaces independent
trajectory sampling:

$$
M_{\cdot a}
\sim
\operatorname{Multinomial}
\left(
n_a,\{\Delta_bZ_{ba}\}_b
\right).
$$

The method is general because no part of these equations depends on a toggle
switch. A model supplies states, a killed generator, physical exit channels,
and a rule for choosing the next subnetwork. ELSE supplies the linear solves
and the state-first sampling.
