# Fundamental Matrix-Analysis Improvements for FESS

## 1. Current mathematical construction

Let \(R\) be the killed generator on a finite FESS subnetwork. The Julia implementation uses the column-oriented convention

$$
\frac{d p(t)}{dt} = R p(t),
$$

and defines the fundamental matrix

$$
Z = (-R)^{-1}.
$$

For an entrance distribution \(p\), the unconditional occupation vector is

$$
u = Zp,
$$

and the mean first-exit time is

$$
\mathbb{E}[T] = \mathbf{1}^{\mathsf T} Zp.
$$

Let \(q_j\) be the total rate of leaving the subnetwork from state \(j\). The probability that state \(j\) is the state immediately before exit is

$$
\mathbb{P}(B=j) = q_j u_j.
$$

The second occupation vector is

$$
u^{(2)} = Zu = Z^2p.
$$

After sampling the pre-exit state \(B=j\), the current method advances the clock by

$$
\Delta t_j
= \frac{u^{(2)}_j}{u_j}
= \mathbb{E}[T \mid B=j].
$$

Thus, the current FESS construction preserves the pre-exit-state distribution and the conditional mean first-exit time, but replaces the random conditional exit time with its mean.

## 1.1 Why the fundamental matrix contains the ELSE quantities

The formulas above follow directly from the killed Markov semigroup. This section derives them and gives a probabilistic interpretation of the entries of \(Z\).

### 1.1.1 The killed semigroup

Let \(S\) be the states currently inside the FESS subnetwork, and define the first-exit time

$$
T = \inf\{t\geq 0:X_t\notin S\}.
$$

The killed generator \(R\) retains all transitions between states in \(S\), while transitions from \(S\) to its complement appear as losses in the column sums of \(R\).

If the initial distribution is \(p\), then

$$
v(t)=e^{Rt}p
$$

is a substochastic vector. Its \(i\)-th entry has the interpretation

$$
v_i(t)
=
\mathbb{P}(X_t=i,\ T>t).
$$

The missing mass is the probability that the trajectory has already left the subnetwork:

$$
1-\mathbf{1}^{\mathsf T}v(t)
=
\mathbb{P}(T\leq t).
$$

Because the killed chain eventually exits, every eigenvalue of \(R\) has negative real part and

$$
\lim_{t\rightarrow\infty}e^{Rt}=0.
$$

Consequently,

$$
R\int_0^\infty e^{Rt}\,dt
=
\left[e^{Rt}\right]_{t=0}^{t=\infty}
=
-I,
$$

so

$$
\int_0^\infty e^{Rt}\,dt
=
-R^{-1}
=
Z.
$$

This integral identity is the source of the occupation, exit-probability, and first-passage-time formulas.

### 1.1.2 Occupation times

For each internal state \(i\), define its random occupation time before exit:

$$
U_i
=
\int_0^T \mathbf{1}_{\{X_t=i\}}\,dt.
$$

Taking expectations and interchanging expectation and integration gives

$$
\begin{aligned}
\mathbb{E}[U_i]
&=
\int_0^\infty
\mathbb{P}(X_t=i,\ T>t)\,dt\\
&=
\int_0^\infty
\left(e^{Rt}p\right)_i\,dt\\
&=
(Zp)_i.
\end{aligned}
$$

Therefore,

$$
u=Zp
$$

is exactly the vector of expected state-occupation times before leaving the subnetwork.

In particular, the matrix entry

$$
Z_{ij}
$$

is the expected time spent in state \(i\) before exit when the chain starts in state \(j\).

This is the Green-function interpretation of \(Z\): its columns describe where time is spent from each possible entrance.

### 1.1.3 Mean exit time

Before exit, the trajectory must be in exactly one internal state. Pathwise,

$$
T=\sum_i U_i.
$$

Taking expectations yields

$$
\mathbb{E}[T]
=
\sum_i \mathbb{E}[U_i]
=
\mathbf{1}^{\mathsf T}Zp.
$$

Thus the mean exit time is the total expected occupation accumulated over all internal states.

### 1.1.4 Exit probabilities as integrated hazards

Let \(q_j\) be the total rate of transitions from internal state \(j\) to states outside the subnetwork. Conditional on being in state \(j\) at time \(t\), the probability of exiting during the interval \([t,t+dt]\) is

$$
q_j\,dt+o(dt).
$$

The joint density of exiting at time \(t\) with \(j\) as the pre-exit state is therefore

$$
f_j(t)
=
q_j\left(e^{Rt}p\right)_j.
$$

Integrating over all possible exit times gives

$$
\begin{aligned}
\mathbb{P}(B=j)
&=
\int_0^\infty f_j(t)\,dt\\
&=
q_j\int_0^\infty
\left(e^{Rt}p\right)_j\,dt\\
&=
q_j(Zp)_j.
\end{aligned}
$$

The intuition is simple: the probability of exiting from state \(j\) is its exit hazard multiplied by the expected amount of time exposed to that hazard.

More generally, let \(C_{cj}\) be the rate of external exit channel \(c\) from state \(j\). Then the matrix

$$
C Z
$$

maps entrance distributions to eventual exit-channel probabilities:

$$
\mathbb{P}(\text{exit through channel }c\mid p)
=
(CZp)_c.
$$

### 1.1.5 Occupation conditioned on the pre-exit state

Suppose first that the chain starts from one entrance state \(a\). The expected occupation of state \(k\), jointly with eventually exiting from state \(b\), can be factored using the Markov property.

The expected time spent in \(k\) before exit is \(Z_{ka}\). Starting from \(k\), the probability of eventually exiting from \(b\) is

$$
q_b Z_{bk}.
$$

Therefore,

$$
\mathbb{E}
\left[
U_k\mathbf{1}_{\{B=b\}}
\mid X_0=a
\right]
=
Z_{ka}\,q_b Z_{bk}.
$$

Since

$$
\mathbb{P}(B=b\mid X_0=a)
=
q_bZ_{ba},
$$

division gives the conditional occupation time

$$
\mathbb{E}[U_k\mid B=b,X_0=a]
=
\frac{Z_{ka}Z_{bk}}{Z_{ba}}.
$$

The exit rate \(q_b\) cancels: conditioning on the pre-exit state removes the absolute scale of its outward hazard.

For a general entrance distribution \(p\), let \(u=Zp\). The same argument gives

$$
\mathbb{E}[U_k\mid B=b,p]
=
\frac{u_k Z_{bk}}{u_b}.
$$

Summing over \(k\) produces the conditional mean exit time:

$$
\begin{aligned}
\mathbb{E}[T\mid B=b,p]
&=
\sum_k \frac{u_kZ_{bk}}{u_b}\\
&=
\frac{(Zu)_b}{u_b}\\
&=
\frac{(Z^2p)_b}{(Zp)_b}.
\end{aligned}
$$

This explains why the implementation first computes

$$
u=Zp
$$

and then

$$
u^{(2)}=Zu=Z^2p.
$$

The ratio \(u^{(2)}_b/u_b\) is not an arbitrary approximation formula: it is the exact conditional mean first-exit time associated with pre-exit state \(b\). The approximation enters only when that mean is used in place of a random draw from the conditional distribution.

### 1.1.6 Higher exit-time moments

Repeated integration of the killed semigroup gives

$$
\int_0^\infty t^m e^{Rt}\,dt
=
m!\,Z^{m+1}.
$$

For example,

$$
\int_0^\infty t e^{Rt}\,dt=Z^2
$$

and

$$
\int_0^\infty t^2 e^{Rt}\,dt=2Z^3.
$$

Hence the joint \(m\)-th moment of the exit time and the event \(B=b\) is

$$
\mathbb{E}
\left[T^m\mathbf{1}_{\{B=b\}}\right]
=
m!\,q_b\left(Z^{m+1}p\right)_b.
$$

After conditioning,

$$
\mathbb{E}[T^m\mid B=b]
=
m!\,
\frac{\left(Z^{m+1}p\right)_b}{(Zp)_b}.
$$

The conditional variance is therefore

$$
\operatorname{Var}(T\mid B=b)
=
2\frac{(Z^3p)_b}{(Zp)_b}
-
\left(
\frac{(Z^2p)_b}{(Zp)_b}
\right)^2.
$$

The current algorithm retains the first conditional moment but discards this conditional variance and all higher-order information.

### 1.1.7 Derivation of the C++ cut-time loss

The C++ shedding score also has a direct probabilistic interpretation. Consider removing an unprotected state \(j\). Reorder the states so that the retained states are collected in \(K\), and partition the fundamental matrix as

$$
Z=
\begin{bmatrix}
Z_{KK} & Z_{Kj}\\
Z_{jK} & Z_{jj}
\end{bmatrix}.
$$

After state \(j\) is removed, the fundamental matrix of the retained principal generator is

$$
\widetilde{Z}
=
Z_{KK}
-
\frac{Z_{Kj}Z_{jK}}{Z_{jj}}.
$$

This is the rank-one inverse downdate used by the implementation.

Assume that the entrance distribution has no initial mass in the removable state, so \(p_j=0\). Its expected occupation of state \(j\) is

$$
u_j=Z_{jK}p_K.
$$

The difference between the original and reduced mean exit times is

$$
\begin{aligned}
\Delta T_j
&=
\mathbf{1}^{\mathsf T}Zp
-
\mathbf{1}_K^{\mathsf T}\widetilde{Z}p_K\\
&=
\frac{u_j}{Z_{jj}}
\mathbf{1}^{\mathsf T}Z e_j.
\end{aligned}
$$

Equivalently,

$$
\Delta T_j
=
\frac{u_j}{Z_{jj}}
\sum_i Z_{ij}.
$$

This is exactly the Julia column-oriented form of the C++ expression

$$
\texttt{Z.row(i).t() / Z.diag() \% sum(Z,1)}.
$$

The two factors have useful interpretations. First,

$$
h_j=\frac{u_j}{Z_{jj}}
$$

is the probability of hitting state \(j\) before exiting. Indeed, expected occupation in \(j\) equals the probability of reaching \(j\), multiplied by the expected occupation in \(j\) after starting there.

Second,

$$
\tau_j
=
\mathbf{1}^{\mathsf T}Ze_j
=
\sum_i Z_{ij}
$$

is the mean remaining exit time when starting from state \(j\).

Therefore,

$$
\Delta T_j=h_j\tau_j.
$$

In words: the time lost by turning state \(j\) into an immediate exit is the probability that the trajectory would have reached \(j\), multiplied by the mean amount of subnetwork lifetime that would have remained after reaching it.

### 1.1.8 Row versus column conventions

The C++ implementation uses row-oriented probability evolution, whereas the Julia implementation uses column-oriented evolution. Their fundamental matrices are transposes of one another:

$$
Z_{\mathrm{Julia}}=Z_{\mathrm{C++}}^{\mathsf T}.
$$

Thus a C++ row describing occupation from entrance \(a\) becomes Julia column \(a\). The formulas are mathematically identical once this transpose is taken into account.

## 2. Main deficiency: FESS is presently a first-moment closure

The first-exit time of a finite continuous-time Markov chain has a phase-type distribution. Replacing it with its conditional mean removes timing fluctuations even when the exit-state probabilities are exact.

This can distort quantities that depend on more than the first moment:

- first-passage-time distributions;
- temporal variance;
- autocorrelation functions;
- stochastic oscillation periods;
- phase diffusion;
- synchronization and dephasing across an ensemble;
- transient probability distributions at fixed physical times.

The effect generally grows when a FESS step crosses a larger subnetwork, because more internal timing variability is replaced by one deterministic conditional mean.

## 3. Exact phase-type first-exit sampling

The most important improvement is to sample the full joint law of the first exit time and exit channel.

### 3.1 Sample the time first

The probability of remaining inside the subnetwork through time \(t\) is

$$
S(t)
= \mathbb{P}(T>t)
= \mathbf{1}^{\mathsf T} e^{Rt}p.
$$

Draw \(U\sim\operatorname{Uniform}(0,1)\) and solve

$$
S(T)=U.
$$

At the sampled time \(T\), the unnormalized probability of exiting from state \(j\) is

$$
q_j \left(e^{RT}p\right)_j.
$$

Therefore,

$$
\mathbb{P}(B=j\mid T)
=
\frac{q_j \left(e^{RT}p\right)_j}
{\sum_k q_k \left(e^{RT}p\right)_k}.
$$

After sampling \(B\), the physical reaction leaving the subnetwork is sampled in proportion to its propensity at that state.

### 3.2 Sample the exit state first

The current state-first ordering can also be retained. First sample

$$
\mathbb{P}(B=j)=q_j(Zp)_j.
$$

The conditional survival function is

$$
\mathbb{P}(T>t\mid B=j)
=
\frac{\left(e^{Rt}Zp\right)_j}{(Zp)_j}.
$$

After drawing \(U\sim\operatorname{Uniform}(0,1)\), solve

$$
\frac{\left(e^{RT}Zp\right)_j}{(Zp)_j}=U.
$$

This is the exact random analogue of the conditional-mean calculation used by the current implementation.

### 3.3 Uniformization

Uniformization is preferable to a general eigendecomposition. Choose

$$
\nu \geq \max_j(-R_{jj})
$$

and define

$$
P = I + \frac{R}{\nu}.
$$

Then \(P\) is substochastic and

$$
e^{Rt}
=
e^{-\nu t}
\sum_{k=0}^{\infty}
\frac{(\nu t)^k}{k!}P^k.
$$

This representation preserves nonnegativity and avoids the instability of eigenvector decompositions for nonnormal generators.

### 3.4 Exactness consequence

Suppose the selected subnetwork contains the current state and uses the original reaction rates. If FESS samples the exact first-exit time and exit destination, then restarting from that destination reproduces the original continuous-time Markov chain by the strong Markov property.

The subnetwork may be redrawn, persistent, or selected adaptively from the past. Its size changes the number of internal reactions skipped by a FESS step, but it does not change the path law, apart from numerical error.

Under exact first-exit sampling, shedding is therefore primarily a question of how much acceleration the retained network provides, rather than a direct simulation-accuracy approximation.

## 4. Preserve the entrance-to-boundary transfer operator

If states are to be compressed rather than simply reclassified as exits, the correct object to preserve is the matrix-valued boundary transfer function.

Let \(B\) embed all entrance distributions into the internal state space, and let \(C\) map internal probability to every external exit channel. Define

$$
H(s) = C(sI-R)^{-1}B.
$$

This function contains the joint exit-channel and exit-time law.

At zero frequency,

$$
H(0)=CZB
$$

is the matrix of exit probabilities. Its first two time moments are

$$
M_1 = CZ^2B
$$

and

$$
M_2 = 2CZ^3B.
$$

For entrance \(a\) and exit channel \(b\), the conditional mean is

$$
\mathbb{E}[T\mid b,a]
=
\frac{(M_1)_{ba}}{H(0)_{ba}}.
$$

The current shedding criterion controls only one scalar output,

$$
\mathbf{1}^{\mathsf T}Zp,
$$

and does not directly control exit probabilities, conditional moments, or minority entrance states.

A more fundamental reduction should preserve \(H(s)\), or at least \(H(0)\), \(M_1\), and \(M_2\), for every active entrance.

Possible matrix-reduction approaches include:

- block rational Krylov moment matching;
- tangential interpolation of \(H(s)\);
- positive-system balanced truncation;
- structure-preserving phase-type reduction;
- matching shifted resolvents \((s_\ell I-R)^{-1}\) over several time scales.

Any reduced realization should remain a valid killed generator: off-diagonal entries must be nonnegative and column sums must be nonpositive.

## 5. Eliminate historical states without converting them into immediate exits

Directly deleting a state makes every transition into that state an immediate subnetwork exit. This gives a valid smaller first-exit problem, but shortens the FESS macrostep and discards the possibility of returning from that state without leaving the compressed model.

Partition the subnetwork into a retained core \(C\) and historical states \(H\):

$$
R=
\begin{bmatrix}
R_{CC} & R_{CH}\\
R_{HC} & R_{HH}
\end{bmatrix}.
$$

Exact elimination in the resolvent gives

$$
R_{\mathrm{eff}}(s)
=
R_{CC}
+
R_{CH}(sI-R_{HH})^{-1}R_{HC}.
$$

At \(s=0\), this becomes the static stochastic complement

$$
R_{\mathrm{eff}}(0)
=
R_{CC}
-
R_{CH}R_{HH}^{-1}R_{HC}.
$$

The exact reduced operator depends on \(s\), so eliminating states generally introduces memory. A static Schur complement preserves zero-frequency behavior, but not the complete exit-time distribution.

A more faithful reduced Markov model could approximate this memory with a small number of auxiliary phase states, chosen to match moments or selected values of \(H(s)\).

## 6. Multiple entrances require a block error criterion

For multiple trajectories using one shared matrix, the current shedding distribution is the population-weighted mixture

$$
p_{\mathrm{mix}}
=
\sum_a w_a p_a.
$$

This optimizes average behavior. An entrance with a small population weight can experience a large change while contributing very little to the shedding score.

Instead, retain the complete entrance matrix

$$
B = [p_1\;p_2\;\cdots\;p_m]
$$

and control an error such as

$$
\max_a
\frac{
\left\|\widehat{H}(0)e_a-H(0)e_a\right\|_1
}{
\left\|H(0)e_a\right\|_1
},
$$

along with corresponding first- and second-moment errors.

Alternatives include:

- a worst-entrance criterion;
- a weighted criterion with a positive minimum weight;
- an induced matrix norm over all entrance columns;
- clustering entrances into several compatible common generators.

If one capacity-limited network cannot represent every entrance accurately, splitting the trajectories into several compatible groups is preferable to forcing a high-loss cut.

## 7. Replace per-cut thresholds with a global error certificate

The C++-matched rule accepts a state whenever its individual relative loss is below the threshold. Many individually small cuts can accumulate into a large total change. A hard capacity can also force a cut with arbitrarily large loss.

A global reduction budget should monitor the final changes in

$$
H(0)=CZB,
$$

$$
M_1=CZ^2B,
$$

and

$$
M_2=2CZ^3B.
$$

For example, shedding could require

$$
\left\|\widehat{H}(0)-H(0)\right\| \leq \varepsilon_0,
$$

$$
\left\|\widehat{M}_1-M_1\right\| \leq \varepsilon_1,
$$

and

$$
\left\|\widehat{M}_2-M_2\right\| \leq \varepsilon_2.
$$

Capacity should be subordinate to this certificate. If the requested capacity cannot satisfy the tolerances, the algorithm should enlarge the network, split the entrance batch, or explicitly report that the requested reduction is not admissible.

## 8. Select states by expected skipped reactions

Graph depth is not a dynamical measure of the value of a subnetwork. Mean residence time is also not always the correct acceleration objective: a slow state may contribute substantial physical time while containing few reaction events.

Let \(\lambda_{\mathrm{int},j}\) be the total rate of transitions from state \(j\) to other states retained inside the subnetwork. The expected number of internal reactions before exit is

$$
N_{\mathrm{int}}(p)
=
\lambda_{\mathrm{int}}^{\mathsf T}Zp.
$$

This quantity more directly measures how many SSA reactions are skipped by a FESS macrostep.

Expansion should prioritize boundary destinations that produce the greatest increase in \(N_{\mathrm{int}}\), while shedding should remove states producing the smallest decrease, subject to the exit-law error certificate.

For multiple entrances, possible objectives include

$$
\sum_a w_a N_{\mathrm{int}}(p_a)
$$

or the fairness-oriented objective

$$
\min_a N_{\mathrm{int}}(p_a).
$$

## 9. Respect generator and M-matrix structure

For a valid killed generator, \(-R\) is a nonsingular M-matrix and

$$
Z=(-R)^{-1}\geq 0
$$

entrywise. This structure should be used as a numerical invariant.

Important checks include:

- nonnegative off-diagonal rates in \(R\);
- nonpositive column sums in \(R\);
- nonnegative entries in \(Zp\);
- exit probabilities summing to one;
- residual-based error bounds for linear solves;
- condition monitoring for nearly metastable subnetworks.

Silently clamping negative occupations or exit probabilities can conceal a loss of numerical reliability. Uniformization and structure-preserving model reduction are attractive partly because they maintain positivity by construction.

Nonnormality is another reason not to base the algorithm only on eigenvalues. Shifted resolvents and uniformization describe the transient behavior more reliably than a potentially ill-conditioned eigenvector decomposition.

## 10. Recommended research order

### Stage 1: Exact first-exit times

Implement phase-type time sampling using uniformization. Validate against SSA on small systems using:

- exit-destination probabilities;
- conditional exit-time distributions;
- first three exit-time moments;
- transient distributions at fixed times;
- autocorrelation and first-passage distributions.

### Fixed-work horizon experiment

A useful diagnostic is to run each of the four FESS variants for a fixed number \(n\) of steps and measure the physical terminal time reached. The variants are persistent and redrawn FESS, each run with either one trajectory or multiple trajectories. For method \(M\), define

$$
t_f^{(M)}(n)=\sum_{k=1}^{n}\Delta T_k^{(M)}.
$$

For one trajectory, one step is one FESS step. For multiple trajectories, one step advances every active trajectory once using the same network and matrix.

The experiment reports the completed step count, number of trajectories, mean terminal time, minimum terminal time, and maximum terminal time. The multiple-trajectory result therefore exposes both the average temporal reach and the spread caused by independently sampled exits and holding times.

This diagnostic answers a specific question: how much physical time does one FESS step traverse? It is not by itself a wall-clock efficiency comparison, because a single-trajectory step and a multiple-trajectory step have different costs and advance different numbers of trajectories.

The experiment is implemented by `compare_fess_steps(initial_state, n, model, rates; ...)`. In the Lotka--Volterra notebook, `n_steps` controls \(n\), and the four returned rows show the average, earliest, and latest final times reached by the four FESS variants.

For statistical conclusions, repeat the fixed-work experiment over many seeds and compare the distribution of \(t_f^{(M)}(n)\), not only its mean. After exact phase-type time sampling is introduced, this distribution becomes an especially useful check that FESS retains the correct timing variability while crossing more physical time per macrostep.

### Stage 2: Matrix-based network value

Replace graph-depth expansion and mean-time-only shedding with an objective based on the expected number of skipped internal reactions. For a common generator, use all entrance columns rather than only their population mixture.

### Stage 3: Global error control

Introduce global, a posteriori tolerances for exit probabilities and time moments. Do not permit hard capacity cuts to silently violate these tolerances.

### Stage 4: Historical-state reduction

Investigate stochastic complements, block rational Krylov methods, and positive phase-type realizations that preserve the transfer function \(H(s)\) while representing historical excursions with fewer degrees of freedom.

## 11. Main conclusion

The next fundamental advance should not be a cheaper way to compute the current scalar shedding score. It should be a change in what FESS preserves.

The highest-priority step is exact phase-type first-exit sampling. Once that is available, an arbitrary adaptive subnetwork gives an exact coarse-grained CTMC step. Matrix analysis can then focus on the genuinely important reduction problem: selecting or compressing states so that each generator construction skips as many internal reactions as possible while preserving the full entrance-to-boundary first-exit law.
