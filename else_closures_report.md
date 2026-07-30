# The recurrent-subnetwork closure for the CME (for §2.2)

A characterization of the mass-conserving subnetwork operator that fills the recurrent-chain results stubbed in §2.2. Stated in column convention ($A_{ij}$ = rate $j\to i$, columns sum to zero); transpose if the draft is row-stochastic.

## Setup

Truncate the CME generator to a finite interior $J$ with exterior $J'$: $$A = \begin{pmatrix} A_{JJ} & A_{JJ'}\\ A_{J'J} & A_{J'J'}\end{pmatrix}, \qquad \Phi_{\text{out}} = \mathbf 1^\top A_{J'J}\,p_J, \quad \Phi_{\text{in}}  = \mathbf 1^\top A_{JJ'}\,p_{J'}.$$ Two ways to close the boundary — the two subnetwork types in the draft:

- **transient / wired** (classic FSP, cemetery state): evolve $A_{JJ}$.
- **free / recurrent** (reflecting): evolve the compressed generator
  $\tilde A_{JJ} = A_{JJ} + \operatorname{diag}\!\big(\mathbf 1^\top A_{J'J}\big).$

## Transient closure: leakage

$\mathbf 1^\top A_{JJ} = -\mathbf 1^\top A_{J'J} \le 0$, so
$$\tfrac{d}{dt}\,\mathbf 1^\top p_J = \mathbf 1^\top A_{JJ}\,p_J = -\Phi_{\text{out}},$$
and the mass lost over a step is exactly $\int \Phi_{\text{out}}\,ds$. This is the Munsky–Khammash leakage the draft cites, and matches ELSE's wired subnetwork with an absorbing cemetery.

## Compressed generator

$\tilde A_{JJ}$ is a bona fide CTMC generator on $J$: off-diagonals $\ge 0$ and $$\mathbf 1^\top \tilde A_{JJ} = \mathbf 1^\top A_{JJ} + \mathbf 1^\top A_{J'J} = 0.$$ Hence $e^{t\tilde A_{JJ}}$ is a **mass-conserving positive (stochastic) semigroup**: $e^{t\tilde A_{JJ}} \ge 0$ elementwise and $\lVert e^{t\tilde A_{JJ}}\rVert_{1\to 1} = 1$ for all $t\ge 0$. This is exactly the *free / recurrent subnetwork operator* §2.2 leaves open — the escaping rate $\mathbf 1^\top A_{J'J}$ is folded back onto the diagonal, so no probability leaves $J$ and none is created. (A CTMC generator has no self-loop, so "reflect to origin" and "erase the edge to $J'$" give the identical operator.)

**Model error (flux bound).** Against the true restricted dynamics
$\dot p_J = A_{JJ}p_J + A_{JJ'}p_{J'}$, using $\lVert e^{sA_{JJ}}\rVert_{1\to1}\le 1$, $$\big\lVert p_J(t_{n+1}) - \tilde p_J(t_{n+1})\big\rVert_1 \;\le\;
\int_{t_n}^{t_{n+1}} \big(\Phi_{\text{in}}(s) + \Phi_{\text{out}}(s)\big)\,ds .$$
So the recurrent closure buys stability ($\lVert\cdot\rVert_{1\to1}=1$, good for long-time integration and flux-adaptive stepping) at a boundary-flux-bounded error, with conserved mass.


$$Z = -A_{JJ}^{-1} = \int_0^\infty e^{tA_{JJ}},dt.$$

Mark's fundamental matrix is the time-integral of the FSP absorbing-closure propagator. FSP propagates $e^{tA_{JJ}}$ in time; his method integrates it via the resolvent. Same operator $A_{JJ}$ (your transient/leak closure), used at two levels. Everything in the method is a time-integral of an FSP transient:

- occupation time $= \int_0^\infty [e^{tA_{JJ}}p_0]_i,dt = (Zp_0)_i$
- MFPT $= \int_0^\infty S(t),dt = \mathbf 1^\top Z p_0$, where $S(t)=\mathbf 1^\top e^{tA_{JJ}}p_0$ is the FSP survival
- exit distribution $= A_{J'J},Z,p_0 = \int_0^\infty A_{J'J}e^{tA_{JJ}}p_0,dt$

The sharp part: your error term is his observable

Here's the reframing that's genuinely FSP-native and, I think, the paper's hook. The first-passage-time density is exactly the boundary outflux:
$$f(t) = -S'(t) = \mathbf 1^\top A_{J'J},e^{tA_{JJ}}p_0 = \Phi_{\text{out}}(t).$$

So your FSP leaked-mass term $\int_0^t \Phi_{\text{out}},ds$ — the truncation error you grow $J$ to suppress — is, read the other way, the escape CDF: the entire quantity his method computes. Same object, opposite intent. FSP treats the boundary flux as error to be minimized; his method resolves it exactly as the answer. That role-reversal is the reframing, and it's stated entirely in your flux/leak language.

Three things it buys (so it's not just relabeling)

1. The closure taxonomy explains why the method needs the absorbing closure. Your reflecting/compressed $\tilde A_{JJ}$ is mass-conserving ⇒ has a zero eigenvalue ⇒ $\tilde A_{JJ}^{-1}$ doesn't exist (no fundamental matrix; instead a stationary law). The leak closure $A_{JJ}$ is invertible precisely because it leaks. Invertibility ⟺ escape — your §2.2 taxonomy becomes the existence criterion for his $Z$.
2. The resolvent generalizes the method for free. $(sI - A_{JJ})^{-1}p_0$ is the Laplace transform of the FSP propagator. His method is the $s=0$ slice (means, integrated quantities); general $s$ gives $\hat f(s)=\mathbf 1^\top A_{J'J}(sI-A_{JJ})^{-1}p_0$ — the full first-passage-time distribution and all its moments (derivatives at 0), not just the mean escape time. That's a concrete new result his current paper doesn't have.
3. FSP adaptive truncation becomes his subnetwork-selection rule. You grow/prune $J$ to control $\int\Phi_{\text{out}}$; that same flux criterion is a principled way to choose which subnetwork to collapse. Your adaptive-FSP machinery ports directly.

Thesis, one line

▎ His method is the resolvent (Laplace-domain, $s{=}0$) form of the absorbing-closure FSP: FSP resolves when, his method integrates when out. The boundary flux that is FSP's error is his first-passage density.

Honest caveat

The core identity ($Z=\int e^{tA}dt$, resolvent-at-0) is classical Kemeny–Snell — so the contribution isn't the identity, it's the unification plus what it produces. To make it more than elegant relabeling, it should deliver at least one of: the full FPT-distribution extension (point 2, the strongest — genuinely new and FSP-derived), a cleaner error/convergence theorem for subnetwork selection (point 3), or the invertibility dividing line as a clean lemma (point 1). I'd anchor the reframing on point 2 — it's the part that gives Mark something he can't get from the mean-escape formulation, and it's pure FSP (transient $\Phi_{\text{out}}(t)$ = FPT density).
