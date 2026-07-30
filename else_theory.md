# Application of ELSE to FSP 

An FSP reframing of the subnetwork method. ELSE is a trajectory method but the idea to resolve persistent manifolds analytically is general enough to be applied to FSP. The hypothesis is that the boundary flux the FP-FSP bounds as *error* (Thm 4.5) is exactly the method's escape density.


---

## ELSE 

**Subnetworks.** For a subnetwork $G$ of the full state space $G_0$, edges from $G$ to its complement are handled in one of two ways (Fig. 1 of the draft): the **recurrent (free)** subnetwork erases them (mass-conserving), and the **transient (wired)** subnetwork redirects them to a single absorbing cemetery state $\omega_\infty$ at equal total rate. Finite state projection uses the transient subnetwork; the method here also builds on the transient subnetwork and resolves its escape exactly.

**Rate matrix (row convention).** For rates $r(\omega_i\to\omega_j)$,
$$\mathrm R_{i,j} = \begin{cases} r(\omega_i\to\omega_j), & i\neq j\\ -\sum_{k\neq i} r(\omega_i\to\omega_k), & i=j,\end{cases} \qquad \frac{d\vec p}{dt} = \mathrm R^\top \vec p, \quad \vec p(t) = (e^{\mathrm R t})^\top \vec p(0).$$
For the transient subnetwork the row sums are $\mathrm R\mathbf 1 = -\boldsymbol\Delta \le 0$, where $\Delta_i$ is the total escape (killing) rate from $\omega_i$; $\mathrm R$ is negative definite.

**Fundamental matrix (Eq. 8).** The core object is
$$\mathrm Z_{i,j} := \int_0^\infty P_i(X_t=j)\,dt = \int_0^\infty \delta_i^\top e^{\mathrm R t}\delta_j\,dt = -(\mathrm R^{-1})_{i,j}, \qquad \mathrm Z = -\mathrm R^{-1}\ \ (\ge 0),$$
the expected time spent in $j$ before escape when started at $i$. From $\mathrm Z$ the method reads off, exactly:
- **escape-state probabilities (Eq. 10):** $\ p_{\text{escape}}(\omega_a) = -(\mathrm Z^\top\delta_a)\circ(\mathrm R\mathbf 1)$, normalized since $-\mathrm Z\,\mathrm R\mathbf 1=\mathbf 1$;
- **mean escape time from $a$ (Eq. 15):** $\ \mathbb E\,t_{a\to} = \delta_a^\top \mathrm Z\mathbf 1$, and conditional occupation times (Eqs. 11–12);
- **occupation-time (co)variances (Eqs. 13–14)** from $\mathrm Z^2$;
- **cut times (Eqs. 16–18)** for pruning: via a killing-rate limit $\mathrm R'(\gamma)=\mathrm R-\gamma\delta_b\delta_b^\top$ and Sherman–Morrison, the drop in subnetwork occupation on removing state $b$, used to shed the least significant state.

**Algorithm (exact linear system embedding).** The Gillespie jump-and-hold iteration is replaced by: from the current subnetwork $G$, use $\mathrm Z$ to draw the exact escape state (Eq. 10) and the exact time to escape (Eqs. 11/15); jump there; augment $G$ with the new state; and prune via cut time to hold $|G|\le n_{\max}$. Escape from an entire subnetwork is collapsed into one linear solve, exactly, at cost set by $|G|$ rather than by the number of internal transitions.

---

## Flux-Preserving FSP

**Convention bridge.** In `ex_article.tex` the generator is column convention, $\mathbf A_{ij}=$ rate $j\to i$, $\dot{\boldsymbol p}=\mathbf A\boldsymbol p$. On the active set $J$ (complement $J'$, block form Eq. 4.1), the transient/wired subnetwork is the truncated generator $\mathbf A_{JJ}$ and the recurrent/free subnetwork is the compressed generator $\widetilde{\mathbf A}_{JJ}=\mathbf A_{JJ}+\operatorname{diag}(\mathbf 1^\top\mathbf A_{J'J})$. Fornace's rate matrix is the transpose, $\mathrm R = \mathbf A_{JJ}^\top$, so his fundamental matrix is $\mathrm Z = (-\mathbf A_{JJ}^{-1})^\top$. We work with the column-convention fundamental matrix
$$Z := -\mathbf A_{JJ}^{-1} = \int_0^\infty e^{t\mathbf A_{JJ}}\,dt \ \ge 0. \tag{$\star$}$$
Equation $(\star)$ is the reframing in one line: **the method's fundamental matrix is the time-integral of the FSP truncated-generator semigroup $e^{t\mathbf A_{JJ}}$** — the same $\mathbf A_{JJ}$ whose contraction property $\lVert e^{t\mathbf A_{JJ}}\rVert_{1\to1}\le1$ carries Cor. 4.4 and Thm. 4.5.

### FSP Error (or leakage)

The retained mass under the truncated closure is $\Gamma_J(t)=\mathbf 1^\top e^{t\mathbf A_{JJ}}\boldsymbol p_0$ (survival in $J$). Differentiating and using the row-sum identity $\mathbf 1^\top\mathbf A_{JJ}=-\mathbf 1^\top\mathbf A_{J'J}$ (Eq. 4.2),
$$f(t) := -\Gamma_J'(t) = \mathbf 1^\top\mathbf A_{J'J}\,e^{t\mathbf A_{JJ}}\boldsymbol p_0 = \Phi_{\text{out}}(t), \qquad F(t) = \int_0^t\Phi_{\text{out}}(s)\,ds = 1-\Gamma_J(t). \tag{$\dagger$}$$
So the **boundary outflux $\Phi_{\text{out}}$ is the first-passage (escape-time) density**, and the leaked mass $1-\Gamma_J$ that FSP drives to $\le\varepsilon$ is the **escape CDF**. The quantity the flux-adaptive method treats as truncation error to suppress (Thm. 4.5, $\lVert\boldsymbol\varepsilon(t_{n+1})\rVert_1\le\int(\Phi_{\text{in}}+\Phi_{\text{out}})\,ds$) is, read from the transient subnetwork, exactly what ELSE resolves as the answer.

### Resolvent

Remark 4.6 already invokes the resolvent $(\lambda I-\mathbf A)^{-1}$ through Hille–Yosida. It is the Laplace transform of the FSP semigroup,
$$(\lambda I - \mathbf A_{JJ})^{-1} = \int_0^\infty e^{-\lambda t}\,e^{t\mathbf A_{JJ}}\,dt, \tag{$\ddagger$}$$
where $\lambda>0$ is the Laplace variable and the weight $e^{-\lambda t}$ down-weights late times. FSP lives in the time domain ($e^{t\mathbf A_{JJ}}\boldsymbol p_0$, the full transient); the method lives at $\lambda\to0$, where the weight becomes $1$ and the transform reduces to the plain time-integral $(\star)$. Whether that $\lambda=0$ limit is finite depends on the eigenvalues of the closure, and the two FSP closures of `ex_article.tex` differ exactly here:

- **transient / truncated $\mathbf A_{JJ}$** (Cor. 4.4, strict contraction): every state escapes, all eigenvalues have $\operatorname{Re}<0$, so $(\ddagger)$ **stays finite at $\lambda=0$** and $(\star)$ gives the fundamental matrix $Z$;
- **recurrent / compressed $\widetilde{\mathbf A}_{JJ}$** (Prop. 4.3, stochastic semigroup): $\mathbf 1^\top\widetilde{\mathbf A}_{JJ}=0$, so $0$ is an eigenvalue and $(\ddagger)$ **blows up as $\lambda\to0$**, diverging as $\Pi/\lambda$ along the stationary mode $\boldsymbol\pi$.

### ELSE on Truncated Generator 

At $\lambda=0$ the truncated resolvent stays finite, giving $Z=-\mathbf A_{JJ}^{-1}$ and, via $(\dagger)$, the whole escape law rather than its mean. With $\tau$ the escape time (phase-type$(\boldsymbol p_0,\mathbf A_{JJ})$):
$$\mathbb E[\tau^n] = n!\;\mathbf 1^\top Z^{\,n}\boldsymbol p_0, \qquad \boldsymbol\pi_{\text{exit}} = \mathbf A_{J'J}\,Z\,\boldsymbol p_0. \tag{1}$$
$n=1$ is the method's mean escape time $\mathbf 1^\top Z\boldsymbol p_0$ (Fornace Eq. 15, transposed); $n=2$ gives $\operatorname{Var}(\tau)=2\,\mathbf 1^\top Z^2\boldsymbol p_0-(\mathbf 1^\top Z\boldsymbol p_0)^2$ (his Eq. 14). The powers $Z^n$ carry every moment and $\boldsymbol\pi_{\text{exit}}$ the exit-state law (his Eq. 10); the escape density itself is $f(t)=\Phi_{\text{out}}(t)$ of $(\dagger)$. The method reports the $n=1$ slice; the resolvent supplies the full distribution.

### ELSE on Compressed Generator

**Motivation.** FP-FSP is a simulation scheme: it truncates the infinite-dimensional generator $\mathbf A$ into a sequence of finite matrices $\mathbf A_{J_nJ_n}$, evolves each over a window $[t_n,t_{n+1}]$, and stitches the evolutions together. Per window the transient closure $\mathbf A_{J_nJ_n}$ answers one question: how mass leaves $J_n$. Its mass drains, so the only surviving content is the escape law above. The compressed closure $\widetilde{\mathbf A}_{J_nJ_n}$ instead reflects escaping mass back, giving a mass-conserving chain on the window with a stationary distribution $\boldsymbol\pi^{(n)}$ (unique when $J_n$ is connected, which the flux pruning preserves) and the quantities leakage destroys: within-window mean first-passage times, and how fast a time-average $\tfrac1T\int_0^T g(X_t)\,dt$ settles. This is the recurrent (free) versus transient (wired) subnetwork distinction of nuts.pdf, where the transient/FSP closure "destroys many averaged quantities" while the free closure preserves them.

The $\boldsymbol\pi^{(n)}$ is local in time: the equilibrium of the reflected dynamics on the current window, not a truncation of any global CME stationary law (which may sit at infinity or fail to exist). For a spreading distribution it is partly an artifact of the reflecting wall, but on a window that captures a metastable basin it is the within-basin quasi-equilibrium and $D$ its relaxation structure. Worth exploring.

**Resolvent.** Write $\widetilde{\mathbf A}:=\widetilde{\mathbf A}_{JJ}$, with $\widetilde{\mathbf A}\boldsymbol\pi=\mathbf 0$ and $\mathbf 1^\top\boldsymbol\pi=1$. Let $\Pi:=\boldsymbol\pi\mathbf 1^\top$; for any distribution $\boldsymbol v$, $\Pi\boldsymbol v=(\mathbf 1^\top\boldsymbol v)\boldsymbol\pi=\boldsymbol\pi$, so $\Pi$ collapses everything onto $\boldsymbol\pi$. It is the rank-one matrix whose every column is $\boldsymbol\pi$ (outer product of the column $\boldsymbol\pi$ with the row $\mathbf 1^\top$), equal to the infinite-time limit $\lim_{t\to\infty}e^{t\widetilde{\mathbf A}}$ and to the eigenprojector for the zero eigenvalue of $\widetilde{\mathbf A}$. Two facts drive the rest:
$$\widetilde{\mathbf A}\,\Pi=(\widetilde{\mathbf A}\boldsymbol\pi)\mathbf 1^\top=0, \qquad \Pi\,\widetilde{\mathbf A}=\boldsymbol\pi(\mathbf 1^\top\widetilde{\mathbf A})=0,$$
i.e. $\Pi$ singles out the zero-eigenvalue direction and $I-\Pi$ is everything else, where $\widetilde{\mathbf A}$ has all eigenvalues with $\operatorname{Re}<0$. Splitting the semigroup along this decomposition,
$$e^{t\widetilde{\mathbf A}} = \Pi + \underbrace{e^{t\widetilde{\mathbf A}}(I-\Pi)}_{\to\,0\ \text{exponentially}},$$
which is just $e^{t\widetilde{\mathbf A}}\boldsymbol p_0\to\boldsymbol\pi$ (the equilibration of Prop. 4.3) restated. Feed this into the resolvent integral:
$$(\lambda I-\widetilde{\mathbf A})^{-1} = \int_0^\infty e^{-\lambda t}e^{t\widetilde{\mathbf A}}\,dt = \Pi\int_0^\infty e^{-\lambda t}\,dt \;+\; \int_0^\infty e^{-\lambda t}\big(e^{t\widetilde{\mathbf A}}-\Pi\big)\,dt.$$
The first integral is $\Pi/\lambda$ and carries the entire blow-up; the second stays finite at $\lambda=0$ because $e^{t\widetilde{\mathbf A}}-\Pi$ decays. Evaluating it at $\lambda=0$ defines the remainder:
$$D := \int_0^\infty\!\big(e^{t\widetilde{\mathbf A}}-\Pi\big)\,dt, \qquad\text{giving}\qquad (\lambda I-\widetilde{\mathbf A})^{-1} = \frac{\Pi}{\lambda} + D + O(\lambda). \tag{2}$$
This is the exact recurrent mirror of $(\star)$: there we integrated $e^{t\mathbf A_{JJ}}$ directly because it decayed to $\mathbf 0$; here $e^{t\widetilde{\mathbf A}}$ decays only to $\boldsymbol\pi$, so we subtract $\boldsymbol\pi$ first and integrate the remainder. $D$ (the *deviation matrix*) is the accumulated gap between the chain and its stationary law, and satisfies $-\widetilde{\mathbf A}\,D = I-\Pi$: it is $\widetilde{\mathbf A}$ inverted everywhere except along the stationary direction.

**Quantities.** From $D$ one reads the long-time quantities that leakage would have destroyed: the within-$J$ mean first-passage times $m_{\boldsymbol x\boldsymbol y}=(D_{\boldsymbol y\boldsymbol y}-D_{\boldsymbol x\boldsymbol y})/\pi_{\boldsymbol y}$ (the recurrent analog of the mean escape time), and how fast a time-average $\tfrac1T\int_0^T g(X_t)\,dt$ settles onto its stationary value $\sum_{\boldsymbol x}\pi_{\boldsymbol x}g(\boldsymbol x)$.

---


**Next steps.** Demonstrate the escape density $f(t)=\Phi_{\text{out}}(t)$ and variance (1) on the fast/slow example (`examples/fast_slow_else.jl`) against SSA escape-time histograms (the mean is one number off a distribution the resolvent gives in full).