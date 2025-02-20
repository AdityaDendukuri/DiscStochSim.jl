# Discrete Stochastic Simulation Solver

This repository contains a Julia implementation of a discrete stochastic simulation solver. It provides examples based on the Chemical Master Equation (CME) in the `examples/` directory, along with code to reproduce the figures from our paper.

## Overview

The solver is designed for simulating the time evolution of probability distributions over the state space in discrete stochastic systems. In the context of chemical reaction networks, it approximates the **Chemical Master Equation (CME)**:

$$
\frac{d}{dt}\boldsymbol{p}(t) = \mathbf{A}\,\boldsymbol{p}(t),
$$

where \(\boldsymbol{p}(t)\) is the probability vector and \(\mathbf{A}\) is the generator matrix. Due to the typically infinite state space, the solver works on a truncated state space and carefully controls the error introduced by the truncation.

### Quantile-Based Pruning

At each time step, the algorithm prunes a fixed fraction \(\alpha\) of the states with the smallest probabilities by:

1. **Sorting states by probability.**
2. **Computing the cumulative probability:**

   $$
   P_{\mathrm{cumulative}}(k,t) = \sum_{i=1}^{k} p(x_i,t)
   $$

3. **Determining the cutoff index \(k^*\):**  
   \(k^*\) is chosen such that 
   
   $$
   P_{\mathrm{cumulative}}(k^*,t) \ge \alpha \quad \text{and} \quad P_{\mathrm{cumulative}}(k^*-1,t) < \alpha.
   $$
   
   The cutoff probability \(q_\alpha = p(x_{k^*},t)\) is then used to remove all states with probabilities \(\le q_\alpha\).

After pruning, the remaining probability mass is renormalized so that the new probability vector sums to one. Theoretical analysis shows that if the total pruned mass is \(m\), then the error in the \(\ell^1\)-norm is bounded by

$$
\|\boldsymbol{p}(t)-\widetilde{\boldsymbol{p}}(t)\|_1 \le 2\,m.
$$

### Error Control and Adaptive Time Stepping

The solver employs an adaptive time-stepping rule to ensure the local error per step—accounting for both the truncation error and the time-stepping error from approximating the matrix exponential—remains within a prescribed tolerance. The adaptive step size is determined by an inequality of the form

$$
\delta t \le \left(\frac{\epsilon_{\mathrm{tol}} - \epsilon_{\mathrm{trunc}}(t)}{C}\right)^{\!1/p},
$$

where \(C\) and \(p\) are constants derived from the error analysis.

## Usage

Detailed examples are available in the `examples/` directory.


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://AdityaDendukuri.github.io/DiscStochSim.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://AdityaDendukuri.github.io/DiscStochSim.jl/dev/)
[![Build Status](https://github.com/AdityaDendukuri/DiscStochSim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AdityaDendukuri/DiscStochSim.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/AdityaDendukuri/DiscStochSim.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/AdityaDendukuri/DiscStochSim.jl)
