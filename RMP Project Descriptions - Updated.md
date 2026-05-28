## **Project 1: Forward Modeling of Stochastic Biochemical Dynamics**

**Subject Area:** Applied Mathematics, Computational Biology, Scientific Computing

### **Background**

Cells are fundamentally noisy: the same gene can be on in one cell and off in its neighbor, and a single unlucky collision between molecules can trigger a cascade that changes a cell's fate. The Chemical Master Equation (CME) is the gold-standard mathematical framework for capturing this randomness, describing the full probability distribution over every possible molecular count at every moment in time. The catch is that the CME is enormous — a cell with just a few dozen molecular species has a state space larger than the number of atoms in the universe — so solving it directly is impossible without smart algorithms. This project develops and extends a Multigrid finite state projection (FSP) solver that makes CME computation tractable by focusing computational effort where it is actually needed: at the propagating wavefront of probability mass, not the settled interior.

### **Project Description**

This collaborative project is split into two distinct but complementary components, allowing each student to make an independent research contribution while working toward a shared scientific goal. Both will work directly in an established Julia codebase and are expected to write, test, and document their own code.

* **Student 1 — Stochastic Waves on Complex Geometries:** Real biological cells have irregular shapes that profoundly affect how molecules spread through them. This student will extend the multigrid FSP solver from its current rectangular grid to arbitrary spatial geometries. The key insight is that voxels are already stored as a graph — nodes connected by diffusion edges — so the core algorithm needs surprisingly little modification to handle non-rectangular domains. The student will discretize 2D cell outlines (and, time permitting, 3D meshes from imaging data) into voxel graphs, simulate stochastic reaction-diffusion waves propagating across these shapes, and measure how geometry affects wave speed and spatial spread. Along the way they will build visualization tools for wavefronts on irregular domains and compare results against theoretical predictions from reaction-diffusion PDE theory.

* **Student 2 — Optimizing Multigrid Speed and Accuracy:** The efficiency of a multigrid solver rests on two pillars: the smoother, which damps local probability-distribution errors within each grid level, and the coarse-grid operators, which transfer information between fine and coarse representations of the state space. This student will systematically study and improve both. On the smoother side, they will implement and benchmark relaxation strategies — including adaptive Krylov time-stepping and block Gauss-Seidel sweeps — and measure how quickly each one reduces the local error across a single V-cycle. On the coarse-grid side, they will develop and test restriction operators that go beyond the current total-count aggregation, quantify the resulting coarsening error, and tune V-cycle parameters (pre- and post-smoothing steps, coarsening ratio) for optimal convergence. The Schlogl bistable model will serve as the primary benchmark throughout, since its bimodal probability distributions are known to stress-test naive coarsening strategies.

### **Prerequisites**

* **(Required)** Strong mathematical foundation equivalent to AP Calculus.
* **(Preferred)** Experience with programming (Julia or Python); familiarity with linear algebra and basic concepts from differential equations.

### **Tentative Schedule**

| Week | Student 1 (Complex Geometries) | Student 2 (Multigrid Solver) |
| :---- | :---- | :---- |
| Week 1–2 | Literature review: FSP, graph-based spatial discretization, stochastic wave theory. Introduction to the Julia codebase. | Literature review: Multigrid smoothers (Krylov, block methods), restriction and prolongation operators. Introduction to the Julia codebase. |
| Week 3 | Implement general graph-based voxel connectivity; replace Cartesian neighbor function with adjacency-list representation. Verify correctness on the rectangular grid as a sanity check. | Implement adaptive Krylov time-stepping as an alternative smoother; profile performance on the birth-death RDME benchmark. |
| Week 4 | Discretize 2D cell outline geometries into voxel graphs; test wavefront propagation on an L-shaped and a circular domain. | Implement block Gauss-Seidel smoother; compare convergence rate and wall time against the Krylov smoother on the Schlogl model. |
| Week 5 | Measure wave speeds on irregular geometries and compare to theoretical KPP predictions; investigate how boundary curvature affects the wavefront. | Develop improved restriction operators (conditional/marginal-based); quantify coarsening error versus the current total-count aggregation. |
| Week 6 | Begin 3D mesh extension (if ahead of schedule) or deepen analysis on 2D: test multiple cell shapes, study stochastic fluctuations in wave arrival time. | Tune V-cycle parameters (pre-/post-smooth steps, coarsening ratio); run full convergence benchmarks across problem sizes. |
| Week 7–8 | Finalize simulations, document code, generate figures, and prepare the joint presentation. | Finalize performance comparison plots, write up convergence analysis, and prepare the joint presentation. |
