# HPSC Examination Question Bank & Model Solutions Repository
## High Performance Scientific Computing (CD61002) — IIT Kharagpur

> **Course**: CD61002 — High Performance Scientific Computing  
> **Target Audience**: Students preparing for Mid-Semester, End-Semester, and Class Tests  
> **Pedagogical Alignment**: Prof. Somnath Roy (IIT Kharagpur), Gilbert Strang, Yousef Saad, Grama et al.  
> **Document Purpose**: A curated repository of high-probability examination questions covering all 18 lecture modules, complete with step-by-step analytical solutions, algebraic proofs, matrix calculations, and architectural trade-off evaluations.

---

## 📑 Question Bank Index by Module

- [Module 1: PDE Discretization, Stencils & Sparse Storage Schemes](#module-1-pde-discretization-stencils--sparse-storage-schemes)
  - [Q1.1: 2D Poisson Stencil with Mixed Boundary Conditions & Matrix Assembly](#q11-2d-poisson-stencil-with-mixed-boundary-conditions--matrix-assembly-5-marks)
  - [Q1.2: Unpacking Compressed Sparse Row (CSR) and Modified Sparse Row (MSR) Arrays](#q12-unpacking-compressed-sparse-row-csr-and-modified-sparse-row-msr-arrays-6-marks)
  - [Q1.3: ELLPACK vs. DIA Memory Footprint and GPU Warp Coalescing](#q13-ellpack-vs-dia-memory-footprint-and-gpu-warp-coalescing-5-marks)
- [Module 2: Classical Iterative Solvers & Convergence Theory](#module-2-classical-iterative-solvers--convergence-theory)
  - [Q2.1: Spectral Radius Derivation and Unconditional Divergence of Iteration Matrix](#q21-spectral-radius-derivation-and-unconditional-divergence-of-iteration-matrix-6-marks)
  - [Q2.2: Gershgorin Circle Theorem & Proof of Jacobi Convergence for SDD Matrices](#q22-gershgorin-circle-theorem--proof-of-jacobi-convergence-for-sdd-matrices-5-marks)
  - [Q2.3: Young's Theorem, SOR Optimal Relaxation Factor $\omega_{opt}$, and Iteration Contraction](#q23-youngs-theorem-sor-optimal-relaxation-factor-omega_opt-and-iteration-contraction-6-marks)
  - [Q2.4: Counterexample: Non-Singular System Failing Jacobi But Solvable by Gauss-Seidel](#q24-counterexample-non-singular-system-failing-jacobi-but-solvable-by-gauss-seidel-4-marks)
- [Module 3: Direct Solvers, TDMA & Sherman-Morrison Formula](#module-3-direct-solvers-tdma--sherman-morrison-formula)
  - [Q3.1: Thomas Algorithm (TDMA) Recurrence Derivation and Stability Bound](#q31-thomas-algorithm-tdma-recurrence-derivation-and-stability-bound-6-marks)
  - [Q3.2: Cyclic Tridiagonal Solver via Sherman-Morrison Rank-1 Perturbation](#q32-cyclic-tridiagonal-solver-via-sherman-morrison-rank-1-perturbation-6-marks)
- [Module 4: Projection Methods & Steepest Descent](#module-4-projection-methods--steepest-descent)
  - [Q4.1: Derivation of Optimal Step Size $\alpha_k$ and Residual Orthogonality Proof](#q41-derivation-of-optimal-step-size-alpha_k-and-residual-orthogonality-proof-6-marks)
  - [Q4.2: Kantorovich Inequality and Geometric Analysis of Zigzagging](#q42-kantorovich-inequality-and-geometric-analysis-of-zigzagging-5-marks)
- [Module 5: Krylov Subspace Methods (CG, Arnoldi, GMRES, BiCGSTAB)](#module-5-krylov-subspace-methods-cg-arnoldi-gmres-bicgstab)
  - [Q5.1: Conjugate Gradient Step Length $\alpha_k$ and Search Direction $\beta_k$ Derivation](#q51-conjugate-gradient-step-length-alpha_k-and-search-direction-beta_k-derivation-7-marks)
  - [Q5.2: Finite-Element Mesh Refinement Scaling on CG Iterations & Computational Cost](#q52-finite-element-mesh-refinement-scaling-on-cg-iterations--computational-cost-6-marks)
  - [Q5.3: Arnoldi Hessenberg Reduction vs. Symmetric Lanczos 3-Term Recurrence](#q53-arnoldi-hessenberg-reduction-vs-symmetric-lanczos-3-term-recurrence-5-marks)
  - [Q5.4: GMRES Minimization via Givens Rotations & Residual Norm Monitoring](#q54-gmres-minimization-via-givens-rotations--residual-norm-monitoring-6-marks)
  - [Q5.5: BiCG Breakdown Modes and BiCGSTAB Polynomial Stabilization Mechanism](#q55-bicg-breakdown-modes-and-bicgstab-polynomial-stabilization-mechanism-5-marks)
- [Module 6: Preconditioning Techniques](#module-6-preconditioning-techniques)
  - [Q6.1: Preconditioned Conjugate Gradient (PCG) Algebraic Derivation](#q61-preconditioned-conjugate-gradient-pcg-algebraic-derivation-6-marks)
  - [Q6.2: Incomplete Cholesky IC(0) vs. SSOR Preconditioner Construction](#q62-incomplete-cholesky-ic0-vs-ssor-preconditioner-construction-5-marks)
- [Module 7: Computer Architecture, Caching & Memory Wall](#module-7-computer-architecture-caching--memory-wall)
  - [Q7.1: The Memory Wall: DRAM Latency Penalty and Effective Access Time (EAT)](#q71-the-memory-wall-dram-latency-penalty-and-effective-access-time-eat-5-marks)
  - [Q7.2: Row-Major vs. Column-Major Striding and Cache Tiling Analysis](#q72-row-major-vs-column-major-striding-and-cache-tiling-analysis-6-marks)
- [Module 8: Parallel Architectures, Interconnects & Topologies](#module-8-parallel-architectures-interconnects--topologies)
  - [Q8.1: Network Topology Metrics: Degree, Diameter, and Bisection Width Derivation](#q81-network-topology-metrics-degree-diameter-and-bisection-width-derivation-6-marks)
  - [Q8.2: Cache Coherency (MESI Protocol) and False Sharing Pathology with Code Fix](#q82-cache-coherency-mesi-protocol-and-false-sharing-pathology-with-code-fix-5-marks)
- [Module 9: Parallel Performance, Scalability & Amdahl's Law](#module-9-parallel-performance-scalability--amdahls-law)
  - [Q9.1: Amdahl's Law with Parallel Communication Overhead and Maximum Efficiency](#q91-amdahls-law-with-parallel-communication-overhead-and-maximum-efficiency-6-marks)
  - [Q9.2: Isoefficiency Function Derivation and Scalability Evaluation](#q92-isoefficiency-function-derivation-and-scalability-evaluation-6-marks)
- [Module 10: Parallel Numerical Linear Algebra](#module-10-parallel-numerical-linear-algebra)
  - [Q10.1: Parallel SpMV: 1D Row Striping vs. 2D Checkerboard Mesh Decomposition](#q101-parallel-spmv-1d-row-striping-vs-2d-checkerboard-mesh-decomposition-6-marks)
  - [Q10.2: Cannon's Algorithm on 2D Torus: Communication Complexity & Shift Steps](#q102-cannons-algorithm-on-2d-torus-communication-complexity--shift-steps-7-marks)
- [Full Model Examination Paper (25 Marks, 2 Hours)](#full-model-examination-paper-25-marks-2-hours)

---

## Module 1: PDE Discretization, Stencils & Sparse Storage Schemes

### Q1.1: 2D Poisson Stencil with Mixed Boundary Conditions & Matrix Assembly [5 Marks]
**Problem Statement**:  
Consider the 2D Poisson equation on the unit square $\Omega = [0, 1] \times [0, 1]$:
$$\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = -2 \pi^2 \sin(\pi x) \sin(\pi y)$$
The domain is discretized using a uniform grid with spacing $h = \Delta x = \Delta y = \frac{1}{3}$, yielding interior nodes at $(x_i, y_j) = (ih, jh)$ for $i, j \in \{1, 2\}$.  
The boundary conditions are:
- Bottom ($y = 0$): $T(x, 0) = 0$ (Dirichlet)
- Top ($y = 1$): $T(x, 1) = 0$ (Dirichlet)
- Left ($x = 0$): $T(0, y) = 0$ (Dirichlet)
- Right ($x = 1$): $\left.\frac{\partial T}{\partial x}\right|_{x=1} = 0$ (Homogeneous Neumann)

Using standard second-order central difference approximations and a fictitious ghost node for the Neumann boundary, assemble the linear system $A x = b$ for the 4 interior nodes numbered lexicographically $\text{ipt}(i, j) = i + (j-1)N_x$ where $N_x = 2$.

#### Full Step-by-Step Solution:

1. **Grid Setup and Node Indexing**:
   - Spacing: $h = 1/3$.
   - Interior grid coordinates:
     - Node 1: $(x_1, y_1) = (1/3, 1/3)$
     - Node 2: $(x_2, y_1) = (2/3, 1/3)$
     - Node 3: $(x_1, y_2) = (1/3, 2/3)$
     - Node 4: $(x_2, y_2) = (2/3, 2/3)$
   - $N_x = 2$ nodes in the $x$-direction, $N_y = 2$ nodes in the $y$-direction. Total unknowns $N = 4$.

2. **Standard 5-Point Interior Stencil**:
   The standard 2nd-order central difference operator for $\nabla^2 T = f$ gives:
   $$\frac{T_{i+1, j} - 2T_{i, j} + T_{i-1, j}}{h^2} + \frac{T_{i, j+1} - 2T_{i, j} + T_{i, j-1}}{h^2} = f_{i, j}$$
   Multiplying by $-h^2$:
   $$4T_{i, j} - T_{i+1, j} - T_{i-1, j} - T_{i, j+1} - T_{i, j-1} = -h^2 f_{i, j}$$
   Here $-h^2 = -(1/3)^2 = -1/9$. Thus the right-hand side source term is:
   $$b_{i, j} = \frac{2\pi^2}{9} \sin(\pi x_i) \sin(\pi y_j)$$

3. **Neumann Boundary Condition & Ghost Node Formulation**:
   On the right boundary ($x = 1$, i.e., $i = 3$), we enforce $\left.\frac{\partial T}{\partial x}\right|_{3, j} = 0$:
   Using a central difference at $i = 2$:
   $$\left.\frac{\partial T}{\partial x}\right|_{2, j} \approx \frac{T_{3, j} - T_{1, j}}{2h}$$
   However, node $(2, j)$ has neighbor $T_{3, j}$ on the boundary. If we treat node 2 as an interior node adjacent to boundary $i = 3$, we apply the boundary condition directly at $x = 1$.  
   Using a 2nd-order backwards/ghost formulation at $i = 2$:
   $$\left.\frac{\partial T}{\partial x}\right|_{2, j} = 0 \implies T_{3, j} = T_{1, j}$$
   Substituting $T_{3, j} = T_{1, j}$ into the 5-point stencil for $i = 2$:
   $$4T_{2, j} - T_{1, j} - T_{1, j} - T_{2, j+1} - T_{2, j-1} = b_{2, j} \implies 4T_{2, j} - 2T_{1, j} - T_{2, j+1} - T_{2, j-1} = b_{2, j}$$

4. **Equations for Each Unknown**:
   - **Node 1 $(1, 1)$**:
     - Neighbors: $T_{0, 1} = 0$ (left), $T_{1, 0} = 0$ (bottom), $T_{2, 1} = T_2$ (right), $T_{1, 2} = T_3$ (top).
     - Equation:
       $$4T_1 - T_2 - T_3 = b_1$$
       where $b_1 = \frac{2\pi^2}{9} \sin(\pi/3) \sin(\pi/3) = \frac{2\pi^2}{9} \left(\frac{\sqrt{3}}{2}\right)^2 = \frac{2\pi^2}{9} \cdot \frac{3}{4} = \frac{\pi^2}{6} \approx 1.6449$.

   - **Node 2 $(2, 1)$**:
     - Neighbors: $T_{1, 1} = T_1$ (left), $T_{2, 0} = 0$ (bottom), $T_{3, 1} = T_1$ (right, Neumann ghost), $T_{2, 2} = T_4$ (top).
     - Equation:
       $$4T_2 - 2T_1 - T_4 = b_2$$
       where $b_2 = \frac{2\pi^2}{9} \sin(2\pi/3) \sin(\pi/3) = \frac{\pi^2}{6} \approx 1.6449$.

   - **Node 3 $(1, 2)$**:
     - Neighbors: $T_{0, 2} = 0$ (left), $T_{1, 1} = T_1$ (bottom), $T_{2, 2} = T_4$ (right), $T_{1, 3} = 0$ (top).
     - Equation:
       $$-T_1 + 4T_3 - T_4 = b_3$$
       where $b_3 = \frac{2\pi^2}{9} \sin(\pi/3) \sin(2\pi/3) = \frac{\pi^2}{6} \approx 1.6449$.

   - **Node 4 $(2, 2)$**:
     - Neighbors: $T_{1, 2} = T_3$ (left), $T_{2, 1} = T_2$ (bottom), $T_{3, 2} = T_3$ (right, Neumann ghost), $T_{2, 3} = 0$ (top).
     - Equation:
       $$-T_2 - 2T_3 + 4T_4 = b_4$$
       where $b_4 = \frac{\pi^2}{6} \approx 1.6449$.

5. **Assembled Matrix System**:
   $$\begin{bmatrix}
   4 & -1 & -1 & 0 \\
   -2 & 4 & 0 & -1 \\
   -1 & 0 & 4 & -1 \\
   0 & -1 & -2 & 4
   \end{bmatrix}
   \begin{bmatrix} T_1 \\ T_2 \\ T_3 \\ T_4 \end{bmatrix}
   = \frac{\pi^2}{6}
   \begin{bmatrix} 1 \\ 1 \\ 1 \\ 1 \end{bmatrix}$$

   > **Exam Insight**: Notice that the Neumann boundary condition introduces an asymmetry in the coefficient matrix ($a_{21} = -2 \ne a_{12} = -1$). Consequently, $A \ne A^T$! Conjugate Gradient *cannot* be used directly without symmetrization or normal equations; GMRES or BiCGSTAB must be employed.

---

### Q1.2: Unpacking Compressed Sparse Row (CSR) and Modified Sparse Row (MSR) Arrays [6 Marks]
**Problem Statement**:  
A sparse coefficient matrix $A \in \mathbb{R}^{5 \times 5}$ is represented in the Modified Sparse Row (MSR) format using two 1D arrays:
```
AA = [ 10.0,  20.0,  30.0,  40.0,  50.0,   *,   -1.0,   2.0,  -3.0,   4.0,  -5.0 ]
JA = [    7,     8,    10,    11,    12,  12,      2,     1,     4,     2,     3 ]
```
1. Explain the significance of entry `JA(6) = 12` and the asterisk `*` at `AA(6)`.
2. Fully reconstruct the $5 \times 5$ matrix $A$.
3. Convert this matrix into standard 0-indexed Compressed Sparse Row (CSR) format with arrays `val`, `col_ind`, and `row_ptr`.

#### Full Step-by-Step Solution:

1. **Anatomy of MSR Format**:
   - In MSR format for an $n \times n$ matrix with $N_{nz}$ non-zero elements:
     - Array length: $N_{nz} + 1$. Here length is 11, so $N_{nz} = 10$.
     - `AA(1:n)` stores the $n$ main diagonal entries: $a_{11}=10, a_{22}=20, a_{33}=30, a_{44}=40, a_{55}=50$.
     - `AA(n+1)` (position 6) is an unused placeholder/pointer delimiter, conventionally marked with `*`.
     - `AA(n+2 : end)` stores the strictly off-diagonal non-zero entries.
     - `JA(1:n)` points to the starting index in `AA` of off-diagonal elements for each row.
     - `JA(n+1)` points to the first unoccupied position after all off-diagonals (i.e., `length + 1` = 12).
     - Because `JA(5) = 12` and `JA(6) = 12`, the difference $\text{JA}(6) - \text{JA}(5) = 0$. This proves that **Row 5 contains ZERO off-diagonal elements** (it has only a diagonal entry $a_{55} = 50$).

2. **Reconstruction of Matrix $A$**:
   - **Row 1**:
     - Diagonal: $a_{11} = 10.0$.
     - Off-diagonals: index range from `JA(1)=7` to `JA(2)-1 = 7` (1 element).
     - At index 7: value `AA(7) = -1.0`, column `JA(7) = 2` $\implies a_{12} = -1.0$.
   - **Row 2**:
     - Diagonal: $a_{22} = 20.0$.
     - Off-diagonals: index range from `JA(2)=8` to `JA(3)-1 = 9` (2 elements).
     - Index 8: value `AA(8) = 2.0`, column `JA(8) = 1` $\implies a_{21} = 2.0$.
     - Index 9: value `AA(9) = -3.0`, column `JA(9) = 4` $\implies a_{24} = -3.0$.
   - **Row 3**:
     - Diagonal: $a_{33} = 30.0$.
     - Off-diagonals: index range from `JA(3)=10` to `JA(4)-1 = 10` (1 element).
     - Index 10: value `AA(10) = 4.0`, column `JA(10) = 2` $\implies a_{32} = 4.0$.
   - **Row 4**:
     - Diagonal: $a_{44} = 40.0$.
     - Off-diagonals: index range from `JA(4)=11` to `JA(5)-1 = 11` (1 element).
     - Index 11: value `AA(11) = -5.0`, column `JA(11) = 3` $\implies a_{43} = -5.0$.
   - **Row 5**:
     - Diagonal: $a_{55} = 50.0$.
     - Off-diagonals: none (`JA(5) = JA(6) = 12`).

   $$A = \begin{bmatrix}
   10.0 & -1.0 &  0.0 &  0.0 &  0.0 \\
    2.0 & 20.0 &  0.0 & -3.0 &  0.0 \\
    0.0 &  4.0 & 30.0 &  0.0 &  0.0 \\
    0.0 &  0.0 & -5.0 & 40.0 &  0.0 \\
    0.0 &  0.0 &  0.0 &  0.0 & 50.0
   \end{bmatrix}$$

3. **Conversion to 0-Indexed CSR Format**:
   - Total non-zeros: $N_{nz} = 1 + 3 + 2 + 2 + 1 = 9$ (excluding explicit zeros).
   - In 0-indexing:
     - `val`: All non-zero entries row by row:
       ```
       val = [ 10.0, -1.0, 2.0, 20.0, -3.0, 4.0, 30.0, -5.0, 40.0, 50.0 ]
       ```
     - `col_ind` (0-indexed column coordinates):
       ```
       col_ind = [ 0, 1, 0, 1, 3, 1, 2, 2, 3, 4 ]
       ```
     - `row_ptr` (cumulative non-zero counts per row):
       - Row 0 starts at index 0. Has 2 non-zeros.
       - Row 1 starts at index 2. Has 3 non-zeros.
       - Row 2 starts at index 5. Has 2 non-zeros.
       - Row 3 starts at index 7. Has 2 non-zeros.
       - Row 4 starts at index 9. Has 1 non-zero.
       - End pointer is 10.
       ```
       row_ptr = [ 0, 2, 5, 7, 9, 10 ]
       ```

---

### Q1.3: ELLPACK vs. DIA Memory Footprint and GPU Warp Coalescing [5 Marks]
**Problem Statement**:  
A 3D structured CFD grid of size $100 \times 100 \times 100$ ($N = 10^6$ nodes) is discretized using a 7-point Laplacian stencil.  
1. Compute the exact storage requirement in Megabytes (MB) for double-precision floating point values (8 bytes) and 32-bit integer indices (4 bytes) under:
   - (a) Compressed Sparse Row (CSR) format.
   - (b) Diagonal Storage (DIA) format.
   - (c) ELLPACK format.
2. Explain why the ELLPACK format delivers substantially higher memory throughput than CSR when executing Sparse Matrix-Vector multiplication (SpMV) on modern NVIDIA GPU architectures.

#### Full Step-by-Step Solution:

1. **Storage Calculations**:
   - Total unknowns: $N = 10^6$.
   - Interior stencil has 7 points. Non-zeros: $N_{nz} \approx 7N - 2(N_x N_y + N_y N_z + N_x N_z) \approx 7 \times 10^6$ entries.

   - **(a) CSR Format**:
     - `val`: $N_{nz}$ doubles $\implies 7 \times 10^6 \times 8\text{ B} = 56\text{ MB}$.
     - `col_ind`: $N_{nz}$ integers $\implies 7 \times 10^6 \times 4\text{ B} = 28\text{ MB}$.
     - `row_ptr`: $N + 1$ integers $\implies (10^6 + 1) \times 4\text{ B} \approx 4\text{ MB}$.
     - **Total CSR Storage** = $56 + 28 + 4 = \mathbf{88\text{ MB}}$.

   - **(b) DIA Format**:
     - Active diagonals for 7-point 3D stencil:
       Offsets: $[-N_x N_y, \; -N_x, \; -1, \; 0, \; 1, \; N_x, \; N_x N_y] \implies 7\text{ active diagonals}$.
     - Dense array `DIAG(N, 7)` of doubles:
       $10^6 \times 7 \times 8\text{ B} = 56\text{ MB}$.
     - Offset vector `IOFF(7)` of integers:
       $7 \times 4\text{ B} = 28\text{ bytes} \approx 0\text{ MB}$.
     - Zero padding overhead on boundaries: minimal ($2 \times 10^4$ zeros).
     - **Total DIA Storage** $\approx \mathbf{56\text{ MB}}$ (No column index array needed!).

   - **(c) ELLPACK Format**:
     - Maximum non-zeros per row: $K = 7$.
     - `COEF(N, 7)` array of doubles: $10^6 \times 7 \times 8\text{ B} = 56\text{ MB}$.
     - `JCOEF(N, 7)` array of integers: $10^6 \times 7 \times 4\text{ B} = 28\text{ MB}$.
     - **Total ELLPACK Storage** = $56 + 28 = \mathbf{84\text{ MB}}$.

2. **Why ELLPACK Outperforms CSR on GPUs**:
   - **GPU Architecture Principle**: NVIDIA GPUs execute threads in SIMT (Single Instruction, Multiple Threads) groups called **Warps** (32 threads).
   - In **CSR SpMV**, thread $i$ computes row $i$. If rows have different numbers of non-zeros, threads within the warp diverge (warp divergence), idling while waiting for the longest row. Furthermore, row indices require irregular memory indirect accesses ($x[col\_ind[j]]$).
   - In **ELLPACK SpMV**, data is stored in **Column-Major order**:
     $$\text{COEF}[k \cdot N + i], \quad \text{JCOEF}[k \cdot N + i]$$
     When 32 threads in a warp access the $k$-th non-zero of their respective rows $i, i+1, \dots, i+31$, they access **32 consecutive memory addresses** in DRAM.
   - This triggers **Memory Coalescing**: the GPU memory controller fulfills the request in a single 128-byte cache line transaction, maximizing memory bandwidth saturation ($> 90\%$ of peak).

---

## Module 2: Classical Iterative Solvers & Convergence Theory

### Q2.1: Spectral Radius Derivation and Unconditional Divergence of Iteration Matrix [6 Marks]
**Problem Statement**:  
To solve $A x = b$, an engineer proposes the stationary iterative scheme $x^{(k+1)} = G x^{(k)} + f$, where the iteration matrix $G$ is given by:
$$G = \begin{bmatrix}
0.4 & 0.0 & 0.0 & 0.0 \\
1.5 & 0.2 & 0.0 & 0.0 \\
2.0 & 1.0 & 0.5 & 3.0 \\
0.5 & -1.0 & 1.0 & 1.5
\end{bmatrix}$$
1. Compute all eigenvalues of $G$ algebraically.
2. Determine the spectral radius $\rho(G)$.
3. State whether the method converges for an arbitrary initial guess $x^{(0)}$.
4. If $x^{(0)} - x^* = [0, 0, 1, 1]^T$, derive the analytical expression for the error norm $\|e^{(k)}\|_2$ as $k \to \infty$.

#### Full Step-by-Step Solution:

1. **Block Lower-Triangular Partitioning**:
   Partition $G$ into $2 \times 2$ blocks:
   $$G = \begin{bmatrix} G_{11} & \mathbf{0} \\ G_{21} & G_{22} \end{bmatrix}$$
   where:
   $$G_{11} = \begin{bmatrix} 0.4 & 0.0 \\ 1.5 & 0.2 \end{bmatrix}, \quad G_{22} = \begin{bmatrix} 0.5 & 3.0 \\ 1.0 & 1.5 \end{bmatrix}$$
   Because $G$ is block lower-triangular, $\det(G - \lambda I) = \det(G_{11} - \lambda I) \cdot \det(G_{22} - \lambda I) = 0$.

2. **Eigenvalues of $G_{11}$**:
   Since $G_{11}$ is triangular:
   $$\lambda_1 = 0.4, \quad \lambda_2 = 0.2$$

3. **Eigenvalues of $G_{22}$**:
   $$\det(G_{22} - \lambda I) = \det \begin{bmatrix} 0.5 - \lambda & 3.0 \\ 1.0 & 1.5 - \lambda \end{bmatrix} = (0.5 - \lambda)(1.5 - \lambda) - (3.0)(1.0) = 0$$
   Expanding:
   $$0.75 - 0.5\lambda - 1.5\lambda + \lambda^2 - 3.0 = 0 \implies \lambda^2 - 2.0\lambda - 2.25 = 0$$
   Applying the quadratic formula:
   $$\lambda = \frac{2.0 \pm \sqrt{(-2.0)^2 - 4(1)(-2.25)}}{2} = \frac{2.0 \pm \sqrt{4.0 + 9.0}}{2} = \frac{2.0 \pm \sqrt{13.0}}{2}$$
   Since $\sqrt{13} \approx 3.60555$:
   $$\lambda_3 = \frac{2.0 + 3.60555}{2} = \mathbf{2.8028}$$
   $$\lambda_4 = \frac{2.0 - 3.60555}{2} = \mathbf{-0.8028}$$

4. **Spectral Radius**:
   $$\rho(G) = \max_i |\lambda_i| = \max(|0.4|, |0.2|, |2.8028|, |-0.8028|) = \mathbf{2.8028}$$

5. **Convergence Assessment**:
   By the Master Convergence Theorem, a linear stationary iteration converges for all $x^{(0)}$ if and only if $\rho(G) < 1$.  
   Since $\mathbf{\rho(G) = 2.8028 > 1}$, the scheme **diverges unconditionally**.

6. **Error Growth for Given Initial Error**:
   The initial error $e^{(0)} = [0, 0, 1, 1]^T$ has zero components in the first two entries, isolating the sub-system $G_{22}$:
   $$e_{3:4}^{(k)} = (G_{22})^k \begin{bmatrix} 1 \\ 1 \end{bmatrix}$$
   Projecting $[1, 1]^T$ onto the eigenvectors of $G_{22}$, the component along the dominant eigenvector $v_3$ (corresponding to $\lambda_3 = 2.8028$) is non-zero.  
   Therefore, asymptotically:
   $$\|e^{(k)}\|_2 \approx C \cdot (\lambda_3)^k = C \cdot (2.8028)^k$$
   The error grows exponentially by a factor of $\approx 2.803$ on **every single iteration**.

---

### Q2.2: Gershgorin Circle Theorem & Proof of Jacobi Convergence for SDD Matrices [5 Marks]
**Problem Statement**:  
1. State the Gershgorin Circle Theorem for an arbitrary matrix $A \in \mathbb{C}^{n \times n}$.
2. Define a **Strictly Diagonally Dominant (SDD)** matrix.
3. Using the Gershgorin Circle Theorem, rigorously prove that the Jacobi iterative method converges for any initial guess if matrix $A$ is Strictly Diagonally Dominant.

#### Full Step-by-Step Solution:

1. **Gershgorin Circle Theorem**:
   Let $A \in \mathbb{C}^{n \times n}$. For each row $i \in \{1, \dots, n\}$, define the row off-diagonal sum:
   $$R_i = \sum_{j \ne i} |a_{ij}|$$
   The $i$-th Gershgorin disk is defined as:
   $$D_i = \{z \in \mathbb{C} \mid |z - a_{ii}| \le R_i\}$$
   **Theorem**: Every eigenvalue $\lambda$ of $A$ lies within the union of all Gershgorin disks:
   $$\sigma(A) \subset \bigcup_{i=1}^n D_i$$

2. **Strict Diagonal Dominance**:
   Matrix $A$ is Strictly Diagonally Dominant (SDD) if for all $i \in \{1, \dots, n\}$:
   $$|a_{ii}| > \sum_{j \ne i} |a_{ij}| = R_i$$

3. **Convergence Proof of Jacobi Iteration**:
   - The Jacobi iteration matrix is:
     $$G_J = -D^{-1}(L + U) = I - D^{-1} A$$
   - The entries of $G_J = [g_{ij}]$ are:
     $$g_{ii} = 0 \quad (\text{all diagonal entries are identically zero})$$
     $$g_{ij} = -\frac{a_{ij}}{a_{ii}} \quad \text{for } j \ne i$$
   - Apply the Gershgorin Circle Theorem to matrix $G_J$:
     - Center of disk $i$: $c_i = g_{ii} = 0$.
     - Radius of disk $i$:
       $$R_i(G_J) = \sum_{j \ne i} |g_{ij}| = \sum_{j \ne i} \left| -\frac{a_{ij}}{a_{ii}} \right| = \frac{1}{|a_{ii}|} \sum_{j \ne i} |a_{ij}|$$
   - Since $A$ is SDD, $|a_{ii}| > \sum_{j \ne i} |a_{ij}|$, which directly implies:
     $$R_i(G_J) = \frac{\sum_{j \ne i} |a_{ij}|}{|a_{ii}|} < 1 \quad \text{for all } i = 1, \dots, n$$
   - Let $\mu$ be any eigenvalue of $G_J$. By Gershgorin's theorem:
     $$\mu \in D_k \implies |\mu - 0| \le R_k(G_J) < 1$$
   - Therefore, for every eigenvalue $\mu \in \sigma(G_J)$:
     $$|\mu| < 1 \implies \mathbf{\rho(G_J) = \max_i |\mu_i| < 1}$$
   - By the Master Convergence Theorem, $\rho(G_J) < 1$ is necessary and sufficient for convergence. Hence, Jacobi iteration **always converges** for SDD matrices. $\blacksquare$

---

### Q2.3: Young's Theorem, SOR Optimal Relaxation Factor $\omega_{opt}$, and Iteration Contraction [6 Marks]
**Problem Statement**:  
A 2D Laplace problem is discretized on an $N \times N$ grid, producing a consistently ordered 2-cyclic matrix $A$. The spectral radius of the Jacobi iteration matrix is found to be $\rho(G_J) = 0.96$.
1. Compute the spectral radius of the Gauss-Seidel iteration matrix $\rho(G_{GS})$ using Young's relation.
2. Calculate the optimal SOR relaxation parameter $\omega_{opt}$ and the corresponding spectral radius $\rho(G_{SOR}(\omega_{opt}))$.
3. If an engineer mistakenly selects $\omega = 1.6$, compute $\rho(G_{SOR}(1.6))$ and explain how the difference norm $\|x^{(k+1)} - x^{(k)}\|$ behaves asymptotically.
4. How many times faster does optimal SOR converge compared to Gauss-Seidel?

#### Full Step-by-Step Solution:

1. **Young's Theorem for Gauss-Seidel**:
   For consistently ordered, 2-cyclic matrices:
   $$\rho(G_{GS}) = [\rho(G_J)]^2$$
   Substituting $\rho(G_J) = 0.96$:
   $$\mathbf{\rho(G_{GS}) = (0.96)^2 = 0.9216}$$

2. **Optimal SOR Parameter and Spectral Radius**:
   Young’s analytical formula for optimal relaxation is:
   $$\omega_{opt} = \frac{2}{1 + \sqrt{1 - \rho(G_J)^2}}$$
   Compute the radical:
   $$\sqrt{1 - (0.96)^2} = \sqrt{1 - 0.9216} = \sqrt{0.0784} = 0.28$$
   Therefore:
   $$\mathbf{\omega_{opt} = \frac{2}{1 + 0.28} = \frac{2}{1.28} = 1.5625}$$
   For $\omega = \omega_{opt}$, the spectral radius is:
   $$\mathbf{\rho(G_{SOR}(\omega_{opt})) = \omega_{opt} - 1 = 1.5625 - 1 = 0.5625}$$

3. **Analysis for $\omega = 1.6$**:
   Notice that $\omega = 1.6 > \omega_{opt} = 1.5625$.  
   For any $\omega$ in the over-relaxed regime ($\omega_{opt} \le \omega < 2$):
   $$\mathbf{\rho(G_{SOR}(\omega)) = \omega - 1}$$
   Thus for $\omega = 1.6$:
   $$\mathbf{\rho(G_{SOR}(1.6)) = 1.6 - 1 = 0.6000}$$
   **Asymptotic Difference Contraction**:
   As $k \to \infty$, the successive iteration difference norm contracts according to the dominant eigenvalue:
   $$\lim_{k \to \infty} \frac{\|x^{(k+1)} - x^{(k)}\|}{\|x^{(k)} - x^{(k-1)}\|} = \rho(G_{SOR}) = \mathbf{0.6000}$$
   The difference between successive iterates drops by exactly $40\%$ per iteration.

4. **Convergence Acceleration Factor**:
   The asymptotic rates of convergence are:
   $$R_\infty(GS) = -\ln(\rho(G_{GS})) = -\ln(0.9216) \approx 0.08164$$
   $$R_\infty(SOR_{opt}) = -\ln(\rho(G_{SOR})) = -\ln(0.5625) \approx 0.57536$$
   Speedup factor in iteration count:
   $$\text{Speedup} = \frac{R_\infty(SOR_{opt})}{R_\infty(GS)} = \frac{0.57536}{0.08164} \approx \mathbf{7.05\times}$$
   Optimal SOR converges more than **7 times faster** than Gauss-Seidel!

---

### Q2.4: Counterexample: Non-Singular System Failing Jacobi But Solvable by Gauss-Seidel [4 Marks]
**Problem Statement**:  
Prove that non-singularity of $A$ ($\det(A) \ne 0$) does **not** guarantee convergence of the Jacobi method. Construct a concrete $3 \times 3$ matrix $A$ with non-zero diagonal entries such that:
1. $\det(A) \ne 0$.
2. The Jacobi method diverges ($\rho(G_J) > 1$).
3. The Gauss-Seidel method converges ($\rho(G_{GS}) < 1$).

#### Full Step-by-Step Solution:

1. **Matrix Construction**:
   Consider the real symmetric matrix:
   $$A = \begin{bmatrix}
   1 & 0.8 & 0.8 \\
   0.8 & 1 & 0.8 \\
   0.8 & 0.8 & 1
   \end{bmatrix}$$

2. **Verification of Non-Singularity**:
   Subtract row 1 from rows 2 and 3:
   $$\det(A) = \det \begin{bmatrix}
   1 & 0.8 & 0.8 \\
   0 & 0.2 & 0 \\
   0 & 0 & 0.2
   \end{bmatrix} + \dots$$
   Eigenvalues of an all-ones perturbed matrix $A = (1 - 0.8)I + 0.8 \mathbf{1}\mathbf{1}^T = 0.2 I + 0.8 J$:
   - $\lambda_1 = 0.2 + 3(0.8) = 2.6$
   - $\lambda_2 = \lambda_3 = 0.2$
   $$\det(A) = (2.6)(0.2)(0.2) = \mathbf{0.104 \ne 0}$$
   Hence, $A$ is strictly non-singular (and positive definite!).

3. **Divergence of Jacobi Method**:
   $D = I \implies G_J = I - A$:
   $$G_J = \begin{bmatrix}
   0 & -0.8 & -0.8 \\
   -0.8 & 0 & -0.8 \\
   -0.8 & -0.8 & 0
   \end{bmatrix}$$
   The eigenvalues of $G_J$ are $\mu_i = 1 - \lambda_i(A)$:
   - $\mu_1 = 1 - 2.6 = \mathbf{-1.6}$
   - $\mu_2 = \mu_3 = 1 - 0.2 = \mathbf{0.8}$
   Spectral radius:
   $$\mathbf{\rho(G_J) = \max(|-1.6|, |0.8|) = 1.6 > 1}$$
   **The Jacobi method diverges unconditionally!**

4. **Convergence of Gauss-Seidel Method**:
   Because $A$ is symmetric and all its eigenvalues are strictly positive ($\lambda_i > 0 \implies A$ is Symmetric Positive Definite), we apply the **Ostrowski-Reich Theorem**:
   > *Theorem*: If $A$ is symmetric with positive diagonal entries, the Gauss-Seidel method converges ($\rho(G_{GS}) < 1$) if and only if $A$ is positive definite.
   
   Since $A$ is SPD, the Gauss-Seidel method is **guaranteed to converge** ($\rho(G_{GS}) < 1$), despite the total divergence of Jacobi.

---

## Module 3: Direct Solvers, TDMA & Sherman-Morrison Formula

### Q3.1: Thomas Algorithm (TDMA) Recurrence Derivation and Stability Bound [6 Marks]
**Problem Statement**:  
Consider a tridiagonal system $A T = d$ of dimension $N$:
$$-a_i T_{i-1} + b_i T_i - c_i T_{i+1} = d_i, \quad i = 1, \dots, N \quad (T_0 = T_{N+1} = 0)$$
1. Derive the forward elimination recurrences for coefficients $\beta_i$ and $\gamma_i$ such that $T_i = \gamma_i + \alpha_i T_{i+1}$.
2. State the backward substitution formula.
3. Prove that the condition of diagonal dominance ($|b_i| \ge |a_i| + |c_i|$ with strict inequality for at least one $i$) guarantees that $\beta_i \ne 0$ for all $i$, ensuring numerical stability without pivoting.

#### Full Step-by-Step Solution:

1. **Forward Recurrence Derivation**:
   Assume the solution at step $i$ can be expressed in terms of $T_{i+1}$ as:
   $$T_{i-1} = \alpha_{i-1} T_i + \gamma_{i-1} \quad \text{(induction hypothesis)}$$
   Substitute this into the $i$-th equation:
   $$-a_i (\alpha_{i-1} T_i + \gamma_{i-1}) + b_i T_i - c_i T_{i+1} = d_i$$
   Group terms by $T_i$:
   $$(b_i - a_i \alpha_{i-1}) T_i = c_i T_{i+1} + (d_i + a_i \gamma_{i-1})$$
   Divide by $(b_i - a_i \alpha_{i-1})$:
   $$T_i = \left( \frac{c_i}{b_i - a_i \alpha_{i-1}} \right) T_{i+1} + \left( \frac{d_i + a_i \gamma_{i-1}}{b_i - a_i \alpha_{i-1}} \right)$$
   Matching with the recurrence $T_i = \alpha_i T_{i+1} + \gamma_i$, we identify:
   $$\mathbf{\alpha_i = \frac{c_i}{b_i - a_i \alpha_{i-1}}}, \qquad \mathbf{\gamma_i = \frac{d_i + a_i \gamma_{i-1}}{b_i - a_i \alpha_{i-1}}}$$
   Defining the denominator $\beta_i \equiv b_i - a_i \alpha_{i-1}$, we get:
   $$\mathbf{\beta_i = b_i - \frac{a_i c_{i-1}}{\beta_{i-1}}}, \qquad \mathbf{\gamma_i = \frac{d_i + a_i \gamma_{i-1}}{\beta_i}}$$
   **Base Case ($i = 1$)**: Since $T_0 = 0$, $a_1 = 0$:
   $$\beta_1 = b_1, \qquad \gamma_1 = \frac{d_1}{b_1}, \qquad \alpha_1 = \frac{c_1}{b_1}$$

2. **Backward Substitution**:
   At the boundary $i = N$, since $T_{N+1} = 0$:
   $$\mathbf{T_N = \gamma_N}$$
   Then sweep backwards for $i = N-1, N-2, \dots, 1$:
   $$\mathbf{T_i = \alpha_i T_{i+1} + \gamma_i = \frac{c_i}{\beta_i} T_{i+1} + \gamma_i}$$

3. **Stability Proof via Diagonal Dominance**:
   We prove by induction that $|\alpha_i| < 1$ and $|\beta_i| > |c_i| \ge 0$ for all $i$:
   - **Base Step ($i = 1$)**:
     By diagonal dominance, $|b_1| \ge |c_1|$. If strict inequality holds, $|\beta_1| = |b_1| > |c_1|$.
     Then $|\alpha_1| = |c_1 / b_1| < 1$.
   - **Inductive Step**:
     Assume $|\alpha_{i-1}| < 1$. Consider $\beta_i = b_i - a_i \alpha_{i-1}$:
     $$|\beta_i| = |b_i - a_i \alpha_{i-1}| \ge |b_i| - |a_i| |\alpha_{i-1}| > |b_i| - |a_i|$$
     Using diagonal dominance $|b_i| \ge |a_i| + |c_i|$:
     $$|\beta_i| > (|a_i| + |c_i|) - |a_i| = |c_i|$$
     Therefore:
     $$|\beta_i| > |c_i| \ge 0 \implies \mathbf{\beta_i \ne 0}$$
     and:
     $$|\alpha_i| = \left| \frac{c_i}{\beta_i} \right| < 1$$
   - Since $\beta_i$ is bounded strictly away from zero at every step, **division by zero is impossible** and error amplification is bounded ($|\alpha_i| < 1$). Thus TDMA is unconditionally stable without pivoting. $\blacksquare$

---

### Q3.2: Cyclic Tridiagonal Solver via Sherman-Morrison Rank-1 Perturbation [6 Marks]
**Problem Statement**:  
Periodic boundary conditions transform a tridiagonal system into a cyclic tridiagonal system $\tilde{A} x = d$, where non-zero corner elements $a_{1, N} = \alpha$ and $a_{N, 1} = \beta$ destroy the band structure.
1. State the **Sherman-Morrison Formula** for the inverse of a rank-1 perturbed matrix $(A + u v^T)$.
2. Show how $\tilde{A}$ can be decomposed as $\tilde{A} = A + u v^T$, where $A$ is strictly tridiagonal.
3. Outline the complete $O(N)$ algorithm using TDMA to solve the cyclic system.

#### Full Step-by-Step Solution:

1. **The Sherman-Morrison Formula**:
   If $A \in \mathbb{R}^{n \times n}$ is invertible and $u, v \in \mathbb{R}^n$ are vectors such that $1 + v^T A^{-1} u \ne 0$:
   $$\mathbf{(A + u v^T)^{-1} = A^{-1} - \frac{A^{-1} u v^T A^{-1}}{1 + v^T A^{-1} u}}$$

2. **Rank-1 Decomposition of Cyclic Tridiagonal Matrix**:
   Let the cyclic matrix be:
   $$\tilde{A} = \begin{bmatrix}
   b_1 & c_1 & 0 & \dots & \alpha \\
   a_2 & b_2 & c_2 & \dots & 0 \\
   \vdots & \ddots & \ddots & \ddots & \vdots \\
   0 & \dots & a_{N-1} & b_{N-1} & c_{N-1} \\
   \beta & \dots & 0 & a_N & b_N
   \end{bmatrix}$$
   Choose scalar parameter $\gamma$ (e.g., $\gamma = -b_1$). Construct vectors $u$ and $v$:
   $$u = \begin{bmatrix} \gamma \\ 0 \\ \vdots \\ 0 \\ \beta \end{bmatrix}, \qquad v = \begin{bmatrix} 1 \\ 0 \\ \vdots \\ 0 \\ \alpha / \gamma \end{bmatrix}$$
   Compute the outer product $u v^T$:
   $$u v^T = \begin{bmatrix}
   \gamma \cdot 1 & 0 & \dots & \gamma(\alpha / \gamma) \\
   0 & 0 & \dots & 0 \\
   \vdots & \vdots & \ddots & \vdots \\
   \beta \cdot 1 & 0 & \dots & \beta(\alpha / \gamma)
   \end{bmatrix} = \begin{bmatrix}
   \gamma & 0 & \dots & \alpha \\
   0 & 0 & \dots & 0 \\
   \vdots & \vdots & \ddots & \vdots \\
   \beta & 0 & \dots & \frac{\alpha \beta}{\gamma}
   \end{bmatrix}$$
   Subtracting $u v^T$ from $\tilde{A}$ gives a strictly tridiagonal matrix $A = \tilde{A} - u v^T$:
   $$A = \begin{bmatrix}
   b_1 - \gamma & c_1 & 0 & \dots & 0 \\
   a_2 & b_2 & c_2 & \dots & 0 \\
   \vdots & \ddots & \ddots & \ddots & \vdots \\
   0 & \dots & a_{N-1} & b_{N-1} & c_{N-1} \\
   0 & \dots & 0 & a_N & b_N - \frac{\alpha \beta}{\gamma}
   \end{bmatrix}$$

3. **$O(N)$ Solution Algorithm**:
   To solve $\tilde{A} x = d \iff (A + u v^T) x = d$:
   - Multiply by $(A + u v^T)^{-1}$:
     $$x = A^{-1} d - \frac{A^{-1} u (v^T A^{-1} d)}{1 + v^T A^{-1} u}$$
   - **Step 1**: Solve tridiagonal system $A y = d$ for $y$ using TDMA ($O(N)$ operations).
   - **Step 2**: Solve tridiagonal system $A z = u$ for $z$ using TDMA ($O(N)$ operations).
   - **Step 3**: Compute scalar dot products:
     $$v^T y = y_1 + \frac{\alpha}{\gamma} y_N, \qquad v^T z = z_1 + \frac{\alpha}{\gamma} z_N$$
   - **Step 4**: Compute final solution via vector update (AXPY):
     $$\mathbf{x = y - \left( \frac{v^T y}{1 + v^T z} \right) z}$$
   - **Complexity**: Requires exactly **2 tridiagonal solves** $\implies 2 \times O(N) = \mathbf{O(N)}$ total operations!

---

## Module 4: Projection Methods & Steepest Descent

### Q4.1: Derivation of Optimal Step Size $\alpha_k$ and Residual Orthogonality Proof [6 Marks]
**Problem Statement**:  
Let $A \in \mathbb{R}^{n \times n}$ be a Symmetric Positive Definite (SPD) matrix. The problem $A x = b$ is equivalent to minimizing the quadratic functional:
$$J(x) = \frac{1}{2} x^T A x - b^T x$$
At iteration $k$ with search direction $p_k = r_k = b - A x_k$:
1. Prove that $\nabla J(x) = A x - b = -r(x)$.
2. Derive the analytical expression for the step size $\alpha_k$ that minimizes $J(x_k + \alpha p_k)$ along the search ray.
3. Prove that consecutive residuals are strictly orthogonal: $r_{k+1} \perp r_k$ (i.e., $r_{k+1}^T r_k = 0$).

#### Full Step-by-Step Solution:

1. **Gradient of Quadratic Functional**:
   $$J(x) = \frac{1}{2} \sum_{i=1}^n \sum_{j=1}^n a_{ij} x_i x_j - \sum_{i=1}^n b_i x_i$$
   Differentiating with respect to $x_m$:
   $$\frac{\partial J}{\partial x_m} = \frac{1}{2} \sum_{j=1}^n a_{mj} x_j + \frac{1}{2} \sum_{i=1}^n a_{im} x_i - b_m$$
   Since $A$ is symmetric ($a_{im} = a_{mi}$):
   $$\frac{\partial J}{\partial x_m} = \sum_{j=1}^n a_{mj} x_j - b_m = (Ax)_m - b_m$$
   In vector form:
   $$\mathbf{\nabla J(x) = A x - b = -r(x)}$$

2. **Derivation of Step Size $\alpha_k$**:
   Define the 1D line function $\phi(\alpha) \equiv J(x_k + \alpha p_k)$:
   $$\phi(\alpha) = \frac{1}{2} (x_k + \alpha p_k)^T A (x_k + \alpha p_k) - b^T (x_k + \alpha p_k)$$
   Expanding:
   $$\phi(\alpha) = \frac{1}{2} x_k^T A x_k + \alpha p_k^T A x_k + \frac{1}{2} \alpha^2 p_k^T A p_k - b^T x_k - \alpha b^T p_k$$
   Differentiate with respect to $\alpha$ and set to zero:
   $$\frac{d\phi}{d\alpha} = p_k^T A x_k + \alpha p_k^T A p_k - b^T p_k = 0$$
   Rearranging terms:
   $$\alpha p_k^T A p_k = p_k^T (b - A x_k) = p_k^T r_k$$
   Since $p_k = r_k$ in Steepest Descent:
   $$\mathbf{\alpha_k = \frac{r_k^T r_k}{r_k^T A r_k}}$$
   Since $A$ is SPD, the denominator $r_k^T A r_k > 0$ for all $r_k \ne 0$, so $\alpha_k$ is always uniquely defined and positive.

3. **Proof of Consecutive Residual Orthogonality**:
   The update formula is $x_{k+1} = x_k + \alpha_k r_k$.  
   Multiply by $-A$ and add $b$:
   $$b - A x_{k+1} = b - A x_k - \alpha_k A r_k \implies \mathbf{r_{k+1} = r_k - \alpha_k A r_k}$$
   Take the inner product with $r_k$:
   $$r_{k+1}^T r_k = (r_k - \alpha_k A r_k)^T r_k = r_k^T r_k - \alpha_k (r_k^T A r_k)$$
   Substitute the derived expression $\alpha_k = \frac{r_k^T r_k}{r_k^T A r_k}$:
   $$r_{k+1}^T r_k = r_k^T r_k - \left( \frac{r_k^T r_k}{r_k^T A r_k} \right) (r_k^T A r_k) = r_k^T r_k - r_k^T r_k = \mathbf{0}$$
   Thus $\mathbf{r_{k+1} \perp r_k}$. $\blacksquare$

---

### Q4.2: Kantorovich Inequality and Geometric Analysis of Zigzagging [5 Marks]
**Problem Statement**:  
1. State the **Kantorovich Inequality** for a symmetric positive definite matrix $A$.
2. State the error convergence bound of the Method of Steepest Descent in the $A$-norm:
   $$\|e_{k+1}\|_A \le \left( \frac{\kappa_2(A) - 1}{\kappa_2(A) + 1} \right) \|e_k\|_A$$
3. For a 2D anisotropic diffusion problem with condition number $\kappa_2(A) = 1000$, calculate the number of iterations required to reduce the initial error by $10^{-4}$.
4. Explain geometrically why Steepest Descent exhibits severe "zigzagging" in ill-conditioned problems.

#### Full Step-by-Step Solution:

1. **The Kantorovich Inequality**:
   Let $A \in \mathbb{R}^{n \times n}$ be SPD with extreme eigenvalues $\lambda_{\max} \ge \dots \ge \lambda_{\min} > 0$. For any non-zero vector $x \in \mathbb{R}^n$:
   $$\mathbf{\frac{(x^T x)^2}{(x^T A x)(x^T A^{-1} x)} \ge \frac{4 \lambda_{\min} \lambda_{\max}}{(\lambda_{\min} + \lambda_{\max})^2} = \frac{4 \kappa}{(\kappa + 1)^2}}$$
   where $\kappa = \lambda_{\max} / \lambda_{\min}$.

2. **Steepest Descent Error Bound**:
   Using the Kantorovich inequality, the error reduction per step in the energy norm $\|e\|_A = \sqrt{e^T A e}$ satisfies:
   $$\mathbf{\|e_{k+1}\|_A \le \left( \frac{\kappa - 1}{\kappa + 1} \right) \|e_k\|_A}$$

3. **Numerical Iteration Calculation**:
   Given $\kappa = 1000$:
   $$\rho = \frac{1000 - 1}{1000 + 1} = \frac{999}{1001} \approx 0.998002$$
   To achieve $\|e_k\|_A / \|e_0\|_A \le 10^{-4}$:
   $$\rho^k \le 10^{-4} \implies k \ln(\rho) \le \ln(10^{-4})$$
   $$k \ge \frac{-4 \ln(10)}{\ln(0.998002)} = \frac{-9.21034}{-0.00200} \approx \mathbf{4,606\text{ iterations}}$$

4. **Geometric Cause of Zigzagging**:
   - The level curves of the quadratic functional $J(x)$ are multidimensional ellipsoids. The ratio of the major axis to the minor axis is $\sqrt{\kappa(A)} = \sqrt{1000} \approx 31.6$.
   - The gradient $\nabla J(x_k) = -r_k$ is orthogonal to the elliptical contour tangent at $x_k$.
   - Because $r_{k+1} \perp r_k$, consecutive search steps are strictly at **$90^\circ$ right angles** to one another.
   - In a steep, elongated valley ($\kappa \gg 1$), the gradient points almost transversely across the valley rather than down along the bottom toward the minimum $x^*$.
   - The iterates bounce back and forth between the steep walls of the canyon in tiny, perpendicular zig-zag steps, resulting in asymptotic stagnation.

---

## Module 5: Krylov Subspace Methods (CG, Arnoldi, GMRES, BiCGSTAB)

### Q5.1: Conjugate Gradient Step Length $\alpha_k$ and Search Direction $\beta_k$ Derivation [7 Marks]
**Problem Statement**:  
In the Conjugate Gradient (CG) method for an SPD system $A x = b$, the search directions $p_k$ are constructed to be $A$-conjugate ($p_i^T A p_j = 0$ for $i \ne j$).
1. Starting from the iterate update $x_{k+1} = x_k + \alpha_k p_k$, derive $\alpha_k$ enforcing that $r_{k+1}$ is orthogonal to $p_k$.
2. In the direction update $p_k = r_k + \beta_{k-1} p_{k-1}$, derive the Fletcher-Reeves formula for $\beta_{k-1}$ enforcing that $p_k^T A p_{k-1} = 0$.
3. Show that this reduces to:
   $$\beta_{k-1} = \frac{r_k^T r_k}{r_{k-1}^T r_{k-1}}$$

#### Full Step-by-Step Solution:

1. **Derivation of Step Size $\alpha_k$**:
   The residual update is:
   $$r_{k+1} = b - A x_{k+1} = b - A(x_k + \alpha_k p_k) = r_k - \alpha_k A p_k$$
   We require that $r_{k+1} \perp p_k$ (projection condition):
   $$p_k^T r_{k+1} = 0 \implies p_k^T (r_k - \alpha_k A p_k) = 0$$
   $$p_k^T r_k - \alpha_k p_k^T A p_k = 0 \implies \mathbf{\alpha_k = \frac{p_k^T r_k}{p_k^T A p_k}}$$
   Since $p_k = r_k + \beta_{k-1} p_{k-1}$ and $r_k \perp p_{k-1}$, $p_k^T r_k = (r_k + \beta_{k-1} p_{k-1})^T r_k = r_k^T r_k$.  
   Therefore:
   $$\mathbf{\alpha_k = \frac{r_k^T r_k}{p_k^T A p_k}}$$

2. **Derivation of Conjugacy Parameter $\beta_{k-1}$**:
   The search direction update is:
   $$p_k = r_k + \beta_{k-1} p_{k-1}$$
   We require $A$-conjugacy between consecutive search directions:
   $$p_{k-1}^T A p_k = 0$$
   Substitute $p_k$:
   $$p_{k-1}^T A (r_k + \beta_{k-1} p_{k-1}) = 0$$
   $$p_{k-1}^T A r_k + \beta_{k-1} p_{k-1}^T A p_{k-1} = 0$$
   Solving for $\beta_{k-1}$:
   $$\mathbf{\beta_{k-1} = -\frac{p_{k-1}^T A r_k}{p_{k-1}^T A p_{k-1}} = -\frac{r_k^T A p_{k-1}}{p_{k-1}^T A p_{k-1}}}$$

3. **Simplification to Fletcher-Reeves Formula**:
   From the residual recurrence at step $k-1$:
   $$r_k = r_{k-1} - \alpha_{k-1} A p_{k-1} \implies A p_{k-1} = \frac{1}{\alpha_{k-1}} (r_{k-1} - r_k)$$
   Substitute $A p_{k-1}$ into the numerator of $\beta_{k-1}$:
   $$r_k^T A p_{k-1} = r_k^T \left[ \frac{1}{\alpha_{k-1}} (r_{k-1} - r_k) \right] = \frac{1}{\alpha_{k-1}} (r_k^T r_{k-1} - r_k^T r_k)$$
   By the orthogonality of the Krylov subspace residuals, $r_k^T r_{k-1} = 0$. Thus:
   $$r_k^T A p_{k-1} = -\frac{1}{\alpha_{k-1}} r_k^T r_k$$
   From step 1, $\alpha_{k-1} = \frac{r_{k-1}^T r_{k-1}}{p_{k-1}^T A p_{k-1}} \implies p_{k-1}^T A p_{k-1} = \frac{r_{k-1}^T r_{k-1}}{\alpha_{k-1}}$.  
   Substitute both into the expression for $\beta_{k-1}$:
   $$\beta_{k-1} = -\frac{-\frac{1}{\alpha_{k-1}} r_k^T r_k}{\frac{1}{\alpha_{k-1}} r_{k-1}^T r_{k-1}} = \mathbf{\frac{r_k^T r_k}{r_{k-1}^T r_{k-1}}}$$
   This completes the proof. $\blacksquare$

---

### Q5.2: Finite-Element Mesh Refinement Scaling on CG Iterations & Computational Cost [6 Marks]
**Problem Statement**:  
A 2D structural mechanics problem is solved using linear triangular finite elements on an initial mesh of size $h_1 = 0.02$ ($N_1 = 2,500$ nodes), requiring $T_1 = 0.8\text{ seconds}$ to converge with the Conjugate Gradient method.
1. The mesh is refined uniformly by a factor of 4 ($h_2 = h_1 / 4 = 0.005$). Determine the new number of degrees of freedom $N_2$.
2. Knowing that the global stiffness matrix condition number scales as $\kappa([K]) = O(h^{-2})$, by what theoretical factor do the CG iterations increase?
3. In a benchmark on a dual-channel workstation, the solver runtime on mesh 2 is observed to be $T_2 = 28.5\text{ seconds}$ (a $35.6\times$ increase), far exceeding the theoretical iteration scaling. Explain the two architectural hardware factors responsible for this discrepancy.

#### Full Step-by-Step Solution:

1. **Node Count Scaling**:
   In 2D, the number of nodes scales inversely with the square of the grid spacing $h$:
   $$N \propto \frac{1}{h^2} \implies N_2 = N_1 \times \left(\frac{h_1}{h_2}\right)^2 = 2,500 \times (4)^2 = 2,500 \times 16 = \mathbf{40,000\text{ nodes}}$$

2. **Iteration Scaling via Chebyshev Bound**:
   The number of iterations $k$ for Conjugate Gradient is bounded by:
   $$k \propto \sqrt{\kappa([K])}$$
   Since $\kappa([K]) = C \cdot h^{-2}$:
   $$\sqrt{\kappa([K])} \propto \sqrt{h^{-2}} = h^{-1} = \frac{1}{h}$$
   Refining the mesh by a factor of 4 ($h_2 = h_1 / 4$) implies:
   $$\frac{k_2}{k_1} \approx \frac{h_1}{h_2} = \mathbf{4\times}$$
   The number of CG iterations scales by a factor of **4**.

3. **Theoretical Work Scaling vs. Architectural Reality**:
   - **Theoretical FLOP Scaling**:
     Work per iteration is dominated by SpMV, which is $O(N)$ for sparse matrices.
     $$\text{Total Work} W \propto k \times N \propto h^{-1} \times h^{-2} = O(h^{-3}) = O(N^{1.5})$$
     Predicted runtime increase:
     $$\frac{T_2}{T_1} \approx 4 \times 16 = \mathbf{64\times \text{ for FLOPs}} \dots \text{ wait!}$$
     Let's check: $k$ scales by 4, $N$ scales by 16. Total operations scale by $4 \times 16 = 64\times$.  
     Why did $T_2$ increase by $35.6\times$?
   - **Hardware Architectural Factors**:
     1. **Cache Capacity Boundary Crossing**:
        For $N_1 = 2,500$, the sparse matrix and vectors occupy:
        $$2,500 \times 7 \times 12\text{ bytes} \approx 210\text{ KB}$$
        This fits entirely inside the ultra-fast L2 cache ($512\text{ KB}$), achieving near-peak compute throughput.  
        For $N_2 = 40,000$, working memory exceeds $3.36\text{ MB}$, overflowing the L2 cache and spilling into the slower L3 cache and DRAM, dropping memory bandwidth from $200\text{ GB/s}$ to $40\text{ GB/s}$.
     2. **Memory Bandwidth Saturation (The Memory Wall)**:
        SpMV has an extremely low arithmetic intensity ($\approx 0.25\text{ FLOP/byte}$). As the problem size expands beyond cache limits, the CPU cores stall waiting for DRAM word requests, shifting the execution from compute-bound to memory-bandwidth-bound.

---

### Q5.3: Arnoldi Hessenberg Reduction vs. Symmetric Lanczos 3-Term Recurrence [5 Marks]
**Problem Statement**:  
1. State the $m$-step Arnoldi procedure using Modified Gram-Schmidt and write the fundamental matrix relation:
   $$A V_m = V_{m+1} \bar{H}_m$$
2. Prove that if $A$ is real symmetric ($A = A^T$), the upper Hessenberg matrix $H_m = V_m^T A V_m$ is **tridiagonal**.
3. State the resulting 3-term Lanczos recurrence relation and explain why this dramatically reduces memory requirements compared to Arnoldi.

#### Full Step-by-Step Solution:

1. **The Arnoldi Process and Fundamental Relation**:
   Given matrix $A \in \mathbb{R}^{n \times n}$ and initial unit vector $v_1 = r_0 / \|r_0\|_2$:
   For $j = 1, \dots, m$:
   - $w = A v_j$
   - For $i = 1, \dots, j$:
     - $h_{ij} = v_i^T w$
     - $w = w - h_{ij} v_i$
   - $h_{j+1, j} = \|w\|_2$
   - If $h_{j+1, j} = 0$, stop; otherwise $v_{j+1} = w / h_{j+1, j}$.
   In matrix form, defining $V_m = [v_1, \dots, v_m] \in \mathbb{R}^{n \times m}$ and the $(m+1) \times m$ upper Hessenberg matrix $\bar{H}_m$:
   $$\mathbf{A V_m = V_{m+1} \bar{H}_m = V_m H_m + h_{m+1, m} v_{m+1} e_m^T}$$

2. **Proof that $H_m$ is Tridiagonal for Symmetric $A$**:
   Since the columns of $V_{m+1}$ are orthonormal ($V_m^T V_m = I_m$):
   Multiply $A V_m = V_{m+1} \bar{H}_m$ from the left by $V_m^T$:
   $$V_m^T A V_m = V_m^T (V_m H_m + h_{m+1, m} v_{m+1} e_m^T) = I_m H_m + \mathbf{0} = H_m$$
   Now take the matrix transpose of $H_m$:
   $$H_m^T = (V_m^T A V_m)^T = V_m^T A^T (V_m^T)^T = V_m^T A^T V_m$$
   Since $A = A^T$:
   $$H_m^T = V_m^T A V_m = H_m$$
   Thus, $H_m$ is **symmetric**!  
   By definition of the Arnoldi process, $H_m$ is an **upper Hessenberg matrix** ($h_{ij} = 0$ for $i > j + 1$).  
   A matrix that is both symmetric ($h_{ij} = h_{ji}$) and upper Hessenberg must have:
   $$h_{ij} = 0 \quad \text{for all } |i - j| > 1$$
   Therefore, $H_m$ is strictly **tridiagonal**:
   $$H_m = T_m = \begin{bmatrix}
   \alpha_1 & \beta_1 & 0 & \dots & 0 \\
   \beta_1 & \alpha_2 & \beta_2 & \dots & 0 \\
   0 & \beta_2 & \alpha_3 & \dots & 0 \\
   \vdots & \ddots & \ddots & \ddots & \vdots \\
   0 & \dots & 0 & \beta_{m-1} & \alpha_m
   \end{bmatrix} \quad \blacksquare$$

3. **The 3-Term Lanczos Recurrence and Memory Savings**:
   Because $h_{ij} = 0$ for $i < j - 1$, the vector $A v_j$ only needs to be orthogonalized against the **immediately preceding two vectors** $v_j$ and $v_{j-1}$:
   $$\mathbf{\beta_j v_{j+1} = A v_j - \alpha_j v_j - \beta_{j-1} v_{j-1}}$$
   - **Arnoldi Storage**: Must store all $m$ orthonormal vectors $\{v_1, \dots, v_m\}$ in RAM to orthogonalize $v_{j+1}$ ($m \times n$ floats $\implies$ memory explodes as $m$ grows).
   - **Lanczos Storage**: Only requires storing **3 vectors** ($v_{j+1}, v_j, v_{j-1}$), reducing memory footprint from $O(m n)$ to $O(n)$!

---

### Q5.4: GMRES Minimization via Givens Rotations & Residual Norm Monitoring [6 Marks]
**Problem Statement**:  
The Generalized Minimal Residual (GMRES) method finds $x_m \in x_0 + \mathcal{K}_m(A, r_0)$ that minimizes the residual norm $\|b - A x_m\|_2$.
1. Show that this least-squares problem reduces to:
   $$\min_{y \in \mathbb{R}^m} \|\beta e_1 - \bar{H}_m y\|_2$$
   where $\beta = \|r_0\|_2$ and $\bar{H}_m \in \mathbb{R}^{(m+1) \times m}$.
2. Describe how a sequence of $m$ Givens plane rotations $\Omega_1, \Omega_2, \dots, \Omega_m$ transforms $\bar{H}_m$ into an upper triangular matrix $R_m$.
3. Explain how the residual norm $\|r_m\|_2$ can be monitored at every step **without** computing the solution vector $x_m$ or performing matrix-vector multiplication.

#### Full Step-by-Step Solution:

1. **Reduction of GMRES Least-Squares Problem**:
   Any iterate $x_m \in x_0 + \mathcal{K}_m$ can be written as $x_m = x_0 + V_m y$ for some $y \in \mathbb{R}^m$.
   The residual is:
   $$r_m = b - A x_m = b - A(x_0 + V_m y) = r_0 - A V_m y$$
   Using the Arnoldi relation $A V_m = V_{m+1} \bar{H}_m$ and setting $r_0 = \beta v_1 = V_{m+1} (\beta e_1)$:
   $$r_m = V_{m+1} (\beta e_1) - V_{m+1} \bar{H}_m y = V_{m+1} (\beta e_1 - \bar{H}_m y)$$
   Since the columns of $V_{m+1}$ are orthonormal, the Euclidean norm is preserved:
   $$\|r_m\|_2 = \|V_{m+1} (\beta e_1 - \bar{H}_m y)\|_2 = \|\beta e_1 - \bar{H}_m y\|_2$$
   Thus, minimizing $\|r_m\|_2$ over $\mathbb{R}^n$ reduces to:
   $$\mathbf{\min_{y \in \mathbb{R}^m} \|\beta e_1 - \bar{H}_m y\|_2}$$

2. **Givens Rotations QR Transformation**:
   At step $j$, the Hessenberg matrix has non-zeros on the sub-diagonal entry $h_{j+1, j}$.  
   A $2 \times 2$ Givens rotation matrix acting on rows $j$ and $j+1$ is constructed:
   $$\Omega_j = \begin{bmatrix}
   c_j & s_j \\
   -s_j & c_j
   \end{bmatrix}, \quad \text{where } c_j = \frac{h_{jj}}{\sqrt{h_{jj}^2 + h_{j+1, j}^2}}, \; s_j = \frac{h_{j+1, j}}{\sqrt{h_{jj}^2 + h_{j+1, j}^2}}$$
   Multiplying $\bar{H}_m$ by $Q_m = \Omega_m \dots \Omega_2 \Omega_1$ eliminates all sub-diagonal elements:
   $$Q_m \bar{H}_m = \begin{bmatrix} R_m \\ \mathbf{0}^T \end{bmatrix}$$
   where $R_m \in \mathbb{R}^{m \times m}$ is upper triangular.

3. **In-Progress Residual Norm Monitoring**:
   Applying the same orthogonal rotations to the right-hand side vector $\beta e_1$:
   $$Q_m (\beta e_1) = \begin{bmatrix} g_m \\ \gamma_{m+1} \end{bmatrix}, \quad g_m \in \mathbb{R}^m, \; \gamma_{m+1} \in \mathbb{R}$$
   The transformed least-squares problem is:
   $$\|\beta e_1 - \bar{H}_m y\|_2 = \left\| \begin{bmatrix} g_m \\ \gamma_{m+1} \end{bmatrix} - \begin{bmatrix} R_m \\ \mathbf{0}^T \end{bmatrix} y \right\|_2 = \|g_m - R_m y\|_2^2 + |\gamma_{m+1}|^2$$
   The minimum is achieved when $R_m y = g_m$, which drives the first term to zero.  
   The remaining residual norm is **identically equal to the scalar entry $|\gamma_{m+1}|$**:
   $$\mathbf{\|r_m\|_2 = |\gamma_{m+1}|}$$
   - **Significance**: The solver checks whether $|\gamma_{m+1}| < \epsilon$ at negligible cost ($O(1)$ scalar operation). The expensive backward substitution $R_m y = g_m$ and vector reconstruction $x_m = x_0 + V_m y$ ($O(m n)$ operations) are performed **only once** when convergence is achieved!

---

### Q5.5: BiCG Breakdown Modes and BiCGSTAB Polynomial Stabilization Mechanism [5 Marks]
**Problem Statement**:  
1. Explain the difference between **Lanczos Biorthogonalization** (BiCG) and standard Arnoldi/Lanczos methods.
2. Define the two distinct breakdown modes that cause BiCG to fail:
   - (a) Division by zero in vector normalization (Pivot Breakdown).
   - (b) Lanczos breakdown (Serious Breakdown).
3. How does the **BiCGSTAB** (Biconjugate Gradient Stabilized) algorithm eliminate the erratic residual oscillations of BiCG and avoid using $A^T$?

#### Full Step-by-Step Solution:

1. **Lanczos Biorthogonalization Principle**:
   For non-symmetric matrices ($A \ne A^T$), standard Lanczos fails because $H_m$ is not tridiagonal.  
   BiCG resolves this by maintaining **two dual Krylov subspaces**:
   $$\mathcal{K}_m(A, v_1) = \text{span}\{v_1, A v_1, \dots, A^{m-1} v_1\}$$
   $$\mathcal{K}_m(A^T, w_1) = \text{span}\{w_1, A^T w_1, \dots, (A^T)^{m-1} w_1\}$$
   The basis vectors satisfy the **biorthogonality condition**:
   $$w_i^T v_j = \delta_{ij}$$
   This restores short 3-term recurrences, but projections are oblique rather than orthogonal.

2. **The Two Breakdown Modes of BiCG**:
   - **(a) Division by Zero (Pivot breakdown)**:
     In the scalar update $\alpha_j = \frac{\langle r_j, \tilde{r}_j \rangle}{\langle A p_j, \tilde{p}_j \rangle}$, the denominator $\langle A p_j, \tilde{p}_j \rangle = 0$ while $r_j \ne 0$. The step size becomes infinite.
   - **(b) Serious Lanczos Breakdown**:
     The dual inner product vanishes:
     $$\langle r_j, \tilde{r}_j \rangle = 0 \quad \text{while } r_j \ne 0 \text{ and } \tilde{r}_j \ne 0$$
     Because the search space and shadow space become mutually orthogonal, biorthogonal vectors $v_{j+1}$ and $w_{j+1}$ cannot be normalized, and the algorithm halts completely.

3. **BiCGSTAB Stabilization Mechanism**:
   - **Elimination of $A^T$**: Conjugate Gradient Squared (CGS) squared the BiCG residual polynomial $r_j^{\text{CGS}} = [\phi_j(A)]^2 r_0$, which removed $A^T$ but amplified round-off errors and caused violent residual spikes.
   - **BiCGSTAB Smoothing**:
     BiCGSTAB decomposes the residual polynomial into the product of the BiCG polynomial $\phi_j(t)$ and a **stabilizing linear factor**:
     $$r_j = \psi_j(A) \phi_j(A) r_0, \quad \text{where } \psi_j(t) = (1 - \omega_1 t)(1 - \omega_2 t) \dots (1 - \omega_j t)$$
   - At each half-step, BiCGSTAB performs a **1D Steepest Descent / Minimal Residual step** to choose $\omega_j$:
     $$\omega_j = \frac{\langle A s_j, s_j \rangle}{\langle A s_j, A s_j \rangle}$$
     This local minimization exponentially dampens oscillations, prevents wild residual spikes, and guarantees smooth, monotonic-like convergence without ever forming $A^T$.

---

## Module 6: Preconditioning Techniques

### Q6.1: Preconditioned Conjugate Gradient (PCG) Algebraic Derivation [6 Marks]
**Problem Statement**:  
Let $A \in \mathbb{R}^{n \times n}$ be SPD and $M \in \mathbb{R}^{n \times n}$ be an SPD preconditioner factored as $M = L L^T$.
1. Transform $A x = b$ into the split preconditioned symmetric system:
   $$\tilde{A} \tilde{x} = \tilde{b}, \quad \text{where } \tilde{A} = L^{-1} A L^{-T}$$
2. Apply standard Conjugate Gradient to $\tilde{A} \tilde{x} = \tilde{b}$.
3. Show that by substituting variables $x = L^{-T} \tilde{x}$, $p = L^{-T} \tilde{p}$, and defining the preconditioned residual solve:
   $$M z_k = r_k$$
   the explicit Cholesky factors $L$ and $L^T$ completely disappear from the algorithm!

#### Full Step-by-Step Solution:

1. **Split Symmetrization**:
   Given $A x = b$ and $M = L L^T$:
   $$L^{-1} A (L^{-T} L^T) x = L^{-1} b \implies (L^{-1} A L^{-T}) (L^T x) = L^{-1} b$$
   Define:
   $$\tilde{A} = L^{-1} A L^{-T}, \quad \tilde{x} = L^T x, \quad \tilde{b} = L^{-1} b$$
   Since $A$ is SPD, $\tilde{A}^T = (L^{-1} A L^{-T})^T = L^{-1} A^T L^{-T} = \tilde{A}$, and $y^T \tilde{A} y = (L^{-T} y)^T A (L^{-T} y) > 0$. Thus $\tilde{A}$ is SPD.

2. **Standard CG on Split System**:
   $$\tilde{\alpha}_k = \frac{\tilde{r}_k^T \tilde{r}_k}{\tilde{p}_k^T \tilde{A} \tilde{p}_k}$$
   $$\tilde{x}_{k+1} = \tilde{x}_k + \tilde{\alpha}_k \tilde{p}_k$$
   $$\tilde{r}_{k+1} = \tilde{r}_k - \tilde{\alpha}_k \tilde{A} \tilde{p}_k$$
   $$\tilde{\beta}_k = \frac{\tilde{r}_{k+1}^T \tilde{r}_{k+1}}{\tilde{r}_k^T \tilde{r}_k}$$
   $$\tilde{p}_{k+1} = \tilde{r}_{k+1} + \tilde{\beta}_k \tilde{p}_k$$

3. **Elimination of Cholesky Factors**:
   - Relationship between residuals:
     $$\tilde{r}_k = \tilde{b} - \tilde{A} \tilde{x}_k = L^{-1} b - L^{-1} A L^{-T} (L^T x_k) = L^{-1} (b - A x_k) = L^{-1} r_k$$
   - Numerator of $\tilde{\alpha}_k$:
     $$\tilde{r}_k^T \tilde{r}_k = (L^{-1} r_k)^T (L^{-1} r_k) = r_k^T (L L^T)^{-1} r_k = r_k^T M^{-1} r_k$$
     Define $\mathbf{z_k \equiv M^{-1} r_k} \iff \mathbf{M z_k = r_k}$. Then:
     $$\tilde{r}_k^T \tilde{r}_k = \mathbf{r_k^T z_k}$$
   - Denominator of $\tilde{\alpha}_k$:
     Let $p_k \equiv L^{-T} \tilde{p}_k$. Then:
     $$\tilde{p}_k^T \tilde{A} \tilde{p}_k = (L^T p_k)^T (L^{-1} A L^{-T}) (L^T p_k) = p_k^T A p_k$$
     Therefore:
     $$\mathbf{\alpha_k = \frac{r_k^T z_k}{p_k^T A p_k}}$$
   - Update of $x_{k+1}$:
     $$L^T x_{k+1} = L^T x_k + \alpha_k (L^T p_k) \implies \mathbf{x_{k+1} = x_k + \alpha_k p_k}$$
   - Update of $r_{k+1}$:
     $$L^{-1} r_{k+1} = L^{-1} r_k - \alpha_k L^{-1} A p_k \implies \mathbf{r_{k+1} = r_k - \alpha_k A p_k}$$
   - Update of $\beta_k$:
     $$\mathbf{\beta_k = \frac{r_{k+1}^T z_{k+1}}{r_k^T z_k}}$$
   - Direction update:
     $$L^T p_{k+1} = L^{-1} r_{k+1} + \beta_k L^T p_k \implies p_{k+1} = (L L^T)^{-1} r_{k+1} + \beta_k p_k \implies \mathbf{p_{k+1} = z_{k+1} + \beta_k p_k}$$
   - **Conclusion**: The algorithm requires only the solution of the linear system $M z = r$. The explicit matrix factor $L$ is never needed! $\blacksquare$

---

### Q6.2: Incomplete Cholesky IC(0) vs. SSOR Preconditioner Construction [5 Marks]
**Problem Statement**:  
1. Define the **Incomplete Cholesky factorization with zero fill-in (IC(0))** for an SPD sparse matrix $A$.
2. Write the matrix formulation of the **Symmetric Successive Over-Relaxation (SSOR)** preconditioner $M_{\text{SSOR}}$.
3. Compare IC(0) and SSOR in terms of:
   - (a) Setup phase CPU cost.
   - (b) Storage overhead beyond matrix $A$.
   - (c) Parallelizability of the forward/backward solve $M z = r$.

#### Full Step-by-Step Solution:

1. **Incomplete Cholesky IC(0)**:
   - Standard Cholesky computes $A = L L^T$, where $L$ suffers from severe fill-in ($l_{ij} \ne 0$ where $a_{ij} = 0$).
   - IC(0) computes an approximate factor $\tilde{L}$ such that:
     $$A = \tilde{L} \tilde{L}^T - R$$
     enforcing the zero-fill condition:
     $$\tilde{l}_{ij} = 0 \quad \text{whenever } a_{ij} = 0$$
   - The non-zero pattern of $\tilde{L}$ is **identical** to the lower triangular pattern of $A$.

2. **SSOR Preconditioner Formulation**:
   Splitting $A = D - E - E^T$ (where $D$ is diagonal, $-E$ is strictly lower triangular):
   A forward Gauss-Seidel step followed by an adjoint backward Gauss-Seidel step yields:
   $$\mathbf{M_{\text{SSOR}} = \frac{1}{\omega(2 - \omega)} (D - \omega E) D^{-1} (D - \omega E^T)}$$
   For $\omega = 1$ (Symmetric Gauss-Seidel / SGS):
   $$\mathbf{M_{\text{SGS}} = (D - E) D^{-1} (D - E^T)}$$

3. **Comparative Evaluation**:

| Evaluation Metric | Incomplete Cholesky (IC(0)) | SSOR Preconditioner ($\omega = 1$) |
| :--- | :--- | :--- |
| **Setup Cost** | Requires explicit incomplete factorization loop ($O(N_{nz})$ FLOPs) with square roots. Risk of breakdown if pivot becomes non-positive. | **Zero setup cost**. Factors $(D - E)$ and $D^{-1}$ are read directly from existing arrays of $A$. |
| **Memory Overhead** | Requires storing $\tilde{L}$ ($N_{nz}/2$ extra floats if not done in-place). | **Zero memory overhead**. Operates directly on the elements of $A$. |
| **Solve Complexity** | Triangular solve $\tilde{L} y = r$, $\tilde{L}^T z = y$. Sequential dependency along elimination tree. | Triangular solve $(D-E)y = r$, $(D-E^T)z = Dy$. Sequential across rows. |
| **Parallelizability** | Low (requires multicoloring or level scheduling). | Low (requires Red-Black or hyperplane reordering). |

---

## Module 7: Computer Architecture, Caching & Memory Wall

### Q7.1: The Memory Wall: DRAM Latency Penalty and Effective Access Time (EAT) [5 Marks]
**Problem Statement**:  
A high-performance compute core runs at clock frequency $f = 2.5\text{ GHz}$ and has a peak execution rate of 4 floating-point operations per cycle ($10\text{ GFLOPS}$).  
The processor has a two-level cache hierarchy with the following parameters:
- L1 Cache: Hit rate $H_1 = 0.95$, latency $t_1 = 1\text{ cycle}$ ($0.4\text{ ns}$).
- L2 Cache: Hit rate $H_2 = 0.80$ (of L1 misses), latency $t_2 = 10\text{ cycles}$ ($4.0\text{ ns}$).
- Main Memory (DRAM): Latency $t_{\text{mem}} = 100\text{ ns}$ ($250\text{ cycles}$).
1. Calculate the Effective Access Time (EAT) in nanoseconds and processor clock cycles.
2. If an uncached sparse solver accesses DRAM on every memory reference, calculate its effective processing speed. By what factor is the CPU throttled by DRAM latency?

#### Full Step-by-Step Solution:

1. **Calculation of Effective Access Time (EAT)**:
   The formula for multi-level hierarchical EAT is:
   $$\text{EAT} = t_1 + (1 - H_1) \cdot \left[ t_2 + (1 - H_2) \cdot t_{\text{mem}} \right]$$
   - $t_1 = 0.4\text{ ns}$
   - Miss rate L1: $1 - H_1 = 1 - 0.95 = 0.05$
   - Miss rate L2: $1 - H_2 = 1 - 0.80 = 0.20$
   - L2 miss penalty: $(1 - H_2) \cdot t_{\text{mem}} = 0.20 \times 100\text{ ns} = 20.0\text{ ns}$
   - Combined L2 access: $t_2 + 20.0\text{ ns} = 4.0\text{ ns} + 20.0\text{ ns} = 24.0\text{ ns}$
   - EAT in nanoseconds:
     $$\text{EAT} = 0.4\text{ ns} + 0.05 \times 24.0\text{ ns} = 0.4\text{ ns} + 1.2\text{ ns} = \mathbf{1.6\text{ ns}}$$
   - EAT in clock cycles (at $0.4\text{ ns/cycle}$):
     $$\text{Cycles} = \frac{1.6\text{ ns}}{0.4\text{ ns/cycle}} = \mathbf{4\text{ clock cycles}}$$

2. **DRAM Throttling (Uncached Execution)**:
   - Cycle time: $\tau = \frac{1}{2.5\text{ GHz}} = 0.4\text{ ns}$.
   - Peak throughput: 4 FLOPs / $0.4\text{ ns} = \mathbf{10\text{ GFLOPS}}$.
   - If every memory load stalls for main DRAM latency ($t_{\text{mem}} = 100\text{ ns}$):
     $$\text{Effective Speed} = \frac{1\text{ operation}}{100\text{ ns}} = \mathbf{10\text{ MFLOPS}} = 0.01\text{ GFLOPS}$$
   - **Throttling Factor**:
     $$\text{Slowdown} = \frac{\text{Peak Speed}}{\text{Effective Speed}} = \frac{10\text{ GFLOPS}}{0.01\text{ GFLOPS}} = \mathbf{1000\times}$$
   - **Conclusion**: The processor sits idle, stalled on DRAM wait states, operating at a catastrophic **$0.1\%$ of its peak capability**. This is the physical reality of the **Memory Wall**.

---

### Q7.2: Row-Major vs. Column-Major Striding and Cache Tiling Analysis [6 Marks]
**Problem Statement**:  
1. Contrast Row-Major (C/C++) and Column-Major (Fortran/MATLAB) memory layouts for a 2D array `A[N][N]`.
2. A programmer writes the following nested loop in C on an array of size $N = 10,000$ (double precision, 8 bytes):
   ```c
   for (j = 0; j < N; j++) {
       for (i = 0; i < N; i++) {
           A[i][j] = 2.0 * A[i][j];
       }
   }
   ```
   Assuming a 64-byte cache line size and a cold cache, calculate the total number of cache misses. How many misses would occur if the loops were interchanged?
3. Formulate the **Cache Tiling (Loop Blocking)** transformation for dense matrix multiplication $C = A B$ and derive the optimal tile size $B$ for an L1 cache of $32\text{ KB}$.

#### Full Step-by-Step Solution:

1. **Row-Major vs. Column-Major Layout**:
   - **Row-Major (C/C++)**: Consecutive elements of the same row (`A[i][j]` and `A[i][j+1]`) are adjacent in physical memory (stride-1 in index $j$).
   - **Column-Major (Fortran)**: Consecutive elements of the same column (`A(i, j)` and `A(i+1, j)`) are adjacent in physical memory (stride-1 in index $i$).

2. **Cache Miss Analysis in C**:
   - Cache line size: $64\text{ bytes} / 8\text{ bytes/double} = \mathbf{8\text{ doubles per cache line}}$.
   - Total elements: $N^2 = (10^4)^2 = 10^8$ doubles ($800\text{ MB}$).
   - **Given Code (Inner loop over $i$, Column-wise access `A[i][j]`)**:
     - At each step $i \to i+1$, the address jumps by $N \times 8\text{ bytes} = 80,000\text{ bytes}$.
     - Since $80\text{ KB} > \text{L1 cache size}$, consecutive iterations access completely different cache lines!
     - Spatial locality is **completely destroyed**. Every single memory access causes a **cache miss**:
       $$\text{Misses}_{\text{bad}} = N \times N = 10^4 \times 10^4 = \mathbf{10^8\text{ cache misses}}$$
   - **Interchanged Code (Inner loop over $j$, Row-wise access `A[i][j]`)**:
     - At each step $j \to j+1$, the address increments by 8 bytes (consecutive memory).
     - The first access `A[i][0]` loads a 64-byte line containing 8 doubles (`A[i][0] ... A[i][7]`).
     - Accesses 1 through 7 are guaranteed **cache hits** (spatial locality).
     - Cache miss occurs only once every 8 elements:
       $$\text{Misses}_{\text{optimal}} = \frac{N \times N}{8} = \frac{10^8}{8} = \mathbf{1.25 \times 10^7\text{ cache misses}}$$
     - **Performance Impact**: Loop interchange achieves an **$8\times$ reduction in memory traffic**!

3. **Cache Tiling (Loop Blocking) for $C = A B$**:
   - Partition matrices into $B \times B$ submatrices (tiles):
     ```c
     for (ii = 0; ii < N; ii += B)
       for (jj = 0; jj < N; jj += B)
         for (kk = 0; kk < N; kk += B)
           for (i = ii; i < ii + B; i++)
             for (j = jj; j < jj + B; j++)
               for (k = kk; k < kk + B; k++)
                 C[i][j] += A[i][k] * B[k][j];
     ```
   - **Optimal Tile Size Derivation**:
     To prevent cache thrashing, one tile of $A$, one tile of $B$, and one tile of $C$ must simultaneously reside in L1 cache:
     $$3 \times B^2 \times 8\text{ bytes} \le \text{Cache Size} = 32\text{ KB} = 32,768\text{ bytes}$$
     $$24 B^2 \le 32,768 \implies B^2 \le \frac{32,768}{24} \approx 1365.3 \implies B \le \sqrt{1365.3} \approx 36.9$$
     Selecting the nearest power of 2:
     $$\mathbf{B = 32}$$
     This reduces memory traffic between cache and DRAM by a factor of $B = 32\times$!

---

## Module 8: Parallel Architectures, Interconnects & Topologies

### Q8.1: Network Topology Metrics: Degree, Diameter, and Bisection Width Derivation [6 Marks]
**Problem Statement**:  
Complete the following comparative table for $p$ processors interconnected in various network topologies:
1. Linear Array ($p$ nodes).
2. 2D Mesh ($\sqrt{p} \times \sqrt{p}$ nodes, without wraparound).
3. 2D Torus ($\sqrt{p} \times \sqrt{p}$ nodes, with wraparound).
4. $d$-dimensional Hypercube ($p = 2^d$ nodes).
For each topology, evaluate: **Node Degree**, **Network Diameter**, and **Bisection Width**.

#### Full Step-by-Step Solution:

1. **Definitions of Network Metrics**:
   - **Node Degree ($d_{\text{node}}$)**: Maximum number of physical communication links connected to any single processor. Determines switch port complexity.
   - **Network Diameter ($D$)**: Maximum shortest path distance (number of hops) between any pair of nodes. Determines worst-case communication latency.
   - **Bisection Width ($B_w$)**: Minimum number of communication links that must be severed to partition the network into two equal halves of $p/2$ nodes. Determines bandwidth bottleneck for global communication (e.g., all-to-all).

2. **Derivations for Each Topology**:

   - **Linear Array ($p$ nodes)**:
     - End nodes have 1 link, interior nodes have 2 links $\implies \text{Degree} = \mathbf{2}$.
     - Longest path is from node 1 to node $p \implies \text{Diameter} = \mathbf{p - 1}$.
     - Splitting in the middle cuts exactly 1 link $\implies \text{Bisection Width} = \mathbf{1}$.

   - **2D Mesh ($\sqrt{p} \times \sqrt{p}$, no wraparound)**:
     - Interior nodes have 4 links (North, South, East, West) $\implies \text{Degree} = \mathbf{4}$.
     - Longest path is from corner $(1, 1)$ to opposite corner $(\sqrt{p}, \sqrt{p})$:
       $D = (\sqrt{p} - 1) + (\sqrt{p} - 1) = \mathbf{2(\sqrt{p} - 1)}$.
     - Cutting vertically down the center severs 1 link per row across $\sqrt{p}$ rows $\implies B_w = \mathbf{\sqrt{p}}$.

   - **2D Torus ($\sqrt{p} \times \sqrt{p}$, with wraparound)**:
     - Every node has exactly 4 links $\implies \text{Degree} = \mathbf{4}$.
     - Wraparound links halve the maximum distance in each dimension:
       $D = 2 \times \lfloor \sqrt{p} / 2 \rfloor = \mathbf{\sqrt{p}}$ (for even $\sqrt{p}$).
     - Bisecting cuts both the regular links and the wraparound links $\implies B_w = \mathbf{2\sqrt{p}}$.

   - **Hypercube ($p = 2^d$, $d = \log_2 p$)**:
     - Each node connects to $d$ neighbors whose binary representations differ by exactly 1 bit $\implies \text{Degree} = d = \mathbf{\log_2 p}$.
     - Maximum bit flips between any two binary addresses is $d \implies \text{Diameter} = d = \mathbf{\log_2 p}$.
     - Bisecting corresponds to fixing 1 bit position, cutting all links along that dimension:
       $B_w = \frac{p}{2} = \mathbf{2^{d-1}}$.

3. **Master Comparative Summary Table**:

| Network Topology | Node Degree ($d_{\text{node}}$) | Network Diameter ($D$) | Bisection Width ($B_w$) | Architectural Cost |
| :--- | :---: | :---: | :---: | :---: |
| **Linear Array** | $2$ | $p - 1$ | $1$ | $O(p)$ |
| **2D Mesh** | $4$ | $2(\sqrt{p} - 1)$ | $\sqrt{p}$ | $O(p)$ |
| **2D Torus** | $4$ | $\sqrt{p}$ | $2\sqrt{p}$ | $O(p)$ |
| **Hypercube ($d$-cube)** | $\log_2 p$ | $\log_2 p$ | $p/2$ | $O(p \log p)$ |

---

### Q8.2: Cache Coherency (MESI Protocol) and False Sharing Pathology with Code Fix [5 Marks]
**Problem Statement**:  
1. Explain the **Cache Coherency Problem** in Shared-Memory multiprocessors and define the 4 states of the **MESI Protocol** (Modified, Exclusive, Shared, Invalid).
2. The following OpenMP code computes the sum of elements in parallel:
   ```c
   int sum[NUM_THREADS];
   #pragma omp parallel
   {
       int id = omp_get_thread_num();
       sum[id] = 0;
       #pragma omp for
       for (int i = 0; i < N; i++) {
           sum[id] += a[i];
       }
   }
   ```
   Explain why this code suffers from **False Sharing** and exhibits terrible parallel performance on multicore CPUs.
3. Rewrite the code to completely eliminate false sharing.

#### Full Step-by-Step Solution:

1. **Cache Coherency and the MESI Protocol**:
   - In Symmetric Multiprocessors (SMP), each CPU core has its own private L1/L2 cache, while sharing central DRAM. If Core 0 modifies variable $X$ in its local cache without informing Core 1, Core 1 will read stale data from its cache or main memory.
   - **MESI Protocol States**:
     - **M (Modified)**: Line is present only in current cache and is dirty (modified relative to DRAM).
     - **E (Exclusive)**: Line is present only in current cache, but is clean (matches DRAM).
     - **S (Shared)**: Line is present in multiple caches, clean (read-only).
     - **I (Invalid)**: Line data is stale and invalid. Any read triggers a cache miss.

2. **The False Sharing Pathology**:
   - Cache coherency hardware tracks memory status at the granularity of **Cache Lines** (typically 64 bytes), *not individual variables*.
   - In the code above, `sum[NUM_THREADS]` is an array of 4-byte integers.
   - For `NUM_THREADS = 8`, all 8 integers occupy $8 \times 4\text{ B} = 32\text{ bytes}$, which fits **entirely within a single 64-byte cache line**!
   - When Thread 0 writes to `sum[0]`, its core marks the cache line as **Modified (M)** and broadcasts a **Bus Invalidate** signal across the interconnect.
   - The private cache lines of Cores 1 through 7 are immediately forced to the **Invalid (I)** state.
   - When Thread 1 tries to write to `sum[1]`, it experiences a cache miss, loads the line from Core 0, marks it Modified, and invalidates Core 0!
   - The cores spend virtually 100% of their execution time ping-ponging the cache line across the interconnect (**Cache Thrashing**), degrading performance below serial execution!

3. **Corrected Code Implementations**:

   - **Method A: Private Local Accumulator (Best Practice)**:
     ```c
     int total_sum = 0;
     #pragma omp parallel
     {
         int local_sum = 0; // Stored in thread-private CPU register
         #pragma omp for
         for (int i = 0; i < N; i++) {
             local_sum += a[i];
         }
         #pragma omp atomic
         total_sum += local_sum;
     }
     ```

   - **Method B: Structure Padding**:
     If an array must be maintained, pad each element to exceed the 64-byte cache line boundary:
     ```c
     struct ThreadData {
         int sum;
         char pad[60]; // 4 + 60 = 64 bytes! Guarantees independent cache lines.
     };
     struct ThreadData sum[NUM_THREADS];
     ```

---

## Module 9: Parallel Performance, Scalability & Amdahl's Law

### Q9.1: Amdahl's Law with Parallel Communication Overhead and Maximum Efficiency [6 Marks]
**Problem Statement**:  
A scientific simulation code has a sequential fraction $f = 0.05$ ($5\%$ inherently serial).  
When running the parallelizable portion ($95\%$) on $p$ processors, communication and synchronization overhead introduces an additional penalty given by:
$$T_{\text{comm}}(p) = 0.001 \times p\text{ seconds}$$
Normalized serial runtime is $T_s = 1.0\text{ second}$.
1. Formulate the total parallel runtime $T_p(p)$ as a function of processor count $p$.
2. Derive the analytical expression for the parallel speedup $S(p)$ and parallel efficiency $E(p)$.
3. Determine the optimal number of processors $p^*$ that maximizes parallel speedup.
4. Calculate the peak speedup $S(p^*)$ and the efficiency at this optimal point.

#### Full Step-by-Step Solution:

1. **Parallel Execution Time Formulation**:
   $$T_p(p) = T_{\text{serial}} + T_{\text{parallel}} + T_{\text{comm}}$$
   Given $T_s = 1.0\text{ s}$, $f = 0.05$:
   $$T_{\text{serial}} = f \cdot T_s = 0.05\text{ s}$$
   $$T_{\text{parallel}} = \frac{(1 - f) T_s}{p} = \frac{0.95}{p}\text{ s}$$
   $$T_{\text{comm}} = 0.001 \cdot p\text{ s}$$
   Therefore:
   $$\mathbf{T_p(p) = 0.05 + \frac{0.95}{p} + 0.001 p}$$

2. **Speedup and Efficiency Expressions**:
   $$\mathbf{S(p) = \frac{T_s}{T_p(p)} = \frac{1}{0.05 + \frac{0.95}{p} + 0.001 p}}$$
   $$\mathbf{E(p) = \frac{S(p)}{p} = \frac{1}{0.05 p + 0.95 + 0.001 p^2}}$$

3. **Optimal Processor Count $p^*$**:
   Maximizing speedup $S(p)$ is mathematically equivalent to minimizing runtime $T_p(p)$.  
   Differentiate $T_p(p)$ with respect to $p$ and set to zero:
   $$\frac{d T_p}{dp} = -\frac{0.95}{p^2} + 0.001 = 0$$
   $$\frac{0.95}{p^2} = 0.001 \implies p^2 = \frac{0.95}{0.001} = 950$$
   $$p^* = \sqrt{950} \approx \mathbf{30.82 \implies 31\text{ processors}}$$
   Verify with second derivative: $\frac{d^2 T_p}{dp^2} = \frac{1.9}{p^3} > 0$ (strict minimum).

4. **Peak Speedup and Efficiency at $p^* = 31$**:
   Substitute $p = 31$ into $T_p$:
   $$T_p(31) = 0.05 + \frac{0.95}{31} + 0.001(31) = 0.05 + 0.030645 + 0.031000 = \mathbf{0.111645\text{ s}}$$
   - **Peak Speedup**:
     $$\mathbf{S(31) = \frac{1.0}{0.111645} \approx \mathbf{8.957\times}}$$
   - **Parallel Efficiency**:
     $$\mathbf{E(31) = \frac{S(31)}{31} = \frac{8.957}{31} \approx \mathbf{28.89\%}}$$
   > **Exam Note**: Beyond 31 processors, adding more cores actually **slows down the simulation** ($S(p)$ decreases) because communication overhead grows faster than the computation reduction!

---

### Q9.2: Isoefficiency Function Derivation and Scalability Evaluation [6 Marks]
**Problem Statement**:  
1. Define **Scalability** of a parallel system and define the **Isoefficiency Metric** developed by Kumar and Grama.
2. An algorithm solves a problem of size $W$ on $p$ processors with an overhead function:
   $$T_o(W, p) = c_1 p \sqrt{W} + c_2 p^2$$
   Derive the isoefficiency function $W(p)$ to maintain a constant efficiency $E$.
3. If processor count $p$ is increased by a factor of 4, by what factor must the workload $W$ increase to maintain the same efficiency?
4. Is this algorithm considered poorly scalable, moderately scalable, or highly scalable?

#### Full Step-by-Step Solution:

1. **Definitions**:
   - **Scalability**: The capacity of a parallel algorithm-architecture combination to maintain a constant efficiency $E$ as the number of processors $p$ is increased by simultaneously increasing the problem size $W$.
   - **Isoefficiency Function $W(p)$**: The required growth rate of problem size $W$ as a function of $p$ to keep efficiency $E$ invariant. The fundamental relation is:
     $$W = K \cdot T_o(W, p), \quad \text{where } K = \frac{E}{1 - E}$$

2. **Derivation of Isoefficiency Function**:
   Given $T_o(W, p) = c_1 p \sqrt{W} + c_2 p^2$:
   $$W = K (c_1 p \sqrt{W} + c_2 p^2)$$
   We analyze the workload required to balance each overhead term independently:
   - **Term 1 Balance ($c_1 p \sqrt{W}$)**:
     $$W = K c_1 p \sqrt{W} \implies \sqrt{W} = K c_1 p \implies \mathbf{W \propto p^2}$$
   - **Term 2 Balance ($c_2 p^2$)**:
     $$W = K c_2 p^2 \implies \mathbf{W \propto p^2}$$
   Both terms yield the same asymptotic growth rate:
   $$\mathbf{W = O(p^2)}$$

3. **Workload Increase for $4\times$ Processors**:
   Let $p_{\text{new}} = 4 p$. Since $W \propto p^2$:
   $$W_{\text{new}} = (4 p)^2 = 16 p^2 = \mathbf{16 W}$$
   The problem size $W$ must increase by a factor of **16**!

4. **Scalability Classification**:
   - An algorithm is **optimally scalable** if $W = O(p)$ (linear growth).
   - An algorithm is **highly scalable** if $W = O(p \log p)$.
   - An algorithm with $W = O(p^2)$ is classified as **moderately scalable**. (It is substantially more scalable than matrix multiplication on rings ($O(p^3)$), but requires problem memory per node to grow linearly as $W/p = O(p)$).

---

## Module 10: Parallel Numerical Linear Algebra

### Q10.1: Parallel SpMV: 1D Row Striping vs. 2D Checkerboard Mesh Decomposition [6 Marks]
**Problem Statement**:  
Consider parallel sparse matrix-vector multiplication $y = A x$ on a distributed memory cluster of $p$ processors with network latency $t_s$ and per-word transfer time $t_w$. Matrix $A$ is $n \times n$ with $m$ non-zeros per row.
1. In a **1D Row-wise Block Partitioning**, describe the communication phase and write the parallel execution time $T_p^{1D}$.
2. In a **2D Checkerboard Partitioning** on a $\sqrt{p} \times \sqrt{p}$ mesh, describe the required communication steps (Column Broadcast & Row Reduction) and formulate $T_p^{2D}$.
3. For very large processor counts ($p \gg 1000$), explain why 2D decomposition is asymptotically superior.

#### Full Step-by-Step Solution:

1. **1D Row-wise Block Striping**:
   - Each of the $p$ processors owns $n/p$ consecutive rows of $A$ and the local vector segment $x_{\text{loc}}$ of length $n/p$.
   - **Communication Phase**: To compute $y_{\text{loc}} = A_{\text{loc}} x$, each processor needs non-local components of $x$. In general unstructured topologies, this requires an **All-to-All Broadcast (MPI_Allgather)** of the vector $x$ across all $p$ processors:
     $$T_{\text{comm}}^{1D} = t_s \log_2 p + t_w \left(\frac{p - 1}{p}\right) n \approx t_s \log_2 p + t_w n$$
   - **Computation Phase**: Each processor computes $n/p$ rows, each with $m$ non-zeros $\implies 2 m (n/p)$ FLOPs:
     $$T_{\text{comp}}^{1D} = 2 m \frac{n}{p} t_{\text{flop}}$$
   - **Total 1D Parallel Runtime**:
     $$\mathbf{T_p^{1D} = 2 m \frac{n}{p} t_{\text{flop}} + t_s \log_2 p + t_w n}$$

2. **2D Checkerboard Partitioning on $\sqrt{p} \times \sqrt{p}$ Grid**:
   - Processors are arranged in a 2D mesh $P_{i, j}$ for $i, j \in \{1, \dots, \sqrt{p}\}$.
   - Each processor holds a submatrix block of size $\frac{n}{\sqrt{p}} \times \frac{n}{\sqrt{p}}$.
   - **Step 1 (Column Broadcast)**:
     Vector segments $x_j$ (size $n/\sqrt{p}$) are broadcast vertically along processor columns ($\sqrt{p}$ processors per column):
     $$T_{\text{comm, 1}} = (t_s + t_w \frac{n}{\sqrt{p}}) \log_2 \sqrt{p} = \frac{1}{2} (t_s \log_2 p + t_w \frac{n}{\sqrt{p}} \log_2 p)$$
   - **Step 2 (Local Block SpMV)**:
     $$T_{\text{comp}}^{2D} = 2 m \frac{n}{p} t_{\text{flop}}$$
   - **Step 3 (Row Reduction)**:
     Partial results are summed horizontally along processor rows via **MPI_Reduce_scatter**:
     $$T_{\text{comm, 2}} = \frac{1}{2} (t_s \log_2 p + t_w \frac{n}{\sqrt{p}} \log_2 p)$$
   - **Total 2D Parallel Runtime**:
     $$\mathbf{T_p^{2D} = 2 m \frac{n}{p} t_{\text{flop}} + t_s \log_2 p + t_w \frac{n}{\sqrt{p}} \log_2 p}$$

3. **Asymptotic Superiority of 2D Decomposition**:
   - Compare the data transfer volume term ($t_w$):
     - In 1D striping: $t_w \cdot \mathbf{n}$ (Independent of $p$! Communication volume does **not decrease** as processors are added).
     - In 2D striping: $t_w \cdot \mathbf{\frac{n}{\sqrt{p}} \log_2 p}$ ($\to 0$ as $p \to \infty$).
   - For massively parallel supercomputers ($p \ge 10,000$), 1D row striping creates an insurmountable communication wall where all processors saturate the network links with full-length vector broadcasts. 2D decomposition restricts communication to small $\sqrt{p}$ sub-communicators, achieving superior weak and strong scaling.

---

### Q10.2: Cannon's Algorithm on 2D Torus: Communication Complexity & Shift Steps [7 Marks]
**Problem Statement**:  
Cannon's algorithm multiplies two dense $n \times n$ matrices $C = A B$ on a 2D Torus of $p$ processors arranged as a $\sqrt{p} \times \sqrt{p}$ mesh.
1. Describe the **Initial Preskewing Phase** for blocks $A_{i, j}$ and $B_{i, j}$.
2. Describe the **Multiply-Shift Loop** and state the exact number of stages.
3. Formulate the total parallel runtime $T_p$, showing computation and communication components.
4. Derive the **Isoefficiency Function** of Cannon's algorithm.

#### Full Step-by-Step Solution:

1. **Phase 1: Initial Preskewing Alignment**:
   To ensure that each processor $P_{i, j}$ initially holds matching blocks $A_{i, k}$ and $B_{k, j}$ with identical inner indices $k$:
   - Row $i$ of block matrix $A$ is circularly shifted **left by $i$ positions** ($i \in \{0, \dots, \sqrt{p}-1\}$).
     $$A_{i, j} \leftarrow A_{i, (j + i) \bmod \sqrt{p}}$$
   - Column $j$ of block matrix $B$ is circularly shifted **up by $j$ positions** ($j \in \{0, \dots, \sqrt{p}-1\}$).
     $$B_{i, j} \leftarrow B_{(i + j) \bmod \sqrt{p}, j}$$
   - **Preskewing Communication Cost**:
     Block size is $\frac{n}{\sqrt{p}} \times \frac{n}{\sqrt{p}}$, so message size is $m = \frac{n^2}{p}$ words.
     Maximum shift distance is $\sqrt{p}-1$:
     $$T_{\text{preskew}} = 2 \left( t_s + t_w \frac{n^2}{p} \right)$$

2. **Phase 2: Multiply-and-Shift Loop**:
   The main loop executes for exactly **$\sqrt{p}$ stages**:
   - For step $k = 1$ to $\sqrt{p}$:
     1. **Compute**: Local matrix multiply-add on block size $m_b = n / \sqrt{p}$:
        $$C_{i, j} = C_{i, j} + A_{i, j} B_{i, j}$$
        Requires $2 \left(\frac{n}{\sqrt{p}}\right)^3 = \frac{2 n^3}{p \sqrt{p}}$ FLOPs per stage.
     2. **Communicate (Shift)**:
        - Circularly shift $A$ **left by 1 block** along row.
        - Circularly shift $B$ **up by 1 block** along column.
        Each shift involves nearest-neighbor communication of size $n^2 / p$:
        $$T_{\text{shift-step}} = 2 \left( t_s + t_w \frac{n^2}{p} \right)$$

3. **Total Parallel Execution Time Formulation**:
   Summing over all $\sqrt{p}$ stages:
   - **Total Computation**:
     $$T_{\text{comp}} = \sqrt{p} \times \left( \frac{2 n^3}{p \sqrt{p}} \right) t_{\text{flop}} = \mathbf{\frac{2 n^3}{p} t_{\text{flop}}}$$
   - **Total Communication** (Preskewing + $\sqrt{p}$ shifts):
     $$T_{\text{comm}} \approx \sqrt{p} \times 2 \left( t_s + t_w \frac{n^2}{p} \right) = \mathbf{2 \sqrt{p} \, t_s + 2 \frac{n^2}{\sqrt{p}} \, t_w}$$
   - **Master Runtime Equation**:
     $$\mathbf{T_p = \frac{2 n^3}{p} t_{\text{flop}} + 2 \sqrt{p} \, t_s + 2 \frac{n^2}{\sqrt{p}} \, t_w}$$

4. **Derivation of Isoefficiency Function**:
   - Serial work: $W = 2 n^3 \implies n = (W/2)^{1/3}$.
   - Total overhead:
     $$T_o = p T_p - W = 2 p \sqrt{p} \, t_s + 2 n^2 \sqrt{p} \, t_w$$
   - Balancing $W = K \cdot T_o$:
     - For the latency term ($2 p^{1.5} t_s$):
       $$W \propto p^{1.5}$$
     - For the bandwidth term ($2 n^2 \sqrt{p} t_w$):
       Substitute $n^2 = (W/2)^{2/3}$:
       $$W = K' W^{2/3} \sqrt{p} \implies W^{1/3} \propto \sqrt{p} \implies W \propto (\sqrt{p})^3 = p^{1.5}$$
   - Both terms scale identically! Therefore, the isoefficiency function of Cannon's algorithm is:
     $$\mathbf{W = O(p^{1.5}) \quad \left(\text{or } W \propto p^{3/2}\right)}$$
   - **Scalability Significance**: To keep parallel efficiency constant when processors $p$ quadruple ($4\times$), the workload $W$ only needs to increase by $4^{1.5} = 8\times$. This makes Cannon's algorithm **exceptionally scalable** for dense linear algebra!

---

## Full Model Examination Paper (25 Marks, 2 Hours)

> **Instructions**: Answer all questions. You may use class lecture notes and reference textbooks. Calculators are permitted. Part A consists of analytical theory and derivations; Part B involves sparse array reconstruction and algorithmic tracing.

---

### Part A (15 Marks)

#### Question 1 [3 Marks]
A linear stationary iterative method $x^{(k+1)} = G x^{(k)} + f$ has the iteration matrix:
$$G = \begin{bmatrix}
0.5 & 0 & 0 & 0 \\
-1.2 & -0.4 & 0 & 0 \\
3.0 & 2.0 & 0.8 & 1.2 \\
-1.0 & 0 & 0.4 & 0.6
\end{bmatrix}$$
1. Find all eigenvalues of $G$ and determine the spectral radius $\rho(G)$.
2. Does the iteration converge? Justify your answer.

**Solution**:
1. Partition $G$ into block lower triangular form:
   $$G_{11} = \begin{bmatrix} 0.5 & 0 \\ -1.2 & -0.4 \end{bmatrix}, \quad G_{22} = \begin{bmatrix} 0.8 & 1.2 \\ 0.4 & 0.6 \end{bmatrix}$$
   - Eigenvalues of $G_{11}$: $\lambda_1 = 0.5, \; \lambda_2 = -0.4$.
   - Eigenvalues of $G_{22}$:
     $$\det(G_{22} - \lambda I) = (0.8 - \lambda)(0.6 - \lambda) - (1.2)(0.4) = 0.48 - 1.4\lambda + \lambda^2 - 0.48 = 0$$
     $$\lambda^2 - 1.4\lambda = 0 \implies \lambda(\lambda - 1.4) = 0 \implies \lambda_3 = 0, \; \lambda_4 = 1.4$$
2. The spectral radius is $\mathbf{\rho(G) = \max(|0.5|, |-0.4|, |0|, |1.4|) = 1.4}$.  
   Since $\mathbf{\rho(G) = 1.4 > 1}$, the scheme **diverges unconditionally**.

---

#### Question 2 [3 Marks]
In the Conjugate Gradient method for an SPD matrix $A$:
1. State the Chebyshev polynomial error bound for $\|e_k\|_A$ in terms of condition number $\kappa_2(A)$.
2. If $\kappa_2(A) = 10,000$, calculate the theoretical number of iterations required to reduce the initial energy norm of error by $10^{-6}$.
3. How many iterations would the Method of Steepest Descent take for the same problem?

**Solution**:
1. The Chebyshev bound is:
   $$\|e_k\|_A \le 2 \left( \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1} \right)^k \|e_0\|_A$$
2. For $\kappa = 10,000 \implies \sqrt{\kappa} = 100$:
   $$\rho_{CG} = \frac{100 - 1}{100 + 1} = \frac{99}{101} \approx 0.980198$$
   To achieve $2 \rho_{CG}^k \le 10^{-6} \implies \rho_{CG}^k \le 0.5 \times 10^{-6}$:
   $$k \ge \frac{\ln(5 \times 10^{-7})}{\ln(0.980198)} = \frac{-14.50866}{-0.02000} \approx \mathbf{726\text{ iterations}}$$
3. For Steepest Descent:
   $$\rho_{SD} = \frac{\kappa - 1}{\kappa + 1} = \frac{9999}{10001} \approx 0.99980002$$
   $$k \ge \frac{\ln(10^{-6})}{\ln(0.99980002)} = \frac{-13.81551}{-0.00020} \approx \mathbf{69,078\text{ iterations}}$$
   **CG is approximately $95\times$ faster!**

---

#### Question 3 [3 Marks]
Explain the concept of **Preconditioning** for Krylov subspace methods. Why is the condition $M \approx A$ insufficient on its own when selecting a practical preconditioner? State the three criteria for an effective preconditioner.

**Solution**:
1. Preconditioning transforms the original ill-conditioned system $Ax = b$ into an equivalent system $M^{-1} A x = M^{-1} b$ whose coefficient matrix $M^{-1} A$ has clustered eigenvalues and a vastly reduced condition number ($\kappa(M^{-1}A) \ll \kappa(A)$).
2. The condition $M \approx A$ alone is insufficient because taking $M = A$ yields $M^{-1} A = I$ (which converges in 1 iteration), but solving $M z = r \iff A z = r$ requires solving the exact original problem, defeating the purpose!
3. The **Three Golden Criteria** are:
   - **Spectral Clustering**: $M^{-1} A$ must be well-conditioned with eigenvalues tightly clustered near 1.
   - **Inexpensive Inversion**: The linear system $M z = r$ must be computationally cheap to solve in each iteration ($O(N)$ operations, e.g., triangular solves).
   - **Low Memory Overhead**: Storing the preconditioner matrix must require minimal additional RAM ($O(N)$ or zero extra storage).

---

#### Question 4 [3 Marks]
A sequential program executes in $T_s = 100\text{ seconds}$, with $75\%$ of its operations parallelizable.
1. What is the theoretical maximum speedup possible according to Amdahl's Law on an infinite number of processors ($p \to \infty$)?
2. If communication overhead on $p$ processors is $T_{\text{comm}}(p) = 0.5 \times \log_2 p\text{ seconds}$, calculate the actual speedup on $p = 64$ processors.

**Solution**:
1. Serial fraction $f = 1 - 0.75 = 0.25$.
   $$S_{\max} = \frac{1}{f} = \frac{1}{0.25} = \mathbf{4.0\times}$$
2. On $p = 64$ processors:
   - Serial computation time: $T_{\text{serial}} = 0.25 \times 100 = 25\text{ s}$.
   - Parallel computation time: $T_{\text{parallel}} = \frac{75}{64} \approx 1.1719\text{ s}$.
   - Communication overhead: $T_{\text{comm}} = 0.5 \times \log_2(64) = 0.5 \times 6 = 3.0\text{ s}$.
   - Total parallel runtime:
     $$T_p(64) = 25 + 1.1719 + 3.0 = \mathbf{29.1719\text{ s}}$$
   - Actual Speedup:
     $$\mathbf{S(64) = \frac{100}{29.1719} \approx \mathbf{3.428\times}}$$

---

#### Question 5 [3 Marks]
Differentiate between **Uniform Memory Access (UMA)** and **Non-Uniform Memory Access (NUMA)** shared-memory architectures. Why does a naive loop stride cause severe performance degradation on NUMA systems?

**Solution**:
1. **UMA (Symmetric Multiprocessing)**: All CPU cores connect to a centralized DRAM through a shared system bus or crossbar. Every processor has identical latency and bandwidth to any physical memory address.
2. **NUMA**: Memory is physically segmented and attached directly to individual processor sockets. All cores can address all memory locations via an integrated interconnect (e.g., Intel UPI, AMD Infinity Fabric). Accessing **local memory** (same socket) is fast ($40\text{ ns}$), whereas accessing **remote memory** (adjacent socket) incurs significantly higher latency ($120\text{ ns}$) and lower bandwidth.
3. **NUMA Degradation**: If an OpenMP thread on Socket 1 frequently accesses data initialized and allocated on Socket 0 (violating NUMA "first-touch" affinity), memory requests must continuously traverse the inter-socket interconnect, saturating cross-socket bandwidth and inducing massive processor wait states.

---

### Part B (10 Marks)

#### Question 6 [10 Marks]
Consider a sparse matrix $A \in \mathbb{R}^{6 \times 6}$ stored in **Diagonal Storage (DIA)** format with diagonal offset array $\text{IOFF} = [-2, 0, 1]$ and 2D array `DIAG` given below:

$$\text{IOFF} = [-2, \quad 0, \quad +1]$$

$$\text{DIAG} = \begin{bmatrix}
* & 4.0 & 1.0 \\
* & 5.0 & -1.0 \\
2.0 & 6.0 & 2.0 \\
-1.0 & 4.0 & 1.0 \\
1.0 & 5.0 & * \\
3.0 & 7.0 & *
\end{bmatrix}, \qquad
b = \begin{bmatrix}
5.0 \\
3.0 \\
16.0 \\
10.0 \\
11.0 \\
25.0
\end{bmatrix}$$
where `*` represents boundary zero padding.

1. **[4 Marks]** Explicitly reconstruct the full $6 \times 6$ coefficient matrix $A$.
2. **[2 Marks]** Verify whether matrix $A$ is Strictly Diagonally Dominant (SDD).
3. **[4 Marks]** Perform the **first iteration of the Gauss-Seidel method** starting with initial guess $x^{(0)} = [0, 0, 0, 0, 0, 0]^T$. Show all numerical calculations for $x_1^{(1)}$ through $x_6^{(1)}$.

#### Full Step-by-Step Solution:

1. **Reconstruction of Matrix $A$**:
   - $\text{IOFF}(1) = -2 \implies$ second sub-diagonal: $a_{i, i-2} = \text{DIAG}(i, 1)$ for $i = 3, 4, 5, 6$.
   - $\text{IOFF}(2) = 0 \implies$ main diagonal: $a_{i, i} = \text{DIAG}(i, 2)$ for $i = 1, \dots, 6$.
   - $\text{IOFF}(3) = +1 \implies$ first super-diagonal: $a_{i, i+1} = \text{DIAG}(i, 3)$ for $i = 1, \dots, 5$.

   Filling entry by entry:
   - **Row 1**: $a_{11} = 4.0, \; a_{12} = 1.0$.
   - **Row 2**: $a_{22} = 5.0, \; a_{23} = -1.0$.
   - **Row 3**: $a_{31} = 2.0, \; a_{33} = 6.0, \; a_{34} = 2.0$.
   - **Row 4**: $a_{42} = -1.0, \; a_{44} = 4.0, \; a_{45} = 1.0$.
   - **Row 5**: $a_{53} = 1.0, \; a_{55} = 5.0$.
   - **Row 6**: $a_{64} = 3.0, \; a_{66} = 7.0$.

   $$A = \begin{bmatrix}
   4.0 &  1.0 &  0.0 &  0.0 &  0.0 &  0.0 \\
   0.0 &  5.0 & -1.0 &  0.0 &  0.0 &  0.0 \\
   2.0 &  0.0 &  6.0 &  2.0 &  0.0 &  0.0 \\
   0.0 & -1.0 &  0.0 &  4.0 &  1.0 &  0.0 \\
   0.0 &  0.0 &  1.0 &  0.0 &  5.0 &  0.0 \\
   0.0 &  0.0 &  0.0 &  3.0 &  0.0 &  7.0
   \end{bmatrix}$$

2. **Strict Diagonal Dominance Verification**:
   Check $|a_{ii}| > \sum_{j \ne i} |a_{ij}|$ for each row:
   - Row 1: $|4.0| > |1.0| = 1.0$ (True: $4 > 1$)
   - Row 2: $|5.0| > |-1.0| = 1.0$ (True: $5 > 1$)
   - Row 3: $|6.0| > |2.0| + |2.0| = 4.0$ (True: $6 > 4$)
   - Row 4: $|4.0| > |-1.0| + |1.0| = 2.0$ (True: $4 > 2$)
   - Row 5: $|5.0| > |1.0| = 1.0$ (True: $5 > 1$)
   - Row 6: $|7.0| > |3.0| = 3.0$ (True: $7 > 3$)
   **Conclusion**: Every row strictly satisfies the dominance condition. Matrix $A$ is **Strictly Diagonally Dominant (SDD)**. Gauss-Seidel and Jacobi iterations are guaranteed to converge!

3. **Gauss-Seidel Iteration 1 ($x^{(0)} = \mathbf{0}$)**:
   Recall Gauss-Seidel updates in-place using newly computed components immediately:
   $$x_i^{(1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j < i} a_{ij} x_j^{(1)} - \sum_{j > i} a_{ij} x_j^{(0)} \right)$$

   - **$i = 1$**:
     $$x_1^{(1)} = \frac{b_1 - a_{12} x_2^{(0)}}{a_{11}} = \frac{5.0 - 1.0(0)}{4.0} = \mathbf{1.25}$$

   - **$i = 2$**:
     $$x_2^{(1)} = \frac{b_2 - a_{23} x_3^{(0)}}{a_{22}} = \frac{3.0 - (-1.0)(0)}{5.0} = \mathbf{0.60}$$

   - **$i = 3$**:
     $$x_3^{(1)} = \frac{b_3 - a_{31} x_1^{(1)} - a_{34} x_4^{(0)}}{a_{33}} = \frac{16.0 - 2.0(1.25) - 2.0(0)}{6.0} = \frac{16.0 - 2.5}{6.0} = \frac{13.5}{6.0} = \mathbf{2.25}$$

   - **$i = 4$**:
     $$x_4^{(1)} = \frac{b_4 - a_{42} x_2^{(1)} - a_{45} x_5^{(0)}}{a_{44}} = \frac{10.0 - (-1.0)(0.60) - 1.0(0)}{4.0} = \frac{10.0 + 0.60}{4.0} = \frac{10.60}{4.0} = \mathbf{2.65}$$

   - **$i = 5$**:
     $$x_5^{(1)} = \frac{b_5 - a_{53} x_3^{(1)}}{a_{55}} = \frac{11.0 - 1.0(2.25)}{5.0} = \frac{8.75}{5.0} = \mathbf{1.75}$$

   - **$i = 6$**:
     $$x_6^{(1)} = \frac{b_6 - a_{64} x_4^{(1)}}{a_{66}} = \frac{25.0 - 3.0(2.65)}{7.0} = \frac{25.0 - 7.95}{7.0} = \frac{17.05}{7.0} \approx \mathbf{2.4357}$$

   **Result after Iteration 1**:
   $$\mathbf{x^{(1)} = \begin{bmatrix} 1.2500 \\ 0.6000 \\ 2.2500 \\ 2.6500 \\ 1.7500 \\ 2.4357 \end{bmatrix}}$$

---

## 🎯 Master Exam Tips & Scoring Strategy

1. **Watch Matrix Symmetry in Krylov Methods**:
   Always verify if $A = A^T$. If an exam problem introduces boundary conditions (like Neumann or Robin) that cause $a_{ij} \ne a_{ji}$, do NOT apply the Conjugate Gradient algorithm! State clearly that GMRES or BiCGSTAB must be used.
2. **Spectral Radius Shortcuts**:
   When analyzing block triangular iteration matrices:
   $$\det(G - \lambda I) = \prod_k \det(G_{kk} - \lambda I)$$
   Do not waste time computing the determinant of the full $4 \times 4$ or $6 \times 6$ matrix.
3. **Sparse Storage Index Offsets**:
   Pay close attention to whether the problem specifies **0-indexing** (C-style) or **1-indexing** (Fortran-style). In MSR, `JA(1:n)` points to off-diagonal start locations; if `JA(k) == JA(k+1)`, row $k$ contains no off-diagonal entries.
4. **Young's Theorem for SOR**:
   Remember that for $\omega \ge \omega_{opt}$, the spectral radius is given directly by:
   $$\rho(G_{SOR}) = \omega - 1$$
   This simple relation allows instant derivation of error contraction without computing matrix powers!
