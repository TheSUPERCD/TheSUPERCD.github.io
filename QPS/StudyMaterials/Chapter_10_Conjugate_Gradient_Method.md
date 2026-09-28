# Chapter 10: Krylov Subspace Methods II: The Conjugate Gradient Method

---

## 1. Introduction and Foundations

The **Conjugate Gradient (CG) Method**, introduced by Magnus Hestenes and Eduard Stiefel (1952), is regarded as the premier iterative algorithm for solving large, sparse **Symmetric Positive Definite (SPD)** linear systems:
$$Ax = b$$
where $A = A^T \in \mathbb{R}^{n \times n}$ and $x^T A x > 0$ for all $x \ne \mathbf{0}$.

### 1.1 Motivation: Eliminating the Zigzagging of Steepest Descent
Recall from Chapter 08 that the Method of Steepest Descent minimizes the quadratic functional:
$$J(x) = \frac{1}{2} x^T A x - x^T b$$
by searching along the negative gradient $r_k = b - Ax_k$.
Because successive search directions are strictly orthogonal in the Euclidean metric ($r_{k+1} \perp r_k$), the search trajectory zigzags across narrow, elongated valleys when $\kappa_2(A) \gg 1$.
**The Breakthrough of Conjugate Gradient**:
Instead of making each search direction orthogonal to the previous one in the standard Euclidean inner product, CG constructs search directions that are **orthogonal in the inner product defined by matrix $A$ (the energy inner product)**!

---

## 2. Theory of $A$-Conjugate Directions

### 2.1 Definition of $A$-Conjugacy ($A$-Orthogonality)
Two non-zero vectors $u, v \in \mathbb{R}^n$ are said to be **$A$-conjugate** (or **$A$-orthogonal**) if:
$$\mathbf{\langle u, v \rangle_A = u^T A v = 0}$$

A set of non-zero vectors $\{p_0, p_1, \dots, p_{k-1}\}$ is mutually $A$-conjugate if:
$$p_i^T A p_j = 0 \quad \text{for all } i \ne j$$

### 2.2 Theorem: Linear Independence of $A$-Conjugate Vectors
> **Theorem**: If $A$ is Symmetric Positive Definite and the non-zero vectors $\{p_0, p_1, \dots, p_{k-1}\}$ are mutually $A$-conjugate, then they are **linearly independent**.

#### Mathematical Proof
Suppose there exist scalars $c_0, c_1, \dots, c_{k-1}$ such that:
$$\sum_{j=0}^{k-1} c_j p_j = \mathbf{0}$$
Multiply from the left by $p_i^T A$ (for any fixed $i \in \{0, \dots, k-1\}$):
$$p_i^T A \left( \sum_{j=0}^{k-1} c_j p_j \right) = \sum_{j=0}^{k-1} c_j (p_i^T A p_j) = 0$$
Due to mutual $A$-conjugacy, $p_i^T A p_j = 0$ for all $j \ne i$. The summation collapses to a single term:
$$c_i (p_i^T A p_i) = 0$$
Since $A$ is positive definite and $p_i \ne \mathbf{0}$, the quadratic form $p_i^T A p_i > 0$.
Dividing by $p_i^T A p_i$ gives:
$$c_i = 0$$
Since this holds for every $i \in \{0, \dots, k-1\}$, all coefficients must be identically zero.
Thus, $\{p_0, \dots, p_{k-1}\}$ are linearly independent. $\blacksquare$

### 2.3 Optimization along Conjugate Directions
If we minimize $J(x)$ successively along $A$-conjugate directions $p_0, p_1, \dots$:
- The minimization along direction $p_k$ **does not undo or degrade** the minimization achieved along previous directions $p_0, \dots, p_{k-1}$!
- After $m$ steps, $x_m$ is the **exact minimizer of $J(x)$ over the entire $m$-dimensional subspace** spanned by $\{p_0, \dots, p_{m-1}\}$.

---

## 3. Complete Derivation of the Conjugate Gradient Algorithm

The algorithm proceeds by generating sequences of iterates $x_j$, residuals $r_j = b - A x_j$, and conjugate search directions $p_j$.

### 3.1 Update Equations
1. **Solution Update**:
   $$x_{j+1} = x_j + \alpha_j p_j \quad \text{--- (1)}$$
2. **Residual Update**:
   Multiplying Eq. (1) by $A$ and subtracting from $b$:
   $$b - A x_{j+1} = b - A x_j - \alpha_j A p_j$$
   $$\mathbf{r_{j+1} = r_j - \alpha_j A p_j \quad \text{--- (2)}}$$
3. **Search Direction Update**:
   The new search direction $p_{j+1}$ is formed by taking the new steepest descent residual $r_{j+1}$ and adding a scaled multiple of the previous direction:
   $$\mathbf{p_{j+1} = r_{j+1} + \beta_j p_j \quad \text{--- (3)}}$$
   with $p_0 = r_0$.

---

### 3.2 Derivation of the Step Length $\alpha_j$
We demand that $x_{j+1}$ minimizes $J(x_j + \alpha p_j)$ along direction $p_j$.
From exact line search, the directional derivative along $p_j$ must vanish:
$$\nabla J(x_{j+1}) \cdot p_j = 0 \iff -r_{j+1}^T p_j = 0 \implies \mathbf{p_j^T r_{j+1} = 0}$$
Substitute Eq. (2) into this condition:
$$p_j^T (r_j - \alpha_j A p_j) = 0 \implies p_j^T r_j - \alpha_j p_j^T A p_j = 0$$
$$\alpha_j = \frac{p_j^T r_j}{p_j^T A p_j}$$
Now, using Eq. (3) at step $j-1$: $p_j = r_j + \beta_{j-1} p_{j-1}$.
Multiply by $r_j^T$:
$$p_j^T r_j = (r_j + \beta_{j-1} p_{j-1})^T r_j = r_j^T r_j + \beta_{j-1} (p_{j-1}^T r_j)$$
Because $p_{j-1}^T r_j = 0$ from the previous line search, the second term vanishes!
$$p_j^T r_j = r_j^T r_j$$
Therefore:
$$\mathbf{\alpha_j = \frac{r_j^T r_j}{p_j^T A p_j}}$$

---

### 3.3 Derivation of the Conjugacy Coefficient $\beta_j$
We enforce that the new direction $p_{j+1}$ is $A$-conjugate to the previous direction $p_j$:
$$\mathbf{p_j^T A p_{j+1} = 0}$$
Substitute Eq. (3) for $p_{j+1}$:
$$p_j^T A (r_{j+1} + \beta_j p_j) = 0$$
$$p_j^T A r_{j+1} + \beta_j p_j^T A p_j = 0$$
Solving for $\beta_j$:
$$\beta_j = -\frac{p_j^T A r_{j+1}}{p_j^T A p_j} = -\frac{r_{j+1}^T A p_j}{p_j^T A p_j} \quad (\text{since } A = A^T) \quad \text{--- (4)}$$

From the residual update Eq. (2):
$$r_{j+1} = r_j - \alpha_j A p_j \implies A p_j = \frac{r_j - r_{j+1}}{\alpha_j}$$
Substitute $A p_j$ into the numerator of Eq. (4):
$$r_{j+1}^T A p_j = r_{j+1}^T \left( \frac{r_j - r_{j+1}}{\alpha_j} \right) = \frac{r_{j+1}^T r_j - r_{j+1}^T r_{j+1}}{\alpha_j}$$
Because $r_{j+1}$ is orthogonal to $r_j$ ($r_{j+1}^T r_j = 0$):
$$r_{j+1}^T A p_j = -\frac{r_{j+1}^T r_{j+1}}{\alpha_j}$$
Substitute this back into Eq. (4):
$$\beta_j = -\frac{-\frac{r_{j+1}^T r_{j+1}}{\alpha_j}}{p_j^T A p_j} = \frac{r_{j+1}^T r_{j+1}}{\alpha_j (p_j^T A p_j)}$$
Now substitute $\alpha_j = \frac{r_j^T r_j}{p_j^T A p_j}$:
$$\beta_j = \frac{r_{j+1}^T r_{j+1}}{\left( \frac{r_j^T r_j}{p_j^T A p_j} \right) (p_j^T A p_j)} = \mathbf{\frac{r_{j+1}^T r_{j+1}}{r_j^T r_j}}$$
This is the famous **Fletcher-Reeves Formula**:
$$\mathbf{\beta_j = \frac{\|r_{j+1}\|_2^2}{\|r_j\|_2^2}}$$

---

## 4. Fundamental Theorems and Invariances of CG

At each iteration $k$ of the Conjugate Gradient algorithm:
1. **Orthogonality of Residuals**:
   $$\mathbf{r_i^T r_j = 0 \quad \text{for all } i \ne j}$$
2. **Mutual $A$-Conjugacy of Directions**:
   $$\mathbf{p_i^T A p_j = 0 \quad \text{for all } i \ne j}$$
3. **Subspace Equivalence**:
   $$\mathbf{\text{span}\{p_0, p_1, \dots, p_{k-1}\} = \text{span}\{r_0, r_1, \dots, r_{k-1}\} = \mathcal{K}_k(A, r_0)}$$
4. **Energy Norm Minimization (Galerkin Optimality)**:
   The iterate $x_k$ is the unique vector in the affine subspace $x_0 + \mathcal{K}_k(A, r_0)$ that minimizes the error in the $A$-norm:
   $$\mathbf{\|x_k - x^*\|_A = \min_{x \in x_0 + \mathcal{K}_k(A, r_0)} \|x - x^*\|_A}$$
   where $\|e\|_A = \sqrt{e^T A e} = \sqrt{(x - x^*)^T A (x - x^*)}$.

---

## 5. The Complete Conjugate Gradient Algorithm

```
Algorithm: Conjugate Gradient (CG) Method
Input: Matrix A (SPD), vector b, initial guess x_0, tolerance ε
Output: Solution vector x

1. Initialization:
   r_0 = b - A * x_0
   p_0 = r_0
   ρ_0 = r_0^T * r_0
   k = 0

2. While (sqrt(ρ_k) > ε) do:
       w_k = A * p_k                      // ONLY ONE SpMV per iteration!
       gamma = p_k^T * w_k                // Dot product
       
       α_k = ρ_k / gamma                  // Step size
       x_{k+1} = x_k + α_k * p_k          // AXPY update
       r_{k+1} = r_k - α_k * w_k          // AXPY residual update
       
       ρ_{k+1} = r_{k+1}^T * r_{k+1}      // Dot product
       β_k = ρ_{k+1} / ρ_k                // Fletcher-Reeves ratio
       p_{k+1} = r_{k+1} + β_k * p_k      // AXPY search direction update
       
       k = k + 1
   end while
```

### 5.1 Per-Iteration Computational and Storage Breakdown (Slide 221)
| Operation | Expression | FLOP Count |
| :--- | :--- | :--- |
| **Sparse Matrix-Vector Multiply** | $w = A p$ | $2 N_{nz}$ FLOPs |
| **Vector Dot Products (2)** | $p^T w$ and $r_{new}^T r_{new}$ | $4n$ FLOPs |
| **Vector AXPY Updates (3)** | $x \leftarrow x + \alpha p$, $r \leftarrow r - \alpha w$, $p \leftarrow r + \beta p$ | $6n$ FLOPs |
| **Total Work per Step** | | **$\mathbf{2 N_{nz} + 10n}$ FLOPs** |
| **Vector Storage Footprint** | Arrays: `x`, `r`, `p`, `w` | **Only $4n$ floating-point words!** |

---

## 6. Convergence Analysis & Chebyshev Error Bounds

### 6.1 Exact Arithmetic Termination
In theoretical exact arithmetic:
- The search directions $p_0, p_1, \dots, p_{n-1}$ form an $A$-orthogonal basis of $\mathbb{R}^n$.
- Therefore, the algorithm must terminate with $r_m = \mathbf{0}$ in **at most $m \le n$ iterations**!
- In practice, round-off errors destroy finite termination; CG behaves as a truly infinite iterative method governed by the eigenvalue distribution.

### 6.2 The Master Chebyshev Convergence Theorem (Slide 222)
> **Theorem**: Let $x_m$ be the $m$-th iterate of the Conjugate Gradient method applied to an SPD matrix $A$, and let $x^*$ be the exact solution. Let $\kappa = \lambda_{\max}/\lambda_{\min}$ be the spectral condition number. Then:
> $$\mathbf{\|x_m - x^*\|_A \le \frac{\|x_0 - x^*\|_A}{C_m(1 + 2\eta)} \le 2 \left( \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1} \right)^m \|x_0 - x^*\|_A}$$
> where:
> - $\eta = \frac{\lambda_{\min}}{\lambda_{\max} - \lambda_{\min}}$.
> - $C_m(t)$ is the **Chebyshev polynomial of the first kind** of degree $m$:
>   $$C_m(t) = \frac{1}{2} \left[ \left( t + \sqrt{t^2 - 1} \right)^m + \left( t - \sqrt{t^2 - 1} \right)^m \right]$$

### 6.3 Decisive Comparison: Steepest Descent vs. Conjugate Gradient
The asymptotic error reduction factors per iteration are:
$$\text{Steepest Descent}: \quad \rho_{SD} = \frac{\kappa - 1}{\kappa + 1}$$
$$\text{Conjugate Gradient}: \quad \mathbf{\rho_{CG} = \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}}$$

#### Quantitative Benchmark ($\kappa = 100$)
- **Steepest Descent**:
  $$\rho_{SD} = \frac{100 - 1}{100 + 1} = \frac{99}{101} \approx 0.9802$$
  Iterations to reduce error by $10^{-6}$:
  $$k \ge \frac{\ln(10^{-6})}{\ln(0.9802)} \approx \frac{-13.8155}{-0.0200} \approx \mathbf{691\text{ iterations}}$$
- **Conjugate Gradient**:
  $$\sqrt{\kappa} = \sqrt{100} = 10 \implies \rho_{CG} = \frac{10 - 1}{10 + 1} = \frac{9}{11} \approx 0.8182$$
  Iterations to reduce error by $10^{-6}$:
  $$k \ge \frac{\ln(10^{-6})}{\ln(0.8182)} \approx \frac{-13.8155}{-0.2007} \approx \mathbf{69\text{ iterations}}$$
- **Result**: Conjugate Gradient converges **$10\times$ faster**! For $\kappa = 10,000$, CG is **$100\times$ faster**!

### 6.4 Superlinear Convergence and Clustered Eigenvalues
If the eigenvalues of $A$ are clustered into $p$ tight clusters (with a few outliers), CG eliminates the outlying eigenvalues within the first few iterations.
Once the outliers are eliminated, CG converges at a rate governed by the clustered condition number, exhibiting **superlinear convergence** (the error decreases faster than any fixed geometric factor)!

### 6.5 Comprehensive Empirical Multi-Solver Benchmark (Slide 221)

To demonstrate the dramatic supremacy of the Conjugate Gradient method over classical stationary and gradient algorithms on discretized 2D elliptic PDEs, the course presents the following benchmark across mesh sizes $N = 16$ to $N = 256$ to reach convergence tolerances $\epsilon \approx 10^{-10}$ to $10^{-11}$:

| Solver Method | $N = 16$ Iter (Residual) | $N = 32$ Iter (Residual) | $N = 64$ Iter (Residual) | $N = 128$ Iter (Residual) | $N = 256$ Iter (Residual) | Dense FLOPs / Iter |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Jacobi (`jac`)** | 1,253 ($9.89 \times 10^{-11}$) | 4,446 ($9.98 \times 10^{-11}$) | 16,106 ($9.99 \times 10^{-11}$) | 58,828 ($1.00 \times 10^{-10}$) | 215,057 ($1.00 \times 10^{-10}$) | $n^2 + 2n$ |
| **Gauss-Seidel (`gs`)** | 624 ($9.72 \times 10^{-11}$) | 2,216 ($9.98 \times 10^{-11}$) | 8,038 ($9.99 \times 10^{-11}$) | 29,383 ($1.00 \times 10^{-10}$) | 107,466 ($1.00 \times 10^{-10}$) | $n^2 + 2n$ |
| **Minimal Residual (`mr`)** | 1,335 ($9.52 \times 10^{-11}$) | 4,746 ($9.99 \times 10^{-11}$) | 17,562 ($9.99 \times 10^{-11}$) | 65,725 ($1.00 \times 10^{-10}$) | 252,191 ($1.00 \times 10^{-10}$) | $n^2 + 3n$ |
| **Steepest Descent (`sd`)** | 1,313 ($9.93 \times 10^{-11}$) | 4,830 ($1.00 \times 10^{-10}$) | 17,891 ($1.00 \times 10^{-10}$) | 67,021 ($1.00 \times 10^{-10}$) | 247,067 ($1.00 \times 10^{-10}$) | $n^2 + 3n$ |
| **SOR (`sor`)** | 202 ($9.84 \times 10^{-11}$) | 762 ($9.81 \times 10^{-11}$) | 2,815 ($9.93 \times 10^{-11}$) | 10,381 ($9.98 \times 10^{-11}$) | 38,221 ($1.00 \times 10^{-10}$) | $n^2 + 4n$ |
| **Conjugate Gradient (`cg`)** | **32** ($8.50 \times 10^{-15}$) | **63** ($7.19 \times 10^{-11}$) | **124** ($6.99 \times 10^{-11}$) | **247** ($7.74 \times 10^{-11}$) | **484** ($7.44 \times 10^{-11}$) | $\mathbf{n^2 + n}$ |

#### Critical Exam Observations from Benchmark
1. **Gauss-Seidel vs. Jacobi**: Gauss-Seidel consistently requires almost exactly **half the iterations of Jacobi** across all mesh resolutions ($624$ vs $1253$ at $N=16$; $107,466$ vs $215,057$ at $N=256$), validating the theoretical asymptotic convergence rate $R_\infty(G) \approx 2 R_\infty(J)$.
2. **Steepest Descent vs. CG**: While Steepest Descent suffers from severe orthogonal zig-zagging requiring $247,067$ iterations at $N=256$, Conjugate Gradient enforces $A$-conjugacy, collapsing the iteration count to just **$484$ iterations** — a staggering **$510\times$ acceleration**!
3. **Algorithmic Scaling with Mesh Refinement ($h = 1/N$)**:
   - For Jacobi and Gauss-Seidel, doubling $N$ quadruples the iteration count ($O(N^2) = O(h^{-2})$ iterations).
   - For SOR, doubling $N$ roughly doubles the iteration count ($O(N) = O(h^{-1})$ iterations).
   - For Conjugate Gradient, doubling $N$ strictly doubles the iteration count ($O(\sqrt{\kappa(A)}) = O(N) = O(h^{-1})$ iterations), but with a dramatically smaller proportionality constant than SOR ($484$ vs $38,221$ at $N=256$, an **$79\times$ speedup over SOR**)!
4. **Computational Cost per Iteration (Slide 221)**:
   - For general dense operations, the floating-point count per step is:
     - **SOR**: $n^2 + 4n$ operations.
     - **Steepest Descent**: $n^2 + 3n$ operations.
     - **Conjugate Gradient**: $n^2 + n$ operations (least overhead per step among non-stationary methods!).
   - In sparse finite-difference implementations, $n^2$ matrix-vector multiplication is replaced by $5n$ or $7n$ stencil operations, making CG overwhelmingly optimal in total runtime.

---

> [!IMPORTANT]
> **Past Exam Focus on Conjugate Gradient Scaling & Constraints**:
> - **2025 Midsem Q1**: In a finite-element 1D Poisson problem ($[K]\{T\} = \{F\}$), why refining the mesh by $10\times$ scales CG iterations as $O(\sqrt{\kappa}) \approx \sqrt{10} \approx 3.16\text{--}4.5\times$, why CPU time jumps from $1.0\text{ s}$ to $11.62\text{ s}$ due to algorithmic complexity and cache spillover, and why BiCGSTAB takes roughly double the CPU time per iteration (2 SpMVs, 4 dot products vs 1 SpMV, 2 dot products):  
>   👉 [**2025 Exam Q1 Detailed Solution**](file:///home/thesupercd/Documents/HPSC_TutorialsNAssignments/ExamNotes/Previous_Years_Question_Papers_Solutions.md#question-1-6-marks).
> - **2018 Mid-Spring Q4**: Construction of a $3 \times 3$ strictly diagonally dominant non-symmetric matrix that can be solved by Gauss-Seidel but **cannot** be solved by Conjugate Gradient:  
>   👉 [**2018 Exam Q4 Detailed Solution**](file:///home/thesupercd/Documents/HPSC_TutorialsNAssignments/ExamNotes/Previous_Years_Question_Papers_Solutions.md#question-4-2-marks).
