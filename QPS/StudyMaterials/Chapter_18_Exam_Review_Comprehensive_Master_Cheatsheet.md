# Chapter 18: Comprehensive Exam Review, Master Synthesis & Formula Cheatsheet

---

## 1. Master Formula Cheatsheet

### 1.1 Computer Architecture & Memory Performance
| Concept | Mathematical Formula | Key Variables & Notes |
| :--- | :--- | :--- |
| **Peak Processor Speed** | $\text{Peak FLOPS} = \text{Cores} \times f \times \frac{\text{FLOPs}}{\text{Cycle}}$ | $f = \text{frequency (Hz)}$. E.g., $1\text{ GHz} \times 4\text{ ops} = 4\text{ GFLOPS}$. |
| **Cycle Time** | $t_{cycle} = \frac{1}{f}$ | At $1\text{ GHz}$, $t_{cycle} = 10^{-9}\text{ s} = 1\text{ ns}$. |
| **Uncached Effective Speed** | $\text{Effective Speed} = \frac{1}{t_{latency}}$ | At $100\text{ ns}$ DRAM latency, effective speed = $10\text{ MFLOPS}$ ($400\times$ drop!). |
| **Effective Access Time (EAT)**| $\text{EAT} = H \cdot t_{cache} + (1 - H) \cdot t_{DRAM}$ | $H$ = Cache Hit Ratio ($0 \le H \le 1$). |
| **Cache Line Size** | Standard = $64\text{ bytes}$ | Holds exactly $8$ double-precision (64-bit) floats. |

---

### 1.2 Discretization & Numerical PDEs
| Concept | Mathematical Formula | Notes |
| :--- | :--- | :--- |
| **Central Difference (2nd Deriv.)**| $\left.\frac{d^2 T}{dx^2}\right|_i = \frac{T_{i+1} - 2T_i + T_{i-1}}{\Delta x^2} - \frac{\Delta x^2}{12} T^{(4)}(\xi)$ | Truncation error is strictly second-order $O(\Delta x^2)$. |
| **2D Laplace 5-Point Stencil**| $4T_{i, j} - T_{i+1, j} - T_{i-1, j} - T_{i, j+1} - T_{i, j-1} = 0$ | Uniform grid ($\Delta x = \Delta y$). Pentadiagonal matrix. |
| **Lexicographical 1D Index** | $\text{ipt}(i, j) = (j - 1)N_x + i$ | Maps $(i, j) \to \text{ipt}$; neighbors at $\pm 1$ and $\pm N_x$. |
| **3D Poisson 7-Point Stencil**| $6T_{\text{ipt}} - \sum_{\text{neighbors}} T = \Delta x^2 f_{\text{ipt}}$ | Septadiagonal matrix; index $\text{ipt} = i + (j-1)N_x + (k-1)N_x N_y$. |
| **FEM 1D Linear Shape Functions**| $N_1(\xi) = 1 - \xi, \quad N_2(\xi) = \xi$ | $\xi = (x - x_m)/\Delta x \in [0, 1]$. |
| **FEM 1D Element Stiffness** | $k^e = \frac{1}{\Delta x} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$ | Assembles to standard central difference stencil! |
| **CSR Storage Requirement** | $\text{Storage} = 2 N_{nz} + N + 1$ | `val` ($N_{nz}$), `col_ind` ($N_{nz}$), `row_ptr` ($N+1$). |

---

### 1.3 Linear Algebra Foundations & Spectral Theory
| Concept | Mathematical Formula | Notes |
| :--- | :--- | :--- |
| **Orthogonal Complement** | $N(A) = (C(A^T))^\perp, \quad N(A^T) = (C(A))^\perp$ | Fundamental Theorem of Linear Algebra. |
| **Rank-Nullity Theorem** | $\dim C(A^T) + \dim N(A) = n$ | For $A \in \mathbb{R}^{m \times n}$ with rank $r$: $r + (n - r) = n$. |
| **Solvability Condition** | $Ax = b \text{ solvable} \iff b \perp N(A^T)$ | Fredholm alternative ($b \in C(A)$). |
| **Right Inverse (Existence)** | $C = A^T (A A^T)^{-1}$ | Exists iff $A$ has full row rank ($r = m \le n$). $AC = I_m$. |
| **Left Inverse (Uniqueness)** | $B = (A^T A)^{-1} A^T$ | Exists iff $A$ has full column rank ($r = n \le m$). $BA = I_n$. |
| **Rayleigh Quotient** | $R(x) = \frac{x^T A x}{x^T x}$ | $\lambda_{\min} \le R(x) \le \lambda_{\max}$ for symmetric $A$. |
| **Spectral Condition Number**| $\kappa_2(A) = \frac{\sigma_{\max}}{\sigma_{\min}} \quad (\text{or } \frac{\lambda_{\max}}{\lambda_{\min}} \text{ for SPD})$ | $\frac{\|\delta x\|}{\|x\|} \le \kappa(A) \frac{\|\delta b\|}{\|b\|}$. |
| **Sylvester's Criterion** | $\det(A_k) > 0 \quad \text{for all } k = 1, \dots, n$ | Necessary & sufficient test for Positive Definiteness. |

---

### 1.4 Direct Solvers & TDMA
| Method / Theorem | Mathematical Formulation | Complexity |
| :--- | :--- | :--- |
| **Gaussian Elimination / LU**| $A = L U \implies L y = b \text{ then } U x = y$ | $\frac{2}{3} N^3$ FLOPs (Dense); severe fill-in for sparse PDEs. |
| **Cholesky Factorization** | $A = L L^T \implies l_{jj} = \sqrt{a_{jj} - \sum l_{jk}^2}$ | $\frac{1}{3} N^3$ FLOPs (Exactly half of LU!); SPD matrices only. |
| **TDMA Forward Recurrences** | $\beta_i = b_i - \frac{a_i c_{i-1}}{\beta_{i-1}}, \quad \gamma_i = \frac{d_i - a_i \gamma_{i-1}}{\beta_i}$ | $\beta_1 = b_1, \quad \gamma_1 = d_1 / b_1$. |
| **TDMA Backward Substitution**| $T_N = \gamma_N, \quad T_i = \gamma_i - \frac{c_i}{\beta_i} T_{i+1}$ | **$O(N)$ operations!** ($10$ ops for $5 \times 5$ vs $75$ for Gauss). |
| **TDMA Stability Condition** | $|b_i| \ge |a_i| + |c_i|$ with strict inequality for one $i$ | Irreducibly diagonally dominant $\implies$ no pivoting needed. |
| **Sherman-Morrison Formula** | $(A + u v^T)^{-1} = A^{-1} - \frac{A^{-1} u v^T A^{-1}}{1 + v^T A^{-1} u}$ | Solves rank-1 perturbed & cyclic tridiagonal systems in $O(N)$. |

---

### 1.5 Classical Iterative Solvers & Convergence Theory
| Method | Iteration Matrix ($G$) & Step | Convergence Condition & Rate |
| :--- | :--- | :--- |
| **General Splitting** | $x^{(k+1)} = G x^{(k)} + f, \quad G = M^{-1} N$ | Converges $\iff \rho(G) < 1$. |
| **Asymptotic Rate** | $\tau = -\ln \rho(G)$ | Iterations for $\epsilon$: $k \ge \frac{-\ln \epsilon}{\tau}$. |
| **Jacobi Method** | $G_J = D^{-1}(E + F) = I - D^{-1} A$ | Converges if $A$ is Strictly Diagonally Dominant. |
| **Gauss-Seidel Method** | $G_{GS} = (D - E)^{-1} F$ | $\rho(G_{GS}) = [\rho(G_J)]^2 \implies \tau_{GS} = 2 \tau_J$ ($2\times$ faster!). |
| **SOR Method** | $G_{SOR} = (D - \omega E)^{-1} [(1 - \omega) D + \omega F]$ | Ostrowski-Reich: Converges for SPD $A \iff 0 < \omega < 2$. |
| **Optimal SOR Parameter** | $\mathbf{\omega_{opt} = \frac{2}{1 + \sqrt{1 - \rho(G_J)^2}}}$ | Optimal spectral radius: $\rho(G_{SOR, opt}) = \omega_{opt} - 1$. |

---

### 1.6 Projection, Steepest Descent & Krylov Subspace Solvers
| Solver | Key Mathematical Formulas | Optimization / Convergence Property |
| :--- | :--- | :--- |
| **Petrov-Galerkin** | $x_m \in x_0 + \mathcal{K}_m, \quad (b - A x_m) \perp \mathcal{L}_m$ | Fundamental framework for all modern iterative solvers. |
| **Steepest Descent (SD)** | $\alpha_k = \frac{r_k^T r_k}{r_k^T A r_k}, \quad x_{k+1} = x_k + \alpha_k r_k$ | Successive residuals are orthogonal ($r_{k+1} \perp r_k$). |
| **SD Convergence Bound** | $\|e_{k+1}\|_A \le \left( \frac{\kappa - 1}{\kappa + 1} \right) \|e_k\|_A$ | Kantorovich inequality; zigzags when $\kappa \gg 1$. |
| **Arnoldi Relation** | $A V_m = V_m H_m + h_{m+1, m} v_{m+1} e_m^T = V_{m+1} \bar{H}_m$ | $V_m^T A V_m = H_m$ (Upper Hessenberg matrix). |
| **Lanczos 3-Term Recur.**| $\beta_j v_{j+1} = A v_j - \alpha_j v_j - \beta_{j-1} v_{j-1}$ | Symmetric $A \implies H_m = T_m$ is **tridiagonal**! |
| **Conjugate Gradient (CG)**| $\alpha_k = \frac{r_k^T r_k}{p_k^T A p_k}, \quad \beta_k = \frac{r_{k+1}^T r_{k+1}}{r_k^T r_k}$ | Directions are $A$-conjugate ($p_i^T A p_j = 0$). |
| **CG Convergence Bound** | $\|e_m\|_A \le 2 \left( \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1} \right)^m \|e_0\|_A$ | Chebyshev bound; vastly faster than Steepest Descent! |
| **GMRES Least-Squares** | $y_m = \arg\min_y \|\beta e_1 - \bar{H}_m y\|_2$ | Transformed via Givens rotations: $\mathbf{\|r_m\|_2 = \|\gamma_{m+1}\|}$. |
| **BiCGSTAB Step** | $s_j = r_j - \alpha_j A p_j, \quad \omega_j = \frac{\langle A s_j, s_j \rangle}{\|A s_j\|_2^2}$ | $r_{j+1} = s_j - \omega_j A s_j$; 2 SpMVs, no $A^T$, smooth. |
| **PCG Preconditioned Step**| $M z_{k+1} = r_{k+1}, \quad \beta_k = \frac{r_{k+1}^T z_{k+1}}{r_k^T z_k}$ | Solves $M z = r$ once per step; preserves SPD symmetry! |

---

### 1.7 Parallel Computing Architecture & Scalability Laws
| Law / Metric | Mathematical Expression | Asymptotic Limits & Physical Meaning |
| :--- | :--- | :--- |
| **Parallel Speedup** | $S(p) = \frac{T_s}{T_p}$ | Linear: $S = p$. Sublinear: $S < p$. Superlinear: $S > p$. |
| **Parallel Efficiency** | $E(p) = \frac{S(p)}{p} = \frac{T_s}{p T_p}$ | $E = \frac{1}{1 + T_o / T_s} \le 1$. |
| **Total Overhead** | $T_o(W, p) = p T_p - T_s$ | $T_p = \frac{T_s + T_o}{p}$. |
| **Isoefficiency Function**| $\mathbf{W = K \cdot T_o(W, p)}$ | $K = \frac{E}{1 - E}$. Growth rate of $W$ to keep $E$ constant. |
| **Amdahl's Law** | $\mathbf{S(p) \le \frac{1}{f + \frac{1-f}{p}}}$ | Serial fraction $f = \Phi_s / T_s$. Upper bound: $\mathbf{S_{\max} = 1/f}$. |
| **Gustafson's Law** | $\mathbf{S_{\text{scaled}}(p) = p - (p - 1) s}$ | Scaled speedup (weak scaling); linear scaling preserved! |
| **Wormhole Routing Time** | $T_{comm} \approx t_s + m \cdot t_w$ | $t_s$ = startup latency, $m$ = words, $t_w = 1/\text{bandwidth}$. |
| **Cannon's Matrix Multiply**| $T_p = \frac{2n^3}{p} + 2\sqrt{p} \, t_s + 2\frac{n^2}{\sqrt{p}} \, t_w$ | 2D Torus algorithm; Isoefficiency $W = O(p^{1.5})$. |

---

## 2. Master Comparison Matrix of Linear System Solvers

| Method | Target Matrix Restrictions | Subspaces $\mathcal{K}_m, \mathcal{L}_m$ | Minimization Property | Recurrence Length | Work per Iteration | Storage Vectors | Numerical Stability |
| :--- | :--- | :--- | :--- | :---: | :--- | :---: | :--- |
| **Gauss-Seidel** | Diagonally dominant or SPD | Fixed-point splitting | None | 1 | 1 SpMV | $1n$ (in-place) | Stable; poor parallelism |
| **SOR** | SPD ($0 < \omega < 2$) | Accelerated splitting | None | 1 | 1 SpMV | $1n$ (in-place) | Fast at $\omega_{opt}$; serial chain |
| **TDMA (Thomas)** | Tridiagonal, Diag. Dominant | Direct Elimination | Exact solve | 1 | $8N$ ops total | $O(N)$ arrays | Unconditionally stable |
| **Steepest Descent**| Symmetric Positive Definite | $\mathcal{K}_1 = r_k, \mathcal{L}_1 = r_k$ | $J(x)$ along line $r_k$ | 1 | 1 SpMV, 2 dots, 2 AXPYs | $3n$ | Severe zigzagging ($\kappa \gg 1$) |
| **Conjugate Gradient**| Symmetric Positive Definite | $\mathcal{K}_m = \mathcal{K}_m(A, r_0), \mathcal{L}_m = \mathcal{K}_m$ | $\|x - x^*\|_A$ over $\mathcal{K}_m$ | 2 (short) | 1 SpMV, 2 dots, 3 AXPYs | $4n$ | **Optimal for SPD**; Chebyshev |
| **PCG** | Symmetric Positive Definite | Split $\tilde{A} = L^{-1} A L^{-T}$ | $\|x - x^*\|_A$ (preconditioned) | 2 (short) | 1 SpMV, 1 Solve $M z = r$ | $5n$ | Fast; reduces $\kappa \to 1$ |
| **FOM** | Any Non-Singular Square | $\mathcal{K}_m = \mathcal{K}_m(A, r_0), \mathcal{L}_m = \mathcal{K}_m$ | Galerkin ($r_m \perp \mathcal{K}_m$) | $m$ (full) | 1 SpMV, $m$ dots, $m$ AXPYs | $(m+2)n$ | Can break down if $H_m$ singular |
| **GMRES** | Any Non-Singular Square | $\mathcal{K}_m = \mathcal{K}_m(A, r_0), \mathcal{L}_m = A \mathcal{K}_m$ | $\|r_m\|_2$ over $\mathcal{K}_m$ | $m$ (full) | 1 SpMV, $m$ dots, Givens | $(m+2)n$ | **Monotonic decrease**; $O(mn)$ memory |
| **GMRES($m$)** | Any Non-Singular Square | Restarted every $m$ steps | Local $\|r_m\|_2$ minimum | $m$ (capped) | 1 SpMV, $m$ dots | $(m+2)n$ | Fixed memory; risk of stagnation |
| **BiCG** | Any Non-Singular Square | $\mathcal{K}_m(A), \mathcal{L}_m = \mathcal{K}_m(A^T)$ | Oblique Petrov-Galerkin | 2 (short) | 2 SpMVs ($A, A^T$), 2 dots | $6n$ | Erratic; requires $A^T$; breakdowns |
| **CGS** | Any Non-Singular Square | Squared polynomial $\phi_j^2(A)$ | None (squared contraction) | 2 (short) | 2 SpMVs (no $A^T$!), 2 dots | $7n$ | Fast; extreme round-off amplification |
| **BiCGSTAB** | Any Non-Singular Square | Stabilized polynomial $\psi_j \phi_j$ | Local 1D SD residual norm | 2 (short) | 2 SpMVs (no $A^T$!), 4 dots | $7n$ | **Industry standard for non-symmetric** |

---

## 3. High-Yield Exam Questions & Comprehensive Answers

### Question 1: Why does a serial calculation of turbulent combustion take lifetimes, and how does parallel computing overcome this?
- **Answer**: In Direct Numerical Simulation (DNS) of 3D turbulent combustion, spatial grid resolution scales as $(L/\eta)^3 \sim Re^{9/4}$, requiring $N \approx 10^9 - 10^{11}$ nodes. Simulating $10\text{ seconds}$ at $\Delta t = 10^{-3}\text{ s}$ demands $10^4$ time steps. At $O(10^{13})\text{ FLOP/step}$, total work is $10^{17}\text{ FLOP}$. A single-core $10\text{ GFLOPS}$ processor requires $10^7\text{ seconds} \approx 3.17\text{ years}$ (and fine grids of $10^{13}$ floats require up to $317\text{ years}$). Furthermore, storing $10^{13}$ double-precision floats requires $160\text{ Terabytes}$ of RAM, exceeding uniprocessor capacity. Parallel HPC clusters overcome this by:
  1. Aggregating hundreds of Terabytes of distributed memory across thousands of nodes.
  2. Delivering sustained PetaFLOPS throughput ($10^{16}\text{ FLOP/s}$), reducing runtime from lifetimes down to minutes.

---

### Question 2: State the exact differences between Verification and Validation (V&V).
- **Answer**:
  - **Verification ("Solving the equations right")**: The process of determining that a computational model correctly implements and solves the intended mathematical model. It evaluates numerical correctness, consistency, grid convergence ($O(\Delta x^2)$), and absence of programming bugs (e.g., using the Method of Manufactured Solutions).
  - **Validation ("Solving the right equations")**: The process of determining the degree to which a mathematical model accurately represents the real physical world from the perspective of the intended engineering application. It evaluates physics fidelity by comparing numerical results directly against high-precision laboratory experiments.

---

### Question 3: Why does total numerical error increase if the grid spacing $\Delta x$ is made excessively small?
- **Answer**: Total numerical error is the sum of discretization/truncation error and floating-point round-off error:
  $$\text{Error}_{total} = \text{Error}_{truncation} + \text{Error}_{round-off} \approx C_1 (\Delta x)^p + \frac{C_2}{\Delta x^q}$$
  Truncation error vanishes monotonically as $\Delta x \to 0$. However, reducing $\Delta x$ increases the total number of grid points and required floating-point operations. Because each operation on a digital computer incurs finite precision round-off (IEEE 754), round-off errors accumulate with the number of operations. Below an optimal threshold $\Delta x_{optimal}$, round-off error dominates, causing total error to explode.

---

### Question 4: State the Fundamental Theorem of Linear Algebra and prove that $N(A) \perp C(A^T)$.
- **Answer**:
  - **Part 1 (Dimensions)**: For $A \in \mathbb{R}^{m \times n}$ with rank $r$: $\dim C(A) = \dim C(A^T) = r$, $\dim N(A) = n - r$, and $\dim N(A^T) = m - r$.
  - **Part 2 (Orthogonality)**: $N(A) = (C(A^T))^\perp$ in $\mathbb{R}^n$, and $N(A^T) = (C(A))^\perp$ in $\mathbb{R}^m$.
  - **Proof of $N(A) \perp C(A^T)$**:
    Let $x \in N(A) \implies Ax = \mathbf{0}$. Writing $A$ in terms of its rows $r_1^T, \dots, r_m^T$:
    $$Ax = \begin{bmatrix} r_1^T x \\ \vdots \\ r_m^T x \end{bmatrix} = \begin{bmatrix} 0 \\ \vdots \\ 0 \end{bmatrix} \implies r_i \cdot x = 0 \quad \text{for all } i = 1, \dots, m$$
    Any vector $v \in C(A^T)$ is a linear combination of rows: $v = \sum_{i=1}^m c_i r_i$.
    Taking the dot product:
    $$v \cdot x = \left( \sum_{i=1}^m c_i r_i \right) \cdot x = \sum_{i=1}^m c_i (r_i \cdot x) = \sum_{i=1}^m c_i (0) = 0$$
    Since $\dim C(A^T) + \dim N(A) = r + (n - r) = n$, $N(A)$ is the complete orthogonal complement of $C(A^T)$.

---

### Question 5: When does a rectangular matrix have a Right Inverse vs. a Left Inverse? What are their implications for $Ax = b$?
- **Answer**:
  - **Right Inverse ($AC = I_m$)**: Exists if and only if $A \in \mathbb{R}^{m \times n}$ has **full row rank** ($r = m \le n$). Formula: $C = A^T (A A^T)^{-1}$. Implication: A solution to $Ax = b$ **always exists** for any $b$ ($x = Cb$), but is non-unique (infinitely many solutions if $m < n$).
  - **Left Inverse ($BA = I_n$)**: Exists if and only if $A \in \mathbb{R}^{m \times n}$ has **full column rank** ($r = n \le m$). Formula: $B = (A^T A)^{-1} A^T$. Implication: If a solution to $Ax = b$ exists, it is **strictly unique** ($x = Bb$). If $b \notin C(A)$, no exact solution exists, and $x = Bb$ provides the unique Ordinary Least Squares solution.

---

### Question 6: Prove that all eigenvalues of a real symmetric matrix are real.
- **Answer**: Let $A = A^T \in \mathbb{R}^{n \times n}$. Let $Ax = \lambda x$ with $x \ne \mathbf{0}$.
  Take the conjugate transpose: $x^H A^H = \bar{\lambda} x^H$.
  Since $A$ is real and symmetric, $A^H = A^T = A \implies x^H A = \bar{\lambda} x^H$.
  Multiply first equation by $x^H$ on left: $x^H A x = \lambda x^H x$.
  Multiply conjugate equation by $x$ on right: $x^H A x = \bar{\lambda} x^H x$.
  Subtracting yields:
  $$0 = (\lambda - \bar{\lambda}) (x^H x)$$
  Since $x \ne \mathbf{0}$, $x^H x = \|x\|_2^2 > 0$. Dividing by $\|x\|_2^2$ gives $\lambda - \bar{\lambda} = 0 \implies \lambda = \bar{\lambda} \implies \lambda \in \mathbb{R}$.

---

### Question 7: Why is Gaussian Elimination / LU decomposition avoided for large 3D PDE systems?
- **Answer**:
  1. **Arithmetic FLOP Explosion**: Dense elimination scales as $O(N^3)$. For a 3D grid with $N = 10^6$, $\frac{2}{3} N^3 \sim 10^{18}\text{ FLOPs}$ (years of runtime).
  2. **Fill-In and Memory Explosion**: Even though $A$ has only 7 non-zeros per row, Gaussian elimination introduces non-zero entries into zero positions within the band. In 3D (bandwidth $W \sim N^{2/3}$), fill-in requires $O(N^{5/3})$ memory, exceeding physical RAM.
  3. **Round-Off Accumulation**: Accumulation of rounding errors over $O(N^3)$ operations destroys numerical precision unless expensive pivoting is used, which ruins parallel scalability.

---

### Question 8: Why is the Thomas Algorithm (TDMA) unconditionally stable without pivoting for diagonally dominant matrices?
- **Answer**: The forward recurrence computes $\beta_i = b_i - \frac{a_i c_{i-1}}{\beta_{i-1}}$.
  For a strictly diagonally dominant matrix, $|b_1| > |c_1| \implies |\beta_1| > |c_1| \implies \frac{|c_1|}{|\beta_1|} < 1$.
  By induction, if $|\beta_{i-1}| > |c_{i-1}|$, then:
  $$|\beta_i| = \left| b_i - a_i \frac{c_{i-1}}{\beta_{i-1}} \right| \ge |b_i| - |a_i| \frac{|c_{i-1}|}{|\beta_{i-1}|} > |b_i| - |a_i| \ge |c_i|$$
  Thus, $|\beta_i| > |c_i|$ for all $i$, which guarantees:
  1. $\beta_i \ne 0$ (division by zero is impossible).
  2. The backward multiplier $|c_i / \beta_i| < 1$, preventing any amplification of round-off error during back-substitution!

---

### Question 9: State the Sherman-Morrison formula and explain its application to periodic boundary conditions.
- **Answer**:
  $$(A + u v^T)^{-1} = A^{-1} - \frac{A^{-1} u v^T A^{-1}}{1 + v^T A^{-1} u}$$
  Under periodic boundary conditions ($T(0) = T(L)$), the matrix is cyclic tridiagonal, having non-zero corner elements $a_1$ and $c_N$. This breaks standard TDMA.
  By expressing $A_{cyclic} = A_{tridiag} + u v^T$ (where $u = [\alpha, 0, \dots, \beta]^T$ and $v = [1, 0, \dots, \gamma]^T$), the system can be solved by performing **two standard TDMA passes** on $A_{tridiag}$ and updating the result via the vector formula, preserving optimal $O(N)$ complexity!

---

### Question 10: Prove that Gauss-Seidel converges twice as fast as Jacobi for consistently ordered matrices.
- **Answer**: By David Young's Theorem, the eigenvalues $\lambda \in \sigma(G_{SOR})$ and $\mu \in \sigma(G_J)$ satisfy:
  $$(\lambda + \omega - 1)^2 = \lambda \omega^2 \mu^2$$
  Setting $\omega = 1$ for Gauss-Seidel:
  $$\lambda^2 = \lambda \mu^2 \implies \lambda = \mu^2$$
  Taking the maximum modulus across all eigenvalues:
  $$\rho(G_{GS}) = \max |\lambda| = \max |\mu|^2 = \left[ \rho(G_J) \right]^2$$
  The asymptotic rate of convergence is:
  $$\tau_{GS} = -\ln \rho(G_{GS}) = -\ln [\rho(G_J)]^2 = -2 \ln \rho(G_J) = 2 \tau_J$$
  Because $\tau_{GS} = 2 \tau_J$, Gauss-Seidel requires exactly half the iterations of Jacobi to achieve the same error reduction!

---

### Question 11: Write the formula for the optimal SOR parameter $\omega_{opt}$ and state its behavior at the limits.
- **Answer**:
  $$\omega_{opt} = \frac{2}{1 + \sqrt{1 - \rho(G_J)^2}}$$
  - If $\rho(G_J) \to 0$ (well-conditioned/decoupled): $\omega_{opt} \to \frac{2}{1 + 1} = 1$ (SOR reduces to standard Gauss-Seidel).
  - If $\rho(G_J) \to 1$ (infinitely refined mesh): $\omega_{opt} \to \frac{2}{1 + 0} = 2$.
  - For any SPD matrix, if $\omega \ge 2$, SOR unconditionally diverges.

---

### Question 12: Why is Classical Gram-Schmidt (CGS) avoided in numerical computing, and how does Modified Gram-Schmidt (MGS) fix it?
- **Answer**: In floating-point arithmetic, when projecting a vector $a_j$ that is nearly parallel to the existing basis, CGS subtracts projections based on the original raw vector $a_j$, suffering catastrophic subtractive cancellation. As $j$ increases, newly generated vectors lose orthogonality to early basis vectors ($q_1^T q_j \gg \epsilon$).
  MGS fixes this by sequentially projecting the **working intermediate vector** onto each $q_i$ and updating the vector in-place before the next projection. This maintains orthogonality to near machine precision.

---

### Question 13: Prove that successive residuals in Steepest Descent are orthogonal ($r_{k+1} \perp r_k$).
- **Answer**: In Steepest Descent, $x_{k+1} = x_k + \alpha_k r_k$ with $\alpha_k = \frac{r_k^T r_k}{r_k^T A r_k}$.
  The new residual is $r_{k+1} = b - A x_{k+1} = r_k - \alpha_k A r_k$.
  Taking the inner product:
  $$r_{k+1}^T r_k = (r_k - \alpha_k A r_k)^T r_k = r_k^T r_k - \alpha_k (r_k^T A r_k)$$
  Substituting $\alpha_k$:
  $$r_{k+1}^T r_k = r_k^T r_k - \left( \frac{r_k^T r_k}{r_k^T A r_k} \right) (r_k^T A r_k) = r_k^T r_k - r_k^T r_k = 0$$
  Thus, $r_{k+1}$ is strictly orthogonal to $r_k$.

---

### Question 14: Explain the "zigzagging" pathology of Steepest Descent.
- **Answer**: When $\kappa_2(A) = \lambda_{\max}/\lambda_{\min} \gg 1$, the iso-contours of the quadratic functional $J(x)$ are extremely eccentric, needle-thin ellipsoidal valleys. Because consecutive search directions are strictly orthogonal ($r_{k+1} \perp r_k$), the trajectory bounces back and forth across the valley walls at $90^\circ$ angles rather than pointing along the valley floor toward the minimum $x^*$. The error contracts by factor $\frac{\kappa - 1}{\kappa + 1} \approx 1 - \frac{2}{\kappa} \to 1$, causing progress to stall.

---

### Question 15: How does Conjugate Gradient eliminate zigzagging? State its optimality property.
- **Answer**: Conjugate Gradient replaces Euclidean orthogonality with **$A$-conjugacy ($p_i^T A p_j = 0$ for $i \ne j$)**.
  - **Optimality**: In CG, minimizing $J(x)$ along direction $p_k$ does not degrade the minimization achieved along previous directions. At step $k$, $x_k$ is the exact minimizer of the error in the $A$-norm ($\|x - x^*\|_A$) over the entire $k$-dimensional Krylov subspace $x_0 + \mathcal{K}_k(A, r_0)$. It never repeats a direction, eliminating zigzagging and converging in at most $n$ steps in exact arithmetic.

---

### Question 16: Contrast the convergence bounds of Steepest Descent and Conjugate Gradient.
- **Answer**:
  - Steepest Descent: $\|e_k\|_A \le \left( \frac{\kappa - 1}{\kappa + 1} \right)^k \|e_0\|_A$.
  - Conjugate Gradient: $\|e_k\|_A \le 2 \left( \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1} \right)^k \|e_0\|_A$.
  - For $\kappa = 10,000$:
    - SD contraction factor: $\frac{9999}{10001} \approx 0.9998$ (requires $\sim 35,000$ iterations).
    - CG contraction factor: $\frac{100 - 1}{100 + 1} = \frac{99}{101} \approx 0.9802$ (requires $\sim 350$ iterations — **$100\times$ faster!**).

---

### Question 17: How does GMRES monitor convergence without explicitly calculating the solution vector $x_m$?
- **Answer**: In GMRES, the least-squares problem $\min \|\beta e_1 - \bar{H}_m y\|_2$ is transformed via accumulated Givens plane rotations $Q_m$ into:
  $$\|Q_m (\beta e_1 - \bar{H}_m y)\|_2^2 = \|g_m - R_m y\|_2^2 + |\gamma_{m+1}|^2$$
  Setting $R_m y = g_m$ sets the first term to zero.
  The minimum residual norm is **strictly equal to the scalar $|\gamma_{m+1}|$**:
  $$\|r_m\|_2 = |\gamma_{m+1}|$$
  Because $|\gamma_{m+1}|$ is computed as a byproduct of updating the Givens rotation at step $m$, the solver tracks the exact residual norm **WITHOUT solving $R_m y = g_m$ and WITHOUT forming $x_m = x_0 + V_m y$**, saving massive memory and arithmetic bandwidth.

---

### Question 18: Why is BiCGSTAB preferred over CGS for non-symmetric systems?
- **Answer**:
  - CGS squares the BiCG contraction polynomial ($r_j = \phi_j^2(A) r_0$) to eliminate $A^T$. However, squaring also squares the round-off errors and residual oscillations, leading to wild spikes and numerical overflow.
  - BiCGSTAB stabilizes the iteration by setting $r_j = \psi_j(A) \phi_j(A) r_0$, where $\psi_j(t) = (1 - \omega_j t) \psi_{j-1}(t)$. Parameter $\omega_j$ is computed at every step via an exact 1D Steepest Descent line search to minimize the residual norm $\|r_{j+1}\|_2$. This smooths out oscillations, maintains stability, avoids $A^T$, and requires only 2 SpMVs per iteration.

---

### Question 19: Why is Split Preconditioning mandatory for Conjugate Gradient, and how does PCG eliminate explicit factor solves?
- **Answer**: Conjugate Gradient mathematically requires an SPD matrix. If left preconditioning ($M^{-1} A$) is applied, $M^{-1} A$ is non-symmetric even if both $A$ and $M$ are symmetric! Split preconditioning factors $M = L L^T$ and solves $\tilde{A} \tilde{x} = \tilde{b}$ where $\tilde{A} = L^{-1} A L^{-T}$. $\tilde{A}$ is strictly symmetric positive definite.
  PCG eliminates explicit factorizations of $L$ by substituting auxiliary vector $z_k = M^{-1} r_k$. The entire algorithm is expressed using only original matrix $A$ and a single solve $M z = r$ per step.

---

### Question 20: What is the difference between ILU(0) and ILUT?
- **Answer**:
  - **ILU(0)**: Performs incomplete LU factorization where fill-in elements are discarded if $a_{ij} = 0$ in the original matrix $A$. Memory footprint is strictly identical to matrix $A$.
  - **ILUT($p, \tau$)**: A threshold-based incomplete LU. An element is retained only if its magnitude exceeds drop tolerance $\tau \cdot \|a_i\|$, and at most $p$ largest non-zero entries are kept per row. Provides higher accuracy at the cost of variable memory.

---

### Question 21: What is False Sharing in multi-core shared memory, and how is it prevented?
- **Answer**: False sharing occurs when two independent threads modify distinct logical variables that reside on the **same 64-byte physical cache line**. Whenever thread 0 writes to its variable, the cache coherency hardware marks the entire 64-byte line invalid in thread 1's cache, forcing continuous cache line invalidations and bus ping-ponging ("cache thrashing") despite zero logical data sharing.
  **Prevention**: Pad variables to 64 bytes (`alignas(64)`), or allocate thread-private arrays so distinct threads never share a cache line.

---

### Question 22: In Fortran, why is column-wise matrix-vector multiplication $2\times$ faster than row-wise access?
- **Answer**: Fortran stores 2D arrays in **Column-Major order** (adjacent elements of the same column are contiguous in memory).
  In column-wise access ($\sum_j A_{ji} x_j$), stepping along index $i$ accesses memory with **stride-1**. When $A_{1i}$ is fetched, the next 7 floats are automatically brought into the 64-byte L1 cache line, yielding $87.5\%$ cache hits.
  In row-wise access ($\sum_j A_{ij} x_j$), stepping along $j$ jumps by $M$ words ($M \times 8\text{ bytes}$). For large $M$, each access misses the cache, forcing continuous high-latency DRAM stalls.

---

### Question 23: State Amdahl’s Law and derive the maximum speedup on an infinite number of processors.
- **Answer**: Let sequential time be $T_s = \Phi_s + \Psi$, where $\Phi_s = f T_s$ is non-parallelizable and $\Psi = (1 - f) T_s$ is parallelizable.
  Parallel execution time: $T_p \ge f T_s + \frac{(1 - f) T_s}{p}$.
  Speedup:
  $$S(p) = \frac{T_s}{T_p} \le \frac{1}{f + \frac{1 - f}{p}}$$
  As $p \to \infty$:
  $$\mathbf{S_{\max} = \lim_{p \to \infty} S(p) = \frac{1}{f + 0} = \frac{1}{f}}$$
  If $5\%$ of a code is sequential ($f = 0.05$), the maximum possible speedup on any supercomputer is $\frac{1}{0.05} = 20\times$.

---

### Question 24: What is the Isoefficiency Function, and what does $W = O(p \log p)$ signify?
- **Answer**: The Isoefficiency Function $W = K \cdot T_o(W, p)$ (where $K = \frac{E}{1 - E}$) defines the rate at which problem size $W$ must increase as a function of processor count $p$ to keep parallel efficiency $E$ strictly constant.
  $W = O(p \log p)$ signifies a **highly scalable parallel algorithm** (e.g., parallel reduction or FFT). It means that when the processor count $p$ is doubled, the problem size $W$ only needs to increase by slightly more than double ($p \log p$) to maintain the exact same parallel efficiency.

---

### Question 25: Describe Cannon’s Algorithm for parallel matrix multiplication and state its communication complexity.
- **Answer**: Cannon’s algorithm multiplies $C = AB$ on a 2D Torus network of $\sqrt{p} \times \sqrt{p}$ processors:
  1. **Preskewing**: Row $i$ of $A$ is circularly shifted left by $i$; column $j$ of $B$ is circularly shifted up by $j$.
  2. **Multiply-Shift Loop ($\sqrt{p}$ stages)**:
     - Local block multiply: $C_{i, j} \mathrel{+}= A_{i, j} B_{i, j}$.
     - Shift $A$ left by 1 block; shift $B$ up by 1 block.
  - **Total Runtime**:
    $$\mathbf{T_p = \frac{2n^3}{p} + 2\sqrt{p} \, t_s + 2\frac{n^2}{\sqrt{p}} \, t_w}$$
  - **Isoefficiency**: $W = O(p^{1.5})$.

---

### Question 26: Compare the convergence and iteration scaling of Jacobi, Gauss-Seidel, SOR, Steepest Descent, and Conjugate Gradient on 2D Poisson systems (Slide 221).
- **Answer**:
  - **Empirical Iteration Comparison (Slide 221 at $N = 256$, residual $\approx 10^{-10}$)**:
    - **Jacobi**: $215,057$ iterations.
    - **Gauss-Seidel**: $107,466$ iterations ($\approx 2\times$ fewer than Jacobi, verifying $R_\infty(GS) = 2 R_\infty(J)$).
    - **Steepest Descent**: $247,067$ iterations (poor performance due to orthogonal zigzagging).
    - **SOR (at $\omega_{opt}$)**: $38,221$ iterations.
    - **Conjugate Gradient**: **$484$ iterations** ($444\times$ fewer than Jacobi, $79\times$ fewer than SOR!).
  - **Scaling with Mesh Refinement ($h = 1/N$)**:
    - Jacobi & Gauss-Seidel: $O(N^2) = O(h^{-2})$ iterations (doubling grid quadruples iterations).
    - SOR: $O(N) = O(h^{-1})$ iterations (doubling grid doubles iterations).
    - Conjugate Gradient: $O(\sqrt{\kappa(A)}) = O(N) = O(h^{-1})$ iterations, but with a vastly smaller constant factor.
  - **Floating Point Count per Iteration (Dense)**:
    - SOR: $n^2 + 4n$ FLOPs.
    - Steepest Descent: $n^2 + 3n$ FLOPs.
    - Conjugate Gradient: $n^2 + n$ FLOPs.

---

### Question 27: A parallel program with workload $W$ running on $p$ processors has total overhead $T_o = p^{3/4} W^{3/4}$. Find its isoefficiency function. If $W$ increases by factor $m$, by what factor must $p$ increase to preserve efficiency? (Slide 380).
- **Answer**:
  1. The isoefficiency relationship is $W = K \cdot T_o(W, p)$ where $K = \frac{E}{1 - E}$:
     $$W = K p^{3/4} W^{3/4} \implies W^{1/4} = K p^{3/4} \implies \mathbf{W = K^4 p^3 = O(p^3)}$$
  2. To find the processor growth required for a workload increase:
     $$p = \left(\frac{W}{K^4}\right)^{1/3} \propto W^{1/3}$$
     If $W$ increases by a factor of $m$ ($W_{\text{new}} = m W$):
     $$p_{\text{new}} \propto (m W)^{1/3} = m^{1/3} W^{1/3} = \mathbf{m^{1/3} p}$$
     The number of processors must increase by **$m^{1/3}$**.

---

## 3. Examination Archive, Predicted Questions & Model Solutions

For exam practice and review, refer to the two dedicated repositories:

1. 🎯 [**Predicted Exam Questions & Solutions Repository**](file:///home/thesupercd/Documents/HPSC_TutorialsNAssignments/ExamNotes/Likely_Exam_Questions_And_Solutions_Repository.md): A curated high-yield question bank across all 10 core modules with full analytical and numerical derivations, plus a complete 25-mark predicted model test paper.
2. 📝 [**Previous Years Examination Question Papers: Detailed Solutions & Analysis**](file:///home/thesupercd/Documents/HPSC_TutorialsNAssignments/ExamNotes/Previous_Years_Question_Papers_Solutions.md): Official past exam papers (2018 Midsem & 2025 Midsem) with 100% complete step-by-step solutions.

### Solved Past Examination Papers:
1. **[Mid Spring Semester Examination 2018 (22-02-2018)](file:///home/thesupercd/Documents/HPSC_TutorialsNAssignments/ExamNotes/Previous_Years_Question_Papers_Solutions.md#1-mid-spring-semester-examination-2018)**:
   - Part A (15 Marks):
     - Q1: Iteration matrix block lower-triangular decomposition, characteristic polynomial, eigenvalues, and divergence proof.
     - Q2: Convergence factor $\mu \equiv \rho(G)$, asymptotic rate $R_\infty = -\ln(\mu)$, and exact derivation of iteration count $k \ge \frac{-\ln \epsilon}{R_\infty}$.
     - Q3: Iterative scheme dependencies: $\kappa_2(A)$ (Steepest Descent), $\sqrt{\kappa_2(A)}$ (CG), and $\rho(G)$ (Jacobi, GS, SOR).
     - Q4: $3 \times 3$ strictly diagonally dominant non-symmetric matrix solvable by Gauss-Seidel but invalid for Conjugate Gradient.
     - Q5: Oblique Petrov-Galerkin projection ($r_m \perp \mathcal{K}_m(A^T, \tilde{r}_0)$) and residual oscillations / breakdowns in BiCG.
     - Q6: Preconditioning principles ($M^{-1} A \approx I$, spectral clustering, condition reduction, ILU/IC/SSOR).
   - Part B (10 Marks):
     - Q6: Diagonal storage (DIA) format reconstruction, full $6 \times 6$ linear system, analytical solution $x^* = [-8, 6.5, 13.5, 4, -7, 18]^T$, Gauss-Seidel in-place implementation, and relative residual stopping criteria.
2. **[Mid-Semester Examination 2025 (26-07-2025)](file:///home/thesupercd/Documents/HPSC_TutorialsNAssignments/ExamNotes/Previous_Years_Question_Papers_Solutions.md#2-mid-semester-examination-2025)**:
   - Part A (18 Marks):
     - Q1: Poisson FEM grid scaling: $\kappa([K]) = O(h^{-2}) = O(N) \implies$ iterations scale as $O(\sqrt{N}) \approx 3.16\text{--}4.5\times$; CPU time scaling with cache boundary ($1.0\text{ s} \to 11.62\text{ s}$); BiCGSTAB $2\times$ FLOP/SpMV overhead vs CG.
     - Q2: 2D Poisson $\nabla^2 T = \sin x \sin y$ central difference stencil, boundary node handling, and why direct solvers fail due to $O(N^3)$ fill-in vs iterative $O(N)$ sparse storage.
     - Q3: Spectral radius comparison ($\rho(G_1) = 0.87 < \rho(G_2) = 0.89 \implies G_1$ converges faster); non-singular $3 \times 3$ matrix failing Jacobi iteration ($\det(A)=5 \ne 0, \rho(G_J)=4 > 1$).
     - Q4: Sequential time $T_s = 80\text{ s}$, $80\%$ parallelizable: communication overhead extraction ($T_{\text{comm}}(8)=8\text{ s}$), linear communication scaling to $p=16$ ($T_p(16)=36\text{ s}$), parallel efficiency $E(16)=13.89\%$.
     - Q5: Architectural comparison table: Shared-Memory SMP vs Distributed-Memory Multicomputers.
   - Part B (7 Marks):
     - Q6: ELLPACK-ITPACK storage matrix reconstruction ($10 \times 3$), SOR update formula, Young's theorem $\rho(G_{SOR}) = \omega - 1$ for $\omega \ge \omega_{opt} \approx 1.31 \implies \rho(1.4)=0.40$ vs $\rho(1.6)=0.60$, iteration differences tracking $\delta^{(k+1)} \approx \rho \delta^{(k)}$ for steps 9, 10, 11; proof that $\omega = 1.4$ converges substantially faster.
