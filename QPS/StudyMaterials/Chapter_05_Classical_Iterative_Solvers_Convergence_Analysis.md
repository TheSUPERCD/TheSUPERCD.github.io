# Chapter 05: Classical Iterative Solvers & Convergence Theory

---

## 1. General Splitting Framework for Iterative Solvers

An **iterative solver** computes a sequence of successive approximations $\{x^{(0)}, x^{(1)}, x^{(2)}, \dots\}$ that converges toward the true solution $x^* = A^{-1} b$ as the iteration count $k \to \infty$.

### 1.1 The Canonical Matrix Splitting
Consider the linear system $Ax = b$, where $A \in \mathbb{R}^{n \times n}$ is non-singular.
We split matrix $A$ into:
$$A = M - N$$
where $M$ is a non-singular, easily invertible matrix (the **splitting matrix**), and $N = M - A$.

The original system $Ax = b$ can be rearranged as:
$$(M - N) x = b \implies M x = N x + b$$
This induces the fundamental fixed-point iteration:
$$M x^{(k+1)} = N x^{(k)} + b$$
Multiplying by $M^{-1}$:
$$\mathbf{x^{(k+1)} = M^{-1} N x^{(k)} + M^{-1} b \iff x^{(k+1)} = G x^{(k)} + f}$$
where:
- $\mathbf{G = M^{-1} N = M^{-1}(M - A) = I - M^{-1} A}$ is the **Iteration Matrix**.
- $\mathbf{f = M^{-1} b}$ is the modified right-hand side vector.

### 1.2 Additive Decomposition of Matrix $A$
Matrix $A$ is canonically decomposed into its diagonal, strictly lower, and strictly upper triangular parts:
$$A = D - E - F \quad (\text{or } A = D + L + U)$$
where:
- $D = \text{diag}(a_{11}, a_{22}, \dots, a_{nn})$ is the diagonal matrix.
- $-E$ (or $L$) is the strictly lower triangular matrix ($e_{ij} = -a_{ij}$ for $i > j$, zero otherwise).
- $-F$ (or $U$) is the strictly upper triangular matrix ($f_{ij} = -a_{ij}$ for $i < j$, zero otherwise).

---

## 2. Jacobi, Gauss-Seidel, and SOR Methods

### 2.1 The Jacobi Iteration Method
- **Core Principle**: Each unknown $x_i^{(k+1)}$ is updated using **only** values from the previous iteration step $k$. No newly computed values from step $k+1$ are used.
- **Component-Wise Formulation**:
  From the $i$-th equation $\sum_{j=1}^n a_{ij} x_j = b_i$:
  $$a_{ii} x_i + \sum_{j \ne i} a_{ij} x_j = b_i \implies \mathbf{x_i^{(k+1)} = \frac{b_i - \sum_{j \ne i} a_{ij} x_j^{(k)}}{a_{ii}}}$$
- **Matrix Splitting**:
  $$M = D, \quad N = E + F$$
  $$D x^{(k+1)} = (E + F) x^{(k)} + b$$
  $$\mathbf{G_J = D^{-1}(E + F) = I - D^{-1} A, \quad f_J = D^{-1} b}$$
- **HPC Characteristics**:
  - **Embarrassingly Parallel**: Every component $x_i^{(k+1)}$ can be computed completely independently of all other components at step $k+1$.
  - **Memory Overhead**: Requires two distinct storage vectors: `x_old` ($x^{(k)}$) and `x_new` ($x^{(k+1)}$). In-place updates are prohibited.

### 2.2 The Gauss-Seidel (G-S) Iteration Method
- **Core Principle**: As soon as a new value $x_j^{(k+1)}$ is computed (for $j < i$), it **immediately replaces** $x_j^{(k)}$ in the calculation of subsequent components $x_i^{(k+1)}$ within the same sweep.
- **Component-Wise Formulation**:
  $$\sum_{j < i} a_{ij} x_j^{(k+1)} + a_{ii} x_i^{(k+1)} + \sum_{j > i} a_{ij} x_j^{(k)} = b_i$$
  $$\mathbf{x_i^{(k+1)} = \frac{b_i - \sum_{j < i} a_{ij} x_j^{(k+1)} - \sum_{j > i} a_{ij} x_j^{(k)}}{a_{ii}}}$$
- **Matrix Splitting**:
  $$M = D - E, \quad N = F$$
  $$(D - E) x^{(k+1)} = F x^{(k)} + b$$
  $$\mathbf{G_{GS} = (D - E)^{-1} F, \quad f_{GS} = (D - E)^{-1} b}$$
- **HPC Characteristics**:
  - **Memory Efficient**: Updates are performed strictly **in-place** within a single array `x`. Old values are only saved if needed for convergence checking.
  - **Parallelism Bottleneck**: Introducing immediate dependencies creates a sequential chain of computations across rows. Natural parallelism is destroyed unless reordered (e.g., Red-Black ordering).

### 2.3 The Successive Over-Relaxation (SOR) Method
- **Core Motivation**:
  Notice the difference between consecutive iterates in Gauss-Seidel:
  $$x^{(k+1)} - x^{(k)} = G x^{(k)} + f - x^{(k)} = (G - I) x^{(k)} + f$$
  If this directional correction step is accelerated (amplified) by an extrapolation factor $\omega > 1$, convergence can be vastly accelerated!
- **Component-Wise Formulation**:
  First compute the tentative Gauss-Seidel target value $x_i^{(k+1)*}$:
  $$x_i^{(k+1)*} = \frac{b_i - \sum_{j < i} a_{ij} x_j^{(k+1)} - \sum_{j > i} a_{ij} x_j^{(k)}}{a_{ii}}$$
  Then update $x_i^{(k+1)}$ as a weighted linear combination:
  $$\mathbf{x_i^{(k+1)} = x_i^{(k)} + \omega \left( x_i^{(k+1)*} - x_i^{(k)} \right) = (1 - \omega) x_i^{(k)} + \omega x_i^{(k+1)*}}$$
- **Matrix Formulation**:
  Multiplying by $a_{ii}$ and assembling:
  $$a_{ii} x_i^{(k+1)} + \omega \sum_{j < i} a_{ij} x_j^{(k+1)} = (1 - \omega) a_{ii} x_i^{(k)} - \omega \sum_{j > i} a_{ij} x_j^{(k)} + \omega b_i$$
  $$(D - \omega E) x^{(k+1)} = \left[ (1 - \omega) D + \omega F \right] x^{(k)} + \omega b$$
  $$\mathbf{G_{SOR} = (D - \omega E)^{-1} \left[ (1 - \omega) D + \omega F \right]}$$
  $$\mathbf{f_{SOR} = \omega (D - \omega E)^{-1} b}$$
- **Role of the Relaxation Parameter $\omega$**:
  - **$\omega = 1$**: Standard Gauss-Seidel iteration.
  - **$1 < \omega < 2$**: **Over-relaxation** (accelerates convergence for slow, stiff systems).
  - **$0 < \omega < 1$**: **Under-relaxation** (stabilizes oscillating or diverging non-linear systems).
  - **$\omega \ge 2$ or $\omega \le 0$**: The iteration is guaranteed to **diverge**!

---

## 3. Rigorous Convergence Analysis

### 3.1 Error Propagation and the Spectral Radius
Let $x^*$ denote the exact solution ($Ax^* = b \iff x^* = G x^* + f$).
Define the error vector at iteration step $k$:
$$e^{(k)} = x^{(k)} - x^*$$

Subtracting the fixed-point identity $x^* = G x^* + f$ from $x^{(k+1)} = G x^{(k)} + f$:
$$x^{(k+1)} - x^* = G (x^{(k)} - x^*) \implies \mathbf{e^{(k+1)} = G e^{(k)}}$$
By mathematical induction:
$$\mathbf{e^{(k)} = G^k e^{(0)}}$$

#### Definition: Spectral Radius
The **spectral radius** $\rho(G)$ of a matrix $G$ is the maximum absolute value of its eigenvalues:
$$\rho(G) = \max_{1 \le i \le n} |\lambda_i(G)|$$
By Gelfand's formula, for any induced matrix norm $\|\cdot\|$:
$$\lim_{k \to \infty} \|G^k\|^{1/k} = \rho(G)$$

### 3.2 Master Convergence Theorem
> **Theorem**: The iterative method $x^{(k+1)} = G x^{(k)} + f$ converges to the unique solution $x^* = (I - G)^{-1} f$ for **any** initial guess $x^{(0)}$ if and only if:
> $$\mathbf{\rho(G) < 1}$$

#### Complete Mathematical Proof
1. **Sufficiency ($\rho(G) < 1 \implies$ Convergence)**:
   By the Jordan canonical form, $G = P J P^{-1}$, where $J = \text{diag}(J_1, \dots, J_m)$ consists of Jordan blocks:
   $$J_i = \begin{bmatrix} \lambda_i & 1 & 0 & \dots \\ 0 & \lambda_i & 1 & \dots \\ \vdots & \vdots & \ddots & 1 \\ 0 & 0 & \dots & \lambda_i \end{bmatrix}$$
   Raising to power $k$:
   $$J_i^k = \begin{bmatrix} \lambda_i^k & \binom{k}{1} \lambda_i^{k-1} & \dots & \binom{k}{p-1} \lambda_i^{k-p+1} \\ 0 & \lambda_i^k & \dots & \binom{k}{p-2} \lambda_i^{k-p+2} \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & \lambda_i^k \end{bmatrix}$$
   If $|\lambda_i| \le \rho(G) < 1$, then $\lim_{k \to \infty} \binom{k}{r} \lambda_i^{k-r} = 0$ for any fixed $r$.
   Therefore, $\lim_{k \to \infty} J^k = \mathbf{0}$, which implies:
   $$\lim_{k \to \infty} G^k = P \left( \lim_{k \to \infty} J^k \right) P^{-1} = \mathbf{0}$$
   Consequently:
   $$\lim_{k \to \infty} e^{(k)} = \lim_{k \to \infty} G^k e^{(0)} = \mathbf{0} \implies \lim_{k \to \infty} x^{(k)} = x^*$$

2. **Neumann Series Convergence**:
   When $\rho(G) < 1$, the matrix $(I - G)$ is guaranteed to be non-singular, and its inverse is given by the absolutely convergent **Neumann series**:
   $$(I - G)^{-1} = \sum_{j=0}^\infty G^j$$
   Thus, $x^* = (I - G)^{-1} f$ is well-defined and unique.

3. **Necessity (Convergence $\implies \rho(G) < 1$)**:
   Suppose $\rho(G) \ge 1$. Let $v$ be an eigenvector corresponding to eigenvalue $|\lambda| \ge 1$, so $G v = \lambda v$.
   If the initial error is chosen along this eigenvector ($e^{(0)} = v$):
   $$e^{(k)} = G^k v = \lambda^k v$$
   Since $|\lambda|^k \ge 1$, $\|e^{(k)}\| = |\lambda|^k \|v\| \not\to 0$. The iteration fails to converge. $\blacksquare$

### 3.3 Convergence Rate and Computational Iteration Estimates
The asymptotic behavior of the error norm satisfies:
$$\lim_{k \to \infty} \frac{\|e^{(k+1)}\|}{\|e^{(k)}\|} = \rho(G)$$
- **Asymptotic Rate of Convergence ($\tau$ or $R_\infty$)**:
  $$\mathbf{\tau = -\ln \rho(G)}$$
- **Number of Iterations Required ($k_\epsilon$)**:
  To reduce the initial error by a factor of $\epsilon$ (i.e., $\|e^{(k)}\| / \|e^{(0)}\| \le \epsilon$):
  $$\|e^{(k)}\| \approx [\rho(G)]^k \|e^{(0)}\| \le \epsilon \|e^{(0)}\|$$
  $$k \ln \rho(G) \le \ln \epsilon \implies -k \tau \le \ln \epsilon \implies \mathbf{k \ge \frac{-\ln \epsilon}{\tau} = \frac{\ln(1/\epsilon)}{-\ln \rho(G)}}$$
- **Crucial Exam Takeaway**: As $\rho(G) \to 1$, $\ln \rho(G) \approx -(1 - \rho(G)) \implies k \approx \frac{\ln(1/\epsilon)}{1 - \rho(G)}$. A tiny increase in spectral radius causes an explosive increase in iteration count!

> [!IMPORTANT]
> **Past Exam Focus on Spectral Radius & Convergence**:
> - **2018 Mid-Spring Q1**: Given a $4 \times 4$ block triangular iteration matrix $G$, find its eigenvalues, determine $\rho(G) \approx 2.788 > 1$, and prove unconditional divergence:  
>   👉 [**2018 Exam Q1 Solution**](file:///home/thesupercd/Documents/HPSC_TutorialsNAssignments/ExamNotes/Previous_Years_Question_Papers_Solutions.md#question-1-3-marks).
> - **2018 Mid-Spring Q2**: Formal definitions of convergence factor $\mu = \rho(G)$, asymptotic rate $R_\infty$, and exact derivation of iteration count $k$:  
>   👉 [**2018 Exam Q2 Solution**](file:///home/thesupercd/Documents/HPSC_TutorialsNAssignments/ExamNotes/Previous_Years_Question_Papers_Solutions.md#question-2-3-marks).
> - **2025 Midsem Q3**: Comparing convergence speed for $\rho(G_1) = 0.87$ vs $\rho(G_2) = 0.89$, and constructing a non-singular $3 \times 3$ matrix for which Jacobi diverges:  
>   👉 [**2025 Exam Q3 Solution**](file:///home/thesupercd/Documents/HPSC_TutorialsNAssignments/ExamNotes/Previous_Years_Question_Papers_Solutions.md#question-3-3-marks).
> - **2025 Midsem Q6**: Tracing successive SOR difference contractions $\delta^{(k+1)} \approx \rho(G_{SOR}) \delta^{(k)}$ for $\omega = 1.4$ vs $1.6$:  
>   👉 [**2025 Exam Q6 Solution**](file:///home/thesupercd/Documents/HPSC_TutorialsNAssignments/ExamNotes/Previous_Years_Question_Papers_Solutions.md#question-6-part-b-7-marks).

---

## 4. Diagonal Dominance and Convergence Theorems

### 4.1 Classes of Diagonally Dominant Matrices
1. **Weakly Diagonally Dominant**:
   $$|a_{ii}| \ge \sum_{j \ne i} |a_{ij}| \quad \text{for all } i = 1, \dots, n$$
2. **Strictly Diagonally Dominant (SDD)**:
   $$|a_{ii}| > \sum_{j \ne i} |a_{ij}| \quad \text{for all } i = 1, \dots, n$$
3. **Irreducibly Diagonally Dominant (IDD)**:
   - Matrix $A$ is irreducible (its directed graph is strongly connected; no row/column permutations can make it block upper triangular).
   - $A$ is weakly diagonally dominant for all rows.
   - Strict inequality holds for **at least one row** $k$: $|a_{kk}| > \sum_{j \ne k} |a_{kj}|$.

### 4.2 Proof of Jacobi Convergence for SDD Matrices (Gershgorin Circle Theorem)
Let $A$ be Strictly Diagonally Dominant. The Jacobi iteration matrix is:
$$G_J = -D^{-1}(E + F) = [g_{ij}] \quad \text{where } g_{ii} = 0, \text{ and } g_{ij} = -\frac{a_{ij}}{a_{ii}} \text{ for } i \ne j$$
By the Gershgorin Circle Theorem, every eigenvalue $\lambda$ of $G_J$ lies in at least one disk in the complex plane:
$$|\lambda - g_{ii}| \le \sum_{j \ne i} |g_{ij}| \implies |\lambda| \le \sum_{j \ne i} \left| -\frac{a_{ij}}{a_{ii}} \right| = \frac{\sum_{j \ne i} |a_{ij}|}{|a_{ii}|}$$
Since $A$ is SDD:
$$\sum_{j \ne i} |a_{ij}| < |a_{ii}| \implies \frac{\sum_{j \ne i} |a_{ij}|}{|a_{ii}|} < 1 \quad \text{for all } i$$
Therefore:
$$\rho(G_J) = \max_i |\lambda_i| \le \max_i \left( \frac{\sum_{j \ne i} |a_{ij}|}{|a_{ii}|} \right) < 1$$
**Conclusion**: Both Jacobi and Gauss-Seidel iterations **unconditionally converge** for any Strictly or Irreducibly Diagonally Dominant matrix, for any initial guess $x^{(0)}$!

### 4.3 Regular Splitting Theorem (Varga)
- A splitting $A = M - N$ is called a **regular splitting** if $M^{-1} \ge 0$ (all entries non-negative) and $N \ge 0$.
- **Theorem**: If $M, N$ is a regular splitting of matrix $A$, then:
  $$\mathbf{\rho(M^{-1}N) < 1 \iff A^{-1} \ge 0}$$

---

## 5. Successive Over-Relaxation (SOR) Theory

### 5.1 The Ostrowski-Reich Theorem
> **Theorem**: If $A$ is a real Symmetric Positive Definite (SPD) matrix with positive diagonal entries ($a_{ii} > 0$), then the SOR method converges for any initial vector $x^{(0)}$ if and only if:
> $$\mathbf{0 < \omega < 2}$$

### 5.2 Young's Theorem for Consistently Ordered Matrices
Many structured sparse matrices arising from PDEs (such as tridiagonal and pentadiagonal 5-point stencils) possess the property of being **consistently ordered and 2-cyclic** (bipartite coloring).

> **Theorem (David M. Young)**:
> Let $A$ be a consistently ordered 2-cyclic matrix. Let $\mu \in \sigma(G_J)$ be an eigenvalue of the Jacobi iteration matrix, and $\lambda \in \sigma(G_{SOR})$ be an eigenvalue of the SOR iteration matrix. Then $\lambda$ and $\mu$ satisfy the fundamental relation:
> $$\mathbf{(\lambda + \omega - 1)^2 = \lambda \omega^2 \mu^2}$$

#### Consequence for Gauss-Seidel ($\omega = 1$)
Setting $\omega = 1$ in Young's relation:
$$(\lambda + 1 - 1)^2 = \lambda (1)^2 \mu^2 \implies \lambda^2 = \lambda \mu^2 \implies \mathbf{\lambda = \mu^2}$$
Taking the maximum modulus:
$$\mathbf{\rho(G_{GS}) = \left[ \rho(G_J) \right]^2}$$
- **Profound Result**: Since $\rho(G_J) < 1$, squaring it makes it strictly smaller:
  $$\rho(G_{GS}) < \rho(G_J)$$
  $$\tau_{GS} = -\ln \rho(G_{GS}) = -\ln [\rho(G_J)]^2 = -2 \ln \rho(G_J) = 2 \tau_J$$
- **Exam Rule**: **Gauss-Seidel converges exactly TWICE as fast as Jacobi** for consistently ordered matrices (requires half the number of iterations)!

### 5.3 Derivation of the Optimal Relaxation Factor ($\omega_{opt}$)
We seek the value of $\omega \in (1, 2)$ that minimizes the spectral radius $\rho(G_{SOR})$.
From $(\lambda + \omega - 1)^2 = \lambda \omega^2 \mu^2$, solving the quadratic for $\lambda$:
$$\lambda - \omega \mu \sqrt{\lambda} + (\omega - 1) = 0$$
$$\sqrt{\lambda} = \frac{\omega \mu \pm \sqrt{\omega^2 \mu^2 - 4(\omega - 1)}}{2}$$
1. If $\omega^2 \mu^2 - 4(\omega - 1) > 0$, the roots are real and distinct; one root increases as $\omega$ increases.
2. If $\omega^2 \mu^2 - 4(\omega - 1) < 0$, the roots are complex conjugates with magnitude:
   $$|\lambda| = \omega - 1$$
   Notice that $|\lambda| = \omega - 1$ is completely **independent of $\mu$**!
3. The minimum of $\max |\lambda|$ occurs at the discriminant inflection point where the radical vanishes for the largest Jacobi eigenvalue $\mu = \rho(G_J)$:
   $$\omega^2 \rho(G_J)^2 - 4(\omega - 1) = 0$$
   $$\rho(G_J)^2 \omega^2 - 4\omega + 4 = 0$$
   Solving for $\omega$ via the quadratic formula:
   $$\omega = \frac{4 \pm \sqrt{16 - 16 \rho(G_J)^2}}{2 \rho(G_J)^2} = \frac{4 [1 \pm \sqrt{1 - \rho(G_J)^2}]}{2 \rho(G_J)^2} = \frac{2 [1 - \sqrt{1 - \rho(G_J)^2}]}{\rho(G_J)^2}$$
   Multiplying numerator and denominator by $[1 + \sqrt{1 - \rho(G_J)^2}]$:
   $$\omega = \frac{2 [1 - (1 - \rho(G_J)^2)]}{\rho(G_J)^2 [1 + \sqrt{1 - \rho(G_J)^2}]} = \mathbf{\frac{2}{1 + \sqrt{1 - \rho(G_J)^2}}}$$

$$\mathbf{\omega_{opt} = \frac{2}{1 + \sqrt{1 - \rho(G_J)^2}}}$$
At this optimal parameter, the SOR spectral radius is:
$$\mathbf{\rho(G_{SOR, opt}) = \omega_{opt} - 1}$$

---

## 6. Numerical Experiments: 1D Poisson Benchmark Data

The lecture slides document explicit numerical experiments solving the 1D Poisson boundary value problem:
$$\frac{d^2 T}{dx^2} = 0, \quad 0 \le x \le 1, \quad T(0) = 0, T(1) = 1$$
Discretized using central differences into tridiagonal matrices of varying dimensions. For example, for a $7 \times 7$ interior grid (*Slide 71*):
$$\begin{bmatrix}
-2 & 1 & 0 & 0 & 0 & 0 & 0 \\
1 & -2 & 1 & 0 & 0 & 0 & 0 \\
0 & 1 & -2 & 1 & 0 & 0 & 0 \\
0 & 0 & 1 & -2 & 1 & 0 & 0 \\
0 & 0 & 0 & 1 & -2 & 1 & 0 \\
0 & 0 & 0 & 0 & 1 & -2 & 1 \\
0 & 0 & 0 & 0 & 0 & 1 & -2
\end{bmatrix}
\begin{bmatrix}
T_2 \\ T_3 \\ T_4 \\ T_5 \\ T_6 \\ T_7 \\ T_8
\end{bmatrix}
=
\begin{bmatrix}
0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ -1
\end{bmatrix}$$
Target convergence tolerance: $\epsilon = 10^{-6}$.

### 6.1 Quantitative Comparison: Jacobi vs. Gauss-Seidel (Slide 73)
| Matrix Size | Method | Spectral Radius $\rho(G)$ | Convergence Rate $\tau = -\ln \rho(G)$ | Iterations for $\epsilon = 10^{-6}$ |
| :---: | :---: | :---: | :---: | :---: |
| **$10 \times 10$** | **Jacobi** | $0.94632$ | $0.05517$ | $200$ |
| | **Gauss-Seidel** | $0.89533$ | $0.11034$ | $106$ |
| **$20 \times 20$** | **Jacobi** | $0.98645$ | $0.01364$ | $709$ |
| | **Gauss-Seidel** | $0.97309$ | $0.02728$ | $374$ |
| **$40 \times 40$** | **Jacobi** | $0.99661$ | $0.00340$ | $2435$ |
| | **Gauss-Seidel** | $0.99322$ | $0.00680$ | $1291$ |

### 6.2 Complete Parametric Sweep: GS vs. SOR across Relaxation Factors (Slide 78)
| Metric | Equations ($N$) | Jacobi | Gauss-Seidel ($\omega=1$) | SOR ($\omega=1.2$) | SOR ($\omega=1.5$) | SOR ($\omega=1.8$) | SOR ($\omega=1.9$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Iterations** | **10** | 200 | 106 | 72 | **30** ($\omega_{opt}$) | 81 | 173 |
| (to reach $\epsilon = 10^{-6}$) | **20** | 709 | 374 | 258 | 133 | **83** ($\omega_{opt}$) | 170 |
| | **40** | 2435 | 1291 | 899 | 479 | **160** ($\omega_{opt}$) | 174 |
| **Spectral Radius** | **10** | 0.94632 | 0.89553 | 0.84206 | 0.59423 | **0.80000** | **0.90000** |
| $\rho(G)$ | **20** | 0.98645 | 0.97309 | 0.95956 | 0.91676 | **0.80000** | **0.90000** |
| | **40** | 0.99661 | 0.99322 | 0.98983 | 0.97952 | 0.92950 | **0.90000** |
| **Convergence Rate** | **10** | 0.05517 | 0.11034 | 0.17190 | 0.52049 | **0.22314** | **0.10536** |
| $\tau = -\ln \rho(G)$ | **20** | 0.01364 | 0.02728 | 0.04128 | 0.08691 | **0.22314** | **0.10536** |
| | **40** | 0.00340 | 0.00680 | 0.01022 | 0.02068 | 0.07311 | **0.10536** |

### 6.3 Crucial Exam Insights from the SOR Experimental Data
1. **Existence of a Unique Optimal Factor ($\omega_{opt}$)**:
   - For $N = 10$: Minimum iterations occur at $\omega \approx 1.5$ (collapsed from 106 down to **30 iterations**, a $3.5\times$ speedup over GS).
   - For $N = 20$: Minimum iterations occur at $\omega \approx 1.8$ (collapsed from 374 down to **83 iterations**, a $4.5\times$ speedup over GS).
   - For $N = 40$: Minimum iterations occur at $\omega \approx 1.8$ (collapsed from 1291 down to **160 iterations**, an **$8.1\times$ speedup over GS** and **$15.2\times$ speedup over Jacobi**!).
2. **The Constant Spectral Radius Property for $\omega \ge \omega_{opt}$**:
   - Notice that for $\omega = 1.8$ and $\omega = 1.9$, whenever $\omega > \omega_{opt}$, the spectral radius is **strictly independent of the matrix dimension $N$**!
     $$\rho(G_{SOR}) = \omega - 1$$
   - Specifically:
     - At $\omega = 1.8$: $\rho(G) = 1.8 - 1 = \mathbf{0.80000}$ for both $N=10$ and $N=20$.
     - At $\omega = 1.9$: $\rho(G) = 1.9 - 1 = \mathbf{0.90000}$ for $N=10, 20$, and $40$.
   - Consequently, the asymptotic convergence rate is constant: $\tau = -\ln(0.8) = \mathbf{0.22314}$, and $\tau = -\ln(0.9) = \mathbf{0.10536}$!
3. **Over-Relaxation Penalty**: Increasing $\omega$ beyond $\omega_{opt}$ causes the complex conjugate eigenvalue magnitudes ($|\lambda| = \omega - 1$) to swell linearly towards 1, degrading convergence rapidly.


## 7. Practical Implementation Details & Stopping Criteria

### 7.1 Stopping Criteria
Three common stopping criteria are used in production codes:
1. **Absolute Coordinate Difference (Chebyshev $L_\infty$ Norm)**:
   $$\|x^{(k+1)} - x^{(k)}\|_\infty = \max_{1 \le i \le n} |x_i^{(k+1)} - x_i^{(k)}| < \epsilon$$
2. **Euclidean Residual Norm ($L_2$ Norm)**:
   $$\|r^{(k)}\|_2 = \|b - Ax^{(k)}\|_2 < \epsilon$$
3. **Relative Residual Norm (Dimensionless & Scale-Independent)**:
   $$\mathbf{\frac{\|b - Ax^{(k)}\|_2}{\|b - Ax^{(0)}\|_2} < \epsilon_{tol} \quad \text{or} \quad \frac{\|r^{(k)}\|_2}{\|b\|_2} < \epsilon_{tol}}}$$
   *(The relative residual norm is the gold standard in scientific software).*

### 7.2 Implementation Comparison Matrix
| Feature | Jacobi Method | Gauss-Seidel Method | SOR Method |
| :--- | :--- | :--- | :--- |
| **Memory Footprint** | $2n$ vectors (`x_old`, `x_new`) | $n$ vector (`x` in-place) | $n$ vector (`x` in-place) |
| **Parallel Suitability** | High (naturally decoupled) | Low (serial dependency chain) | Low (serial dependency chain) |
| **Convergence Speed** | Slowest ($\rho$) | $2\times$ faster than Jacobi ($\rho^2$) | Extremely fast with $\omega_{opt}$ |
| **Storage Scheme Friendly**| Compatible with CSR/DIA | Sensitive to update order | Sensitive to update order |
