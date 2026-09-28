# Chapter 11: Krylov Subspace Methods III: GMRES, BiCG, CGS & BiCGSTAB

---

## 1. The Challenge of Non-Symmetric Linear Systems

When mathematical models incorporate convective transport, velocity-pressure coupling, or non-reciprocal physics (e.g., Navier-Stokes equations, convection-diffusion, magnetohydrodynamics), the resulting discretization matrix is **non-symmetric**:
$$A \ne A^T$$

### 1.1 Why Standard Conjugate Gradient Fails
1. **Loss of Energy Norm**: If $A$ is non-symmetric, the quadratic functional $J(x) = \frac{1}{2} x^T A x - b^T x$ has gradient $\nabla J(x) = \frac{1}{2}(A + A^T)x - b \ne Ax - b$. Minimizing $J(x)$ no longer solves $Ax = b$.
2. **Loss of Short Recurrences**: For non-symmetric matrices, the Faber-Manteuffel Theorem proves that one cannot simultaneously have an optimal Krylov projection method and short (3-term) vector recurrences!
3. **The Algorithmic Dilemma**:
   - **Approach 1 (Full Orthogonalization / GMRES)**: Retains optimal minimization of the residual ($\|r_m\|_2$), but must maintain long recurrences spanning all previous basis vectors ($O(m)$ work and memory per step).
   - **Approach 2 (Biorthogonalization / BiCG, BiCGSTAB)**: Retains short 3-term recurrences, but sacrifices the guaranteed monotonic minimization property.

---

## 2. Generalized Minimal Residual Method (GMRES)

Introduced by Yousef Saad and Martin Schultz (1986), **GMRES** is the gold standard projection solver for general non-symmetric, non-singular linear systems.

### 2.1 Mathematical Formulation
- **Search Subspace**: $\mathcal{K}_m(A, r_0) = \text{span}\{r_0, A r_0, \dots, A^{m-1} r_0\}$.
- **Constraint Subspace**: $\mathcal{L}_m = A \mathcal{K}_m(A, r_0)$.
- **Minimization Objective**:
  Find $x_m \in x_0 + \mathcal{K}_m$ such that the Euclidean norm of the residual is minimized:
  $$\mathbf{\|r_m\|_2 = \|b - A x_m\|_2 = \min_{x \in x_0 + \mathcal{K}_m} \|b - Ax\|_2}$$

### 2.2 Transformation to a Hessenberg Least-Squares Problem
From Arnoldi’s method (Chapter 09):
- Orthonormal basis $V_m = [v_1, \dots, v_m] \in \mathbb{R}^{n \times m}$ for $\mathcal{K}_m$.
- Extended upper Hessenberg matrix $\bar{H}_m \in \mathbb{R}^{(m+1) \times m}$.
- First Arnoldi relation:
  $$A V_m = V_{m+1} \bar{H}_m$$
Any vector in the affine space $x \in x_0 + \mathcal{K}_m$ can be parameterized by $y \in \mathbb{R}^m$:
$$x = x_0 + V_m y$$
The residual vector is:
$$r = b - A (x_0 + V_m y) = (b - A x_0) - A V_m y = r_0 - V_{m+1} \bar{H}_m y$$
Recall that $v_1 = r_0 / \beta$, where $\beta = \|r_0\|_2$, so $r_0 = \beta v_1 = V_{m+1} (\beta e_1)$, where $e_1 = [1, 0, \dots, 0]^T \in \mathbb{R}^{m+1}$.
$$r = V_{m+1} (\beta e_1) - V_{m+1} \bar{H}_m y = V_{m+1} (\beta e_1 - \bar{H}_m y)$$

Now evaluate the Euclidean norm $\|r\|_2$.
Because the columns of $V_{m+1}$ are strictly orthonormal ($V_{m+1}^T V_{m+1} = I_{m+1}$):
$$\|r\|_2^2 = r^T r = (\beta e_1 - \bar{H}_m y)^T (V_{m+1}^T V_{m+1}) (\beta e_1 - \bar{H}_m y) = \|\beta e_1 - \bar{H}_m y\|_2^2$$
$$\mathbf{\|b - Ax\|_2 = \|\beta e_1 - \bar{H}_m y\|_2}$$
**The Core GMRES Reduction**:
The original massive $n$-dimensional minimization is transformed into a small **$(m+1) \times m$ linear least-squares problem**:
$$\mathbf{y_m = \arg\min_{y \in \mathbb{R}^m} \|\beta e_1 - \bar{H}_m y\|_2}$$
where $m \ll n$ (typically $m \in [20, 50]$).

---

### 2.3 Solving the Least-Squares Problem via Givens Plane Rotations
To solve $\min_y \|\beta e_1 - \bar{H}_m y\|_2$ progressively, we transform $\bar{H}_m$ into upper triangular form using a sequence of **Givens plane rotations** $\Omega_1, \Omega_2, \dots, \Omega_m$.

#### 1. Definition of the Givens Rotation Matrix ($\Omega_i$)
At step $i$, to eliminate the subdiagonal entry $h_{i+1, i}$ using diagonal entry $h_{ii}^{(i-1)}$:
$$\Omega_i = \begin{bmatrix} 1 & & & & & \\ & \ddots & & & & \\ & & c_i & s_i & & \\ & & -s_i & c_i & & \\ & & & & \ddots & \\ & & & & & 1 \end{bmatrix} \quad \text{at rows } i \text{ and } i+1$$
where:
$$c_i = \frac{h_{ii}^{(i-1)}}{\sqrt{(h_{ii}^{(i-1)})^2 + (h_{i+1, i})^2}}, \quad s_i = \frac{h_{i+1, i}}{\sqrt{(h_{ii}^{(i-1)})^2 + (h_{i+1, i})^2}}, \quad c_i^2 + s_i^2 = 1$$

#### 2. Orthogonal Transformation of the System
Define the accumulated orthogonal matrix $Q_m = \Omega_m \Omega_{m-1} \dots \Omega_1 \in \mathbb{R}^{(m+1) \times (m+1)}$.
Multiplying $\bar{H}_m$ and $\beta e_1$:
$$\bar{R}_m = Q_m \bar{H}_m = \begin{bmatrix} R_m \\ \mathbf{0}^T \end{bmatrix} = \begin{bmatrix} r_{11} & r_{12} & \dots & r_{1m} \\ 0 & r_{22} & \dots & r_{2m} \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & r_{mm} \\ 0 & 0 & \dots & 0 \end{bmatrix}$$
$$\bar{g}_m = Q_m (\beta e_1) = \begin{bmatrix} \gamma_1 \\ \gamma_2 \\ \vdots \\ \gamma_m \\ \gamma_{m+1} \end{bmatrix} = \begin{bmatrix} g_m \\ \gamma_{m+1} \end{bmatrix}$$
where $R_m \in \mathbb{R}^{m \times m}$ is square upper triangular, and $g_m \in \mathbb{R}^m$.

#### 3. Complete Solution and Residual Monitoring (Slides 233–234)
Because $Q_m$ is an orthogonal matrix, it preserves the Euclidean norm:
$$\|\beta e_1 - \bar{H}_m y\|_2^2 = \|Q_m (\beta e_1 - \bar{H}_m y)\|_2^2 = \|\bar{g}_m - \bar{R}_m y\|_2^2$$
Partitioning into upper $m$ rows and the final row:
$$\|\bar{g}_m - \bar{R}_m y\|_2^2 = \left\| \begin{bmatrix} g_m - R_m y \\ \gamma_{m+1} - 0 \end{bmatrix} \right\|_2^2 = \mathbf{\|g_m - R_m y\|_2^2 + |\gamma_{m+1}|^2}$$
Since $R_m$ is non-singular, we can set the first term to zero by solving:
$$\mathbf{R_m y_m = g_m \implies y_m = R_m^{-1} g_m}$$
Then the residual norm is:
$$\mathbf{\min_{y \in \mathbb{R}^m} \|\beta e_1 - \bar{H}_m y\|_2 = |\gamma_{m+1}|}$$

#### Decisive Exam Takeaway: In-Progress Residual Monitoring
The Euclidean norm of the residual at step $m$ is **strictly equal to $|\gamma_{m+1}|$**:
$$\mathbf{\|r_m\|_2 = |\gamma_{m+1}|}$$
**We do NOT need to compute $y_m$ or form $x_m$ to monitor convergence!** We simply check if $|\gamma_{m+1}| < \epsilon$. Only after convergence is reached do we perform the single triangular back-solve $R_m y_m = g_m$ and compute $x_m = x_0 + V_m y_m$.

---

### 2.4 Complete GMRES Algorithm
```
Algorithm: Full GMRES (Generalized Minimal Residual)
Input: Matrix A, vector b, initial guess x_0, dimension limit m, tolerance ε
Output: Solution vector x_m

1. Initialize:
   r_0 = b - A * x_0
   β = norm(r_0, 2)
   v_1 = r_0 / β
   g = [β, 0, ..., 0]^T of length m+1

2. For j = 1 to m do:
       // Arnoldi Modified Gram-Schmidt step
       w = A * v_j
       for i = 1 to j do
           h[i, j] = v_i^T * w
           w = w - h[i, j] * v_i
       end for
       h[j+1, j] = norm(w, 2)
       v_{j+1} = w / h[j+1, j]

       // Apply existing Givens rotations to column j
       for i = 1 to j-1 do
           temp = c[i] * h[i, j] + s[i] * h[i+1, j]
           h[i+1, j] = -s[i] * h[i, j] + c[i] * h[i+1, j]
           h[i, j] = temp
       end for

       // Compute new Givens rotation to zero out h[j+1, j]
       c[j] = h[j, j] / sqrt(h[j, j]^2 + h[j+1, j]^2)
       s[j] = h[j+1, j] / sqrt(h[j, j]^2 + h[j+1, j]^2)

       // Eliminate subdiagonal entry
       h[j, j] = c[j] * h[j, j] + s[j] * h[j+1, j]
       h[j+1, j] = 0

       // Apply Givens rotation to RHS vector g
       g[j+1] = -s[j] * g[j]
       g[j] = c[j] * g[j]

       // Check convergence WITHOUT forming x!
       if |g[j+1]| <= ε then
           m = j
           break
       end if
   end for

3. Form approximate solution:
   Solve upper triangular system R_m * y = g[1:m] via back-substitution
   x_m = x_0 + V_m * y
```

### 2.5 Restarted GMRES: GMRES($m$)
- **The Memory / Work Wall**:
  At iteration step $m$, GMRES must store $m$ vectors of length $n$. The memory cost is $O(m n)$ words, and the orthogonalization work scales as $O(m^2 n)$ FLOPs.
  For $m = 1000, n = 10^7$, storing $V_m$ requires $80\text{ Gigabytes}$ of RAM!
- **Restart Strategy**:
  Limit the subspace dimension to a fixed size $m$ (e.g., $m = 30$). Run $m$ steps of GMRES to compute $x_m$. If not converged, set $x_0 \leftarrow x_m$, wipe the basis $V_m$, and restart.
- **Drawback (Stagnation)**:
  Restarting discards the Krylov history. If the residual is orthogonal to the newly formed subspace, GMRES($m$) can **stagnate** (fail to make any further progress).

---

## 3. Lanczos Biorthogonalization for Non-Symmetric Systems

To restore short 3-term recurrences without storing all basis vectors, Lanczos biorthogonalization constructs two mutually orthogonal Krylov subspaces.

### 3.1 Dual Krylov Subspaces
Starting from vectors $v_1$ and $w_1$ with $v_1^T w_1 = 1$:
$$\mathcal{K}_m(A, v_1) = \text{span}\{v_1, A v_1, \dots, A^{m-1} v_1\}$$
$$\mathcal{K}_m(A^T, w_1) = \text{span}\{w_1, A^T w_1, \dots, (A^T)^{m-1} w_1\}$$

### 3.2 The Biorthogonality Conditions
The algorithm constructs bases $V_m = [v_1, \dots, v_m]$ and $W_m = [w_1, \dots, w_m]$ such that:
$$\mathbf{W_m^T V_m = I_m \iff w_i^T v_j = \delta_{ij}}$$
$$\mathbf{W_m^T A V_m = T_m = \begin{bmatrix} \alpha_1 & \gamma_1 & & \\ \beta_1 & \alpha_2 & \gamma_2 & \\ & \ddots & \ddots & \ddots \end{bmatrix} \quad (\text{Tridiagonal Matrix!})}$$

The coupled 3-term recurrences are:
$$\mathbf{\beta_j v_{j+1} = A v_j - \alpha_j v_j - \gamma_{j-1} v_{j-1}}$$
$$\mathbf{\gamma_j w_{j+1} = A^T w_j - \alpha_j w_j - \beta_{j-1} w_{j-1}}$$
where $\alpha_j = w_j^T A v_j$, and $\beta_j, \gamma_j$ are chosen such that $w_{j+1}^T v_{j+1} = 1$.

---

## 4. The Bi-Conjugate Gradient (BiCG) Algorithm

BiCG applies an oblique Petrov-Galerkin projection:
- Search space: $\mathcal{K}_m(A, r_0)$
- Constraint space: $\mathcal{L}_m = \mathcal{K}_m(A^T, r_0^*)$ (where $r_0^*$ is an arbitrary shadow residual, typically $r_0^* = r_0$).

### 4.1 Severe Practical Limitations of BiCG (Slide 91)
1. **Requires $A^T$ Operations**: Every step requires multiplying by the transpose matrix $A^T p^*$. In modern "matrix-free" simulation codes where $A x$ is evaluated via subroutine calls, $A^T$ is often impossible to compute.
2. **Double Risk of Breakdown**:
   - **Division by Zero in Step Size**: $p_j^T A p_j = 0$ while $p_j \ne \mathbf{0}$.
   - **Biorthogonality Breakdown**: $w_j^T v_j = 0$ while neither vector is zero.
3. **Wildly Erratic Residual History**: Unlike GMRES, the residual norm $\|r_k\|_2$ oscillates wildly, causing severe numerical round-off corruption.

---

## 5. Transpose-Free Variants: CGS and BiCGSTAB

To eliminate the problematic matrix transpose $A^T$, researchers reformulated the Krylov residuals in polynomial form:
$$r_j = \phi_j(A) r_0, \quad r_j^* = \phi_j(A^T) r_0^*$$
where $\phi_j(t)$ is a polynomial with $\phi_j(0) = 1$.

### 5.1 Conjugate Gradient Squared (CGS)
Peter Sonneveld (1989) noted that the inner products in BiCG satisfy:
$$\langle r_j, r_j^* \rangle = \langle \phi_j(A) r_0, \phi_j(A^T) r_0^* \rangle = \langle \phi_j^2(A) r_0, r_0^* \rangle$$
By defining the new residual as the **squared polynomial**:
$$\mathbf{r_j^{CGS} = \phi_j^2(A) r_0}$$
- **Advantages**: Completely eliminates $A^T$! If BiCG converges, CGS contracts error with polynomial $\phi_j^2(t)$, converging **twice as fast**!
- **Fatal Pathology (Slide 101)**:
  Squaring the polynomial also **squares the round-off errors and residual oscillations**! Near the solution, CGS often suffers catastrophic numerical overflow or severe floating-point degradation.

---

### 5.2 Bi-Conjugate Gradient Stabilized (BiCGSTAB)
H. A. van der Vorst (1992) resolved the instability of CGS by replacing the second oscillating polynomial $\phi_j(A)$ with a **smooth, monotonic 1D Steepest Descent operator**:
$$\mathbf{r_j = \psi_j(A) \phi_j(A) r_0}$$
where $\psi_j(t)$ is generated recursively by linear contraction factors:
$$\psi_j(t) = (1 - \omega_j t) \psi_{j-1}(t)$$
At each step, parameter $\omega_j$ is chosen via **exact local line search (Steepest Descent)** to minimize the 2-norm of the residual!

#### Complete Derivation of $\omega_j$ (Slide 103)
The intermediate residual vector after the BiCG step is:
$$s_j = r_j - \alpha_j A p_j$$
The updated residual after the stabilizing factor is:
$$r_{j+1} = (I - \omega_j A) s_j = s_j - \omega_j A s_j$$
We choose $\omega_j$ to minimize $\|r_{j+1}\|_2^2 = \|s_j - \omega_j A s_j\|_2^2$:
$$f(\omega) = s_j^T s_j - 2\omega (A s_j)^T s_j + \omega^2 (A s_j)^T (A s_j)$$
Setting derivative $\frac{df}{d\omega} = 0$:
$$-2 (A s_j)^T s_j + 2\omega (A s_j)^T (A s_j) = 0$$
$$\mathbf{\omega_j = \frac{(A s_j)^T s_j}{(A s_j)^T (A s_j)} = \frac{\langle A s_j, s_j \rangle}{\|A s_j\|_2^2}}$$

### 5.3 Complete BiCGSTAB Algorithm
```
Algorithm: BiCGSTAB (Bi-Conjugate Gradient Stabilized)
Input: Matrix A (non-symmetric, non-singular), vector b, initial guess x_0, tolerance ε
Output: Solution vector x

1. r_0 = b - A * x_0
   Choose r_0^* such that (r_0, r_0^*) != 0 (typically r_0^* = r_0)
   p_0 = r_0
   ρ_0 = α = ω_0 = 1, v_0 = 0

2. While (norm(r_k) > ε) do:
       ρ_k = (r_k, r_0^*)
       β = (ρ_k / ρ_{k-1}) * (α / ω_{k-1})
       p_k = r_k + β * (p_{k-1} - ω_{k-1} * v_{k-1})
       
       v_k = A * p_k                        // SpMV 1
       α = ρ_k / (v_k, r_0^*)
       
       s_k = r_k - α * v_k                  // Intermediate residual
       if norm(s_k) < ε then
           x_k = x_k + α * p_k
           break
       end if
       
       t_k = A * s_k                        // SpMV 2
       ω_k = (t_k, s_k) / (t_k, t_k)        // Steepest descent parameter
       
       x_{k+1} = x_k + α * p_k + ω_k * s_k  // Update solution
       r_{k+1} = s_k - ω_k * t_k            // Update residual
       k = k + 1
   end while
### 5.4 Numerical Benchmark: Convergence History of SOR vs. BCG vs. BiCGSTAB (Slide 106)

To demonstrate the real-world convergence properties on an unsymmetric PDE discretization, the lecture presents an exact iteration-by-iteration residual norm history for target tolerance $\epsilon \approx 10^{-9}$:

| Iteration ($k$) | SOR Residual Norm | BCG Residual Norm | BiCGSTAB Residual Norm |
| :---: | :---: | :---: | :---: |
| **1** | $2.184 \times 10^{-1}$ | $9.790$ | $3.726$ |
| **2** | $1.576 \times 10^{-1}$ | $3.372$ | $1.852$ |
| **3** | $1.062 \times 10^{-1}$ | $2.869$ | $2.239$ |
| **5** | $4.976 \times 10^{-3}$ | $1.972$ | $2.294 \times 10^{-1}$ |
| **7** | $1.946 \times 10^{-3}$ | $6.202 \times 10^{-1}$ | $2.010 \times 10^{-2}$ |
| **10** | $4.890 \times 10^{-4}$ | $1.139 \times 10^{-1}$ | $3.177 \times 10^{-3}$ |
| **12** | $1.305 \times 10^{-4}$ | $7.859 \times 10^{-2}$ | $1.148 \times 10^{-3}$ |
| **15** | $3.129 \times 10^{-5}$ | $2.756 \times 10^{-4}$ | $3.629 \times 10^{-4}$ |
| **17** | $7.726 \times 10^{-6}$ | $1.305 \times 10^{-2}$ *(oscillation!)* | $4.467 \times 10^{-6}$ |
| **18** | $4.868 \times 10^{-6}$ | $2.220 \times 10^{-4}$ | $2.046 \times 10^{-8}$ |
| **20** | $1.923 \times 10^{-6}$ | $2.230 \times 10^{-5}$ | $\mathbf{2.489 \times 10^{-9}}$ **(CONVERGED!)** |
| **22** | $7.603 \times 10^{-7}$ | $1.716 \times 10^{-7}$ | — |
| **25** | $1.889 \times 10^{-7}$ | $\mathbf{8.023 \times 10^{-9}}$ **(CONVERGED!)** | — |
| **30** | $4.695 \times 10^{-8}$ | — | — |
| **35** | $4.611 \times 10^{-9}$ | — | — |
| **44** | $\mathbf{7.076 \times 10^{-9}}$ **(CONVERGED!)** | — | — |

#### Critical Exam Observations:
1. **Convergence Speed**: BiCGSTAB reaches tolerance in **20 iterations**, outperforming BCG (25 iterations) and requiring **less than half the iterations of SOR (44 iterations)**.
2. **Wild BCG Oscillations**: Notice that at iteration 17, the BCG residual surges from $2.75 \times 10^{-4}$ back up to $1.305 \times 10^{-2}$ (an increase of two orders of magnitude!) due to near-breakdown in the biorthogonality condition.
3. **Smooth BiCGSTAB Stabilization**: In contrast, BiCGSTAB's intermediate steepest descent projection successfully damps these oscillations, maintaining steady quadratic-like residual reduction.

---

## 6. Master Comparative Synthesis of Iterative Solvers

The table below provides a comprehensive comparison across all iterative solvers studied in Chapters 05–11:


| Solver | Matrix Restrictions | Minimization Property | Recurrence Length | SpMVs / Iteration | Transpose $A^T$ Req.? | Stability & Convergence |
| :--- | :--- | :--- | :---: | :---: | :---: | :--- |
| **Jacobi** | Strictly Diagonally Dominant | None (fixed point) | 1 | 1 | No | Slow ($\rho$); parallel-friendly |
| **Gauss-Seidel** | Diagonally Dominant / SPD | None (fixed point) | 1 | 1 | No | $2\times$ faster than Jacobi ($\rho^2$) |
| **SOR** | SPD ($0 < \omega < 2$) | None (accelerated) | 1 | 1 | No | Fast at $\omega_{opt}$; serial chain |
| **Steepest Descent**| Symmetric Positive Definite | $J(x)$ along line $r_k$ | 1 | 1 | No | Slow (zigzags when $\kappa \gg 1$) |
| **Conjugate Gradient**| Symmetric Positive Definite | $\|e_k\|_A$ over $\mathcal{K}_k$ | 2 (short) | 1 | No | **Optimal for SPD**; Chebyshev bound |
| **FOM** | Any Non-Singular Matrix | Galerkin ($r_m \perp \mathcal{K}_m$) | $m$ (full) | 1 | No | Non-monotonic; can divide by zero |
| **GMRES** | Any Non-Singular Matrix | $\|r_m\|_2$ over $\mathcal{K}_m$ | $m$ (full) | 1 | No | **Monotonically decreasing**; $O(m n)$ memory |
| **BiCG** | Any Non-Singular Matrix | Oblique Petrov-Galerkin | 2 (short) | 2 ($A$ and $A^T$) | **YES** | Erratic oscillations; breakdowns |
| **CGS** | Any Non-Singular Matrix | Squaring of BiCG | 2 (short) | 2 | No | Fast when converging; highly unstable |
| **BiCGSTAB** | Any Non-Singular Matrix | Stabilized by SD step | 2 (short) | 2 | No | **Fast, smooth, robust standard for non-symmetric** |
