# Chapter 12: Preconditioning Techniques

---

## 1. Fundamentals and Motivation

The convergence rate of any Krylov subspace method (such as Conjugate Gradient or GMRES) is fundamentally dictated by two spectral properties of matrix $A$:
1. The **condition number** $\kappa_2(A) = \frac{\lambda_{\max}}{\lambda_{\min}}$ (or $\frac{\sigma_{\max}}{\sigma_{\min}}$).
2. The **clustering of eigenvalues** in the complex plane.

### 1.1 The Stiff PDE Dilemma
When solving elliptic or parabolic boundary value problems (e.g., Laplace, Poisson, diffusion equations) on a grid with spacing $h$:
$$\kappa_2(A) \sim O\left(\frac{1}{h^2}\right) = O(N)$$
As the grid is refined to improve scientific accuracy ($h \to 0$):
- For a $100 \times 100$ grid ($h = 0.01$): $\kappa_2(A) \sim 10^4$.
- For a $1000 \times 1000$ grid ($h = 0.001$): $\kappa_2(A) \sim 10^6$.
The theoretical iteration count of Conjugate Gradient scales as $O(\sqrt{\kappa}) \sim O(1/h) = O(\sqrt{N})$. The solver slows down drastically, and round-off error destroys orthogonality.

### 1.2 The Concept of Preconditioning
**Preconditioning** is the process of transforming the ill-conditioned system $Ax = b$ into an equivalent linear system that possesses:
- A condition number close to 1: $\kappa_2(M^{-1} A) \approx 1$.
- Highly clustered eigenvalues around 1.

### 1.3 The Three Golden Rules of Preconditioner Design
A practical preconditioner $M$ must satisfy three competing criteria:
1. **Spectral Closeness**: $M$ must be a close mathematical approximation to $A$ so that $M^{-1} A \approx I$.
2. **Computational Inversion Cost**: Solving the auxiliary system:
   $$\mathbf{M z = r}$$
   must be computationally very inexpensive ($O(N)$ operations).
3. **Storage Footprint**: The memory required to represent $M$ must be minimal, maintaining the sparse footprint of the original problem.

---

## 2. Three Classical Preconditioned Forms

Depending on how the preconditioning operator $M$ is applied to the linear system, there are three distinct formulations:

### 2.1 Left Preconditioning
Apply the inverse operator $M^{-1}$ from the left:
$$\mathbf{M^{-1} A x = M^{-1} b}$$
- **Operator**: $\tilde{A} = M^{-1} A$, $\tilde{b} = M^{-1} b$.
- **Residual Evaluated**:
  The solver operates on the **preconditioned residual**:
  $$\tilde{r} = M^{-1} b - M^{-1} A x = M^{-1} (b - Ax) = \mathbf{M^{-1} r}$$
- **Exam Consideration**: Stopping criteria testing $\|\tilde{r}\|_2 < \epsilon$ test the norm of $M^{-1} r$, which can differ substantially from the true physical residual $\|r\|_2 = \|b - Ax\|_2$ if $\|M^{-1}\|$ is very large or small.

### 2.2 Right Preconditioning
Transform the solution variable by setting $x = M^{-1} y$, leading to:
$$\mathbf{A M^{-1} y = b, \quad \text{with } x = M^{-1} y}$$
- **Operator**: $\tilde{A} = A M^{-1}$, unknown $y$.
- **Residual Evaluated**:
  $$\tilde{r} = b - (A M^{-1}) y = b - A (M^{-1} y) = \mathbf{b - Ax = r}$$
- **Decisive Advantage**: The solver tests the **TRUE physical residual** $r = b - Ax$ at every step!
- **Post-Processing**: Requires one final back-solve $x = M^{-1} y$ after convergence.

### 2.3 Split (Symmetric) Preconditioning
When $A$ and $M$ are Symmetric Positive Definite (SPD), factoring $M = L L^T$ (e.g., via Cholesky) enables symmetric split preconditioning:
$$L^{-1} A L^{-T} (L^T x) = L^{-1} b$$
$$\mathbf{\tilde{A} \tilde{x} = \tilde{b}}$$
where:
$$\mathbf{\tilde{A} = L^{-1} A L^{-T}, \quad \tilde{x} = L^T x, \quad \tilde{b} = L^{-1} b}$$
- **Preservation of Symmetry**:
  $$\tilde{A}^T = (L^{-1} A L^{-T})^T = (L^{-T})^T A^T (L^{-1})^T = L^{-1} A L^{-T} = \mathbf{\tilde{A}}$$
  Because $\tilde{A}$ is strictly symmetric positive definite, the **Conjugate Gradient method can be applied directly**!

---

## 3. Preconditioned Conjugate Gradient (PCG)

Direct application of split preconditioning would seem to require computing the Cholesky factor $L$ and solving systems with $L$ and $L^T$.
**The Master Stroke of PCG**:
Through algebraic transformation, all occurrences of $L$ and $L^T$ can be eliminated, expressing the entire algorithm in terms of the original matrix $A$ and a **single preconditioning solve $M z = r$** per iteration!

### 3.1 Step-by-Step Derivation of PCG
Apply standard CG to the symmetric split system $\tilde{A} \tilde{x} = \tilde{b}$:
1. **Define Transformed Variables**:
   $$\tilde{x}_k = L^T x_k, \quad \tilde{r}_k = L^{-1} r_k, \quad \tilde{p}_k = L^T p_k$$
2. **Transform Search Direction and Residual**:
   $$\tilde{x}_{k+1} = \tilde{x}_k + \alpha_k \tilde{p}_k \implies L^T x_{k+1} = L^T x_k + \alpha_k L^T p_k \implies \mathbf{x_{k+1} = x_k + \alpha_k p_k}$$
   $$\tilde{r}_{k+1} = \tilde{r}_k - \alpha_k \tilde{A} \tilde{p}_k \implies L^{-1} r_{k+1} = L^{-1} r_k - \alpha_k (L^{-1} A L^{-T}) (L^T p_k)$$
   Multiplying by $L$:
   $$\mathbf{r_{k+1} = r_k - \alpha_k A p_k}$$

3. **Introduce Auxiliary Vector $z_k$**:
   Define:
   $$\mathbf{z_k = M^{-1} r_k \iff M z_k = r_k}$$
   Notice the scalar inner product in the numerator of $\alpha_k$:
   $$\tilde{r}_k^T \tilde{r}_k = (L^{-1} r_k)^T (L^{-1} r_k) = r_k^T (L^{-T} L^{-1}) r_k = r_k^T (L L^T)^{-1} r_k = r_k^T M^{-1} r_k = \mathbf{r_k^T z_k}$$
4. **Denominator of $\alpha_k$**:
   $$\tilde{p}_k^T \tilde{A} \tilde{p}_k = (L^T p_k)^T (L^{-1} A L^{-T}) (L^T p_k) = p_k^T (L L^{-1}) A (L^{-T} L^T) p_k = \mathbf{p_k^T A p_k}$$
   Thus:
   $$\mathbf{\alpha_k = \frac{r_k^T z_k}{p_k^T A p_k}}$$
5. **Fletcher-Reeves Coefficient $\beta_k$**:
   $$\beta_k = \frac{\tilde{r}_{k+1}^T \tilde{r}_{k+1}}{\tilde{r}_k^T \tilde{r}_k} = \mathbf{\frac{r_{k+1}^T z_{k+1}}{r_k^T z_k}}$$
6. **Search Direction Update**:
   $$\tilde{p}_{k+1} = \tilde{r}_{k+1} + \beta_k \tilde{p}_k \implies L^T p_{k+1} = L^{-1} r_{k+1} + \beta_k L^T p_k$$
   Multiplying by $L^{-T}$:
   $$p_{k+1} = (L L^T)^{-1} r_{k+1} + \beta_k p_k = M^{-1} r_{k+1} + \beta_k p_k \implies \mathbf{p_{k+1} = z_{k+1} + \beta_k p_k}$$

### 3.2 The Complete PCG Algorithm
```
Algorithm: Preconditioned Conjugate Gradient (PCG)
Input: Matrix A (SPD), Preconditioner M (SPD), vector b, guess x_0, tolerance ε
Output: Solution vector x

1. r_0 = b - A * x_0
   Solve M * z_0 = r_0              // Preconditioner solve 0
   p_0 = z_0
   γ_0 = r_0^T * z_0
   k = 0

2. While (norm(r_k) > ε) do:
       w_k = A * p_k                // SpMV
       α_k = γ_k / (p_k^T * w_k)    // Step size
       
       x_{k+1} = x_k + α_k * p_k    // Update solution
       r_{k+1} = r_k - α_k * w_k    // Update residual
       
       Solve M * z_{k+1} = r_{k+1}  // Preconditioner solve (ONLY 1 per step!)
       
       γ_{k+1} = r_{k+1}^T * z_{k+1}
       β_k = γ_{k+1} / γ_k
       p_{k+1} = z_{k+1} + β_k * p_k
       
       k = k + 1
   end while
```

- **Per-Iteration Work**:
  - $1$ SpMV ($A p_k$).
  - **$1$ Preconditioner Solve** ($M z_{k+1} = r_{k+1}$).
  - $2$ Dot Products ($p_k^T A p_k$ and $r_{k+1}^T z_{k+1}$).
  - $3$ AXPY vector updates.

### 3.3 Alternative Derivation: Left Preconditioned CG via the $M$-Inner Product (Slides 406–410, 415)

In lecture (Slides 406–410), PCG is also derived directly from the left-preconditioned system:
$$M^{-1} A x = M^{-1} b$$
where both $A$ and $M$ are SPD ($M = L L^T$).
Although the product matrix $M^{-1} A$ is not symmetric under the standard Euclidean dot product ($u^T v$), it **is strictly self-adjoint under the $M$-inner product**:
$$\langle u, v \rangle_M \equiv u^T M v$$

#### Proof of Symmetry under $\langle \cdot, \cdot \rangle_M$
$$\langle M^{-1} A u, v \rangle_M = (M^{-1} A u)^T M v = u^T A^T M^{-T} M v = u^T A v = u^T M (M^{-1} A v) = \langle u, M^{-1} A v \rangle_M \quad \checkmark$$

#### Left PCG Derivation Steps (Slides 407–410)
1. **Residual of Preconditioned System**:
   $$z_0 = M^{-1} b - M^{-1} A x_0 = M^{-1} r_0$$
   Search direction initialized to $p_0 = z_0$.
2. **Residual Update**:
   $$z_{j+1} = z_j - \alpha_j M^{-1} A p_j$$
3. **Imposing $M$-Orthogonality**:
   $$\langle z_{j+1}, z_j \rangle_M = 0 \implies \langle z_j - \alpha_j M^{-1} A p_j, z_j \rangle_M = 0$$
   $$\alpha_j = \frac{\langle z_j, z_j \rangle_M}{\langle M^{-1} A p_j, z_j \rangle_M} = \frac{z_j^T M z_j}{(M^{-1} A p_j)^T M z_j}$$
   - Numerator: $z_j^T M z_j = z_j^T M (M^{-1} r_j) = z_j^T r_j = (r_j, z_j)$.
   - Denominator: Since $z_j = p_j - \beta_j p_{j-1}$ and $p_j$ is $A$-conjugate to previous directions:
     $$(M^{-1} A p_j)^T M z_j = (A p_j)^T M^{-1} M (p_j - \beta_j p_{j-1}) = p_j^T A p_j = (A p_j, p_j)$$
   Therefore:
   $$\mathbf{\alpha_j = \frac{r_j^T z_j}{p_j^T A p_j}}$$
4. **Direction Parameter $\beta_j$**:
   $$\beta_j = \frac{\langle z_{j+1}, z_{j+1} \rangle_M}{\langle z_j, z_j \rangle_M} = \mathbf{\frac{r_{j+1}^T z_{j+1}}{r_j^T z_j}}$$
5. **Search Direction Update**:
   $$p_{j+1} = z_{j+1} + \beta_j p_j$$

#### Master Equivalence Theorem (Slide 415)
> **Theorem**:
> 1. All variants of Preconditioned CG (Left PCG, Right PCG, and Split PCG) strictly preserve the symmetry of the underlying linear system (via the $M$-inner product in Left PCG or similarity/congruence in Split PCG).
> 2. **All variants generate the exact same mathematical iterates**: Starting with the same initial guess $x_0$, the computed solution $x_k$ at every single iteration step $k$ is identical across all three formulations!
> 3. Although preconditioning fundamentally accelerates convergence by compressing the spectrum, each iteration introduces solving an auxiliary linear system ($M z_j = r_j$). Preconditioning is only computationally profitable if $M z = r$ can be solved rapidly (e.g., triangular solves via Incomplete Cholesky or diagonal scaling).

---

## 4. Preconditioning for GMRES

### 4.1 Left-Preconditioned GMRES
Applies Arnoldi's process to the preconditioned operator $\tilde{A} = M^{-1} A$ starting from initial vector $\tilde{r}_0 = M^{-1} (b - A x_0)$:
- Solves:
  $$\min_{y \in \mathbb{R}^m} \|M^{-1} r_0 - \bar{H}_m y\|_2$$
- At each Arnoldi step, computing $w = \tilde{A} v_j$ requires:
  1. SpMV: $u = A v_j$.
  2. Preconditioner solve: $w = M^{-1} u \iff M w = u$.

### 4.2 Right-Preconditioned GMRES
Solves $A M^{-1} y = b$, with $x = M^{-1} y$:
- Generates orthonormal basis for $\mathcal{K}_m(A M^{-1}, r_0)$.
- Minimizes the **true Euclidean residual norm** $\|b - Ax_m\|_2$ at every step.
- At each Arnoldi step:
  1. Preconditioner solve: $u = M^{-1} v_j \iff M u = v_j$.
  2. SpMV: $w = A u$.
- Solution recovered at termination:
  $$x_m = x_0 + M^{-1} (V_m y_m)$$

### 4.3 Flexible GMRES (FGMRES) (Slide 422)
- **Problem**: In standard GMRES, $M$ must be a fixed, constant linear operator throughout all $m$ iterations. If $M$ is non-linear, or is itself an inner iterative solver with varying tolerances, standard GMRES fails.
- **Solution (FGMRES)**:
  Stores the preconditioned basis vectors $z_j = M_j^{-1} v_j$ in a separate matrix $Z_m = [z_1, \dots, z_m]$.
  The Arnoldi relation becomes:
  $$A Z_m = V_{m+1} \bar{H}_m$$
  The solution is updated as:
  $$x_m = x_0 + Z_m y_m$$
  Allows the preconditioner to change arbitrarily at each step!

---

## 5. Standard Classes of Preconditioners

### 5.1 Jacobi (Diagonal) Preconditioner
- **Definition**: $M = D = \text{diag}(a_{11}, a_{22}, \dots, a_{nn})$.
- **Solving $M z = r$**:
  $$z_i = \frac{r_i}{a_{ii}}$$
- **Pros**: Zero memory overhead; embarrassingly parallel.
- **Cons**: Only accounts for row scales; ineffective for highly coupled grid stencils.

### 5.2 Symmetric Gauss-Seidel (SGS) / SSOR Preconditioner (Slides 426–428)
Combining forward and backward Gauss-Seidel sweeps yields a symmetric preconditioner:
$$\mathbf{M_{SGS} = (D + L) D^{-1} (D + U)}$$
where $A = D + L + U$.
- **Solving $M_{SGS} z = r$** in three stages:
  1. Forward triangular solve: $(D + L) y = r$.
  2. Diagonal scaling: $w = D y$.
  3. Backward triangular solve: $(D + U) z = w$.
- **Cost**: Exactly $O(N_{nz})$ FLOPs.
- **Advantage**: Guaranteed to be symmetric positive definite whenever $A$ is SPD!

### 5.3 Incomplete LU Factorization (ILU) (Slide 429)
The ideal direct preconditioner is $M = A = LU$, which would solve the system in 1 step. But full LU suffers from memory-destroying fill-in.
**Incomplete LU (ILU)** computes approximate triangular factors $\tilde{L}$ and $\tilde{U}$ while dropping fill-in elements according to a sparsity rule:
$$A = \tilde{L} \tilde{U} - R = M - R$$
where $R$ is the residual error matrix.

#### 1. ILU with Zero Fill-In: ILU(0)
- **Zero Fill-In Rule**: During Gaussian elimination, if an entry in position $(i, j)$ was originally zero in matrix $A$ ($a_{ij} = 0$), any newly generated fill-in is **immediately discarded (set to zero)**!
- **Sparsity Pattern**:
  $$\text{NonzeroPattern}(\tilde{L} + \tilde{U}) = \text{NonzeroPattern}(A)$$
- **Memory Footprint**: Exactly equal to the memory of the original sparse matrix $A$!
- **Solving $M z = r$**: Two triangular sparse solves ($\tilde{L} y = r$, then $\tilde{U} z = y$).

#### 2. Threshold-Based Incomplete LU: ILUT($p, \tau$)
- Controls fill-in using two parameters:
  - **Threshold $\tau$**: Any element $l_{ij}$ or $u_{ij}$ with absolute value smaller than $\tau \cdot \|a_i\|$ is dropped.
  - **Fill Limit $p$**: Only the $p$ largest entries in each row are retained.
- Bridges the gap between cheap ILU(0) and exact direct factorization.

#### 3. Incomplete Cholesky: IC(0)
- For SPD matrices, $A \approx \tilde{L} \tilde{L}^T = M$. Preserves symmetry and positive definiteness.

---

## 6. Preconditioner Performance Trade-Off

The total wall-clock time of a preconditioned iterative solve is:
$$\mathbf{T_{total} = T_{setup}(M) + N_{iter}(M) \times \left[ T_{SpMV} + T_{solve}(M z = r) \right]}$$

```
Compute Time
  ▲
  │ \                                           /  Expensive Preconditioner
  │  \                                         /   (High setup & solve time,
  │   \                                       /     e.g., ILUT with high fill)
  │    \             Total Wall-Clock         /
  │     \                 Time               /
  │      \                 │                /
  │       \                ▼               /
  │        \_        Optimal Point        /
  │          \______       •       ______/
  │                 \____     ____/
  │                      \___/
  │ Naive / No Precond.   │
  │ (Excessive iterations)│
  └───────────────────────┼────────────────────────► Preconditioner Complexity
                        ILU(0) / SGS
```

- **Poor Preconditioner (e.g., None or Naive Jacobi)**: $T_{solve}$ is zero, but $N_{iter}$ is massive $\implies$ high total time.
- **Overly Dense Preconditioner (e.g., Near-Exact ILU)**: $N_{iter}$ is tiny (1 or 2), but $T_{setup}$ and $T_{solve}$ are huge ($O(N^2)$) $\implies$ defeats the purpose of iterative methods.
- **Optimal Sweet Spot**: **ILU(0)** or **SGS**, which typically reduce $N_{iter}$ by $5\times$ to $50\times$ while keeping $T_{solve}(M z = r)$ strictly $O(N)$!
