# Chapter 08: Projection Methods & the Method of Steepest Descent

---

## 1. General Petrov-Galerkin Projection Framework

Most modern iterative solvers for large sparse linear systems $Ax = b$ can be understood within the unified mathematical framework of **projection methods**.

```
                           R^n Vector Space
                ┌───────────────────────────────────────┐
                │                                       │
                │        Initial Guess x₀               │
                │              •                        │
                │               \                       │
                │                \   Search Subspace K  │
                │                 ▼  (dimension m)      │
                │      Approximate Solution x_m         │
                │              •                        │
                │              │                        │
                │              │ Residual r_m = b - Ax_m │
                │              ▼                        │
                │   ────────────────────────            │
                │     Constraint Subspace L             │
                │       (r_m is orthogonal to L)        │
                └───────────────────────────────────────┘
```

### 1.1 Formal Definition of a Projection Process
Let $A \in \mathbb{R}^{n \times n}$ be a real square matrix. Let $\mathcal{K}$ and $\mathcal{L}$ be two $m$-dimensional subspaces of $\mathbb{R}^n$:
- $\mathcal{K}$ is the **Search Subspace** (from which candidate corrections are drawn).
- $\mathcal{L}$ is the **Constraint Subspace** (or Test Subspace, against which orthogonality is enforced).

Given an initial guess $x_0 \in \mathbb{R}^n$ with initial residual $r_0 = b - A x_0$:
A **projection technique onto $\mathcal{K}$ orthogonal to $\mathcal{L}$** finds an approximate solution $x_m \in x_0 + \mathcal{K}$ such that the new residual $r_m = b - A x_m$ satisfies the **Petrov-Galerkin Condition**:
$$\mathbf{r_m \perp \mathcal{L} \iff \langle w, b - A x_m \rangle = 0 \quad \text{for all } w \in \mathcal{L}}$$

- **Orthogonal Projection (Galerkin Method)**: When $\mathcal{L} = \mathcal{K}$ (search and constraint spaces are identical).
- **Oblique Projection (Petrov-Galerkin Method)**: When $\mathcal{L} \ne \mathcal{K}$.

### 1.2 Matrix Formulation of the General Projection Method
Let the columns of $V_m \in \mathbb{R}^{n \times m}$ form a basis for search space $\mathcal{K}$:
$$V_m = [v_1, v_2, \dots, v_m]$$
Let the columns of $W_m \in \mathbb{R}^{n \times m}$ form a basis for constraint space $\mathcal{L}$:
$$W_m = [w_1, w_2, \dots, w_m]$$

Any candidate vector $x_m \in x_0 + \mathcal{K}$ can be uniquely written as:
$$\mathbf{x_m = x_0 + V_m y_m}$$
where $y_m \in \mathbb{R}^m$ is a coordinate vector to be determined.

The new residual vector is:
$$r_m = b - A x_m = b - A (x_0 + V_m y_m) = (b - A x_0) - A V_m y_m = r_0 - A V_m y_m$$

Enforcing the Petrov-Galerkin condition $r_m \perp \mathcal{L}$:
Every basis vector $w_i$ of $\mathcal{L}$ must be orthogonal to $r_m$:
$$w_i^T r_m = 0 \quad \text{for all } i = 1, \dots, m \iff W_m^T r_m = \mathbf{0}$$
Substitute the expression for $r_m$:
$$W_m^T (r_0 - A V_m y_m) = \mathbf{0} \implies \mathbf{(W_m^T A V_m) y_m = W_m^T r_0}$$

Assuming the projected $m \times m$ matrix $W_m^T A V_m$ is non-singular:
$$\mathbf{y_m = (W_m^T A V_m)^{-1} W_m^T r_0}$$
The updated solution is:
$$\mathbf{x_m = x_0 + V_m (W_m^T A V_m)^{-1} W_m^T r_0}$$

---

## 2. One-Dimensional Projection Processes ($m = 1$)

When $m = 1$, the subspaces $\mathcal{K}$ and $\mathcal{L}$ are 1-dimensional lines spanned by vectors $v \in \mathbb{R}^n$ and $w \in \mathbb{R}^n$, respectively.
The update takes the form of a line search:
$$x^{(k+1)} = x^{(k)} + \alpha_k v_k$$
where $\alpha_k$ is a scalar step length.
The residual is:
$$r_{k+1} = r_k - \alpha_k A v_k$$
Enforcing $w_k^T r_{k+1} = 0$:
$$w_k^T (r_k - \alpha_k A v_k) = 0 \implies w_k^T r_k - \alpha_k w_k^T A v_k = 0$$
$$\mathbf{\alpha_k = \frac{w_k^T r_k}{w_k^T A v_k}}$$

Depending on the choices of $v_k$ and $w_k$, we obtain three fundamental 1D iterative algorithms:

| Method | Search Direction ($v_k$) | Constraint Direction ($w_k$) | Step Length ($\alpha_k$) | Matrix Applicability |
| :--- | :---: | :---: | :---: | :--- |
| **Steepest Descent (SD)** | $r_k$ | $r_k$ | $\frac{r_k^T r_k}{r_k^T A r_k}$ | Symmetric Positive Definite |
| **Minimum Residual (MR)**| $r_k$ | $A r_k$ | $\frac{r_k^T A r_k}{\|A r_k\|_2^2}$ | Positive Definite ($\frac{A+A^T}{2} > 0$) |
| **Residue Norm SD** | $A^T r_k$ | $A v_k = A A^T r_k$ | $\frac{\|v_k\|_2^2}{\|A v_k\|_2^2}$ | Any Non-Singular Square Matrix |

### 2.1 The Minimum Residual (MR) Method
- **Formulation**: $v_k = r_k$ and $w_k = A r_k$.
- **Minimization Property**: At each step, $\alpha_k$ minimizes the Euclidean norm of the new residual:
  $$\alpha_k = \arg\min_\alpha \|b - A(x_k + \alpha r_k)\|_2^2$$
- **Convergence Theorem**:
  If the symmetric part of $A$, $M_s = \frac{A + A^T}{2}$, is positive definite with minimum eigenvalue $\mu = \lambda_{\min}(M_s) > 0$, and $\sigma = \|A\|_2$:
  $$\mathbf{\|r_{k+1}\|_2 \le \left( 1 - \frac{\mu^2}{\sigma^2} \right)^{1/2} \|r_k\|_2}$$
  The method converges monotonically for any initial guess $x_0$!

### 2.2 The Residue Norm Steepest Descent Method
- **Formulation**: $v_k = A^T r_k$ and $w_k = A v_k$.
- **Equivalent System**: This is mathematically identical to applying the Steepest Descent algorithm to the **Normal Equations**:
  $$A^T A x = A^T b$$
  Since $A^T A$ is always symmetric positive definite for any non-singular matrix $A$, this method converges for **any square non-singular matrix**!

---

## 3. The Method of Steepest Descent in Depth

The Steepest Descent method is the foundational gradient-based iterative solver for Symmetric Positive Definite (SPD) systems.

### 3.1 The Quadratic Energy Functional
Let $A \in \mathbb{R}^{n \times n}$ be an SPD matrix ($A = A^T$ and $x^T A x > 0$ for $x \ne \mathbf{0}$), and $b \in \mathbb{R}^n$.
Define the **quadratic functional** $J: \mathbb{R}^n \to \mathbb{R}$:
$$\mathbf{J(x) = \frac{1}{2} x^T A x - x^T b = \frac{1}{2} \sum_{i=1}^n \sum_{j=1}^n a_{ij} x_i x_j - \sum_{i=1}^n b_i x_i}$$

### 3.2 Equivalence of Minimizing $J(x)$ and Solving $Ax = b$
> **Theorem**: If $A$ is Symmetric Positive Definite, then $x^*$ is the unique global minimizer of $J(x)$ if and only if $Ax^* = b$.

#### Complete Mathematical Proof (From Slides 476–480)
1. **Compute the Gradient $\nabla J(x)$**:
   Evaluate the partial derivative with respect to component $x_k$:
   $$\frac{\partial J}{\partial x_k} = \frac{\partial}{\partial x_k} \left[ \frac{1}{2} \sum_{i=1}^n \sum_{j=1}^n a_{ij} x_i x_j - \sum_{i=1}^n b_i x_i \right]$$
   The double sum contains $x_k$ when $i = k$, when $j = k$, and when $i = j = k$:
   $$\frac{\partial J}{\partial x_k} = \frac{1}{2} \left[ \sum_{j=1}^n a_{kj} x_j + \sum_{i=1}^n a_{ik} x_i \right] - b_k$$
   Because $A$ is symmetric ($a_{ik} = a_{ki}$), the two sums are identical:
   $$\frac{\partial J}{\partial x_k} = \frac{1}{2} \left[ 2 \sum_{j=1}^n a_{kj} x_j \right] - b_k = \sum_{j=1}^n a_{kj} x_j - b_k = (Ax)_k - b_k$$
   In vector form:
   $$\mathbf{\nabla J(x) = Ax - b = -r(x)}$$

2. **Compute the Hessian $\nabla^2 J(x)$**:
   $$\frac{\partial^2 J}{\partial x_k \partial x_l} = \frac{\partial}{\partial x_l} \left( \sum_{j=1}^n a_{kj} x_j - b_k \right) = a_{kl} \implies \mathbf{\nabla^2 J(x) = A}$$
   Since $A$ is positive definite ($A > 0$), the Hessian is strictly positive definite everywhere!
   Therefore, $J(x)$ is **strictly convex**, possessing no local maxima or saddle points, but exactly **one unique global minimum** where $\nabla J(x) = \mathbf{0}$:
   $$\nabla J(x^*) = \mathbf{0} \iff A x^* - b = \mathbf{0} \iff \mathbf{Ax^* = b} \quad \blacksquare$$

### 3.3 Exact Line Search along the Negative Gradient
At any point $x_k$, the direction of steepest descent (fastest decrease of $J$) is the negative gradient:
$$v_k = -\nabla J(x_k) = -(Ax_k - b) = b - Ax_k = \mathbf{r_k}$$
We search along the ray $x(\alpha) = x_k + \alpha r_k$ to find the optimal step length $\alpha_k$ that minimizes $J$:
$$\alpha_k = \arg\min_{\alpha \in \mathbb{R}} J(x_k + \alpha r_k)$$

#### Complete Algebraic Derivation of $\alpha_k$ (From Slides 486–487)
Define $f(\alpha) = J(x_k + \alpha r_k)$:
$$f(\alpha) = \frac{1}{2} (x_k + \alpha r_k)^T A (x_k + \alpha r_k) - (x_k + \alpha r_k)^T b$$
Expanding:
$$f(\alpha) = \frac{1}{2} \left[ x_k^T A x_k + 2\alpha r_k^T A x_k + \alpha^2 r_k^T A r_k \right] - x_k^T b - \alpha r_k^T b$$
$$f(\alpha) = \left[ \frac{1}{2} x_k^T A x_k - x_k^T b \right] + \alpha \left[ r_k^T A x_k - r_k^T b \right] + \frac{\alpha^2}{2} r_k^T A r_k$$
$$f(\alpha) = J(x_k) - \alpha r_k^T (b - Ax_k) + \frac{\alpha^2}{2} r_k^T A r_k = J(x_k) - \alpha r_k^T r_k + \frac{\alpha^2}{2} r_k^T A r_k$$

To find the minimum, set the derivative $\frac{df}{d\alpha} = 0$:
$$\frac{df}{d\alpha} = -r_k^T r_k + \alpha r_k^T A r_k = 0$$
$$\mathbf{\alpha_k = \frac{r_k^T r_k}{r_k^T A r_k}}$$

### 3.4 Crucial Geometric Property: Orthogonality of Successive Residuals
> **Theorem**: In the Steepest Descent method, consecutive residual vectors are **strictly orthogonal**:
> $$\mathbf{r_{k+1}^T r_k = 0 \quad (\text{i.e., } r_{k+1} \perp r_k)}$$

#### Mathematical Proof
The residual at step $k+1$ is:
$$r_{k+1} = b - A x_{k+1} = b - A (x_k + \alpha_k r_k) = (b - A x_k) - \alpha_k A r_k = \mathbf{r_k - \alpha_k A r_k}$$
Take the dot product of $r_{k+1}$ with $r_k$:
$$r_{k+1}^T r_k = (r_k - \alpha_k A r_k)^T r_k = r_k^T r_k - \alpha_k (r_k^T A r_k)$$
Substitute $\alpha_k = \frac{r_k^T r_k}{r_k^T A r_k}$:
$$r_{k+1}^T r_k = r_k^T r_k - \left( \frac{r_k^T r_k}{r_k^T A r_k} \right) (r_k^T A r_k) = r_k^T r_k - r_k^T r_k = 0 \quad \blacksquare$$

- **Geometric Interpretation**: At the minimum along a line search, the search trajectory is exactly **tangent to an iso-contour of $J(x)$**. The gradient at that point ($r_{k+1}$) is perpendicular to the tangent ($r_k$). Thus, the search must turn $90^\circ$ at every step!

---

## 4. Operational Algorithm and Work Reduction

A naive implementation computes $r_{k+1} = b - A x_{k+1}$ via matrix multiplication at each step, requiring two SpMVs per iteration.
By reusing $p_k = A r_k$, we reduce this to **only ONE matrix-vector product per iteration**:

```
Algorithm: Optimized Method of Steepest Descent
Input: Matrix A (SPD), vector b, initial guess x_0, tolerance ε
Output: Converged solution x

1. Compute initial state:
   r = b - A * x_0
   p = A * r

2. While (norm(r) > ε) do:
   alpha = (r^T * r) / (r^T * p)     // 2 dot products
   x = x + alpha * r                 // Vector AXPY
   r = r - alpha * p                 // Vector AXPY (cheap residual update!)
   p = A * r                         // ONLY 1 SpMV per iteration!
end while
```

- **Per-Iteration Computational Cost**:
  - $1$ Sparse Matrix-Vector product ($p = A r$).
  - $2$ Vector Dot Products ($r^T r$ and $r^T p$).
  - $2$ Vector Updates (AXPY).
  - **Total Work**: $\approx 2 N_{nz} + 6n$ FLOPs.

---

## 5. Convergence Rate and the "Zigzagging" Pathology

### 5.1 Convergence Bound in the $A$-Norm
Let $e_k = x_k - x^*$ be the error vector. The $A$-norm (energy norm) is $\|e_k\|_A = \sqrt{e_k^T A e_k}$.
> **Theorem (Kantorovich)**: For an SPD matrix $A$, the error in Steepest Descent satisfies:
> $$\mathbf{\|e_{k+1}\|_A \le \left( \frac{\lambda_{\max} - \lambda_{\min}}{\lambda_{\max} + \lambda_{\min}} \right) \|e_k\|_A = \left( \frac{\kappa_2(A) - 1}{\kappa_2(A) + 1} \right) \|e_k\|_A}$$
> where $\kappa_2(A) = \lambda_{\max}/\lambda_{\min}$ is the spectral condition number.

### 5.2 The "Zigzagging" Pathology
- **Case 1: Well-Conditioned Matrix ($\kappa_2(A) \approx 1 \implies \lambda_{\max} \approx \lambda_{\min}$)**:
  - The iso-contours of $J(x)$ are spherical (circular).
  - The negative gradient $-\nabla J(x_0)$ points directly toward the center (solution $x^*$).
  - The method reaches the exact solution in **a single iteration**!
- **Case 2: Ill-Conditioned Matrix ($\kappa_2(A) \gg 1$)**:
  - The iso-contours are extremely elongated, needle-thin ellipsoidal valleys.
  - Because each step is strictly orthogonal to the previous one ($r_{k+1} \perp r_k$), the trajectory bounces back and forth across the valley walls at $90^\circ$ angles (**zigzagging**).
  - The convergence factor $\frac{\kappa - 1}{\kappa + 1} \approx 1 - \frac{2}{\kappa} \to 1$.
  - Progress along the bottom of the valley toward $x^*$ becomes infinitesimally slow!

```
Iso-contours of J(x) for ill-conditioned A (κ >> 1):
        ┌────────────────────────────────────────────────────────┐
        │                 .-'""'-.                               │
        │             .-'          '-.                           │
        │         .-'        x₀        '-.                       │
        │      .-'            \           '-.                    │
        │    .'                \             '.                  │
        │   /     x*            \ 90°          \                 │
        │  │      • ◄──┐         \              │  Zigzagging    │
        │   \          │ 90°    .'             /   trajectory!   │
        │    '.        └───┘  .'             .'                  │
        │      '-.          .'            .-'                    │
        │         '-.    .-'           .-'                       │
        │            '--'          .-'                           │
        └────────────────────────────────────────────────────────┘
```

### 5.3 Algebraic Emergence of the Krylov Subspace (Slide 196)
Examining the algebraic progression of Steepest Descent iterates reveals a fundamental mathematical pattern:

1. **Step 0 to 1**:
   $$r_0 = b - A x_0$$
   $$x_1 = x_0 + \alpha_0 r_0$$
2. **Step 1 to 2**:
   $$r_1 = b - A x_1 = b - A(x_0 + \alpha_0 r_0) = r_0 - \alpha_0 A r_0$$
   $$x_2 = x_1 + \alpha_1 r_1 = x_0 + \alpha_0 r_0 + \alpha_1 (r_0 - \alpha_0 A r_0) = x_0 + (\alpha_0 + \alpha_1) r_0 - \alpha_0 \alpha_1 A r_0$$
3. **Step 2 to 3**:
   $$r_2 = r_1 - \alpha_1 A r_1 = r_0 - (\alpha_0 + \alpha_1) A r_0 + \alpha_0 \alpha_1 A^2 r_0$$
   $$x_3 = x_2 + \alpha_2 r_2 = x_0 + (\alpha_0 + \alpha_1 + \alpha_2) r_0 - (\alpha_1 + \alpha_2) A r_0 + \alpha_2 A^2 r_0$$
4. **General Iterate ($p$-th Iteration)**:
   By induction, expanding in matrix powers of $A$:
   $$\mathbf{x_p = x_0 + \beta_0 r_0 + \beta_1 A r_0 + \beta_2 A^2 r_0 + \dots + \beta_{p-1} A^{p-1} r_0}$$
   $$\mathbf{x_p \in x_0 + \text{span}\{r_0, A r_0, A^2 r_0, \dots, A^{p-1} r_0\} = x_0 + \mathcal{K}_p(A, r_0)}$$

> [!IMPORTANT]
> **Key Conceptual Revelation**: Every iterate of Steepest Descent naturally belongs to the Krylov subspace $\mathcal{K}_p(A, r_0)$. However, because Steepest Descent chooses the scalar coefficients $\alpha_k$ greedily one step at a time, it falls victim to severe zigzagging. If we instead search for the **globally optimal** linear combination of all $p$ vectors $\{r_0, A r_0, \dots, A^{p-1} r_0\}$ simultaneously, we arrive directly at **Krylov Subspace Projection Methods (Arnoldi, FOM, Lanczos, and Conjugate Gradients)**!

