# Chapter 06: Spectral Theory, SVD, Positive Definiteness & Condition Numbers

---

## 1. Eigenvalues and Eigenvectors

### 1.1 Definitions and Fundamental Concepts
For a square matrix $A \in \mathbb{R}^{n \times n}$, a non-zero vector $x \in \mathbb{C}^n$ is an **eigenvector** of $A$ corresponding to scalar **eigenvalue** $\lambda \in \mathbb{C}$ if:
$$A x = \lambda x \iff (A - \lambda I) x = \mathbf{0}$$
Since $x \ne \mathbf{0}$, the operator $(A - \lambda I)$ must be singular:
$$\det(A - \lambda I) = 0$$
This determinant expands into the **characteristic polynomial** of degree $n$:
$$P_n(\lambda) = (-1)^n \lambda^n + c_{n-1} \lambda^{n-1} + \dots + c_1 \lambda + \det(A) = 0$$

- **Algebraic Multiplicity ($m_a$)**: The number of times $\lambda$ appears as a root of $P_n(\lambda) = 0$.
- **Geometric Multiplicity ($m_g$)**: The dimension of the eigenspace $N(A - \lambda I)$ (number of linearly independent eigenvectors associated with $\lambda$).
- **Defective Matrix**: A matrix where $m_g < m_a$ for at least one eigenvalue (cannot be diagonalized).

---

## 2. Spectral Properties of Real Symmetric Matrices

Real symmetric matrices ($A = A^T$) govern almost all physical equilibrium and energy problems (diffusion, elasticity, structural dynamics). They possess extraordinary mathematical properties that underpin scientific computing.

### 2.1 Theorem 1: All Eigenvalues of a Real Symmetric Matrix are Real
> **Theorem**: If $A \in \mathbb{R}^{n \times n}$ and $A = A^T$, then all $n$ eigenvalues $\lambda_1, \dots, \lambda_n$ are **strictly real numbers** ($\lambda_i \in \mathbb{R}$).

#### Complete Mathematical Proof (From Slides 126–131)
Let $\lambda \in \mathbb{C}$ be an eigenvalue of $A$ with non-zero eigenvector $x \in \mathbb{C}^n$:
$$A x = \lambda x \quad \text{--- (1)}$$
Take the complex conjugate transpose (Hermitian conjugate, $H$) of equation (1):
$$(A x)^H = (\lambda x)^H \implies x^H A^H = \bar{\lambda} x^H \quad \text{--- (2)}$$
where $\bar{\lambda}$ denotes the complex conjugate of $\lambda$, and $x^H = (\bar{x})^T$.

Because $A$ has real entries ($A = \bar{A}$) and is symmetric ($A = A^T$):
$$A^H = (\bar{A})^T = A^T = A$$
Substitute $A^H = A$ into equation (2):
$$x^H A = \bar{\lambda} x^H \quad \text{--- (3)}$$

Now multiply equation (1) from the left by $x^H$:
$$x^H (A x) = x^H (\lambda x) = \lambda (x^H x) \quad \text{--- (4)}$$
Multiply equation (3) from the right by $x$:
$$(x^H A) x = (\bar{\lambda} x^H) x = \bar{\lambda} (x^H x) \quad \text{--- (5)}$$

Subtract equation (5) from equation (4):
$$x^H A x - x^H A x = \lambda (x^H x) - \bar{\lambda} (x^H x)$$
$$0 = (\lambda - \bar{\lambda}) (x^H x)$$
Since $x$ is a non-zero eigenvector:
$$x^H x = \sum_{i=1}^n \bar{x}_i x_i = \sum_{i=1}^n |x_i|^2 = \|x\|_2^2 > 0$$
Dividing by $\|x\|_2^2$:
$$\lambda - \bar{\lambda} = 0 \implies \mathbf{\lambda = \bar{\lambda}}$$
A complex number that equals its conjugate must have zero imaginary part ($\text{Im}(\lambda) = 0$).
Therefore, **$\lambda \in \mathbb{R}$**. $\blacksquare$

### 2.2 Theorem 2: Orthogonality of Eigenvectors
> **Theorem**: For a real symmetric matrix $A = A^T$, eigenvectors corresponding to **distinct eigenvalues** are **mutually orthogonal**.

#### Complete Mathematical Proof
Let $x_1$ and $x_2$ be eigenvectors corresponding to distinct eigenvalues $\lambda_1 \ne \lambda_2$:
$$A x_1 = \lambda_1 x_1 \quad \text{--- (A)}$$
$$A x_2 = \lambda_2 x_2 \quad \text{--- (B)}$$
Multiply Eq. (A) on the left by $x_2^T$:
$$x_2^T A x_1 = x_2^T (\lambda_1 x_1) = \lambda_1 (x_2^T x_1) \quad \text{--- (C)}$$
Multiply Eq. (B) on the left by $x_1^T$:
$$x_1^T A x_2 = x_1^T (\lambda_2 x_2) = \lambda_2 (x_1^T x_2) \quad \text{--- (D)}$$
Take the transpose of both sides of Eq. (C):
$$(x_2^T A x_1)^T = x_1^T A^T (x_2^T)^T = x_1^T A x_2$$
Substitute this into Eq. (D):
$$x_2^T A x_1 = \lambda_2 (x_1^T x_2) = \lambda_2 (x_2^T x_1) \quad \text{--- (E)}$$
Subtract Eq. (E) from Eq. (C):
$$\lambda_1 (x_2^T x_1) - \lambda_2 (x_2^T x_1) = 0 \implies (\lambda_1 - \lambda_2) (x_2^T x_1) = 0$$
Since the eigenvalues are distinct, $\lambda_1 - \lambda_2 \ne 0$. Dividing gives:
$$\mathbf{x_2^T x_1 = x_1 \cdot x_2 = 0} \quad \blacksquare$$

### 2.3 The Spectral Theorem (Orthogonal Diagonalization)
> **Spectral Theorem**: Every real symmetric matrix $A \in \mathbb{R}^{n \times n}$ has $n$ real eigenvalues and can be orthogonally diagonalized:
> $$\mathbf{A = Q \Lambda Q^T}$$
> where:
> - $Q = [q_1, q_2, \dots, q_n]$ is an **orthogonal matrix** ($Q^T Q = Q Q^T = I_n$) whose columns are normalized eigenvectors ($q_i^T q_j = \delta_{ij}$).
> - $\Lambda = \text{diag}(\lambda_1, \lambda_2, \dots, \lambda_n)$ is the diagonal matrix of real eigenvalues.

Spectral decomposition allows expressing $A$ as a sum of rank-1 projections:
$$A = \sum_{i=1}^n \lambda_i q_i q_i^T$$

### 2.4 Rayleigh Quotient & Min-Max Theorem
For any non-zero vector $x \in \mathbb{R}^n$, the **Rayleigh Quotient** $R(x)$ is defined as:
$$R(x) = \frac{x^T A x}{x^T x}$$
- **Rayleigh-Ritz Theorem**:
  $$\mathbf{\lambda_{\min} \le \frac{x^T A x}{x^T x} \le \lambda_{\max}}$$
  $$\lambda_{\min} = \min_{x \ne \mathbf{0}} R(x), \quad \lambda_{\max} = \max_{x \ne \mathbf{0}} R(x)$$

---

## 3. Positive Definite Matrices & Quadratic Forms

### 3.1 Calculus Motivation: Minima of a Function (Slide 387)
Consider a general quadratic surface in two variables:
$$f(x, y) = ax^2 + bxy + cy^2$$
Positive definiteness physically means that $f(x, y) > 0$ for all $(x, y) \ne (0, 0)$, meaning the origin $(0, 0)$ is a **strict global minimum** (a bowl-shaped paraboloid opening upwards).

At the critical point $(x = 0, y = 0)$, the gradients vanish:
$$\frac{\partial f}{\partial x} = 2ax + by = 0, \quad \frac{\partial f}{\partial y} = bx + 2cy = 0$$
From multivariable calculus, for the critical point to be a strict local minimum:
1. $\frac{\partial^2 f}{\partial x^2} = 2a > 0 \implies a > 0$
2. $\frac{\partial^2 f}{\partial y^2} = 2c > 0 \implies c > 0$
3. The determinant of the Hessian matrix must be strictly positive:
   $$\left[ \frac{\partial^2 f}{\partial x^2} \right] \left[ \frac{\partial^2 f}{\partial y^2} \right] > \left[ \frac{\partial^2 f}{\partial x \partial y} \right]^2 \implies (2a)(2c) > b^2 \implies \mathbf{4ac > b^2}$$

In matrix form (Slides 388):
$$f(x, y) = \begin{bmatrix} x & y \end{bmatrix} \begin{bmatrix} a & b/2 \\ b/2 & c \end{bmatrix} \begin{bmatrix} x \\ y \end{bmatrix} > 0$$
In higher dimensions $\mathbb{R}^n$:
$$X^T A X = \sum_{i=1}^n \sum_{j=1}^n a_{ij} x_i x_j > 0 \quad \text{for all non-zero } X \in \mathbb{R}^n$$

### 3.2 Classification of Real Symmetric Matrices
1. **Positive Definite (SPD)**:
   $$x^T A x > 0 \quad \text{for all } x \ne \mathbf{0}$$
2. **Positive Semi-Definite (SPSD)**:
   $$x^T A x \ge 0 \quad \text{for all } x$$
3. **Negative Definite**:
   $$x^T A x < 0 \quad \text{for all } x \ne \mathbf{0}$$
4. **Indefinite**:
   $x^T A x > 0$ for some vectors, and $x^T A x < 0$ for others (saddle-point behavior).

### 3.3 Five Equivalent Characterizations of an SPD Matrix (Slide 389)
In an exam, you can prove a real symmetric matrix $A = A^T$ is Positive Definite using any of the following necessary and sufficient tests:
1. **Quadratic Form Test**: $X^T A X > 0$ for all non-zero real vectors $X$.
2. **Eigenvalue Test**: All eigenvalues are strictly positive: $\lambda_i > 0$ for all $i = 1, \dots, n$.
3. **Sylvester’s Criterion (Leading Principal Minors)**: All upper-left sub-matrices $A_k$ have positive determinants:
   $$\det(A_1) = a_{11} > 0, \quad \det(A_2) = \det \begin{bmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{bmatrix} > 0, \quad \dots, \quad \det(A_n) = \det(A) > 0$$
4. **Pivots in Gaussian Elimination**: All $n$ pivots (without row exchanges) satisfy $d_k > 0$.
5. **Factorization Characterization**: A real symmetric matrix $A$ is positive definite if and only if there exists a matrix $R$ with **linearly independent columns** such that:
   $$\mathbf{A = R^T R}$$
   (When $R$ is chosen upper triangular, this corresponds to the Cholesky factorization $A = R^T R = L L^T$).
6. **Energy Norm**: $A$ defines a valid inner product $\langle u, v \rangle_A = u^T A v$ and norm $\|x\|_A = \sqrt{x^T A x}$.

### 3.4 Conditions for Positive Semi-Definite Matrices (Slide 390)
A real symmetric matrix is **Positive Semi-Definite (SPSD)** if:
1. $X^T A X \ge 0$ for all real vectors $X$.
2. All eigenvalues satisfy $\lambda_i \ge 0$.
3. No principal sub-matrices have negative determinants (all principal minors $\ge 0$).
4. No pivots are negative ($d_k \ge 0$).
5. There exists a matrix $R$ (possibly with dependent columns) such that $A = R^T R$.

---

## 4. Singular Value Decomposition (SVD)

While eigenvalue decomposition is strictly defined for square matrices, SVD applies to **any arbitrary rectangular matrix**.

### 4.1 Theorem of SVD
> **Theorem**: Let $A \in \mathbb{R}^{m \times n}$ be an arbitrary real matrix with rank $r$. Then $A$ can be uniquely factored as:
> $$\mathbf{A = U \Sigma V^T}$$
> where:
> - $U \in \mathbb{R}^{m \times m}$ is an orthogonal matrix ($U^T U = I_m$); its columns $u_i$ are the **left singular vectors**.
> - $V \in \mathbb{R}^{n \times n}$ is an orthogonal matrix ($V^T V = I_n$); its columns $v_i$ are the **right singular vectors**.
> - $\Sigma \in \mathbb{R}^{m \times n}$ is a rectangular diagonal matrix containing the **singular values** $\sigma_i$ ordered non-increasingly:
>   $$\sigma_1 \ge \sigma_2 \ge \dots \ge \sigma_r > \sigma_{r+1} = \dots = \sigma_{\min(m, n)} = 0$$

```
              A             =          U           ×       Σ       ×       V^T
        ┌───────────┐           ┌─────────────┐        ┌───────┐       ┌───────────┐
      m │           │         m │             │      m │σ₁     │     n │           │
        │           │           │             │        │  σ₂   │       └───────────┘
        └───────────┘           └─────────────┘        │    0  │             n
              n                        m               └───────┘
                                                           n
```

### 4.2 Relationship Between SVD and Symmetric Eigensystems
Multiplying $A$ by its transpose reveals the origin of singular values and vectors:
1. **Computing $A^T A$ (an $n \times n$ symmetric matrix)**:
   $$A^T A = (U \Sigma V^T)^T (U \Sigma V^T) = V \Sigma^T U^T U \Sigma V^T = V (\Sigma^T \Sigma) V^T$$
   - Columns of $V$ ($v_i$) are the **eigenvectors of $A^T A$**.
   - The non-zero eigenvalues of $A^T A$ are the **squares of the singular values**:
     $$\lambda_i(A^T A) = \sigma_i^2 \implies \mathbf{\sigma_i = \sqrt{\lambda_i(A^T A)}}$$

2. **Computing $A A^T$ (an $m \times m$ symmetric matrix)**:
   $$A A^T = (U \Sigma V^T) (U \Sigma V^T)^T = U \Sigma V^T V \Sigma^T U^T = U (\Sigma \Sigma^T) U^T$$
   - Columns of $U$ ($u_i$) are the **eigenvectors of $A A^T$**.

3. **Connection for Symmetric Positive Definite Matrices**:
   If $A = A^T$ and $A$ is positive definite, then $U = V = Q$, and the singular values are **identical to the eigenvalues**:
   $$\sigma_i = \lambda_i$$

### 4.3 Geometric Meaning of SVD
Matrix multiplication $y = A x$ can be decomposed into three geometric steps:
1. **$V^T x$**: An orthogonal coordinate rotation / reflection in $\mathbb{R}^n$.
2. **$\Sigma (V^T x)$**: Axis-aligned stretching/compression along coordinate directions by scale factors $\sigma_i$.
3. **$U (\Sigma V^T x)$**: An orthogonal coordinate rotation / reflection in $\mathbb{R}^m$.
**Geometric Summary**: SVD transforms the unit sphere $\|x\|_2 = 1$ in $\mathbb{R}^n$ into a **hyper-ellipsoid** in $\mathbb{R}^m$, where the semi-axes have lengths equal to $\sigma_i$ pointing along directions $u_i$.

### 4.4 Fundamental Remarks on SVD (Slide 392)
1. **Identity for Symmetric Positive Definite Matrices**:
   For any real SPD matrix $A$:
   $$\mathbf{A = U \Sigma V^T \quad \text{is strictly identical to} \quad A = Q \Lambda Q^T}$$
   where $U = V = Q$ and $\sigma_i = \lambda_i > 0$.
2. **Orthonormal Bases for the Four Fundamental Subspaces**:
   The singular vectors $U$ and $V$ furnish explicit orthonormal bases for all four fundamental subspaces of $A \in \mathbb{R}^{m \times n}$ with rank $r$:
   - **First $r$ columns of $U$** ($u_1, \dots, u_r$): Orthonormal basis for the **Column Space $C(A)$**.
   - **Last $m - r$ columns of $U$** ($u_{r+1}, \dots, u_m$): Orthonormal basis for the **Left Nullspace $N(A^T)$**.
   - **First $r$ columns of $V$** ($v_1, \dots, v_r$): Orthonormal basis for the **Row Space $C(A^T)$**.
   - **Last $n - r$ columns of $V$** ($v_{r+1}, \dots, v_n$): Orthonormal basis for the **Nullspace $N(A)$**.
3. **Action on Basis Vectors**:
   Matrix $A$ maps the $j$-th orthonormal right singular vector $v_j$ directly to $\sigma_j$ times the $j$-th orthonormal left singular vector $u_j$:
   $$A V = U \Sigma \implies \mathbf{A v_j = \sigma_j u_j \quad (1 \le j \le r)}$$
   For any vector in the nullspace ($j > r$), $A v_j = \mathbf{0}$.

---

## 5. Vector and Matrix Norms

### 5.1 Vector $p$-Norms
For $x \in \mathbb{R}^n$:
$$\|x\|_p = \left( \sum_{i=1}^n |x_i|^p \right)^{1/p}$$
- **$L_1$ Norm (Manhattan / Taxicab)**: $\|x\|_1 = \sum_{i=1}^n |x_i|$.
- **$L_2$ Norm (Euclidean Norm)**: $\|x\|_2 = \sqrt{\sum_{i=1}^n x_i^2} = \sqrt{x^T x}$.
- **$L_\infty$ Norm (Chebyshev / Maximum Norm)**: $\|x\|_\infty = \max_{1 \le i \le n} |x_i|$.

### 5.2 Induced (Subordinate) Matrix Norms
The induced matrix norm measures the maximum amplification factor that $A$ imparts to any vector $x$:
$$\|A\|_p = \sup_{x \ne \mathbf{0}} \frac{\|Ax\|_p}{\|x\|_p} = \max_{\|x\|_p = 1} \|Ax\|_p$$
- **1-Norm (Maximum Column-Sum Norm)**:
  $$\mathbf{\|A\|_1 = \max_{1 \le j \le n} \sum_{i=1}^m |a_{ij}|}$$
- **$\infty$-Norm (Maximum Row-Sum Norm)**:
  $$\mathbf{\|A\|_\infty = \max_{1 \le i \le m} \sum_{j=1}^n |a_{ij}|}$$
- **2-Norm (Spectral Norm)**:
  $$\mathbf{\|A\|_2 = \sigma_{\max}(A) = \sqrt{\lambda_{\max}(A^T A)}}$$
  *(For real symmetric matrices, $\|A\|_2 = \max_i |\lambda_i| = \rho(A)$).*
- **Submultiplicative Property**: For any compatible matrices and vectors:
  $$\|A x\| \le \|A\| \|x\|, \quad \|A B\| \le \|A\| \|B\|$$

---

## 6. Matrix Condition Number & Error Bounds

### 6.1 Definition of Condition Number
The **condition number** $\kappa(A)$ measures the sensitivity of the solution $x$ of $Ax = b$ to perturbations in the input data ($b$ or $A$).
$$\mathbf{\kappa(A) = \|A\| \|A^{-1}\|}$$
Under the matrix 2-norm:
$$\|A\|_2 = \sigma_{\max}, \quad \|A^{-1}\|_2 = \frac{1}{\sigma_{\min}}$$
$$\mathbf{\kappa_2(A) = \frac{\sigma_{\max}(A)}{\sigma_{\min}(A)}}$$

For Symmetric Positive Definite matrices:
$$\mathbf{\kappa_2(A) = \frac{\lambda_{\max}(A)}{\lambda_{\min}(A)}}$$

- **Range**: $\kappa(A) \ge 1$ always.
- **Well-conditioned**: $\kappa(A) \approx 1$ (orthogonal matrices have $\kappa_2(Q) = 1$, the theoretical minimum).
- **Ill-conditioned**: $\kappa(A) \gg 1$ ($\kappa \sim 10^4$ to $10^{12}$).
- **Singular Matrix**: $\sigma_{\min} = 0 \implies \kappa(A) = \infty$.

### 6.2 Rigorous Derivation of the Error Propagation Inequality
Consider the unperturbed system:
$$A x = b \quad \text{--- (1)}$$
Suppose the right-hand side vector $b$ is perturbed by a small error $\delta b$ (measurement error, round-off). The resulting solution becomes $x + \delta x$:
$$A (x + \delta x) = b + \delta b$$
By linearity:
$$A x + A \delta x = b + \delta b \implies A \delta x = \delta b$$
Since $A$ is non-singular:
$$\delta x = A^{-1} \delta b \quad \text{--- (2)}$$

Take norms on both sides of Eq. (2) and apply the submultiplicative property:
$$\|\delta x\| = \|A^{-1} \delta b\| \le \|A^{-1}\| \|\delta b\| \quad \text{--- (3)}$$
Now, from Eq. (1), take norms:
$$\|b\| = \|A x\| \le \|A\| \|x\| \implies \frac{1}{\|x\|} \le \frac{\|A\|}{\|b\|} \quad \text{--- (4)}$$
Multiply inequality (3) by inequality (4):
$$\frac{\|\delta x\|}{\|x\|} \le \|A^{-1}\| \|\delta b\| \cdot \frac{\|A\|}{\|b\|} = \left( \|A\| \|A^{-1}\| \right) \frac{\|\delta b\|}{\|b\|}$$
$$\mathbf{\frac{\|\delta x\|}{\|x\|} \le \kappa(A) \frac{\|\delta b\|}{\|b\|}}$$

### 6.3 General Perturbation Bound (Both $A$ and $b$ Perturbed)
If both the operator $A$ and the RHS $b$ suffer perturbations $(A + \delta A)(x + \delta x) = b + \delta b$, with $\|\delta A\| < 1/\|A^{-1}\|$:
$$\mathbf{\frac{\|\delta x\|}{\|x\|} \le \frac{\kappa(A)}{1 - \kappa(A) \frac{\|\delta A\|}{\|A\|}} \left( \frac{\|\delta b\|}{\|b\|} + \frac{\|\delta A\|}{\|A\|} \right)}$$

### 6.4 Practical Significance in Numerical Computation
- **Rule of Thumb (Loss of Significance)**:
  If a problem has condition number $\kappa(A) \approx 10^k$, solving the system in IEEE 754 double precision ($16\text{ decimal digits}$) can lose up to $k$ significant digits:
  $$\text{Digits of precision remaining} \approx 16 - \log_{10}(\kappa(A))$$
- **Ill-Conditioning in Discretized PDEs**:
  For the central difference Laplacian on a grid with spacing $h$:
  $$\lambda_{\min} \sim \pi^2, \quad \lambda_{\max} \sim \frac{4}{h^2} \implies \mathbf{\kappa_2(A) \sim O\left(\frac{1}{h^2}\right) = O(N)}$$
  As the mesh is refined to improve accuracy ($h \to 0$), the condition number blows up quadratically! This makes iterative solvers stall, directly necessitating **preconditioning**.
