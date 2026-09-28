# Chapter 04: Direct Solvers, Factorization Techniques & the Thomas Algorithm

---

## 1. Limitations of Direct Solvers in Large-Scale Scientific Computing

A **direct solver** computes the exact solution to a linear system $Ax = b$ in a finite, predetermined sequence of arithmetic operations (in the absence of round-off error).

### 1.1 The Classical Direct Solvers
- **Gaussian Elimination**: Converts $A$ into an upper triangular system $U x = c$ via elementary row operations, followed by backward substitution.
- **LU Factorization**: Decomposes $A = L U$, where $L$ is unit lower triangular and $U$ is upper triangular.
- **Cholesky Factorization**: Decomposes symmetric positive definite $A = L L^T$.
- **Cramer’s Rule**: Explicit determinant-based formula.

### 1.2 The Three Fundamental Failures of Direct Solvers for Large PDEs
1. **Prohibitive Computational Complexity**:
   - For a general dense $N \times N$ matrix, Gaussian elimination and LU decomposition require:
     $$\text{FLOPs} = \frac{2}{3} N^3 + O(N^2)$$
   - In a 3D physical simulation on an $M \times M \times M$ grid, the number of unknowns is $N = M^3$.
   - For a modest grid of $M = 100 \implies N = 10^6$ unknowns:
     $$\text{FLOPs} = \frac{2}{3} (10^6)^3 = \frac{2}{3} \times 10^{18}\text{ FLOPs} \approx 6.7 \times 10^{17}\text{ FLOPs}$$
   - On a fast workstation ($10\text{ GFLOPS}$), this single solve takes $> 2\text{ years}$!

2. **The "Fill-In" Phenomenon and Memory Explosion**:
   - The original PDE coefficient matrix $A$ is **sparse** (only $O(N)$ non-zero entries).
   - During Gaussian elimination, subtracting row multiples introduces non-zero entries into positions that were originally zero. This is known as **fill-in**.
   - For a 2D grid ($N = M^2$, bandwidth $W \sim M = \sqrt{N}$), the banded LU factors suffer fill-in within the band, requiring $O(N^{3/2})$ storage and $O(N^2)$ operations.
   - For a 3D grid ($N = M^3$, bandwidth $W \sim M^2 = N^{2/3}$), fill-in requires $O(N^{5/3})$ memory and $O(N^{7/3})$ operations!
   - This destroys the $O(N)$ memory advantage of sparse discretization.

3. **Accumulation of Finite-Precision Round-Off Errors**:
   - As $N$ grows, error in each elimination step accumulates over $O(N^3)$ operations, leading to severe numerical loss of significance unless expensive pivoting strategies are employed.

---

## 2. Gaussian Elimination & LU Decomposition

### 2.1 Standard LU Decomposition ($A = LU$)
Any non-singular matrix $A$ (admitting non-zero pivots without row interchanges) can be uniquely factored into:
$$A = L U$$
where:
$$L = \begin{bmatrix} 1 & 0 & \dots & 0 \\ l_{21} & 1 & \dots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ l_{n1} & l_{n2} & \dots & 1 \end{bmatrix}, \quad U = \begin{bmatrix} u_{11} & u_{12} & \dots & u_{1n} \\ 0 & u_{22} & \dots & u_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & u_{nn} \end{bmatrix}$$
- **Doolittle Algorithm**: $l_{ii} = 1$ (unit lower triangular).
- **Crout Algorithm**: $u_{ii} = 1$ (unit upper triangular).

### 2.2 Solving $Ax = b$ via Two-Stage Triangular Substitution
Once $A = LU$ is established:
1. **Forward Substitution** ($L y = b$):
   $$y_1 = b_1$$
   $$y_i = b_i - \sum_{j=1}^{i-1} l_{ij} y_j, \quad i = 2, 3, \dots, n \quad \left[\text{Work: } \frac{n^2}{2}\text{ FLOPs}\right]$$
2. **Backward Substitution** ($U x = y$):
   $$x_n = \frac{y_n}{u_{nn}}$$
   $$x_i = \frac{1}{u_{ii}} \left( y_i - \sum_{j=i+1}^n u_{ij} x_j \right), \quad i = n-1, n-2, \dots, 1 \quad \left[\text{Work: } \frac{n^2}{2}\text{ FLOPs}\right]$$

### 2.3 Multiple Right-Hand Sides
- In time-dependent or parametric problems, the operator $A$ remains fixed while the source term $b(t)$ changes.
- The expensive $LU$ factorization ($\frac{2}{3} N^3$) is performed **once**.
- Subsequent solutions require only forward and backward substitution ($O(N^2)$ FLOPs per time-step).

### 2.4 Partial Pivoting ($P A = L U$)
- If a pivot $|a_{kk}^{(k-1)}| \approx 0$, division by near-zero causes catastrophic round-off error amplification.
- **Partial Pivoting**: At step $k$, search the current column $k$ for the element with maximum absolute magnitude:
  $$p = \arg\max_{i \ge k} |a_{ik}^{(k-1)}|$$
- Swap row $k$ with row $p$. Encoded as permutation matrix $P$:
  $$P A = L U$$
- **HPC Drawback of Pivoting**: Row interchanges introduce global dependencies, hinder fine-grained parallelization, and ruin banded matrix storage.

---

## 3. Cholesky Factorization for Symmetric Positive Definite Matrices

When matrix $A$ is **Symmetric Positive Definite (SPD)** ($A = A^T$ and $x^T A x > 0$ for all $x \ne \mathbf{0}$):
- Pivoting is **never required**! The diagonal entries remain strictly positive throughout elimination ($a_{kk}^{(k-1)} > 0$).
- $A$ can be factored into:
  $$\mathbf{A = L L^T}$$
  where $L$ is a real lower triangular matrix with strictly positive diagonal elements.

### 3.1 Derivation of the Cholesky Algorithm
Equating entries of $A = L L^T$:
$$a_{ij} = \sum_{k=1}^n l_{ik} l_{jk} = \sum_{k=1}^{\min(i, j)} l_{ik} l_{jk}$$

1. **For Diagonal Elements ($i = j$)**:
   $$a_{jj} = \sum_{k=1}^{j-1} l_{jk}^2 + l_{jj}^2 \implies \mathbf{l_{jj} = \sqrt{a_{jj} - \sum_{k=1}^{j-1} l_{jk}^2}}$$
2. **For Off-Diagonal Elements ($i > j$)**:
   $$a_{ij} = \sum_{k=1}^{j-1} l_{ik} l_{jk} + l_{ij} l_{jj} \implies \mathbf{l_{ij} = \frac{1}{l_{jj}} \left( a_{ij} - \sum_{k=1}^{j-1} l_{ik} l_{jk} \right)}$$

### 3.2 Computational Savings
- **FLOP Count**: $\frac{1}{3} N^3$ FLOPs — **exactly half** the operations of general LU factorization!
- **Memory Footprint**: Only the lower triangular part of $L$ ($N(N+1)/2$ words) needs to be computed and stored.

---

## 4. Tridiagonal Matrix Algorithm (TDMA / Thomas Algorithm)

For 1D boundary value problems (or directional sweeps in multidimensional ADI methods), the discretization matrix is **tridiagonal**. The Thomas Algorithm is a streamlined, specialized Gaussian elimination taking advantage of zero entries.

### 4.1 Problem Formulation
Consider the $N \times N$ tridiagonal linear system:
$$\begin{bmatrix} b_1 & c_1 & 0 & 0 & \dots & 0 \\ a_2 & b_2 & c_2 & 0 & \dots & 0 \\ 0 & a_3 & b_3 & c_3 & \dots & 0 \\ \vdots & \ddots & \ddots & \ddots & \ddots & \vdots \\ 0 & \dots & 0 & a_{N-1} & b_{N-1} & c_{N-1} \\ 0 & \dots & 0 & 0 & a_N & b_N \end{bmatrix} \begin{bmatrix} T_1 \\ T_2 \\ T_3 \\ \vdots \\ T_{N-1} \\ T_N \end{bmatrix} = \begin{bmatrix} d_1 \\ d_2 \\ d_3 \\ \vdots \\ d_{N-1} \\ d_N \end{bmatrix}$$
where:
- $b_i$ are the main diagonal elements.
- $c_i$ are the super-diagonal elements ($c_N = 0$).
- $a_i$ are the sub-diagonal elements ($a_1 = 0$).
- $d_i$ are the right-hand side source values.

Component-wise equation $i$:
$$a_i T_{i-1} + b_i T_i + c_i T_{i+1} = d_i, \quad i = 1, 2, \dots, N$$

### 4.2 Complete Mathematical Derivation
1. **First Row ($i = 1$)**:
   $$b_1 T_1 + c_1 T_2 = d_1 \implies T_1 = \frac{d_1}{b_1} - \frac{c_1}{b_1} T_2$$
   Express this in the canonical form:
   $$\mathbf{T_1 = \gamma_1 - \frac{c_1}{\beta_1} T_2}$$
   where:
   $$\mathbf{\beta_1 = b_1, \quad \gamma_1 = \frac{d_1}{b_1} = \frac{d_1}{\beta_1}}$$

2. **Inductive Hypothesis**:
   Assume that for row $i-1$, $T_{i-1}$ has already been expressed in terms of $T_i$:
   $$T_{i-1} = \gamma_{i-1} - \frac{c_{i-1}}{\beta_{i-1}} T_i \quad \text{--- (Eq. A)}$$

3. **Elimination at Row $i$**:
   Substitute Eq. (A) into the $i$-th equation:
   $$a_i \left( \gamma_{i-1} - \frac{c_{i-1}}{\beta_{i-1}} T_i \right) + b_i T_i + c_i T_{i+1} = d_i$$
   Grouping terms containing $T_i$:
   $$\left( b_i - \frac{a_i c_{i-1}}{\beta_{i-1}} \right) T_i + c_i T_{i+1} = d_i - a_i \gamma_{i-1}$$
   Isolating $T_i$:
   $$T_i = \frac{d_i - a_i \gamma_{i-1}}{b_i - \frac{a_i c_{i-1}}{\beta_{i-1}}} - \frac{c_i}{b_i - \frac{a_i c_{i-1}}{\beta_{i-1}}} T_{i+1} \quad \text{--- (Eq. B)}$$

4. **Derivation of the Recursive Relations**:
   Matching Eq. (B) with the canonical recurrence $T_i = \gamma_i - \frac{c_i}{\beta_i} T_{i+1}$:
   $$\mathbf{\beta_i = b_i - \frac{a_i c_{i-1}}{\beta_{i-1}}, \quad i = 2, 3, \dots, N}$$
   $$\mathbf{\gamma_i = \frac{d_i - a_i \gamma_{i-1}}{\beta_i}, \quad i = 2, 3, \dots, N}$$

5. **Terminal Row ($i = N$)**:
   At the $N$-th row, $c_N = 0$:
   $$T_N = \gamma_N - \frac{c_N}{\beta_N} T_{N+1} = \gamma_N - 0 \implies \mathbf{T_N = \gamma_N}$$

6. **Backward Substitution Sweep**:
   Knowing $T_N$, compute preceding unknowns backwards from $N-1$ down to 1:
   $$\mathbf{T_i = \gamma_i - \left(\frac{c_i}{\beta_i}\right) T_{i+1}, \quad i = N-1, N-2, \dots, 1}$$

### 4.3 TDMA Algorithm Summary
```
Algorithm: Tridiagonal Matrix Algorithm (TDMA / Thomas)
Input: Vectors a (sub-diag), b (main-diag), c (super-diag), d (RHS) of length N
Output: Solution vector T

1. Initialize forward elimination:
   β[1] = b[1]
   γ[1] = d[1] / β[1]

2. Forward sweep (elimination of sub-diagonal):
   for i = 2 to N do
       β[i] = b[i] - (a[i] * c[i-1]) / β[i-1]
       γ[i] = (d[i] - a[i] * γ[i-1]) / β[i]
   end for

3. Terminal boundary condition:
   T[N] = γ[N]

4. Backward sweep (back substitution):
   for i = N-1 down to 1 do
       T[i] = γ[i] - (c[i] / β[i]) * T[i+1]
   end for
```

### 4.4 Complexity and Performance Comparison
- **Operation Count**:
  - Forward sweep: $2$ multiplications, $2$ divisions, $2$ subtractions per node $\approx 5N$ operations.
  - Backward sweep: $1$ multiplication, $1$ subtraction per node $\approx 2N$ operations.
  - **Total Work**: $\mathbf{\approx 8N \text{ operations} = O(N)}$!
- **Concrete Example Comparison (From Slide 285)**:
  - For an $N = 5$ matrix:
    - **TDMA**: Takes $\mathbf{10\text{ operations}}$!
    - **Gaussian Elimination**: Takes $\mathbf{75\text{ operations}}$!
    - Speedup factor for $N = 5$ is $7.5\times$. For $N = 1000$, speedup is $> 80,000\times$!

### 4.5 Condition for Numerical Stability
The Thomas algorithm does not use pivoting. It is guaranteed to be unconditionally stable against round-off error accumulation if the matrix is **strictly or irreducibly diagonally dominant**:
$$|b_i| \ge |a_i| + |c_i| \quad \text{for all } i, \text{ with } |b_i| > |a_i| + |c_i| \text{ for at least one } i$$
Under this condition, $\beta_i \ne 0$ for all $i$, guaranteeing that division by zero never occurs.

---

## 5. The Sherman-Morrison Formula for Perturbed Systems

### 5.1 Theorem Statement
Let $A \in \mathbb{R}^{n \times n}$ be an invertible square matrix, and let $u, v \in \mathbb{R}^n$ be column vectors. The matrix $A + u v^T$ is a **rank-1 update (perturbation)** of $A$.
If $1 + v^T A^{-1} u \ne 0$, then $A + u v^T$ is invertible, and its inverse is given exactly by:
$$\mathbf{(A + u v^T)^{-1} = A^{-1} - \frac{A^{-1} u v^T A^{-1}}{1 + v^T A^{-1} u}}$$

### 5.2 Mathematical Proof
We multiply $(A + u v^T)$ by the proposed inverse and verify that it equals the identity matrix $I$:
$$(A + u v^T) \left[ A^{-1} - \frac{A^{-1} u v^T A^{-1}}{1 + v^T A^{-1} u} \right]$$
Expanding:
$$= A A^{-1} - \frac{A A^{-1} u v^T A^{-1}}{1 + v^T A^{-1} u} + u v^T A^{-1} - \frac{u v^T A^{-1} u v^T A^{-1}}{1 + v^T A^{-1} u}$$
Notice that $v^T A^{-1} u$ is a scalar quantity, denoted $\alpha = v^T A^{-1} u$.
$$= I - \frac{u v^T A^{-1}}{1 + \alpha} + u v^T A^{-1} - \frac{u \alpha (v^T A^{-1})}{1 + \alpha}$$
Factoring out $u v^T A^{-1}$:
$$= I + u v^T A^{-1} \left[ 1 - \frac{1}{1 + \alpha} - \frac{\alpha}{1 + \alpha} \right] = I + u v^T A^{-1} \left[ 1 - \frac{1 + \alpha}{1 + \alpha} \right] = I + u v^T A^{-1} (0) = I \quad \blacksquare$$

### 5.3 Solving Perturbed Systems Without Re-Factoring
To solve $(A + u v^T) x = b$:
1. Solve $A y = b$ (cost of 1 solver call).
2. Solve $A z = u$ (cost of 1 solver call).
3. Compute scalars: $v^T y$ and $v^T z$.
4. Then $x$ is obtained by a simple vector update:
   $$\mathbf{x = y - \left( \frac{v^T y}{1 + v^T z} \right) z}$$
This requires **zero modification** to the original solver of $A$!

### 5.4 Application to Cyclic Tridiagonal Matrices
In physical systems with **periodic boundary conditions** ($T(0) = T(L)$):
- The boundary equations connect node 1 to node $N$, placing non-zero entries in the top-right ($c_N$) and bottom-left ($a_1$) corners.
- This breaks the standard tridiagonal band:
  $$A_{cyclic} = \begin{bmatrix} b_1 & c_1 & 0 & \dots & a_1 \\ a_2 & b_2 & c_2 & \dots & 0 \\ \vdots & \ddots & \ddots & \ddots & \vdots \\ c_N & \dots & 0 & a_N & b_N \end{bmatrix} = A_{tridiag} + u v^T$$
- Using the Sherman-Morrison formula, $A_{cyclic} x = b$ can be solved via **two standard TDMA passes** on $A_{tridiag}$, preserving $O(N)$ computational complexity without expensive general elimination!

---

## 6. Alternating Direction Implicit (ADI) Method

For 2D elliptic PDEs ($\nabla^2 T = 0$), direct TDMA cannot be applied directly because the 5-point stencil produces a pentadiagonal matrix with wide band offsets $\pm N_x$. The **ADI method** circumvents this by splitting each multidimensional iteration into two alternating 1D directional steps.

```
       Step 1: Horizontal Sweep (x-sweep)         Step 2: Vertical Sweep (y-sweep)
         Implicit along x, Explicit in y            Implicit along y, Explicit in x
              j=3 ──► TDMA along row 3                    i=1     i=2     i=3
              j=2 ──► TDMA along row 2                     ▲       ▲       ▲
              j=1 ──► TDMA along row 1                     │       │       │
                                                         TDMA    TDMA    TDMA
```

### 6.1 Mathematical Formulation of ADI Sweeps
Discretizing the 2D equation:
$$\frac{T_{i-1, j} - 2 T_{i, j} + T_{i+1, j}}{\Delta x^2} + \frac{T_{i, j-1} - 2 T_{i, j} + T_{i, j+1}}{\Delta y^2} = 0$$

1. **Sweep 1 ($x$-Direction Implicit Sweep)**:
   Treat the $x$-derivatives implicitly at intermediate state $T^{(k+1/2)}$ while evaluating $y$-derivatives explicitly using available values $T^{(k)}$:
   $$\frac{T_{i-1, j}^{(k+1/2)} - 2 T_{i, j}^{(k+1/2)} + T_{i+1, j}^{(k+1/2)}}{\Delta x^2} = -\left[ \frac{T_{i, j-1}^{(k)} - 2 T_{i, j}^{(k)} + T_{i, j+1}^{(k)}}{\Delta y^2} \right]$$
   - Along each horizontal grid line $j$, this forms an independent **1D tridiagonal system** in unknowns $\{T_{i, j}^{(k+1/2)}\}_{i=1}^{N_x}$.
   - Solved in $O(N_x)$ operations per row using TDMA!

2. **Sweep 2 ($y$-Direction Implicit Sweep)**:
   Treat the $y$-derivatives implicitly at final state $T^{(k+1)}$ using updated intermediate values $T^{(k+1/2)}$:
   $$\frac{T_{i, j-1}^{(k+1)} - 2 T_{i, j}^{(k+1)} + T_{i, j+1}^{(k+1)}}{\Delta y^2} = -\left[ \frac{T_{i-1, j}^{(k+1/2)} - 2 T_{i, j}^{(k+1/2)} + T_{i+1, j}^{(k+1/2)}}{\Delta x^2} \right]$$
   - Along each vertical grid line $i$, this forms an independent **1D tridiagonal system** in unknowns $\{T_{i, j}^{(k+1)}\}_{j=1}^{N_y}$.
   - Solved in $O(N_y)$ operations per column using TDMA!

### 6.2 Key Advantages for HPC
- **Optimal Complexity**: Each full ADI iteration requires $O(N_x N_y) = O(N)$ arithmetic operations.
- **Inherent Parallelism**: All horizontal rows in Sweep 1 are completely decoupled and can be solved simultaneously on separate processor cores! Likewise, all vertical columns in Sweep 2 are completely decoupled.
