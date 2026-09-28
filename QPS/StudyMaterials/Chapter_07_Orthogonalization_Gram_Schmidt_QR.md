# Chapter 07: Orthogonalization, the Gram-Schmidt Process & QR Factorization

---

## 1. Orthogonality and Orthonormality

In numerical linear algebra, algorithms operating on orthogonal vectors are the gold standard because they do not amplify round-off errors and preserve vector norms.

### 1.1 Definitions
Let $u, v \in \mathbb{R}^n$.
- **Inner Product (Dot Product)**:
  $$\langle u, v \rangle = u^T v = \sum_{i=1}^n u_i v_i$$
- **Orthogonality**: Two vectors $u$ and $v$ are **orthogonal** ($u \perp v$) if:
  $$u^T v = 0$$
- **Orthonormality**: A set of vectors $\{q_1, q_2, \dots, q_k\}$ is **orthonormal** if:
  $$q_i^T q_j = \delta_{ij} = \begin{cases} 1 & \text{if } i = j \\ 0 & \text{if } i \ne j \end{cases}$$
  where $\delta_{ij}$ is the Kronecker delta.

### 1.2 Properties of Orthogonal Matrices
A square matrix $Q \in \mathbb{R}^{n \times n}$ whose columns form an orthonormal basis is called an **orthogonal matrix**:
$$Q = [q_1, q_2, \dots, q_n] \implies Q^T Q = I_n \iff \mathbf{Q^{-1} = Q^T}$$

#### Fundamental Invariances of Orthogonal Matrices
1. **Preservation of the 2-Norm (Isometry)**:
   $$\|Q x\|_2^2 = (Q x)^T (Q x) = x^T (Q^T Q) x = x^T I x = x^T x = \|x\|_2^2 \implies \mathbf{\|Q x\|_2 = \|x\|_2}$$
   Multiplying a vector by an orthogonal matrix causes zero distortion of its Euclidean length.
2. **Preservation of Angles and Inner Products**:
   $$\langle Q x, Q y \rangle = (Q x)^T (Q y) = x^T Q^T Q y = x^T y = \langle x, y \rangle$$
3. **Optimal Condition Number**:
   All singular values of $Q$ are identically 1 ($\sigma_1 = \dots = \sigma_n = 1$). Therefore:
   $$\mathbf{\kappa_2(Q) = \frac{\sigma_{\max}}{\sigma_{\min}} = \frac{1}{1} = 1}$$
   Orthogonal transformations are immune to numerical ill-conditioning!

---

## 2. Classical Gram-Schmidt (CGS) Algorithm

Given a set of linearly independent vectors $\{a_1, a_2, \dots, a_n\} \subset \mathbb{R}^m$, the Gram-Schmidt process constructs an orthonormal basis $\{q_1, q_2, \dots, q_n\}$ that spans the exact same hierarchy of subspaces:
$$\text{span}\{q_1, \dots, q_k\} = \text{span}\{a_1, \dots, a_k\}, \quad \text{for all } k = 1, \dots, n$$

### 2.1 Derivation of CGS
1. **Step 1**:
   Normalize $a_1$:
   $$v_1 = a_1, \quad r_{11} = \|v_1\|_2, \quad \mathbf{q_1 = \frac{v_1}{r_{11}}}$$
2. **Step 2**:
   Subtract the orthogonal projection of $a_2$ onto $q_1$:
   $$\text{proj}_{q_1}(a_2) = (q_1^T a_2) q_1$$
   $$v_2 = a_2 - (q_1^T a_2) q_1$$
   $$r_{22} = \|v_2\|_2, \quad \mathbf{q_2 = \frac{v_2}{r_{22}}}$$
3. **Step $j$ ($j = 2, 3, \dots, n$)**:
   Subtract the projections of $a_j$ onto all previously determined basis vectors $q_1, \dots, q_{j-1}$:
   $$r_{ij} = q_i^T a_j \quad (i = 1, 2, \dots, j-1)$$
   $$\mathbf{v_j = a_j - \sum_{i=1}^{j-1} r_{ij} q_i}$$
   $$r_{jj} = \|v_j\|_2, \quad \mathbf{q_j = \frac{v_j}{r_{jj}}}$$

### 2.2 Classical Gram-Schmidt Pseudocode
```
Algorithm: Classical Gram-Schmidt (CGS)
Input: Linearly independent vectors a_1, a_2, ..., a_n
Output: Orthonormal vectors q_1, q_2, ..., q_n and upper triangular entries r_ij

for j = 1 to n do
    v_j = a_j
    for i = 1 to j-1 do
        r[i, j] = q_i^T * a_j       // Dot product with original a_j
        v_j = v_j - r[i, j] * q_i   // Subtract projection
    end for
    r[j, j] = norm(v_j)
    if r[j, j] == 0 then
        error("Vectors are linearly dependent!")
    end if
    q_j = v_j / r[j, j]
end for
```

---

## 3. Loss of Orthogonality in CGS & The Modified Gram-Schmidt (MGS) Method

### 3.1 The Catastrophic Failure of CGS in Floating-Point Arithmetic
In theoretical exact arithmetic, CGS is mathematically flawless. However, on finite-precision digital computers (IEEE 754):
- When vector $a_j$ is nearly collinear with the subspace spanned by $\{q_1, \dots, q_{j-1}\}$, the vector difference:
  $$v_j = a_j - \sum_{i=1}^{j-1} (q_i^T a_j) q_i$$
  incurs severe **subtractive cancellation**.
- Small round-off errors committed during early projections are amplified.
- As $j$ increases, the newly generated vector $q_j$ catastrophically loses orthogonality with respect to early vectors:
  $$|q_1^T q_j| \gg \epsilon_{\text{machine}}$$
  In pathological cases, $q_j$ can become completely non-orthogonal ($q_1^T q_j \sim 0.1$ or higher!).

```
Loss of Orthogonality Progression:
  Exact Math:           q_1 ┴ q_2 ┴ q_3 ┴ q_4 ... (perfect 90° angles)
  CGS with Round-off:   q_1 ┴ q_2, but q_4 drifts, eventually angle(q_1, q_4) << 90°!
```

### 3.2 The Modified Gram-Schmidt (MGS) Formulation
To prevent cancellation errors from corrupting the projections, MGS rearranges the order of projections.
Instead of projecting the original vector $a_j$ against all $q_i$, MGS takes the vector and **sequentially projects it onto each $q_i$, immediately replacing the vector with its updated orthogonalized residual**:

1. Initialize working vectors:
   $$v_j^{(1)} = a_j \quad \text{for } j = 1, 2, \dots, n$$
2. At step $i$ ($i = 1, \dots, n$):
   Normalize the pivot vector:
   $$r_{ii} = \|v_i^{(i)}\|_2, \quad \mathbf{q_i = \frac{v_i^{(i)}}{r_{ii}}}$$
   Immediately project **all remaining unprocessed vectors** $v_j$ ($j = i+1, \dots, n$) onto the newly completed $q_i$:
   $$r_{ij} = q_i^T v_j^{(i)}$$
   $$\mathbf{v_j^{(i+1)} = v_j^{(i)} - r_{ij} q_i, \quad j = i+1, \dots, n}$$

### 3.3 Modified Gram-Schmidt Pseudocode
```
Algorithm: Modified Gram-Schmidt (MGS)
Input: Vectors a_1, a_2, ..., a_n
Output: Orthonormal basis q_1, q_2, ..., q_n and upper triangular entries r_ij

for i = 1 to n do
    v_i = a_i                     // Working vector initialization
end for

for i = 1 to n do
    r[i, i] = norm(v_i)
    q_i = v_i / r[i, i]
    for j = i+1 to n do
        r[i, j] = q_i^T * v_j     // Dot product with UPDATED v_j!
        v_j = v_j - r[i, j] * q_i // Update v_j immediately
    end for
end for
```

### 3.4 Rigorous Comparison: CGS vs. MGS
| Property | Classical Gram-Schmidt (CGS) | Modified Gram-Schmidt (MGS) |
| :--- | :--- | :--- |
| **Exact Arithmetic** | Identical output to MGS | Identical output to CGS |
| **Floating-Point Stability**| **Numerically Unstable** (severe loss of orthogonality) | **Numerically Stable** ($\|I - Q^T Q\|_2 \approx O(\kappa(A) \epsilon)$) |
| **Projection Vector** | Uses raw vector $a_j$ for all dot products | Uses intermediate orthogonalized vector $v_j^{(i)}$ |
| **Parallel Efficiency** | Higher (dot products can be computed as matrix-vector $Q^T a_j$) | Lower (sequential dependencies between outer steps) |
| **Role in Course** | Pedagogical baseline | **Engine of Arnoldi & GMRES algorithms** |

---

## 4. The QR Factorization ($A = QR$)

The Gram-Schmidt orthogonalization can be expressed cleanly as a matrix factorization.

### 4.1 Matrix Representation
From the Gram-Schmidt relations:
$$a_1 = r_{11} q_1$$
$$a_2 = r_{12} q_1 + r_{22} q_2$$
$$a_3 = r_{13} q_1 + r_{23} q_2 + r_{33} q_3$$
$$\vdots$$
$$a_n = \sum_{i=1}^n r_{in} q_i$$

In matrix form:
$$[a_1, a_2, \dots, a_n] = [q_1, q_2, \dots, q_n] \begin{bmatrix} r_{11} & r_{12} & \dots & r_{1n} \\ 0 & r_{22} & \dots & r_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & r_{nn} \end{bmatrix}$$
$$\mathbf{A = Q R}$$
where:
- $Q \in \mathbb{R}^{m \times n}$ has orthonormal columns ($Q^T Q = I_n$).
- $R \in \mathbb{R}^{n \times n}$ is **upper triangular** with strictly positive diagonal entries ($r_{ii} > 0$).

### 4.2 Solving $Ax = b$ via QR Factorization
For a square invertible matrix $A \in \mathbb{R}^{n \times n}$:
$$A x = b \implies (Q R) x = b$$
Multiply by $Q^T$ on the left (since $Q^{-1} = Q^T$):
$$Q^T Q R x = Q^T b \implies \mathbf{R x = Q^T b}$$
Since $R$ is upper triangular, $x$ is solved immediately by backward substitution in $O(n^2)$ FLOPs without any matrix inversion!

### 4.3 Solving Linear Least-Squares Problems
For an overdetermined rectangular system $A \in \mathbb{R}^{m \times n}$ ($m > n$):
$$\min_{x \in \mathbb{R}^n} \|Ax - b\|_2$$
- **Standard Normal Equations**:
  $$A^T A x = A^T b$$
  *Severe Hazard*: The condition number of $A^T A$ is squared:
  $$\kappa_2(A^T A) = [\kappa_2(A)]^2$$
  If $\kappa(A) = 10^5$, then $\kappa(A^T A) = 10^{10}$, causing severe loss of accuracy!
- **The QR Approach**:
  $$\|Ax - b\|_2 = \|Q R x - b\|_2 = \|Q (R x - Q^T b)\|_2 = \|R x - Q^T b\|_2$$
  The least-squares solution is obtained by solving the triangular system:
  $$\mathbf{R x = Q^T b}$$
  *Advantage*: Solves the least-squares problem with condition number $\kappa_2(A)$, **completely avoiding condition number squaring**!

---

## 5. Alternative Orthogonalization Methods

While Modified Gram-Schmidt computes $Q$ column-by-column, other methods compute QR factorization via orthogonal matrix transformations:

### 5.1 Householder Reflections
A **Householder reflection** reflects any vector $x \in \mathbb{R}^m$ across a hyperplane orthogonal to unit vector $v$:
$$H = I - 2 \frac{v v^T}{v^T v}$$
- $H$ is symmetric ($H = H^T$) and orthogonal ($H^T H = I$).
- **Action**: Can transform an arbitrary vector $x$ so that all entries below the first entry are zeroed out simultaneously:
  $$H x = [\pm \|x\|_2, 0, 0, \dots, 0]^T$$
- Applied repeatedly, Householder reflections reduce $A$ to upper triangular form $R$:
  $$H_n \dots H_2 H_1 A = R \implies Q = H_1 H_2 \dots H_n$$

### 5.2 Givens Rotations
A **Givens plane rotation** $\Omega(i, j, \theta)$ rotates vectors in the 2D plane spanned by coordinate axes $i$ and $j$:
$$\Omega = \begin{bmatrix} 1 & & & & & \\ & \ddots & & & & \\ & & c & s & & \\ & & -s & c & & \\ & & & & \ddots & \\ & & & & & 1 \end{bmatrix}, \quad \text{where } c = \cos\theta, \, s = \sin\theta, \, c^2 + s^2 = 1$$
- **Action**: Selects $c$ and $s$ to zero out **one single entry** at a specific row and column $(j, i)$:
  $$c = \frac{x_i}{\sqrt{x_i^2 + x_j^2}}, \quad s = \frac{x_j}{\sqrt{x_i^2 + x_j^2}}$$
- **Significance in GMRES**: In GMRES, an upper Hessenberg matrix has non-zeros only on the subdiagonal $h_{i+1, i}$. Applying a sequence of $m$ Givens rotations zeroes out each subdiagonal element sequentially, transforming the Hessenberg least-squares problem into an upper triangular system with minimal work ($O(m)$ per step)!
