# Chapter 02: Mathematical Foundations of Linear Algebra & Vector Spaces

---

## 1. Vector Spaces and Subspaces

### 1.1 Formal Definition of a Vector Space
A **vector space** $(V, +, \cdot)$ over a field $\mathbb{F}$ (typically $\mathbb{R}$) is a non-empty set of objects (vectors) equipped with two closed operations:
1. **Vector Addition** ($+$): $u + v \in V$ for all $u, v \in V$.
2. **Scalar Multiplication** ($\cdot$): $c \cdot u \in V$ for all $c \in \mathbb{F}, u \in V$.

These operations must strictly satisfy the **eight fundamental axioms**:
- **Commutativity of Addition**: $u + v = v + u$.
- **Associativity of Addition**: $(u + v) + w = u + (v + w)$.
- **Additive Identity**: There exists a unique zero vector $\mathbf{0} \in V$ such that $u + \mathbf{0} = u$ for all $u \in V$.
- **Additive Inverse**: For every $u \in V$, there exists a vector $-u \in V$ such that $u + (-u) = \mathbf{0}$.
- **Associativity of Scalar Multiplication**: $a(bu) = (ab)u$ for all scalars $a, b \in \mathbb{F}$.
- **Distributivity over Vector Addition**: $a(u + v) = au + av$.
- **Distributivity over Scalar Addition**: $(a + b)u = au + bu$.
- **Scalar Identity**: $1 \cdot u = u$.

### 1.2 Coordinate Spaces vs. Infinite-Dimensional Function Spaces
1. **Finite-Dimensional Coordinate Space ($\mathbb{R}^n$)**:
   - Elements are ordered $n$-tuples of real numbers: $\mathbf{x} = [x_1, x_2, \dots, x_n]^T$.
2. **Continuous Function Spaces ($C[a, b]$)**:
   - Elements are continuous real-valued functions $f(x)$ defined on an interval $[a, b]$.
   - Addition: $(f + g)(x) = f(x) + g(x)$.
   - Scalar multiplication: $(c \cdot f)(x) = c f(x)$.
   - **Crucial Exam Concept**: The space of continuous functions is an **infinite-dimensional vector space** ($\mathbb{R}^\infty$). Each value $f(x)$ at an uncountably infinite number of points $x \in [a, b]$ behaves like an independent coordinate!

### 1.3 Subspaces of a Vector Space
A subset $S \subseteq V$ is a **subspace** of $V$ if and only if $S$ is itself a vector space under the inherited operations of $V$.

#### The Subspace Test (Necessary & Sufficient Conditions)
To verify if a subset $S$ is a valid subspace, verify three criteria:
1. **Non-emptiness / Zero Vector Inclusion**: The zero vector $\mathbf{0} \in S$.
2. **Closure under Addition**: If $u, v \in S$, then $u + v \in S$.
3. **Closure under Scalar Multiplication**: If $u \in S$ and $c \in \mathbb{F}$, then $c u \in S$.

Combining criteria 2 and 3: For all $u, v \in S$ and $a, b \in \mathbb{F}$, the linear combination $a u + b v \in S$.

#### Common Exam Traps & Counterexamples
- **First Quadrant in $\mathbb{R}^2$ ($\{(x, y) \mid x \ge 0, y \ge 0\}$)**:
  - *Is it a subspace?* **NO!** Fails closure under scalar multiplication. If $[1, 2]^T \in S$, multiplying by $c = -1$ yields $[-1, -2]^T \notin S$.
- **Parabola in $\mathbb{R}^2$ ($\{(x, y) \mid y = x^2\}$)**:
  - *Is it a subspace?* **NO!** Contains $[1, 1]^T$ and $[2, 4]^T$. Their sum is $[3, 5]^T$. But $3^2 = 9 \ne 5$. Fails closure under addition.
- **Straight Line passing through origin in $\mathbb{R}^2$ ($y = mx$)**:
  - *Is it a subspace?* **YES!** A 1D subspace of $\mathbb{R}^2$.
- **Straight Line not passing through origin ($y = mx + c, c \ne 0$)**:
  - *Is it a subspace?* **NO!** Fails zero-vector test ($\mathbf{0} \notin S$).
- **Symmetric $3 \times 3$ Matrices in $\mathbb{R}^{3 \times 3}$ ($\{A \mid A = A^T\}$)**:
  - *Is it a subspace?* **YES!** If $A = A^T$ and $B = B^T$, then $(A + B)^T = A^T + B^T = A + B$, and $(cA)^T = cA^T = cA$. Forms a 6-dimensional subspace of $\mathbb{R}^9$.

---

## 2. Linear Independence, Span, Basis, and Dimension

### 2.1 Linear Combination and Span
Given a set of vectors $S = \{v_1, v_2, \dots, v_k\} \subset V$, a **linear combination** is any vector of the form:
$$v = c_1 v_1 + c_2 v_2 + \dots + c_k v_k, \quad c_i \in \mathbb{F}$$
The **span** of $S$, denoted $\text{span}(S)$, is the set of all possible linear combinations of vectors in $S$. The span of any subset is guaranteed to be a subspace of $V$.

### 2.2 Linear Independence
A set of vectors $\{v_1, v_2, \dots, v_k\}$ is **linearly independent** if the only solution to the vector equation:
$$c_1 v_1 + c_2 v_2 + \dots + c_k v_k = \mathbf{0}$$
is the trivial solution:
$$c_1 = c_2 = \dots = c_k = 0$$
If there exist non-zero scalars $c_i$ satisfying the equation, the set is **linearly dependent**. In that case, at least one vector can be expressed as a linear combination of the others.

### 2.3 Fundamental Dependency Theorem
> **Theorem**: If a set $\{v_1, v_2, \dots, v_n\}$ spans a vector space $V$, and $\{w_1, w_2, \dots, w_m\}$ is any set of vectors in $V$ with $m > n$, then $\{w_1, \dots, w_m\}$ is **linearly dependent**.

#### Mathematical Proof
Each vector $w_j$ can be written as a linear combination of the spanning set $\{v_i\}$:
$$w_j = \sum_{i=1}^n a_{ij} v_i \quad \text{for } j = 1, 2, \dots, m$$
In matrix form:
$$W = V A$$
where $V = [v_1, \dots, v_n]$ is $k \times n$, $W = [w_1, \dots, w_m]$ is $k \times m$, and $A = [a_{ij}]$ is an $n \times m$ coefficient matrix.
Consider the linear combination:
$$\sum_{j=1}^m c_j w_j = W \mathbf{c} = V A \mathbf{c} = \mathbf{0}$$
The matrix $A$ has dimensions $n \times m$ with $m > n$ (more columns than rows). By the rank theorem, $A$ has at most rank $n$.
Therefore, the homogeneous system $A \mathbf{c} = \mathbf{0}$ has at least $m - n > 0$ free variables, guaranteeing the existence of a non-trivial solution $\mathbf{c} \ne \mathbf{0}$.
Thus:
$$W \mathbf{c} = V (A \mathbf{c}) = V (\mathbf{0}) = \mathbf{0}$$
with $\mathbf{c} \ne \mathbf{0}$, proving that the set $\{w_1, \dots, w_m\}$ is linearly dependent. $\blacksquare$

### 2.4 Basis and Dimension
- **Basis**: A sequence of vectors $\mathcal{B} = \{v_1, v_2, \dots, v_n\}$ is a basis for $V$ if:
  1. The vectors are linearly independent.
  2. The vectors span $V$ ($\text{span}(\mathcal{B}) = V$).
- **Uniqueness of Representation**: Every vector $v \in V$ can be expressed in a **unique** way as a linear combination of basis vectors:
  $$v = \sum_{i=1}^n c_i v_i$$
- **Dimension ($\dim V$)**: The number of vectors in any basis of $V$. All bases of a given finite-dimensional space $V$ possess the exact same number of vectors.

### 2.5 Row vs. Column Geometric Representations of Linear Systems (Slide 187)

Consider the general $3 \times 3$ linear system:
$$\begin{aligned}
a_1 x + a_2 y + a_3 z &= d_1 \\
b_1 x + b_2 y + b_3 z &= d_2 \\
c_1 x + c_2 y + c_3 z &= d_3
\end{aligned} \iff \begin{bmatrix} a_1 & a_2 & a_3 \\ b_1 & b_2 & b_3 \\ c_1 & c_2 & c_3 \end{bmatrix} \begin{bmatrix} x \\ y \\ z \end{bmatrix} = \begin{bmatrix} d_1 \\ d_2 \\ d_3 \end{bmatrix}$$

```
                Row Picture                                    Column Picture
    (Intersection of 3 Planes in R^3)               (Linear Combination of Vectors in R^3)

               Plane 1 (a·x = d₁)                                  z c
                    /                                             .´
            ───────/───────                                      /
           /  Plane 2 (b·x = d₂)                                /     .´ y b
          /     /                                              /   .´
         ──────/────────                                      / .´
        /  Plane 3 (c·x = d₃)                                ┌─────────► d = x a + y b + z c
       /                                                    / .´
      /                                                    / .´
     Point of Intersection: x* = (x, y, z)                └─────────────► x a
```

#### 1. The Row Picture (Intersection of Hyperplanes)
Each equation represents a 2D plane in 3D Euclidean space ($\mathbb{R}^3$):
- **Unique Solution**: The three planes intersect at exactly one single point $(x, y, z)$.
- **No Solution (Inconsistent)**:
  - Two or more planes are parallel and non-coincident.
  - Two planes intersect along a line, but the third plane is parallel to that line (the line of intersection never penetrates the third plane).
- **Infinite Solutions**:
  - All three planes are coplanar (completely coincide).
  - Two planes are identical/coplanar, and the third plane cuts them along an entire line of intersection.

#### 2. The Column Picture (Linear Combination of Vectors)
The system is rewritten as a vector combination:
$$x \begin{bmatrix} a_1 \\ a_2 \\ a_3 \end{bmatrix} + y \begin{bmatrix} b_1 \\ b_2 \\ b_3 \end{bmatrix} + z \begin{bmatrix} c_1 \\ c_2 \\ c_3 \end{bmatrix} = \begin{bmatrix} d_1 \\ d_2 \\ d_3 \end{bmatrix} \iff x \mathbf{a} + y \mathbf{b} + z \mathbf{c} = \mathbf{d}$$
- **Unique Solution**: The column vectors $\mathbf{a}, \mathbf{b}, \mathbf{c}$ span all of $\mathbb{R}^3$ (they are linearly independent; $\det(A) \ne 0$). Any right-hand side vector $\mathbf{d}$ has unique coordinates $(x, y, z)$.
- **Linearly Dependent Columns**: If any column vector is a linear combination of the other two, all three columns lie in a single 2D plane (or line):
  - **Infinite Solutions**: If $\mathbf{d}$ also lies in that same plane ($\mathbf{d} \in \operatorname{span}\{\mathbf{a}, \mathbf{b}\}$).
  - **No Solution**: If $\mathbf{d}$ points out of that plane ($\mathbf{d} \notin \operatorname{span}\{\mathbf{a}, \mathbf{b}\}$).

---

## 3. The Four Fundamental Subspaces of a Matrix

For any real matrix $A \in \mathbb{R}^{m \times n}$ with rank $r$:

```
                        The Big Picture of Linear Algebra
        Row Space C(A^T)                          Column Space C(A)
       ┌──────────────────┐                      ┌──────────────────┐
       │   dim = r        │     x ──► A x        │   dim = r        │
       │   in R^n         │ ───────────────────► │   in R^m         │
       └──────────────────┘                      └──────────────────┘
                 ▲                                         ▲
                 │ (orthogonal)                            │ (orthogonal)
                 ▼                                         ▼
       ┌──────────────────┐                      ┌──────────────────┐
       │   dim = n - r    │     x ──► 0          │   dim = m - r    │
       │   in R^n         │ ───────────────────► │   in R^m         │
       │   Null Space     │                      │  Left Null Space │
       │      N(A)        │                      │     N(A^T)       │
       └──────────────────┘                      └──────────────────┘
```

### 3.1 Definitions and Properties
1. **Column Space $C(A)$**:
   - The subspace of $\mathbb{R}^m$ spanned by the columns of $A$.
   - Represents all possible vectors $b$ for which $Ax = b$ is solvable.
   - $\dim C(A) = r$.
2. **Row Space $C(A^T)$**:
   - The subspace of $\mathbb{R}^n$ spanned by the rows of $A$ (columns of $A^T$).
   - $\dim C(A^T) = r$.
3. **Null Space $N(A)$**:
   - The set of all vectors $x \in \mathbb{R}^n$ such that $Ax = \mathbf{0}$.
   - $\dim N(A) = n - r$ (number of free variables).
4. **Left Null Space $N(A^T)$**:
   - The set of all vectors $y \in \mathbb{R}^m$ such that $A^T y = \mathbf{0}$ (or $y^T A = \mathbf{0}^T$).
   - $\dim N(A^T) = m - r$.

---

## 4. Fundamental Theorem of Linear Algebra

### 4.1 Part 1: Equality of Row Rank and Column Rank
> **Theorem**: For any matrix $A \in \mathbb{R}^{m \times n}$:
> $$\text{rank}(A) = \dim C(A) = \dim C(A^T) = r$$
> Furthermore:
> $$\dim C(A^T) + \dim N(A) = n \quad (\text{Rank-Nullity Theorem})$$
> $$\dim C(A) + \dim N(A^T) = m$$

#### Proof using Reduced Row Echelon Form (RREF)
Row operations do not alter the linear combinations among rows; hence, the row space is preserved: $C(A^T) = C(R^T)$, where $R = \text{rref}(A)$.
In $R$, there are exactly $r$ non-zero pivot rows, which are clearly linearly independent. Thus, $\dim C(A^T) = r$.
Simultaneously, the pivot columns of $A$ (corresponding to pivot columns of $R$) form a basis for $C(A)$. Thus, $\dim C(A) = r$.
The remaining $n - r$ non-pivot columns correspond to free variables in the solution of $Ax = \mathbf{0}$. Each free variable determines a unique special solution in $N(A)$. Thus, $\dim N(A) = n - r$. $\blacksquare$

### 4.2 Part 2: Orthogonal Complements
> **Theorem**: The fundamental subspaces are mutually orthogonal complements in their respective ambient spaces:
> 1. In $\mathbb{R}^n$: $N(A) = (C(A^T))^\perp$. The null space is orthogonal to the row space.
> 2. In $\mathbb{R}^m$: $N(A^T) = (C(A))^\perp$. The left null space is orthogonal to the column space.

#### Mathematical Proof of Orthogonality
1. **Proof that $N(A) \perp C(A^T)$**:
   Let $x \in N(A)$. By definition:
   $$A x = \mathbf{0}$$
   Writing matrix $A$ in terms of its rows $r_1^T, r_2^T, \dots, r_m^T$:
   $$A x = \begin{bmatrix} r_1^T \\ r_2^T \\ \vdots \\ r_m^T \end{bmatrix} x = \begin{bmatrix} r_1^T x \\ r_2^T x \\ \vdots \\ r_m^T x \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \\ \vdots \\ 0 \end{bmatrix}$$
   This explicitly shows that:
   $$r_i^T x = r_i \cdot x = 0 \quad \text{for every row } i = 1, \dots, m$$
   Any arbitrary vector $v \in C(A^T)$ is a linear combination of rows: $v = \sum_{i=1}^m c_i r_i$.
   Taking the inner product:
   $$v \cdot x = \left(\sum_{i=1}^m c_i r_i\right) \cdot x = \sum_{i=1}^m c_i (r_i \cdot x) = \sum_{i=1}^m c_i (0) = 0$$
   Since $\dim C(A^T) + \dim N(A) = r + (n - r) = n$, the null space is the complete **orthogonal complement** of the row space:
   $$C(A^T) \oplus N(A) = \mathbb{R}^n$$

2. **Proof that $N(A^T) \perp C(A)$**:
   Let $y \in N(A^T)$, so $A^T y = \mathbf{0}$.
   Taking the transpose:
   $$y^T A = \mathbf{0}^T$$
   Writing $A$ in terms of its columns $c_1, c_2, \dots, c_n$:
   $$y^T [c_1, c_2, \dots, c_n] = [y^T c_1, y^T c_2, \dots, y^T c_n] = [0, 0, \dots, 0]$$
   Thus $y \cdot c_j = 0$ for every column $j$.
   Any vector $u \in C(A)$ is $u = \sum d_j c_j$. Hence $y \cdot u = 0$.
   Since $\dim C(A) + \dim N(A^T) = r + (m - r) = m$:
   $$C(A) \oplus N(A^T) = \mathbb{R}^m \quad \blacksquare$$

---

## 5. Comprehensive Worked Example (From Lecture Slides)

Consider the $4 \times 3$ matrix $A$ ($m = 4, n = 3$):
$$A = \begin{bmatrix} 1 & 1 & 2 \\ 2 & 0 & 2 \\ 0 & 1 & 1 \\ 1 & 0 & 1 \end{bmatrix}$$

### Step 1: Determine Rank and Column Space $C(A)$
- Column 1: $c_1 = [1, 2, 0, 1]^T$
- Column 2: $c_2 = [1, 0, 1, 0]^T$
- Column 3: $c_3 = [2, 2, 1, 1]^T$
Notice that $c_3 = c_1 + c_2$:
$$\begin{bmatrix} 1 \\ 2 \\ 0 \\ 1 \end{bmatrix} + \begin{bmatrix} 1 \\ 0 \\ 1 \\ 0 \end{bmatrix} = \begin{bmatrix} 2 \\ 2 \\ 1 \\ 1 \end{bmatrix}$$
Since $c_1$ and $c_2$ are clearly linearly independent, the rank $r = 2$.
$$\dim C(A) = 2, \quad \text{Basis for } C(A) = \left\{ \begin{bmatrix} 1 \\ 2 \\ 0 \\ 1 \end{bmatrix}, \begin{bmatrix} 1 \\ 0 \\ 1 \\ 0 \end{bmatrix} \right\}$$
$C(A)$ is a 2D plane embedded in $\mathbb{R}^4$.

### Step 2: Determine Row Space $C(A^T)$
The rows are:
- $r_1 = [1, 1, 2]^T$
- $r_2 = [2, 0, 2]^T$
- $r_3 = [0, 1, 1]^T$
- $r_4 = [1, 0, 1]^T$
Notice that $r_4 = \frac{1}{2} r_2$, and $r_3 = r_1 - r_4$.
There are only 2 independent rows:
$$\dim C(A^T) = 2, \quad \text{Basis for } C(A^T) = \left\{ \begin{bmatrix} 1 \\ 1 \\ 2 \end{bmatrix}, \begin{bmatrix} 2 \\ 0 \\ 2 \end{bmatrix} \right\} \subset \mathbb{R}^3$$

### Step 3: Solve for Null Space $N(A)$
Solve $Ax = \mathbf{0}$ for $x = [u, v, w]^T \in \mathbb{R}^3$:
$$\begin{bmatrix} 1 & 1 & 2 \\ 2 & 0 & 2 \\ 0 & 1 & 1 \\ 1 & 0 & 1 \end{bmatrix} \begin{bmatrix} u \\ v \\ w \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \\ 0 \\ 0 \end{bmatrix}$$
From row 4: $u + w = 0 \implies u = -w$.
From row 3: $v + w = 0 \implies v = -w$.
Setting free variable $w = -1 \implies u = 1, v = 1$.
$$N(A) = \kappa \begin{bmatrix} 1 \\ 1 \\ -1 \end{bmatrix}, \quad \kappa \in \mathbb{R}$$
$$\dim N(A) = n - r = 3 - 2 = 1 \quad (\text{a line passing through the origin in } \mathbb{R}^3)$$

#### Orthogonality Verification ($N(A) \perp C(A^T)$)
Take the dot product of the basis of $N(A)$ with each basis vector of $C(A^T)$:
$$\begin{bmatrix} 1 \\ 1 \\ -1 \end{bmatrix} \cdot \begin{bmatrix} 1 \\ 1 \\ 2 \end{bmatrix} = (1)(1) + (1)(1) + (-1)(2) = 1 + 1 - 2 = 0 \quad \checkmark$$
$$\begin{bmatrix} 1 \\ 1 \\ -1 \end{bmatrix} \cdot \begin{bmatrix} 2 \\ 0 \\ 2 \end{bmatrix} = (1)(2) + (1)(0) + (-1)(2) = 2 + 0 - 2 = 0 \quad \checkmark$$
Orthogonality holds unconditionally!

### Step 4: Solve for Left Null Space $N(A^T)$
Solve $A^T y = \mathbf{0}$ for $y = [p, q, r, s]^T \in \mathbb{R}^4$:
$$A^T y = \begin{bmatrix} 1 & 2 & 0 & 1 \\ 1 & 0 & 1 & 0 \\ 2 & 2 & 1 & 1 \end{bmatrix} \begin{bmatrix} p \\ q \\ r \\ s \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \\ 0 \end{bmatrix}$$
1. Equation (2): $p + r = 0 \implies r = -p$.
2. Equation (1): $p + 2q + s = 0 \implies s = -p - 2q$.
3. Equation (3): $2p + 2q + r + s = 2p + 2q + (-p) + (-p - 2q) = 0 = 0$ (redundant).
Free variables: $p$ and $q$.
- Case 1: Set $p = 0, q = 1 \implies r = 0, s = -2 \implies y_1 = [0, 1, 0, -2]^T$.
- Case 2: Set $p = -1, q = 0 \implies r = 1, s = 1 \implies y_2 = [-1, 0, 1, 1]^T$.
$$N(A^T) = \omega \begin{bmatrix} 0 \\ 1 \\ 0 \\ -2 \end{bmatrix} + \zeta \begin{bmatrix} -1 \\ 0 \\ 1 \\ 1 \end{bmatrix}$$
$$\dim N(A^T) = m - r = 4 - 2 = 2 \quad (\text{a 2D plane in } \mathbb{R}^4)$$

#### Orthogonality Verification ($N(A^T) \perp C(A)$)
Check dot products with $C(A)$ basis:
$$\begin{bmatrix} 0 \\ 1 \\ 0 \\ -2 \end{bmatrix} \cdot \begin{bmatrix} 1 \\ 2 \\ 0 \\ 1 \end{bmatrix} = 0(1) + 1(2) + 0(0) + (-2)(1) = 2 - 2 = 0 \quad \checkmark$$
$$\begin{bmatrix} -1 \\ 0 \\ 1 \\ 1 \end{bmatrix} \cdot \begin{bmatrix} 1 \\ 2 \\ 0 \\ 1 \end{bmatrix} = (-1)(1) + 0(2) + 1(0) + 1(1) = -1 + 1 = 0 \quad \checkmark$$
Both basis vectors of $N(A^T)$ are strictly orthogonal to $C(A)$!

---

## 6. Solvability of Linear Systems: $Ax = b$

### 6.1 The Consistency Condition (Fredholm Alternative)
For a system $Ax = b$ ($A \in \mathbb{R}^{m \times n}$):
- A solution $x$ exists if and only if $b$ lies in the column space:
  $$b \in C(A)$$
- Since $C(A) = (N(A^T))^\perp$, this is mathematically identical to stating:
  $$b \perp N(A^T)$$
- **Crucial Exam Rule**: If $b$ has any non-zero projection/component along the left null space $N(A^T)$, then $Ax = b$ has **NO SOLUTION** (the system is inconsistent)!

```
            Geometric Solvability in R^m (Slide 552)
                     N(A^T)  (Left Null Space)
                       ▲
                       │     .´ b (unsolvable if b has a component along N(A^T))
                       │   .´
                       │ .´
                       └───────────────► C(A) (Column Space)
                         (b must lie entirely in C(A) for Ax = b to have a solution)
```
In our $4 \times 3$ worked example, $C(A)$ is a 2D plane and $N(A^T)$ is an orthogonal 2D plane in $\mathbb{R}^4$. If $b$ has any projection along either $y_1 = [0, 1, 0, -2]^T$ or $y_2 = [-1, 0, 1, 1]^T$, $Ax = b$ is immediately inconsistent.

### 6.2 The Four Cases of Solvability (Rank vs Dimensions)
| Case | Dimensions & Rank | Inverses | Number of Solutions to $Ax = b$ |
| :--- | :--- | :--- | :--- |
| **$r = m = n$** | Square, full rank | Both left and right ($A^{-1}$ exists) | Exactly **1 unique solution** for any $b$. |
| **$r = m < n$** | Short & wide (more columns) | **Right inverse** exists ($A C = I_m$) | **Infinitely many solutions** for any $b$. |
| **$r = n < m$** | Tall & thin (more rows) | **Left inverse** exists ($B A = I_n$) | **0 or 1 unique solution** (1 if $b \in C(A)$, 0 if $b \notin C(A)$). |
| **$r < m$ and $r < n$** | Rank-deficient | Neither inverse exists | **0 or infinitely many solutions**. |

---

## 7. Left and Right Inverses

### 7.1 Right Inverse (Existence of Solutions)
- A matrix $C \in \mathbb{R}^{n \times m}$ is a **right inverse** of $A \in \mathbb{R}^{m \times n}$ if:
  $$A C = I_{m \times m}$$
- **Condition for Existence**: Matrix $A$ must have **full row rank**:
  $$r = m \le n$$
- **Formula**:
  $$C = A^T (A A^T)^{-1}$$
- *Derivation*:
  Since $A$ has independent rows, the $m \times m$ matrix $A A^T$ is symmetric, positive definite, and invertible!
  Multiplying:
  $$A C = A \left[ A^T (A A^T)^{-1} \right] = (A A^T) (A A^T)^{-1} = I_m \quad \checkmark$$
- **Implication**: If a right inverse exists, a solution to $Ax = b$ **always exists** for any right-hand side vector $b$:
  $$x = C b \implies A x = A (C b) = (A C) b = I_m b = b$$

### 7.2 Left Inverse (Uniqueness of Solutions)
- A matrix $B \in \mathbb{R}^{n \times m}$ is a **left inverse** of $A \in \mathbb{R}^{m \times n}$ if:
  $$B A = I_{n \times n}$$
- **Condition for Existence**: Matrix $A$ must have **full column rank**:
  $$r = n \le m$$
- **Formula**:
  $$B = (A^T A)^{-1} A^T$$
- *Derivation*:
  Since $A$ has independent columns, the $n \times n$ matrix $A^T A$ is symmetric, positive definite, and invertible!
  Multiplying:
  $$B A = \left[ (A^T A)^{-1} A^T \right] A = (A^T A)^{-1} (A^T A) = I_n \quad \checkmark$$
- **Implication**: If a solution to $Ax = b$ exists, it is **strictly unique**:
  $$A x = b \implies B (A x) = B b \implies (B A) x = B b \implies x = B b$$
  (Note: $B = (A^T A)^{-1} A^T$ is the classic Moore-Penrose pseudoinverse used in Ordinary Least Squares!).

### 7.3 Two-Sided Inverse of a Square Matrix
If $A$ is square ($m = n$) and full rank ($r = n$):
- Both left inverse $B$ and right inverse $C$ exist.
- Proof that $B = C$:
  $$B = B I = B (A C) = (B A) C = I C = C = A^{-1}$$

---

## 8. Determinants and Their Geometric Properties

### 8.1 Geometric Meaning of the Determinant
- **In $\mathbb{R}^2$**: The absolute value $|\det([u, v])|$ equals the **area of the parallelogram** formed by vectors $u$ and $v$.
- **In $\mathbb{R}^3$**: The absolute value $|\det([u, v, w])|$ equals the **volume of the parallelepiped** formed by column vectors $a = [a_1, a_2, a_3]^T$, $b = [b_1, b_2, b_3]^T$, $c = [c_1, c_2, c_3]^T$ (Slide 169):
  $$A = \begin{bmatrix} a_1 & b_1 & c_1 \\ a_2 & b_2 & c_2 \\ a_3 & b_3 & c_3 \end{bmatrix}$$
  $$\mathbf{\text{Volume} = \det(A) = a_1(b_2 c_3 - b_3 c_2) + b_1(c_2 a_3 - a_2 c_3) + c_1(a_2 b_3 - a_3 b_2)}$$
  Because $\det(A^T) = \det(A)$, **this volume is identical whether computed from the columns or from the rows of $A$**!
- **Orientation / Sign**:
  - Positive determinant $\implies$ right-handed coordinate system.
  - Negative determinant $\implies$ left-handed coordinate system (reflection).
  - Zero determinant ($\det(A) = 0$) $\implies$ vectors are linearly dependent; the parallelepiped collapses to a flat plane or line (zero volume)!

### 8.2 Ten Key Properties of Determinants
1. $\det(I) = 1$.
2. Swapping two rows reverses the sign: $\det(P A) = -\det(A)$.
3. Linearity in each individual row:
   $$\det \begin{bmatrix} c a_1 & c a_2 \\ b_1 & b_2 \end{bmatrix} = c \det \begin{bmatrix} a_1 & a_2 \\ b_1 & b_2 \end{bmatrix}$$
   $$\det(c A) = c^n \det(A) \quad (\text{for an } n \times n \text{ matrix!})$$
4. If two rows of $A$ are equal, $\det(A) = 0$.
5. Subtracting a multiple of one row from another leaves $\det(A)$ unchanged (fundamental basis of Gaussian elimination!).
6. A row of zeros implies $\det(A) = 0$.
7. For a triangular matrix, the determinant is the product of diagonal pivot elements:
   $$\det(A) = d_1 d_2 \dots d_n$$
8. $\det(A) \ne 0$ if and only if $A$ is non-singular (invertible).
9. $\det(A B) = \det(A) \det(B)$. Consequently:
   $$\det(A^{-1}) = \frac{1}{\det(A)}$$
10. Transposition preserves the determinant: $\det(A^T) = \det(A)$.

### 8.3 Cramer's Rule and Its Computational Limitations
For an $n \times n$ system $Ax = b$ with $\det(A) \ne 0$:
$$x_i = \frac{\det(B_i)}{\det(A)}$$
where $B_i$ is the matrix formed by replacing the $i$-th column of $A$ with vector $b$.
- **Computational Cost**:
  - Evaluating an $n \times n$ determinant via Laplace expansion requires $O(n!)$ operations.
  - Evaluating via Gaussian elimination requires $O(n^3)$ operations.
  - Solving $n$ unknowns requires computing $n + 1$ determinants:
    $$\text{Work} = (n + 1) \times \frac{2}{3} n^3 \approx O(n^4)\text{ operations!}$$
- **Exam Takeaway**: Cramer's rule is an invaluable theoretical formula for algebraic proofs and small ($2 \times 2, 3 \times 3$) systems, but computationally disastrous for scientific computing ($n = 10^6$).
