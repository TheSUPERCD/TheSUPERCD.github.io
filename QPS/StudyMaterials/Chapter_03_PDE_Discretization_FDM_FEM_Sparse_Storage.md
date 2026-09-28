# Chapter 03: Discretization of PDEs & Sparse Matrix Assembly

---

## 1. Physical Governing Equations and Problem Formulation

In scientific computing, physical conservation laws (mass, momentum, energy) are expressed as Partial Differential Equations (PDEs):

### 1.1 One-Dimensional Steady Heat Conduction
For a 1D heat conducting rod of thermal conductivity $k$:
$$\frac{d}{dx} \left( k \frac{dT}{dx} \right) = 0 \quad \xrightarrow{\text{constant } k} \quad \frac{d^2 T}{dx^2} = 0, \quad 0 \le x \le L$$
Subject to boundary conditions:
$$T(0) = T_A, \quad T(L) = T_B$$

### 1.2 Two-Dimensional Laplace and Poisson Equations
In a 2D planar domain $\Omega \subset \mathbb{R}^2$:
$$\nabla^2 T = \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = f(x, y)$$
- If $f(x, y) = 0$: **Laplace Equation** (steady potential flow, electrostatics, source-free heat conduction).
- If $f(x, y) \ne 0$: **Poisson Equation** (heat conduction with internal volumetric heat generation, pressure Poisson equation in incompressible CFD).

### 1.3 Classification of Linear Boundary Conditions
Any linear boundary condition on $\Gamma = \partial \Omega$ belongs to one of three types:
1. **Dirichlet Boundary Condition (Type I / Essential)**:
   - Prescribes the value of the variable itself:
     $$T = T_0 \quad \text{on } \Gamma_1$$
2. **Neumann Boundary Condition (Type II / Natural)**:
   - Prescribes the normal gradient (heat flux):
     $$-k \frac{\partial T}{\partial n} = -k (\nabla T \cdot \mathbf{n}) = q_0 \quad \text{on } \Gamma_2$$
   - Insulated / adiabatic wall: $\frac{\partial T}{\partial n} = 0$.
3. **Robin Boundary Condition (Type III / Mixed)**:
   - Prescribes a linear combination of the variable and its normal flux (e.g., convective heat transfer into surrounding fluid):
     $$-k \frac{\partial T}{\partial n} = h (T - T_\infty) \quad \text{on } \Gamma_3$$

---

## 2. Analytical Solution of the 2D Laplace Equation (Separation of Variables)

Understanding the exact analytical solution provides the baseline for verifying numerical solvers.

### 2.1 Boundary Value Problem Formulation
Consider a rectangular domain $\Omega = [0, L] \times [0, H]$ governed by:
$$\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0$$
Subject to Dirichlet boundary conditions:
- Left wall: $T(0, y) = 0$
- Right wall: $T(L, y) = 0$
- Bottom wall: $T(x, 0) = 0$
- Top wall: $T(x, H) = f(x)$ (e.g., $f(x) = T_0 = \text{constant}$)

```
                      T(x, H) = f(x)
               (0, H) ┌──────────────┐ (L, H)
                      │              │
       T(0, y) = 0    │   ∇²T = 0    │   T(L, y) = 0
                      │              │
               (0, 0) └──────────────┘ (L, 0)
                      T(x, 0) = 0
```

### 2.2 Complete Step-by-Step Derivation
1. **Separation of Variables Hypothesis**:
   Assume the solution decomposes into a product of single-variable functions:
   $$T(x, y) = F(x) G(y)$$
   Differentiating:
   $$\frac{\partial^2 T}{\partial x^2} = F''(x) G(y), \quad \frac{\partial^2 T}{\partial y^2} = F(x) G''(y)$$
   Substituting into Laplace's equation:
   $$F''(x) G(y) + F(x) G''(y) = 0$$
   Dividing across by $F(x) G(y)$:
   $$\frac{F''(x)}{F(x)} = -\frac{G''(y)}{G(y)} = k \quad (\text{Separation Constant})$$

2. **Choice of Separation Constant $k$**:
   Since $F(x)$ must satisfy homogeneous boundary conditions at both ends ($x = 0$ and $x = L$), it must possess an oscillatory (trigonometric) solution.
   - If $k > 0$: $F(x) = c_1 e^{\sqrt{k}x} + c_2 e^{-\sqrt{k}x}$, which cannot satisfy $F(0) = F(L) = 0$ without being identically zero.
   - If $k = 0$: $F(x) = c_1 x + c_2 \implies c_1 = c_2 = 0$ (trivial).
   - Therefore, $k$ must be strictly negative:
     $$k = -\lambda^2 \quad (\lambda > 0)$$

3. **Solving for $F(x)$ (The $x$-Direction Eigenvalue Problem)**:
   $$F''(x) + \lambda^2 F(x) = 0$$
   General solution:
   $$F(x) = A \cos(\lambda x) + B \sin(\lambda x)$$
   Applying boundary conditions:
   - At $x = 0$: $F(0) = A (1) + B (0) = 0 \implies A = 0$.
   - At $x = L$: $F(L) = B \sin(\lambda L) = 0$.
   For a non-trivial solution ($B \ne 0$), we obtain the **quantized eigenvalues**:
   $$\sin(\lambda L) = 0 \implies \lambda_n = \frac{n \pi}{L}, \quad n = 1, 2, 3, \dots$$
   The corresponding eigenfunctions are:
   $$F_n(x) = \sin\left(\frac{n \pi x}{L}\right)$$

4. **Solving for $G(y)$ (The $y$-Direction ODE)**:
   $$-\frac{G''(y)}{G(y)} = -\lambda_n^2 \implies G''(y) - \lambda_n^2 G(y) = 0$$
   General solution expressed in hyperbolic functions:
   $$G(y) = C \cosh(\lambda_n y) + D \sinh(\lambda_n y)$$
   Applying boundary condition at $y = 0$:
   $$G(0) = C \cosh(0) + D \sinh(0) = C (1) + 0 = 0 \implies C = 0$$
   Therefore:
   $$G_n(y) = D_n \sinh\left(\frac{n \pi y}{L}\right)$$

5. **Principle of Linear Superposition**:
   By linearity, the general solution is the infinite series:
   $$T(x, y) = \sum_{n=1}^\infty B_n \sin\left(\frac{n \pi x}{L}\right) \sinh\left(\frac{n \pi y}{L}\right)$$
   where $B_n$ represents the combined constant.

6. **Determining Fourier Coefficients via Top Boundary Condition**:
   At $y = H$:
   $$T(x, H) = \sum_{n=1}^\infty \left[ B_n \sinh\left(\frac{n \pi H}{L}\right) \right] \sin\left(\frac{n \pi x}{L}\right) = f(x)$$
   This is a classical **Fourier Sine Series**. Exploiting the orthogonality of sine functions over $[0, L]$:
   $$\int_0^L \sin\left(\frac{n \pi x}{L}\right) \sin\left(\frac{m \pi x}{L}\right) dx = \begin{cases} 0 & n \ne m \\ \frac{L}{2} & n = m \end{cases}$$
   Multiplying both sides by $\sin\left(\frac{m \pi x}{L}\right)$ and integrating:
   $$B_n \sinh\left(\frac{n \pi H}{L}\right) \frac{L}{2} = \int_0^L f(x) \sin\left(\frac{n \pi x}{L}\right) dx$$
   $$\mathbf{B_n = \frac{2}{L \sinh\left(\frac{n \pi H}{L}\right)} \int_0^L f(x) \sin\left(\frac{n \pi x}{L}\right) dx}$$

7. **Special Case: Uniform Temperature at Top Wall ($f(x) = T_0 = \text{const}$)**:
   $$\int_0^L T_0 \sin\left(\frac{n \pi x}{L}\right) dx = T_0 \left[ -\frac{L}{n \pi} \cos\left(\frac{n \pi x}{L}\right) \right]_0^L = \frac{L T_0}{n \pi} [1 - (-1)^n]$$
   - For even $n$ ($n = 2, 4, 6, \dots$): $1 - (-1)^n = 0 \implies B_n = 0$.
   - For odd $n$ ($n = 1, 3, 5, \dots$): $1 - (-1)^n = 2 \implies \int = \frac{2 L T_0}{n \pi}$.
   Thus:
   $$B_n = \frac{4 T_0}{n \pi \sinh\left(\frac{n \pi H}{L}\right)} \quad (\text{for odd } n)$$
   $$\mathbf{T(x, y) = \frac{4 T_0}{\pi} \sum_{n=1, 3, 5, \dots}^\infty \frac{1}{n} \frac{\sinh\left(\frac{n \pi y}{L}\right)}{\sinh\left(\frac{n \pi H}{L}\right)} \sin\left(\frac{n \pi x}{L}\right)}$$

---

## 3. Finite Difference Method (FDM) Discretization

### 3.1 Taylor Series Derivation of Central Differences
Let a 1D grid have uniform spacing $\Delta x$. Expanding $T(x + \Delta x)$ and $T(x - \Delta x)$ around $x$:
$$T(x + \Delta x) = T(x) + \Delta x T'(x) + \frac{\Delta x^2}{2} T''(x) + \frac{\Delta x^3}{6} T'''(x) + \frac{\Delta x^4}{24} T^{(4)}(x) + O(\Delta x^5) \quad \text{--- (A)}$$
$$T(x - \Delta x) = T(x) - \Delta x T'(x) + \frac{\Delta x^2}{2} T''(x) - \frac{\Delta x^3}{6} T'''(x) + \frac{\Delta x^4}{24} T^{(4)}(x) - O(\Delta x^5) \quad \text{--- (B)}$$

Adding (A) and (B):
$$T(x + \Delta x) + T(x - \Delta x) = 2 T(x) + \Delta x^2 T''(x) + \frac{\Delta x^4}{12} T^{(4)}(x) + O(\Delta x^6)$$
Rearranging for $T''(x)$:
$$\left.\frac{d^2 T}{dx^2}\right|_i = \frac{T_{i+1} - 2 T_i + T_{i-1}}{\Delta x^2} - \frac{\Delta x^2}{12} T^{(4)}(\xi) = \frac{T_{i+1} - 2 T_i + T_{i-1}}{\Delta x^2} + O(\Delta x^2)$$
The truncation error is strictly of order $O(\Delta x^2)$ (second-order accurate).

### 3.2 2D Laplace Discretization (The 5-Point Stencil)
On a 2D Cartesian grid with spacings $\Delta x$ and $\Delta y$:
$$\left.\nabla^2 T\right|_{i, j} = \frac{T_{i+1, j} - 2 T_{i, j} + T_{i-1, j}}{\Delta x^2} + \frac{T_{i, j+1} - 2 T_{i, j} + T_{i, j-1}}{\Delta y^2} = 0$$
Assuming an isotropic grid ($\Delta x = \Delta y$):
$$T_{i+1, j} + T_{i-1, j} + T_{i, j+1} + T_{i, j-1} - 4 T_{i, j} = 0$$

```
                           (i, j+1)  [North]
                              │
                              │
       (i-1, j) [West] ───────┼─────── (i+1, j) [East]
                              │  (i, j)
                              │  [Center, weight = -4]
                           (i, j-1)  [South]
```

### 3.3 Lexicographical Ordering and Matrix Formulation
To transform the 2D grid values $T(i, j)$ into a 1D column vector $\mathbf{T}$, we use lexicographical (row-by-row) indexing:
$$\text{ipt}(i, j) = (j - 1) N_x + i$$
where $i \in \{1, 2, \dots, N_x\}$ and $j \in \{1, 2, \dots, N_y\}$.
Under this linear indexing:
- North neighbor $(i, j+1) \implies \text{ipt} + N_x$
- South neighbor $(i, j-1) \implies \text{ipt} - N_x$
- East neighbor $(i+1, j) \implies \text{ipt} + 1$
- West neighbor $(i-1, j) \implies \text{ipt} - 1$

The algebraic equation becomes:
$$4 T_{\text{ipt}} - T_{\text{ipt}-1} - T_{\text{ipt}+1} - T_{\text{ipt}-N_x} - T_{\text{ipt}+N_x} = 0$$
When assembled into $A \mathbf{T} = \mathbf{b}$, matrix $A$ is:
- **Pentadiagonal**: Non-zeros exist only on 5 diagonals: the main diagonal ($+4$), sub- and super-diagonals ($-1$), and far diagonals shifted by $\pm N_x$ ($-1$).
- **Symmetric**: $A = A^T$.
- **Strictly / Irreducibly Diagonally Dominant**: On interior nodes: $|a_{ii}| = 4 = \sum_{j \ne i} |a_{ij}|$. On boundary-adjacent nodes, Dirichlet values move to the right-hand side, making $|a_{ii}| > \sum |a_{ij}|$.
- **Positive Definite**: All eigenvalues $\lambda_k > 0$.

### 3.4 3D Poisson Discretization (The 7-Point Stencil)
In 3D with isotropic grid spacing $\Delta x = \Delta y = \Delta z$:
$$\text{ipt}(i, j, k) = i + (j - 1) N_x + (k - 1) N_x N_y$$
The resulting stencil connects 7 points (Center, East, West, North, South, Top, Bottom):
$$6 T_{\text{ipt}} - T_{\text{ipt}-1} - T_{\text{ipt}+1} - T_{\text{ipt}-N_x} - T_{\text{ipt}+N_x} - T_{\text{ipt}-N_x N_y} - T_{\text{ipt}+N_x N_y} = \Delta x^2 f_{\text{ipt}}$$
Assembled matrix $A$ is **septadiagonal**.

---

## 4. Finite Volume Method (FVM)

In FVM, the domain is subdivided into contiguous control volumes. The PDE is integrated over control volume $V_P$:
$$\int_{V_P} \nabla \cdot (k \nabla T) \, dV = \int_{\partial V_P} (k \nabla T) \cdot \mathbf{n} \, dA = \sum_{f \in \{e, w, n, s\}} Q_f = 0$$
Evaluating conductive heat fluxes across faces:
$$Q_e = k \left( \frac{T_E - T_P}{\Delta x} \right) \Delta y, \quad Q_w = k \left( \frac{T_W - T_P}{\Delta x} \right) \Delta y$$
$$Q_n = k \left( \frac{T_N - T_P}{\Delta y} \right) \Delta x, \quad Q_s = k \left( \frac{T_S - T_P}{\Delta y} \right) \Delta x$$
For uniform Cartesian grids ($\Delta x = \Delta y$), summing fluxes directly yields:
$$T_E + T_W + T_N + T_S - 4 T_P = 0$$
**Key Invariance**: On uniform orthogonal meshes, FVM and FDM produce mathematically identical algebraic stencils. However, FVM guarantees exact local and global conservation on arbitrary meshes.

---

## 5. Finite Element Method (FEM) & Stiffness Matrix Assembly

### 5.1 Weak Formulation (1D Conduction Example)
Consider:
$$-\frac{d^2 \phi}{dx^2} = 0, \quad 0 \le x \le 1, \quad \phi(0) = 0, \phi(1) = 1$$
Multiply by an arbitrary test (weighting) function $w(x)$ with $w(0) = w(1) = 0$ and integrate:
$$-\int_0^1 w \frac{d^2 \phi}{dx^2} \, dx = 0$$
Applying integration by parts:
$$\left[ -w \frac{d\phi}{dx} \right]_0^1 + \int_0^1 \frac{dw}{dx} \frac{d\phi}{dx} \, dx = 0$$
Since $w$ vanishes at boundaries, the **weak form** is:
$$\int_0^1 \frac{dw}{dx} \frac{d\phi}{dx} \, dx = 0$$

### 5.2 Galerkin Approximation with 1D Linear Elements
Discretize $[0, 1]$ into elements of size $\Delta x$. Within element $m$ with local coordinate $\xi = \frac{x - x_m}{\Delta x} \in [0, 1]$:
$$\phi(\xi) = \phi_m N_1(\xi) + \phi_{m+1} N_2(\xi)$$
where the **linear shape functions** are:
$$N_1(\xi) = 1 - \xi, \quad N_2(\xi) = \xi$$
Derivatives:
$$\frac{dN_1}{dx} = \frac{1}{\Delta x} \frac{dN_1}{d\xi} = -\frac{1}{\Delta x}, \quad \frac{dN_2}{dx} = \frac{1}{\Delta x} \frac{dN_2}{d\xi} = \frac{1}{\Delta x}$$

### 5.3 Derivation of Element Stiffness Matrix ($k^e$)
The element integral is:
$$k_{ij}^e = \int_{x_m}^{x_{m+1}} \frac{dN_i}{dx} \frac{dN_j}{dx} \, dx = \int_0^1 \left(\frac{dN_i}{dx}\right) \left(\frac{dN_j}{dx}\right) \Delta x \, d\xi$$
- $k_{11}^e = \int_0^1 \left(-\frac{1}{\Delta x}\right) \left(-\frac{1}{\Delta x}\right) \Delta x \, d\xi = \frac{1}{\Delta x}$
- $k_{12}^e = \int_0^1 \left(-\frac{1}{\Delta x}\right) \left(\frac{1}{\Delta x}\right) \Delta x \, d\xi = -\frac{1}{\Delta x}$
- $k_{21}^e = -\frac{1}{\Delta x}$
- $k_{22}^e = \frac{1}{\Delta x}$

$$\mathbf{k^e = \frac{1}{\Delta x} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}}$$

### 5.4 Global Assembly
Assembling adjacent elements $m-1$ and $m$ at shared node $m$:
$$\left( k_{22}^{e(m-1)} + k_{11}^{e(m)} \right) \phi_m + k_{21}^{e(m-1)} \phi_{m-1} + k_{12}^{e(m)} \phi_{m+1} = 0$$
$$\left( \frac{1}{\Delta x} + \frac{1}{\Delta x} \right) \phi_m - \frac{1}{\Delta x} \phi_{m-1} - \frac{1}{\Delta x} \phi_{m+1} = 0$$
$$\frac{-\phi_{m-1} + 2\phi_m - \phi_{m+1}}{\Delta x} = 0 \iff \mathbf{\frac{\phi_{m-1} - 2\phi_m + \phi_{m+1}}{\Delta x^2} = 0}$$
FEM recovers the exact second-order central difference stencil!

### 5.5 2D Triangular Elements (Unstructured Meshes)
For 2D complex geometries, triangular elements with 3 nodes $(x_1, y_1), (x_2, y_2), (x_3, y_3)$ are used:
- Area of triangle:
  $$A = \frac{1}{2} \det \begin{bmatrix} 1 & x_1 & y_1 \\ 1 & x_2 & y_2 \\ 1 & x_3 & y_3 \end{bmatrix}$$
- Linear shape functions $N_i(x, y) = \frac{1}{2A} (a_i + b_i x + c_i y)$, where:
  $$a_1 = x_2 y_3 - x_3 y_2, \quad b_1 = y_2 - y_3, \quad c_1 = x_3 - x_2$$
- Gradient matrix $B$:
  $$B = \begin{bmatrix} \frac{\partial N_1}{\partial x} & \frac{\partial N_2}{\partial x} & \frac{\partial N_3}{\partial x} \\ \frac{\partial N_1}{\partial y} & \frac{\partial N_2}{\partial y} & \frac{\partial N_3}{\partial y} \end{bmatrix} = \frac{1}{2A} \begin{bmatrix} b_1 & b_2 & b_3 \\ c_1 & c_2 & c_3 \end{bmatrix}$$
- Element stiffness matrix:
  $$A^e = \iint_A B^T B \, dx dy = A (B^T B)$$
- **Crucial Distinction**: The resulting global assembled matrix is **sparse, symmetric, positive definite, but UNSTRUCTURED** (non-zeros do not follow constant diagonal bands). This necessitates dedicated sparse storage formats.

---

## 6. Sparse Matrix Storage Schemes (Saad's Canonical Formats)

In an $N \times N$ matrix derived from a PDE ($N = 10^6$):
- Storing full dense matrix: $N^2 \times 8\text{ bytes} = 10^{12} \times 8\text{ B} = 8\text{ Terabytes}$ of memory!
- Number of non-zero entries ($N_{nz}$): Each row has only 5 non-zeros $\implies N_{nz} \approx 5 \times 10^6$ entries.
- Storing only non-zeros: $5 \times 10^6 \times 8\text{ B} = 40\text{ Megabytes}$ (a $200,000\times$ reduction in RAM!).

To illustrate the sparse storage schemes, we use the **canonical $5 \times 5$ unstructured matrix ($N = 5, N_z = 12$)** presented in the course lectures (*Yousef Saad, Iterative Methods for Sparse Linear Systems*):

$$A = \begin{bmatrix} 
1. & 0. & 0. & 2. & 0. \\
3. & 4. & 0. & 5. & 0. \\
6. & 0. & 7. & 8. & 9. \\
0. & 0. & 10. & 11. & 0. \\
0. & 0. & 0. & 0. & 12. 
\end{bmatrix}$$

### 6.1 Coordinate Format (COO)
Stores triplets $(i, j, A_{ij})$ in arbitrary order:
- `AA`: Real array containing all non-zero entries (size $N_z$).
- `JR`: Integer array containing 1-based row indices (size $N_z$).
- `JC`: Integer array containing 1-based column indices (size $N_z$).
- **Total Storage**: $3 N_z$ words ($3 \times 12 = 36$ elements).

```
Index:   1    2    3    4    5    6    7    8    9   10   11   12
───────────────────────────────────────────────────────────────────
AA:   [12.   9.   7.   5.   1.   2.  11.   3.   6.   4.   8.  10.]
JR:   [ 5    3    3    2    1    1    4    2    3    2    3    4 ]
JC:   [ 5    5    3    4    1    4    4    1    1    2    4    3 ]
```

### 6.2 Compressed Sparse Row (CSR / CRS)
Eliminates redundant row storage by using row pointers:
- `AA`: Real values of non-zeros, stored row-by-row, left-to-right (size $N_z = 12$).
- `JA`: Integer array containing column indices corresponding to `AA` (size $N_z = 12$).
- `IA`: Integer pointers to the beginning of each row in `AA` and `JA` (size $n + 1 = 6$).
  - `IA(i)` gives the starting position of row $i$.
  - Number of non-zeros in row $i$ = `IA(i+1) - IA(i)`.
  - `IA(n+1) = IA(1) + N_z = 13` (points to fictitious row $n+1$).
- **Total Storage**: $2 N_z + n + 1 = 2(12) + 5 + 1 = 30$ words.

```
Index:   1    2    3    4    5    6    7    8    9   10   11   12
───────────────────────────────────────────────────────────────────
AA:   [ 1.   2.   3.   4.   5.   6.   7.   8.   9.  10.  11.  12.]
JA:   [ 1    4    1    2    4    1    3    4    5    3    4    5 ]

IA:   [ 1    3    6   10   12   13]  (length n + 1 = 6)
```

### 6.3 Compressed Sparse Column (CSC / CCS)
Column-oriented counterpart of CSR:
- `AA`: Non-zeros stored column-by-column.
- `JA`: Row indices for each entry in `AA`.
- `IA`: Column pointers (size $n + 1$).
- Native storage standard for MATLAB and Fortran sparse routines.

### 6.4 Modified Sparse Row (MSR)
Exploits the fact that diagonal elements are always non-zero and accessed most frequently. Reduces storage to **exactly two arrays** (`AA` and `JA`), each of length $N_z + 1 = 13$:
- **`AA` Array**:
  - Positions $1 \dots n$: Main diagonal elements in order (`AA(1:5) = [1., 4., 7., 11., 12.]`).
  - Position $n + 1$: Unused dummy position (marked `*`).
  - Positions $n + 2 \dots N_z + 1$: Strictly off-diagonal non-zeros stored row-by-row (`[2., 3., 5., 6., 8., 9., 10.]`).
- **`JA` Array**:
  - Positions $1 \dots n + 1$: Pointers to the start of each row's off-diagonals in `AA`.
    - `JA(1:6) = [7, 8, 10, 13, 14, 14]`.
    - *Exam Note*: $JA(5) = JA(6) = 14$ indicates that row 5 has NO off-diagonal entries!
  - Positions $n + 2 \dots N_z + 1$: Column indices of the off-diagonal entries (`[4, 1, 4, 1, 4, 5, 3]`).

```
Index:   1    2    3    4    5    6    7    8    9   10   11   12   13
──────────────────────────────────────────────────────────────────────────
AA:   [ 1.   4.   7.  11.  12.   *    2.   3.   5.   6.   8.   9.  10.]
       └──────────┬──────────┘        └─────────────────┬────────────────┘
          Diagonal Elements                 Off-Diagonal Elements

JA:   [ 7    8   10   13   14   14    4    1    4    1    4    5    3 ]
       └──────────┬──────────┘        └─────────────────┬────────────────┘
       Row Off-Diagonal Pointers            Off-Diagonal Column Indices
```

---

## 7. Storage Schemes for Structured Matrices (DIA & ELLPACK)

For banded matrices originating from structured stencils, specialized formats eliminate pointer lookups:

Consider the $5 \times 5$ tridiagonal/banded matrix (*Slide 306*):
$$A = \begin{bmatrix}
1. & 0. & 2. & 0. & 0. \\
3. & 4. & 0. & 5. & 0. \\
0. & 6. & 7. & 0. & 8. \\
0. & 0. & 9. & 10. & 0. \\
0. & 0. & 0. & 11. & 12.
\end{bmatrix}$$

### 7.1 Diagonal Storage Format (DIA)
Non-zero diagonals are packed into a 2D dense array `DIAG(n, N_diag)`:
- Mapping: $\text{DIAG}(i, j) \leftarrow a_{i, i + \text{IOFF}(j)}$
- `IOFF`: Array of diagonal offsets relative to main diagonal ($0$). Here, active diagonals are sub-diagonal ($-1$), main diagonal ($0$), and second super-diagonal ($+2$):
  $$\text{IOFF} = [-1, \quad 0, \quad +2]$$
- `DIAG` Array (dimension $5 \times 3$, where `*` denotes zero padding outside matrix boundaries):
  $$\text{DIAG} = \begin{bmatrix}
  * & 1. & 2. \\
  3. & 4. & 5. \\
  6. & 7. & 8. \\
  9. & 10. & * \\
  11. & 12. & *
  \end{bmatrix}$$
- **Advantage**: Stride-1 sequential memory streaming with zero indirection. Reduces 2D Laplacian Jacobi solve time from $11.57\text{ s}$ to $7.12\text{ s}$!

> [!IMPORTANT]
> **Past Exam Focus (2018 Mid-Spring Q6 Part B)**:
> Exam problems frequently test unpacking DIA storage arrays back into the explicit linear system $Ax = b$ and executing in-place relaxation sweeps (e.g., Gauss-Seidel). See complete step-by-step solution in:  
> 👉 [**2018 Exam Q6 Part B: DIA Storage Reconstruction & Gauss-Seidel Analysis**](file:///home/thesupercd/Documents/HPSC_TutorialsNAssignments/ExamNotes/Previous_Years_Question_Papers_Solutions.md#question-6-part-b-10-marks).

### 7.2 ELLPACK-ITPACK Format
Designed for SIMD vector processors and GPUs where each row has at most $Nd$ non-zeros ($Nd = 3$ in this example):
- Stores two rectangular arrays of dimension $n \times Nd$ ($5 \times 3$):
  - `COEF`: Non-zero values of each row, padded with zeros as necessary.
  - `JCOEF`: Column indices of each entry in `COEF` (zero-padded positions duplicate valid column indices):

$$\text{COEF} = \begin{bmatrix}
1. & 2. & 0. \\
3. & 4. & 5. \\
6. & 7. & 8. \\
9. & 10. & 0. \\
11. & 12. & 0.
\end{bmatrix}, \qquad
\text{JCOEF} = \begin{bmatrix}
1 & 3 & 1 \\
1 & 2 & 4 \\
2 & 3 & 5 \\
3 & 4 & 4 \\
4 & 5 & 5
\end{bmatrix}$$

- **Advantage**: Column-major layout allows GPU warps to access non-zero elements across parallel threads simultaneously without thread divergence or memory uncoalescing.

> [!IMPORTANT]
> **Past Exam Focus (2025 Midsem Q6 Part B)**:
> In the 2025 exam, students were given a $10 \times 3$ ELLPACK-ITPACK representation (`COEF`, `JCOEF`, and vector `b`) and asked to solve with SOR at $\omega = 1.4$ and $1.6$, tracing the successive difference norm $\delta^{(k)}$ for iterations 9, 10, and 11. See full derivation and analytical proof using Young's theorem:  
> 👉 [**2025 Exam Q6 Part B: ELLPACK-ITPACK SOR Convergence Analysis**](file:///home/thesupercd/Documents/HPSC_TutorialsNAssignments/ExamNotes/Previous_Years_Question_Papers_Solutions.md#question-6-part-b-7-marks).

---

## 8. Comparative Summary of Sparse Storage Formats

| Format | Storage Footprint | SpMV Memory Efficiency | Best Hardware Target |
| :--- | :--- | :--- | :--- |
| **COO** | $3 N_z$ | Poor (random indirect reads) | Matrix generation, FEM assembly, Matrix Market I/O |
| **CSR** | $2 N_z + n + 1$ | High (row-contiguous streaming) | General unstructured sparse solvers (Multicore CPU) |
| **CSC** | $2 N_z + n + 1$ | High (column-contiguous streaming) | Direct column-oriented solvers (MATLAB, Fortran) |
| **MSR** | $2 N_z + 1$ (2 arrays) | Very High (direct diagonal access) | Relaxation solvers (Jacobi, GS, SSOR) needing $a_{ii}$ |
| **DIA** | $n \times N_{diag}$ | Optimal (zero indirection, stride-1) | Structured FDM grids on CPUs/Vector systems |
| **ELLPACK** | $2 \cdot n \times Nd$ | Optimal (fully coalesced warps) | GPU accelerators, CUDA SIMT architecture |

