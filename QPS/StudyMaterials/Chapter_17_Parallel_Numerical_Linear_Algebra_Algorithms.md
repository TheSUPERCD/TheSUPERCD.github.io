# Chapter 17: Parallel Numerical Linear Algebra & Memory Management

---

## 1. Fundamentals of Parallel Matrix Computations

Linear algebra computations—vector dot products, matrix-vector products ($y = Ax$), and matrix-matrix multiplications ($C = AB$)—constitute $> 90\%$ of the floating-point execution time in scientific PDE simulations.
Achieving near-peak performance on parallel supercomputers requires:
1. Decomposing data arrays across processors to balance arithmetic loads.
2. Minimizing communication volume and frequency across the interconnect network.
3. Structuring local memory accesses to maximize cache line reuse and avoid false sharing.

---

## 2. Parallel Vector Operations

### 2.1 Parallel Vector Inner Product (Dot Product)
Compute scalar $\alpha = x^T y = \sum_{i=1}^n x_i y_i$ for $x, y \in \mathbb{R}^n$.

```
       Processor 0:   [x_0..x_{n/p-1}] · [y_0..y_{n/p-1}]  ──► Local Sum s_0
       Processor 1:   [...           ] · [...           ]  ──► Local Sum s_1
       ...
       Processor p-1: [...           ] · [...           ]  ──► Local Sum s_{p-1}
                             │
                             ▼  Binary Tree Reduction (MPI_Allreduce)
                      Global Sum α = Σ s_k (broadcast to all)
```

#### Distributed Memory Implementation
1. **Data Distribution**: Vectors $x$ and $y$ are partitioned into contiguous blocks of size $n/p$ across $p$ processors (block distribution).
2. **Local Computation**: Each processor computes the dot product of its local vectors:
   $$s_k = \sum_{i=1}^{n/p} x_{k, i} y_{k, i} \quad \left[\text{Work: } \frac{2n}{p}\text{ FLOPs}\right]$$
3. **Global Communication**: A global reduction tree (`MPI_Allreduce`) aggregates the local partial sums $s_k$ into the global scalar $\alpha$ and broadcasts it to all processors.
   - Using a binary reduction tree, communication requires $\log_2 p$ steps.
   - Message size is a single floating-point word ($m = 1$).
   - Communication time:
     $$T_{comm} = (t_s + t_w) \log_2 p$$
4. **Total Parallel Execution Time**:
   $$\mathbf{T_p = \frac{2n}{p} + (t_s + t_w) \log_2 p}$$

#### Shared Memory Implementation Pitfalls & Solutions (Slides 467–468)
Consider $np$ processors, each assigned $m = n/np$ elements of vectors $a$ and $b$:

1. **The Naive Shared-Memory Flaw (Slide 467)**:
   If all processors directly accumulate into a single shared scalar variable $c$:
   ```fortran
   ! Executed concurrently by all np processors
   do i = 1, m
       c = c + a(proc_id * m + i) * b(proc_id * m + i)
   end do
   ```
   - **Hazard**: Every write to $c$ forces the hardware cache coherency controller to broadcast an invalidation across the shared bus.
   - All other caches holding $c$ are invalidated, causing extreme **bus contention, false sharing, and serialization**! The program runs significantly *slower* than a single-threaded loop.

2. **The Correct Memory Management Scheme (Slide 468)**:
   - Each processor computes its partial dot product into a **private local variable $\hat{c}$** residing in its local registers/cache:
     $$\hat{c}_{procid} = \sum_{k=1}^m \hat{a}_k \hat{b}_k$$
   - Once local accumulations are complete, processors participate in a structured **Binary Tree Reduction (`MPI_Allreduce` / OpenMP reduction)** to accumulate the global scalar:
     $$c = \sum_{procid=0}^{np-1} \hat{c}_{procid}$$
   - This completely eliminates bus thrashing during the $O(n/p)$ arithmetic loop, restricting synchronization to a single $O(\log_2 p)$ step at the end.

---

### 2.2 Vector Update (AXPY)
Compute $y \leftarrow y + \alpha x$ where $\alpha$ is a scalar and $x, y \in \mathbb{R}^n$.
- **Communication Cost**: **ZERO** (Embarrassingly parallel).
- Each processor updates its local $n/p$ elements in-place:
  $$\mathbf{T_p = \frac{2n}{p}}$$

---

## 3. Parallel Matrix-Vector Multiplication ($y = Ax$)

Let $A \in \mathbb{R}^{n \times n}$, $x \in \mathbb{R}^n$, and $y \in \mathbb{R}^n$. Total sequential work is $2n^2$ FLOPs.

### 3.1 1D Row-Wise Block Striping (Slide 470)
- **Data Distribution**:
  - Matrix $A$ is sliced into horizontal bands: each processor owns $n/p$ rows of $A$.
  - Vector $x$ is initially partitioned so that processor $k$ owns local subvector $x_k$ of length $n/p$.

```
       Processor 0:   [  Row Band 0 (n/p × n)  ]   ×   [ Full Vector x (n × 1) ] ──► y_0 (n/p)
       Processor 1:   [  Row Band 1 (n/p × n)  ]   ×   [ Full Vector x (n × 1) ] ──► y_1 (n/p)
       ...
       Processor p-1: [  Row Band p-1          ]   ×   [ Full Vector x (n × 1) ] ──► y_{p-1} (n/p)
```

#### Execution Steps and Complexity
1. **Communication Phase (All-Gather)**:
   To compute the row dot product $y_i = \sum_{j=1}^n a_{ij} x_j$, processor $k$ requires **all $n$ elements of vector $x$**!
   Therefore, each processor must broadcast its local $n/p$ chunk of $x$ to all other $p-1$ processors using an **All-Gather (all-to-all broadcast)** collective call:
   $$T_{comm} = t_s \log_2 p + t_w \left( \frac{p - 1}{p} \right) n \approx \mathbf{t_s \log_2 p + t_w n}$$
2. **Computation Phase**:
   Each processor multiplies its $n/p$ rows by the fully assembled vector $x$:
   $$T_{comp} = \frac{n}{p} \times (2n) = \mathbf{\frac{2n^2}{p}}$$
3. **Total Parallel Time**:
   $$\mathbf{T_p = \frac{2n^2}{p} + t_s \log_2 p + t_w n}$$

---

### 3.2 1D Column-Wise Block Striping
- **Data Distribution**: Each processor owns $n/p$ vertical columns of $A$ and local subvector $x_k$ (size $n/p$).
- **Computation Phase (Zero Initial Communication)**:
  Each processor computes a partial contribution to the entire output vector:
  $$y^{(k)} = A^{(k)} x^{(k)} \quad \left[\text{Work: } \frac{2n^2}{p}\text{ FLOPs}\right]$$
- **Communication Phase (Reduce-Scatter / All-to-One)**:
  The global vector $y = \sum_{k=0}^{p-1} y^{(k)}$ is accumulated via a global vector reduction across the network.
  $$T_{comm} \approx t_s \log_2 p + t_w n$$

---

### 3.3 2D Block Decomposition on a $\sqrt{p} \times \sqrt{p}$ Mesh
For massive processor counts ($p > 10,000$), 1D row striping causes communication time to dominate because message size scales as $O(n)$.
**2D Block Partitioning** drastically reduces communication by dividing both rows and columns:
- Processors are arranged as a 2D logical grid of size $\sqrt{p} \times \sqrt{p}$.
- Processor $(i, j)$ holds a submatrix block $A_{i, j}$ of size $\frac{n}{\sqrt{p}} \times \frac{n}{\sqrt{p}}$.

```
       Processor Grid (√p × √p):
       ┌──────────┬──────────┬──────────┐
       │ A_00     │ A_01     │ A_02     │   ◄── Row Group 0
       ├──────────┼──────────┼──────────┤
       │ A_10     │ A_11     │ A_12     │   ◄── Row Group 1
       ├──────────┼──────────┼──────────┤
       │ A_20     │ A_21     │ A_22     │   ◄── Row Group 2
       └──────────┴──────────┴──────────┘
            ▲          ▲          ▲
       Col Group 0 Col Group 1 Col Group 2
```

#### Three-Step Algorithm
1. **Column Broadcast**:
   Processors on the diagonal broadcast their vector chunk $x_j$ (size $n/\sqrt{p}$) vertically along their column group of $\sqrt{p}$ processors:
   $$T_{comm, 1} = t_s \log_2 \sqrt{p} + t_w \left(\frac{n}{\sqrt{p}}\right) \log_2 \sqrt{p}$$
2. **Local Matrix Multiply**:
   Each processor multiplies its $(n/\sqrt{p}) \times (n/\sqrt{p})$ block by the column-replicated vector:
   $$T_{comp} = 2 \left(\frac{n}{\sqrt{p}}\right)^2 = \mathbf{\frac{2n^2}{p}}$$
3. **Row Reduction**:
   Partial product vectors are summed horizontally along each row group of $\sqrt{p}$ processors:
   $$T_{comm, 2} = t_s \log_2 \sqrt{p} + t_w \left(\frac{n}{\sqrt{p}}\right) \log_2 \sqrt{p}$$
4. **Decisive Scalability Advantage**:
   Total communication time scales as $\mathbf{O\left(\frac{n}{\sqrt{p}} \log p\right)}$, which is **substantially smaller** than the $O(n)$ communication time of 1D striping!

---

## 4. Parallel Matrix-Matrix Multiplication ($C = AB$)

Multiplying two dense matrices $A, B \in \mathbb{R}^{n \times n}$ requires $2n^3$ FLOPs.

### 4.1 Cannon's Algorithm on a 2D Torus
Cannon’s Algorithm (1969) is the classic memory-optimal parallel algorithm for dense matrix multiplication on distributed-memory machines arranged as a 2D Torus network ($\sqrt{p} \times \sqrt{p}$).

#### Data Layout
- Let $q = \sqrt{p}$.
- Matrices $A, B, C$ are partitioned into $q \times q$ square blocks:
  $$A_{i, j}, B_{i, j}, C_{i, j} \in \mathbb{R}^{(n/q) \times (n/q)}$$
  assigned to processor $P_{i, j}$ where $i, j \in \{0, 1, \dots, q-1\}$.

```
       Phase 1: Initial Preskewing
         Row i of A shifted LEFT by i steps.
         Column j of B shifted UP by j steps.
                               │
                               ▼
       Phase 2: Shift-Multiply Loop (repeated √p times)
         1. Local Multiply: C_ij += A_ij * B_ij
         2. Shift A left by 1 step (circular)
         3. Shift B up by 1 step (circular)
```

---

### 4.2 Complete Step-by-Step Derivation of Cannon's Algorithm
1. **Initial Preskewing (Alignment Phase)**:
   Initially, processor $P_{i, j}$ holds $A_{i, j}$ and $B_{i, j}$. But to compute $C_{i, j}$, it needs $A_{i, k} B_{k, j}$.
   To ensure every processor begins with compatible blocks:
   - Circularly shift row $i$ of matrix $A$ **to the left by $i$ positions**:
     $$A_{i, j} \leftarrow A_{i, (j + i) \pmod q}$$
   - Circularly shift column $j$ of matrix $B$ **upwards by $j$ positions**:
     $$B_{i, j} \leftarrow B_{(i + j) \pmod q, j}$$
   Now, processor $P_{i, j}$ holds blocks whose inner indices match ($k = (i + j) \pmod q$)!

2. **The Iterative Multiply-and-Shift Loop**:
   Execute the following two steps $q = \sqrt{p}$ times:
   - **Local Block Multiply**:
     Each processor computes the matrix product of its currently resident blocks and accumulates into $C_{i, j}$:
     $$C_{i, j} \leftarrow C_{i, j} + A_{i, j} B_{i, j}$$
     Computational cost per stage:
     $$\text{Work} = 2 \left(\frac{n}{\sqrt{p}}\right)^3 = \frac{2n^3}{p^{3/2}}\text{ FLOPs}$$
   - **Circular Shift**:
     - Circularly shift each block of $A$ **left by 1 processor** along its row.
     - Circularly shift each block of $B$ **up by 1 processor** along its column.

3. **Restoration (Optional)**:
   After $\sqrt{p}$ stages, $A$ and $B$ are shifted back to their original positions.

---

### 4.3 Total Execution Time Derivation of Cannon's Algorithm
1. **Computational Time**:
   Across all $\sqrt{p}$ stages:
   $$T_{comp} = \sqrt{p} \times \left[ 2 \left(\frac{n}{\sqrt{p}}\right)^3 \right] = \mathbf{\frac{2n^3}{p}}$$
2. **Communication Time**:
   - Initial preskewing requires at most $\sqrt{p}-1$ shifts of size $(n/\sqrt{p})^2 = n^2/p$.
   - The main loop requires $\sqrt{p}-1$ single-hop shifts of size $n^2/p$.
   - Combining both phases (sending 1 block of $A$ and 1 block of $B$ per shift):
     $$T_{comm} \approx 2 \sqrt{p} \cdot t_s + 2 \sqrt{p} \cdot \left(\frac{n^2}{p}\right) t_w = \mathbf{2 \sqrt{p} \, t_s + 2 \frac{n^2}{\sqrt{p}} \, t_w}$$
3. **Total Parallel Time**:
   $$\mathbf{T_p = \frac{2n^3}{p} + 2 \sqrt{p} \, t_s + 2 \frac{n^2}{\sqrt{p}} \, t_w}$$

- **Isoefficiency of Cannon’s Algorithm**:
  $$\frac{T_o}{W} = \frac{p T_p - W}{W} = \frac{2 p^{1.5} t_s + 2 n^2 \sqrt{p} t_w}{2n^3} = \text{const} \implies \mathbf{W = O(p^{1.5})}$$
  Cannon’s algorithm exhibits excellent scalability ($O(p^{1.5})$), vastly outperforming naive 1D striping.

---

## 5. Shared Memory Matrix Multiplication & Cache Tiling

On multi-core shared memory CPUs and GPUs, all threads share the physical memory space (Slides 471–472).

```
       Naive 3-Loop (Cache thrashing):       Tiled / Blocked (Cache-friendly):
       for i = 1 to N                         for ii = 1 to N step B
         for j = 1 to N                         for jj = 1 to N step B
           for k = 1 to N                         for kk = 1 to N step B
             C[i,j] += A[i,k] * B[k,j]              // Sub-block BxB fits entirely in L1!
                                                    Multiply_Block(A, B, C)
```

### 5.1 Pitfalls of Naive Shared Memory Multiplication
1. **Cache Thrashing in Column Access**:
   When evaluating $C_{ij} = \sum_k A_{ik} B_{kj}$, accessing $B_{kj}$ traverses down columns. In C (row-major), this causes a cache miss on every access!
2. **Thread Contention and False Sharing**:
   If multiple threads concurrently update adjacent entries of $C$ residing on the same 64-byte line, false sharing destroys performance.

### 5.2 Cache Tiling (Loop Blocking)
- Matrix multiplication is reordered into $B \times B$ sub-blocks (tiles), where $B$ is chosen such that three $B \times B$ tiles ($A_{tile}, B_{tile}, C_{tile}$) fit entirely within the fast on-chip L1 data cache:
  $$3 \times B^2 \times 8\text{ bytes} \le \text{Size of L1 Cache (e.g., 32 KB)} \implies B \approx 32 \text{ to } 64$$
- Data is brought into L1 cache once, reused $B$ times, and written back to DRAM only once.
- Reduces DRAM memory bandwidth traffic by a factor of $B$ ($32\times$ reduction!), allowing the CPU to execute at near-peak arithmetic speed.

### 5.3 Staggered / Red-Black Scheduling to Prevent Bus Contention (Slide 472)
- In shared-memory systems, even with tiling, if multiple cores attempt to write to contiguous block tiles $C(i_1-i_2, j_1-j_2)$ simultaneously, severe interconnect and cache coherency traffic occurs.
- **Red-Black / Checkerboard Scheduling**:
  - The matrix tiles are partitioned into alternating "Red" and "Black" blocks (like a chessboard).
  - Processors update only "Red" tiles during Phase 1, and "Black" tiles during Phase 2.
  - This ensures that adjacent memory blocks are never modified concurrently, eliminating false sharing, reducing bus contention, and maximizing memory bandwidth utilization.

