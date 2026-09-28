# Chapter 13: Computer Architecture, Memory Hierarchy & Cache Optimization

---

## 1. Classical Von Neumann Sequential Computer Architecture

The foundational architecture of sequential digital computers is the **Von Neumann Model**, consisting of three fundamental subsystems:
1. **Central Processing Unit (CPU)**:
   - **Control Unit (CU)**: Decodes instructions, manages instruction pipelining, and controls data pathways.
   - **Arithmetic Logic Unit (ALU)**: Executes floating-point arithmetic (multiplication, addition) and integer logic.
   - **Registers**: Ultra-high-speed internal storage holding active instruction operands.
2. **Main Memory (DRAM)**:
   - Random Access Memory storing both program instructions and numerical data arrays.
3. **System Bus**:
   - The shared data, address, and control interconnect linking CPU and Memory.

```
       ┌───────────────────────────────┐
       │             CPU               │
       │  ┌────────────┐ ┌──────────┐  │
       │  │  Control   │ │   ALU    │  │
       │  │    Unit    │ │ (FLOPS)  │  │
       │  └────────────┘ └──────────┘  │
       │         ▲             ▲       │
       │         └──────┬──────┘       │
       │                ▼              │
       │          [ Registers ]        │
       └────────────────┬──────────────┘
                        │
                        ▼  Processor-Memory Bus (Latency & Bandwidth Bottleneck)
       ┌────────────────────────────────┐
       │          Main Memory           │
       │             (RAM)              │
       └────────────────────────────────┘
```

---

## 2. Processor Speed vs. Memory Latency: The Memory Wall

A primary theme of High Performance Scientific Computing is that arithmetic is cheap, while data movement is expensive.

### 2.1 Computer Performance Metrics
- **Clock Frequency ($f$)**: The rate at which the CPU executes clock ticks (e.g., $1.0\text{ GHz} = 10^9\text{ cycles/second}$).
- **Cycle Time ($t_{cycle}$)**: Duration of one clock tick:
  $$t_{cycle} = \frac{1}{f} = \frac{1}{10^9\text{ Hz}} = 10^{-9}\text{ seconds} = \mathbf{1.0\text{ nanosecond (ns)}}$$
- **Ideal Peak Processor Speed**:
  Assuming the ALU executes 4 floating-point operations per cycle (e.g., Fused Multiply-Add over SIMD lanes):
  $$\text{Ideal Speed} = 10^9\text{ cycles/s} \times 4\text{ FLOP/cycle} = \mathbf{4.0\text{ GigaFLOPS (GFLOPS)}}$$
- **RAM Bandwidth**: The rate at which continuous data streams can be fetched from RAM (bytes or words per second).
- **RAM Latency**: The elapsed time between the CPU requesting data at a memory address and the data arriving at CPU registers.

### 2.2 The Quantitative Memory Wall Derivation (Slides 2 & 344–345)
Consider what happens when a processor executes a linear algebra operation (e.g., $y_i = y_i + a_{ij} x_j$) without caching:
- Standard DRAM latency is approximately:
  $$t_{latency} \approx \mathbf{100\text{ ns}} = 10^{-7}\text{ seconds}$$
- Because $t_{cycle} = 1\text{ ns}$, a $100\text{ ns}$ latency corresponds to **100 wasted CPU clock cycles** during which the ALU is completely starved of data and forced to idle (CPU stall)!
- If the CPU must fetch operands directly from main RAM for every operation:
  $$\text{Effective Execution Frequency} = \frac{1}{\text{Latency}} = \frac{1}{100\text{ ns}} = \frac{1}{10^{-7}\text{ s}} = \mathbf{10\text{ MegaFLOPS (MFLOPS)}}$$
- **Catastrophic Degradation**:
  $$\frac{\text{Ideal Peak Speed}}{\text{Uncached Effective Speed}} = \frac{4\text{ GFLOPS}}{10\text{ MFLOPS}} = \mathbf{400\times \text{ Slowdown!}}$$
  The CPU delivers less than $0.25\%$ of its rated floating-point potential!
- **Governing Law of Sequential Architecture**:
  $$\mathbf{\text{Computer Performance} = \min(\text{Processor Speed}, \, \text{Memory Latency / Bandwidth})}$$

---

## 3. The Memory Hierarchy & Locality Principles

To mitigate the 400-fold processor-memory performance gap, computer architects introduce a tiered hierarchy of cache memories between the CPU and DRAM.

```
       Fastest, Smallest,   ▲   [ CPU Registers ]      (< 1 ns, ~1 KB)
       Most Expensive       │   [ L1 Data / Inst Cache] (1 - 2 ns, 32 - 64 KB/core)
                            │   [ L2 Cache ]            (3 - 10 ns, 512 KB - 1 MB/core)
                            │   [ L3 Cache (Shared) ]   (10 - 20 ns, 16 - 64 MB/socket)
       Slowest, Largest,    │   [ Main Memory (DRAM) ]  (50 - 100 ns, 16 - 512 GB)
       Cheapest             ▼   [ Secondary Storage ]   (Microseconds to Milliseconds, TBs)
```

### 3.1 Principles of Locality
Caches function effectively only because scientific algorithms exhibit two predictable access patterns:
1. **Temporal Locality**:
   - If a memory location is referenced at time $t$, it is highly likely to be referenced again in the near future.
   - *Examples*: Loop indices, scalar multipliers $\alpha$ and $\beta$ in Conjugate Gradient, reused pivot elements.
2. **Spatial Locality**:
   - If a memory location at address $A$ is referenced, memory locations with adjacent addresses $A + 1, A + 2, \dots$ are highly likely to be referenced soon.
   - *Examples*: Stepping sequentially through elements of a contiguous vector $x[i]$, streaming matrix rows.

### 3.2 Cache Lines and the Cache Hit Ratio
- **Cache Line (Cache Block)**: Data is never transferred between DRAM and cache as single words; it is moved in contiguous chunks called **cache lines** (typically **$64\text{ bytes} = 8\text{ double-precision (64-bit) floats}$**).
- **Cache Hit**: The requested memory address is already resident in cache. Access latency is $1 - 4$ cycles.
- **Cache Miss**: The address is absent from cache. Execution halts while a full 64-byte line is fetched from DRAM ($100 - 200$ cycles).
- **Cache Hit Ratio ($H$)**:
  $$H = \frac{\text{Number of Cache Hits}}{\text{Total Memory References}}$$
  $$\text{Effective Access Time (EAT)} = H \cdot t_{cache} + (1 - H) \cdot t_{DRAM}$$
  Even a tiny $1\%$ drop in hit ratio ($H$ falling from $99\%$ to $98\%$) can double the effective memory access latency!

### 3.3 The Three C's of Cache Misses
1. **Compulsory (Cold) Misses**: The very first time a memory block is referenced; unavoidable.
2. **Capacity Misses**: The working set of the problem exceeds the physical capacity of the cache.
3. **Conflict Misses**: Multiple active memory locations map to the exact same cache set (in direct-mapped or set-associative caches).

---

## 4. Memory Layout in High-Level Languages: Row-Major vs. Column-Major

Matrix entries $A_{ij}$ exist mathematically on a 2D grid, but physical RAM is an unrolled 1D linear array of contiguous byte addresses.

```
       Row-Major (C/C++):      Row 0 [0,0][0,1][0,2] ──► Row 1 [1,0][1,1][1,2] ──► ...
       Column-Major (Fortran): Col 0 [0,0][1,0][2,0] ──► Col 1 [0,1][1,1][2,1] ──► ...
```

### 4.1 Row-Major Order (C, C++, Python/NumPy)
- Elements of the same row are placed in contiguous memory addresses.
- Physical address of entry $A[i][j]$ for an $M \times N$ matrix:
  $$\text{Address}(A[i][j]) = \text{Base} + (i \cdot N + j) \times \text{sizeof(element)}$$
- **Cache-Friendly Loop (Stride-1)**:
  Looping over $j$ in the inner loop:
  ```c
  for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
          sum += A[i][j] * x[j];  // Stride-1 access: A[i][j] and A[i][j+1] are adjacent!
      }
  }
  ```

### 4.2 Column-Major Order (Fortran, MATLAB)
- Elements of the same column are placed in contiguous memory addresses.
- Physical address of entry $A(i, j)$ for an $M \times N$ matrix:
  $$\text{Address}(A(i, j)) = \text{Base} + ((j - 1) \cdot M + (i - 1)) \times \text{sizeof(element)}$$
- **Cache-Friendly Loop (Stride-1)**:
  Looping over $i$ in the inner loop:
  ```fortran
  do j = 1, N
      do i = 1, M
          y(i) = y(i) + A(i, j) * x(j)  ! Stride-1 access: A(i, j) and A(i+1, j) are adjacent!
      end do
  end do
  ```

---

### 2.2 Concrete Parameters of the Memory Hierarchy (Slide 117)
The lecture presents typical latency and capacity values across tiers:

| Tier | Component | Typical Capacity | Access Latency | Bus / Interface Speed |
| :--- | :--- | :--- | :--- | :--- |
| **Level 0** | CPU Registers | $1\text{ KB}$ | $\mathbf{300\text{ ps}}$ ($0.3\text{ ns}$) | Internal CPU cycle |
| **Level 1** | L1 Cache | $64\text{ KB}$ | $\mathbf{1\text{ ns}}$ ($1 - 4$ cycles) | Dedicated L1 bus |
| **Level 2** | L2 Cache | $256\text{ KB}$ | $\mathbf{\sim 10\text{ ns}}$ | Dedicated L2 bus |
| **Level 3** | L3 Cache (Shared) | $2 - 4\text{ MB}$ | $\mathbf{\sim 20\text{ ns}}$ | On-die interconnect |
| **Level 4** | Main Memory (DRAM) | Gigabytes ($16 - 128\text{ GB}$) | $\mathbf{\sim 100\text{ ns}}$ ($200\times$ slower than L0!) | Off-chip Memory Bus |

---

## 5. Cache Miss Issues in Matrix-Vector Multiplication ($y = Ax$)

### 5.1 Naive Access in Fortran vs. Cache Hits (Slides 116–118)
Consider evaluating the matrix-vector dot product $y = Ax$ in Fortran:
$$\sum_{j=1}^n a_{ij} x_j$$
If written naively with the row index $j$ varying in the inner loop:
- In Fortran, matrices are stored in **column-major order**. Therefore, elements in the same row $a_{ij}$ and $a_{i, j+1}$ are separated in memory by $M$ elements ($M \times 8\text{ bytes}$).
- For a large matrix ($M = 10,000$), $M \times 8\text{ B} = 80\text{ KB}$, which exceeds the $64\text{ KB}$ L1 cache size!
- Every single iteration of the inner loop jumps across cache lines, triggering a **cache miss on almost every memory operation**.
- Conversely, storing the matrix transposed ($A^T$) or iterating down columns ensures contiguous memory access ($a_{1i}, a_{2i}, a_{3i}, \dots$), achieving stride-1 access. When $a_{1i}$ is fetched, the next 7 double-precision floats are pulled into the 64-byte cache line automatically, yielding an **$87.5\%$ cache hit rate**!

### 5.2 Quantitative Experimental Benchmark (Slide 116)
Executing matrix-vector multiplication with $n = 10,000$ on a physical test system:

```fortran
! Program 1: Row-wise access in Fortran (Cache Unfriendly)
do i = 1, n
    tmp = 0.0
    do j = 1, n
        tmp = tmp + Mat(i, j) * B(j)
    end do
    C(i) = tmp
end do
! CPU Time: 0.92885900 seconds

! Program 2: Column-wise access in Fortran (Cache Friendly, symmetric Mat)
do i = 1, n
    tmp = 0.0
    do j = 1, n
        tmp = tmp + Mat(j, i) * B(j)
    end do
    C(i) = tmp
end do
! CPU Time: 0.45792997 seconds
```

> [!NOTE]
> **Key Exam Observation**: Column-wise memory access is **more than $2.03\times$ faster** than row-wise access ($0.9288\text{ s} \to 0.4579\text{ s}$) solely due to spatial cache locality, with zero difference in floating-point operations!

### 5.3 C Implementation with Pointer Arithmetic (Slide 115)
In C, 2D matrices flattened as 1D arrays require row-major indexing `*(A + (i*Ndim + j))`:
```c
for (i = 0; i < Ndim; i++) {
    tmp = 0.0;
    for (j = 0; j < Ndim; j++) {
        tmp += *(A + (i*Ndim + j)) * *(B + j);
    }
    *(C + i) = tmp;
}
```
Serial time for dense $N = 1000$ was measured at $0.007596\text{ seconds}$. Operations scale as $O(N^2)$.

---

## 6. Cache Optimization in Sparse Iterative Solvers (Slide 121)

### 6.1 The 2D Laplacian Jacobi Solver Bottleneck
Discretizing the 2D Laplace equation $\nabla^2 T = 0$ on an $N_x \times N_y$ mesh ($N = N_x N_y = 10,000$) produces a sparse pentadiagonal matrix where each row has at most 5 non-zeros:
$$x_i^{(k+1)} = \frac{b_i - a_{i, i-N_x} x_{i-N_x}^{(k)} - a_{i, i-1} x_{i-1}^{(k)} - a_{i, i+1} x_{i+1}^{(k)} - a_{i, i+N_x} x_{i+N_x}^{(k)}}{a_{ii}}$$

### 6.2 Quantitative Timing Comparison Across 3 Implementations
The lecture presents three distinct algorithmic approaches for $N = 10,000$:

| Implementation Scheme | Storage Format | Memory Access Pattern | Measured CPU Time (Slide 121) | Relative Speedup |
| :--- | :--- | :--- | :--- | :--- |
| **1. Full Dense Row Multiplication** | Dense $N \times N$ matrix | Scans all $10,000$ columns per row (mostly zeros) | **$710\text{ seconds}$** | $1.0\times$ (Baseline) |
| **2. Non-Zero Sparse Access** | Naive COO/CSR sparse arrays | Computes only the 5 non-zeros per row, but indirect lookups | **$11.57\text{ seconds}$** | **$61.4\times$ faster** |
| **3. Diagonal Storage (DIA) Access** | DIA format ($N \times 5$ array) | Contiguous stride-1 streaming of diagonal vectors | **$7.12\text{ seconds}$** | **$100\times$ faster** ($1.62\times$ over naive sparse) |

> [!IMPORTANT]
> **Exam Takeaway**: Transitioning from full row to non-zero sparse multiplication eliminates redundant arithmetic ($O(N^2) \to O(N)$), saving $698.4\text{ seconds}$. Transitioning from indirect sparse access to DIA diagonal streaming eliminates cache misses, saving an additional $38.5\%$ execution time!


## 6. Case Study: Optimization of the 2D Laplacian Jacobi Solver

The lecture slides detail a concrete numerical experiment solving the 2D Laplace equation via Jacobi iteration on an $N = 10,000$ matrix (Slide 121):
$$x_i^{(k+1)} = \frac{b_i - \sum_{j \ne i} a_{ij} x_j^{(k)}}{a_{ii}}$$

### 6.1 Naive Sparse Implementation vs. Diagonal Storage Scheme
1. **Naive Non-Zero Element Access**:
   - Implemented using general row-sparse pointers with indirect memory addressing.
   - Non-zero multiplications trigger repeated cache misses due to irregular strides across far diagonals ($\pm N_x$).
   - **Sequential Execution Time**: **$11.57\text{ seconds}$**.
2. **Cache-Friendly Diagonal Storage (DIA) Implementation**:
   - The pentadiagonal matrix is stored as 5 contiguous linear column vectors (`offsets` = $[-N_x, -1, 0, 1, N_x]$).
   - Traversing the diagonals involves strictly sequential, unit-stride memory streams.
   - Cache prefetchers automatically stream subsequent cache lines into L1/L2 caches ahead of time.
   - **Sequential Execution Time**: **$7.12\text{ seconds}$**.

### 6.2 Key Takeaways for HPC
- **Runtime Reduction**: $\frac{11.57 - 7.12}{11.57} \times 100\% = \mathbf{38.46\% \text{ Faster}}$ on the exact same CPU hardware, with **zero changes to the mathematical algorithm**!
- Algorithmic FLOP count alone does **not** determine execution time; data locality, cache line utilization, and stride alignment dominate real-world HPC performance.
