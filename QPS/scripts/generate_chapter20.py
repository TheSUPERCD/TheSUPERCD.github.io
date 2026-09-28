#!/usr/bin/env python3
"""
generate_chapter20.py
Generates chapters/chapter-20.html and chapters/data/chapter-20.js
for Chapter 20: Sparse Matrix Storage Formats & High-Performance Kernels:
COO, CSR, CSC, DIA, ELLPACK, ELLPACK-ITPACK (Ellpack-ltpack), and MSR.
"""
import os
import json

def get_chapter20_body_html():
    return r"""
<h1 id="chapter-20-sparse-matrix-storage-formats--high-performance-kernels">Chapter 20: Sparse Matrix Storage Formats &amp; High-Performance Kernels: CSR, CSC, DIA, ELLPACK, ELLPACK-ITPACK, and COO</h1>
<p class="lead" style="font-size: 1.15rem; color: var(--text-muted); margin-bottom: 2rem; border-bottom: 1px dashed var(--header-border-color); padding-bottom: 1rem;">
<strong>Master Engineering Guide &amp; Theoretical Synthesis:</strong> Exhaustive mathematical foundations, storage architectures, step-by-step conversion hand calculations, memory footprint derivations, GPU coalescing mechanics, and production-grade C/OpenMP/CUDA SpMV kernels spanning <strong>CSR</strong>, <strong>CSC</strong>, <strong>DIA</strong>, <strong>ELLPACK</strong>, <strong>ELLPACK-ITPACK</strong> (Ellpack-ltpack), <strong>COO</strong>, and <strong>MSR</strong>.
</p>

<!-- SECTION 1 -->
<section id="ch20-sec1">
<h2 id="1-foundations-of-sparse-representations-in-scientific-computing">1. Foundations of Sparse Representations in Scientific Computing</h2>

<div class="cram-box intuition-box" style="border-left: 4px solid var(--primary); background: var(--bg-alt); padding: 1.25rem; border-radius: 0 8px 8px 0; margin: 1.5rem 0;">
<h4 style="margin-top: 0; color: var(--primary);">💡 Why Sparse Formats Dictate the Limits of Scientific Computing</h4>
<p>
In continuum physics—whether computing turbulent Navier-Stokes flow around a hypersonic re-entry vehicle, simulating mantle convection across geological epochs, or solving the Schrödinger wave equation—physical laws are inherently <em>local</em>. A molecule of air or a patch of electric potential only directly exchanges momentum or flux with its immediate spatial neighbors. When we discretize these continuous differential equations onto a discrete mesh of $N$ nodes, each node interacts with only a tiny constant number of adjacent nodes ($k \ll N$).
<br/><br/>
If we store the resulting coefficient matrix $A \in \mathbb{R}^{N \times N}$ in a conventional dense 2D array, we commit two fatal computational sins:
<br/>
1. <strong>Memory Suicide:</strong> We store millions of billions of literal zeroes, exhausting RAM and forcing the operating system to thrash or crash.
<br/>
2. <strong>Flop Waste:</strong> In matrix-vector multiplication $y = Ax$, our CPUs and GPUs execute trillions of pointless multiplications by zero ($0.0 \times x_j = 0.0$), stalling high-performance supercomputing clusters on useless work.
<br/><br/>
<strong>Sparse matrix storage formats</strong> are the fundamental data structures invented to store <em>only</em> the non-zero values and their topological coordinates. Selecting the correct format is often the difference between a simulation completing in 10 minutes on a single GPU or running out of memory on a 1,000-node cluster.
</p>
</div>

<h3 id="11-definition--geometry-of-matrix-sparsity">1.1 Definition &amp; Geometry of Matrix Sparsity</h3>
<p>
Let $A \in \mathbb{R}^{M \times N}$ be a matrix, and let $N_{nz}$ (or $\text{nnz}$) denote the total count of non-zero entries:
</p>
$$N_{nz} = \# \left\{ (i, j) \mid A_{ij} \ne 0 \right\}$$

<p>We quantify the sparsity of matrix $A$ using two dual metrics:</p>
<ul>
  <li><strong>Density ($\rho$):</strong> The ratio of non-zero elements to total matrix capacity:
  $$\rho = \frac{N_{nz}}{M \cdot N}, \quad 0 \le \rho \le 1$$
  </li>
  <li><strong>Sparsity Ratio ($S$):</strong> The fraction of elements that are strictly zero:
  $$S = 1 - \rho = 1 - \frac{N_{nz}}{M \cdot N} = \frac{M \cdot N - N_{nz}}{M \cdot N}$$
  </li>
</ul>

<div class="exam-highlight-box" style="border: 1px solid var(--primary-border); background: var(--primary-subtle); padding: 1rem 1.25rem; border-radius: 8px; margin: 1rem 0;">
  <strong>Exam Rule of Thumb:</strong> In High Performance Scientific Computing, a matrix is formally classified as <strong>sparse</strong> if:
  $$N_{nz} = O(N) \quad \text{or} \quad N_{nz} = O(N \log N) \quad \text{as } N \to \infty$$
  Under this asymptotic scaling, the matrix density vanishes towards zero as the mesh is refined ($\lim_{N \to \infty} \rho = 0$). For a million-node grid ($N = 10^6$), a typical PDE matrix has $\rho \approx 0.000005$ ($99.9995\%$ zeroes!).
</div>

<h3 id="12-physical-and-computational-origins-of-sparsity">1.2 Physical and Computational Origins of Sparsity</h3>
<p>
Sparsity is not an arbitrary mathematical artifact; it arises directly from the local discretization of continuous partial differential equations (PDEs):
</p>

<table>
<thead>
<tr>
  <th style="text-align: left;">Physical Discretization Scheme</th>
  <th style="text-align: center;">Spatial Dimension</th>
  <th style="text-align: center;">Local Stencil Connectivity</th>
  <th style="text-align: center;">Non-Zeros per Row ($k$)</th>
  <th style="text-align: left;">Matrix Structure &amp; Topology</th>
</tr>
</thead>
<tbody>
<tr>
  <td><strong>1D Heat Conduction (FDM)</strong></td>
  <td style="text-align: center;">1D</td>
  <td style="text-align: center;">$[-1, \; 2, \; -1]$ (3-point)</td>
  <td style="text-align: center;">$k = 3$</td>
  <td>Tridiagonal, symmetric positive definite (SPD)</td>
</tr>
<tr>
  <td><strong>2D Laplace / Poisson (FDM / FVM)</strong></td>
  <td style="text-align: center;">2D</td>
  <td style="text-align: center;">5-Point Star Stencil</td>
  <td style="text-align: center;">$k = 5$</td>
  <td>Pentadiagonal, block-tridiagonal with bandwidth $N_x$</td>
</tr>
<tr>
  <td><strong>2D Convection-Diffusion (FDM)</strong></td>
  <td style="text-align: center;">2D</td>
  <td style="text-align: center;">9-Point Compact Stencil</td>
  <td style="text-align: center;">$k = 9$</td>
  <td>Non-symmetric, banded block-tridiagonal</td>
</tr>
<tr>
  <td><strong>3D Poisson Equation (FDM / FVM)</strong></td>
  <td style="text-align: center;">3D</td>
  <td style="text-align: center;">7-Point Cross Stencil</td>
  <td style="text-align: center;">$k = 7$</td>
  <td>Septadiagonal, block-tridiagonal with bandwidth $N_x N_y$</td>
</tr>
<tr>
  <td><strong>3D Unstructured Mesh (FEM / FVM)</strong></td>
  <td style="text-align: center;">3D</td>
  <td style="text-align: center;">Tetrahedral / Hexahedral elements</td>
  <td style="text-align: center;">$k \approx 14 \text{ to } 27$</td>
  <td>Unstructured sparsity pattern (random row non-zeros)</td>
</tr>
<tr>
  <td><strong>Network Graphs / Circuit Sim</strong></td>
  <td style="text-align: center;">Arbitrary</td>
  <td style="text-align: center;">Kirchhoff nodal laws</td>
  <td style="text-align: center;">Variable ($k = 2 \text{ to } 10^4$)</td>
  <td>Power-law degree distribution, highly irregular hubs</td>
</tr>
</tbody>
</table>

<h3 id="13-the-catastrophic-costs-of-dense-storage">1.3 The Catastrophic Costs of Dense Storage</h3>
<p>
To appreciate the sheer necessity of sparse matrix representations, consider a modest 3D grid of size $100 \times 100 \times 100$, producing $N = 10^6$ linear algebraic equations with double-precision floating-point coefficients (IEEE 754 64-bit float, $8\text{ bytes}$ per entry):
</p>

<table>
<thead>
<tr>
  <th style="text-align: left;">Resource / Metric</th>
  <th style="text-align: left;">Dense Storage ($N \times N$ 2D Array)</th>
  <th style="text-align: left;">Sparse 7-Point Stencil Storage</th>
  <th style="text-align: center;">Savings Factor</th>
</tr>
</thead>
<tbody>
<tr>
  <td><strong>RAM Memory Footprint</strong></td>
  <td>$N^2 \times 8\text{ B} = 10^{12} \times 8\text{ B} = \mathbf{8{,}000\text{ Gigabytes (8 Terabytes!)}}$</td>
  <td>$7 \times 10^6 \times (8 + 4)\text{ B} \approx \mathbf{84\text{ Megabytes}}$</td>
  <td style="text-align: center;"><strong style="color: #10b981;">95,238× Less RAM!</strong></td>
</tr>
<tr>
  <td><strong>Matrix-Vector Product (FLOPs)</strong></td>
  <td>$2 N^2 = 2 \times 10^{12}\text{ FLOPs (2 TeraFLOPs)}$</td>
  <td>$2 N_{nz} = 2 \times 7 \times 10^6 = \mathbf{1.4 \times 10^7\text{ FLOPs (14 MFLOPs)}}$</td>
  <td style="text-align: center;"><strong style="color: #10b981;">142,857× Fewer FLOPs!</strong></td>
</tr>
<tr>
  <td><strong>Compute Time at 10 GFLOPS</strong></td>
  <td>$T = \frac{2 \times 10^{12}}{10^{10}} = \mathbf{200\text{ seconds (3.3 minutes)}}$ per iteration</td>
  <td>$T = \frac{1.4 \times 10^7}{10^{10}} = \mathbf{0.0014\text{ seconds (1.4 ms)}}$ per iteration</td>
  <td style="text-align: center;"><strong style="color: #10b981;">142,857× Faster!</strong></td>
</tr>
<tr>
  <td><strong>CG Solver (1,000 Iterations)</strong></td>
  <td>$\sim 55\text{ hours}$ (excluding swapping death!)</td>
  <td>$\mathbf{\sim 1.4\text{ seconds}}$ total wall-clock time</td>
  <td style="text-align: center;"><strong style="color: #10b981;">Real-Time Feasibility</strong></td>
</tr>
</tbody>
</table>

<h3 id="14-the-sparse-computational-dilemma-the-memory-wall">1.4 The Sparse Computational Dilemma: The Memory Wall</h3>
<p>
While sparse formats save massive memory and arithmetic operations, they introduce severe architectural bottlenecks on modern CPU and GPU microarchitectures:
</p>

<ol>
  <li><strong>The Indirection Penalty (Memory Gather):</strong>
  In dense matrix-vector multiply $y_i = \sum_{j} A_{ij} x_j$, elements of $x$ are read in continuous sequential order. In sparse multiplication:
  $$y_i = \sum_{k} \text{val}[k] \cdot x[\text{col}[k]]$$
  The column indices $\text{col}[k]$ force the CPU/GPU memory controller to execute an <em>indirect gather</em> ($x[\text{col}[k]]$). If non-zero columns are far apart, successive loads pull distinct 64-byte cache lines from slow DRAM, causing frequent cache misses and TLB misses.
  </li>
  <li><strong>Abysmally Low Arithmetic Intensity ($I$):</strong>
  Recall the Roofline Model: Arithmetic Intensity is $I = \frac{\text{FLOPs}}{\text{Bytes Transferred}}$.
  In double-precision Sparse Matrix-Vector Multiply (SpMV):
  $$\text{FLOPs per non-zero} = 2 \quad (\text{1 multiplication} + 1 \text{ addition})$$
  $$\text{Bytes transferred per non-zero} \ge 8\text{ B (double val)} + 4\text{ B (int col)} = 12\text{ bytes}$$
  $$\mathbf{I = \frac{2\text{ FLOPs}}{12\text{ Bytes}} \approx 0.167\text{ FLOPs/Byte}}$$
  On an NVIDIA H100 GPU ($P_{\text{peak}} \approx 60\text{ TFLOPS}$, $BW_{\text{peak}} \approx 3,350\text{ GB/s}$), the machine balance is $I^* = 60{,}000 / 3{,}350 \approx 17.9\text{ FLOPs/Byte}$.
  Since $0.167 \ll 17.9$, <strong>SpMV is 100% memory-bandwidth bound</strong>! The compute cores are starved of data $> 98\%$ of the time.
  </li>
  <li><strong>Hardware Divergence on Parallel Systems:</strong>
  When parallelizing across threads or GPU warps (groups of 32 threads):
  <ul>
    <li>If Row $i$ has 2 non-zeros and Row $i+1$ has 50 non-zeros, threads assigned to short rows finish immediately and idle, suffering from <strong>warp divergence</strong> and <strong>load imbalance</strong>.</li>
    <li>If threads access scattered column indices, memory requests cannot be merged, destroying <strong>memory coalescing</strong>.</li>
  </ul>
  </li>
</ol>

<h3 id="15-the-two-canonical-benchmark-matrices">1.5 The Two Canonical Benchmark Matrices Used Throughout This Chapter</h3>
<p>
To ensure absolute clarity, cross-comparability, and pedagogical consistency, we use <strong>two canonical matrices</strong> across every single format in this chapter:
</p>

<div class="cram-box notation-box" style="border-left: 4px solid #8b5cf6; background: var(--bg-alt); padding: 1.25rem; border-radius: 0 8px 8px 0; margin: 1.5rem 0;">
<h4 style="margin-top: 0; color: #8b5cf6;">📐 Canonical Matrix A: General Unstructured $5 \times 5$ Matrix ($N = 5, N_{nz} = 12$)</h4>
<p>From Yousef Saad's classic textbook (<em>Iterative Methods for Sparse Linear Systems</em>), representing general unstructured meshes with variable row lengths (1 to 4 non-zeros per row):</p>
$$A = \begin{bmatrix}
1.0 & 0.0 & 0.0 & 2.0 & 0.0 \\
3.0 & 4.0 & 0.0 & 5.0 & 0.0 \\
6.0 & 0.0 & 7.0 & 8.0 & 9.0 \\
0.0 & 0.0 & 10.0 & 11.0 & 0.0 \\
0.0 & 0.0 & 0.0 & 0.0 & 12.0
\end{bmatrix}$$
<ul>
  <li>Row 0: 2 non-zeros at columns 0, 3</li>
  <li>Row 1: 3 non-zeros at columns 0, 1, 3</li>
  <li>Row 2: 4 non-zeros at columns 0, 2, 3, 4</li>
  <li>Row 3: 2 non-zeros at columns 2, 3</li>
  <li>Row 4: 1 non-zero at column 4</li>
</ul>
</div>

<div class="cram-box notation-box" style="border-left: 4px solid #06b6d4; background: var(--bg-alt); padding: 1.25rem; border-radius: 0 8px 8px 0; margin: 1.5rem 0;">
<h4 style="margin-top: 0; color: #06b6d4;">📐 Canonical Matrix B: Structured Banded $5 \times 5$ Stencil Matrix ($N = 5, N_{nz} = 11$)</h4>
<p>From course lecture slides, representing structured PDE finite-difference stencils where non-zeros lie strictly along diagonals:</p>
$$B = \begin{bmatrix}
1.0 & 0.0 & 2.0 & 0.0 & 0.0 \\
3.0 & 4.0 & 0.0 & 5.0 & 0.0 \\
0.0 & 6.0 & 7.0 & 0.0 & 8.0 \\
0.0 & 0.0 & 9.0 & 10.0 & 0.0 \\
0.0 & 0.0 & 0.0 & 11.0 & 12.0
\end{bmatrix}$$
<ul>
  <li>Sub-diagonal (offset $k = -1$): $[3.0, 6.0, 9.0, 11.0]$</li>
  <li>Main diagonal (offset $k = 0$): $[1.0, 4.0, 7.0, 10.0, 12.0]$</li>
  <li>Super-diagonal (offset $k = +2$): $[2.0, 5.0, 8.0]$</li>
</ul>
</div>

</section>

<!-- SECTION 2 -->
<section id="ch20-sec2">
<h2 id="2-coordinate-format-coo--triplet-format">2. Coordinate Format (COO / Triplet Format)</h2>

<h3 id="21-coo-concept--data-structure">2.1 Concept &amp; Data Structure</h3>
<p>
The <strong>Coordinate Format (COO)</strong>, also universally known as the <strong>Triplet format</strong> or <strong>ijv format</strong>, is the most intuitive and conceptually primitive sparse matrix representation.
It represents a matrix purely as an unordered list of non-zero coordinates and their corresponding numerical values.
</p>
<p>It utilizes <strong>three parallel 1D arrays</strong>, each of length exactly equal to $N_{nz}$:</p>
<ul>
  <li><code>val</code> (or <code>AA</code> / <code>DATA</code>): Floating-point array containing the numerical values of the non-zeros.</li>
  <li><code>row</code> (or <code>IR</code> / <code>I</code>): Integer array containing the row index of each entry.</li>
  <li><code>col</code> (or <code>JC</code> / <code>J</code>): Integer array containing the column index of each entry.</li>
</ul>

<div class="exam-highlight-box" style="border: 1px solid var(--primary-border); background: var(--primary-subtle); padding: 1rem 1.25rem; border-radius: 8px; margin: 1rem 0;">
  <strong>Key Properties of COO:</strong>
  <ul>
    <li><strong>Order Invariance:</strong> Triplets $(i, j, A_{ij})$ may appear in <em>any arbitrary order</em>. They do not need to be sorted by row or column.</li>
    <li><strong>Additive Duplicate Entries:</strong> Most numerical assembly libraries (e.g. SciPy <code>scipy.sparse.coo_matrix</code>, PETSc, SuiteSparse) define the coordinate format such that if the same coordinate $(i, j)$ appears multiple times, their values are <strong>summed together</strong>:
    $$A_{ij} = \sum_{k: \text{row}[k]=i, \, \text{col}[k]=j} \text{val}[k]$$
    This makes COO the supreme, undisputed format for <strong>Finite Element Method (FEM) global matrix assembly</strong>!
    </li>
  </ul>
</div>

<h3 id="22-worked-example-converting-canonical-matrix-a-to-coo">2.2 Worked Example: Converting Canonical Matrix A to COO</h3>
<p>Let us write out the COO representation for Canonical Matrix $A$ ($N = 5, N_{nz} = 12$):</p>

<p><strong>Option 1: Row-Major Sorted COO (0-Indexed, C/Python Convention):</strong></p>
<table>
<thead>
<tr>
  <th style="text-align: center;">Array Index $k$</th>
  <th style="text-align: center;">0</th>
  <th style="text-align: center;">1</th>
  <th style="text-align: center;">2</th>
  <th style="text-align: center;">3</th>
  <th style="text-align: center;">4</th>
  <th style="text-align: center;">5</th>
  <th style="text-align: center;">6</th>
  <th style="text-align: center;">7</th>
  <th style="text-align: center;">8</th>
  <th style="text-align: center;">9</th>
  <th style="text-align: center;">10</th>
  <th style="text-align: center;">11</th>
</tr>
</thead>
<tbody>
<tr>
  <td><code>row[k]</code></td>
  <td style="text-align: center;">0</td>
  <td style="text-align: center;">0</td>
  <td style="text-align: center;">1</td>
  <td style="text-align: center;">1</td>
  <td style="text-align: center;">1</td>
  <td style="text-align: center;">2</td>
  <td style="text-align: center;">2</td>
  <td style="text-align: center;">2</td>
  <td style="text-align: center;">2</td>
  <td style="text-align: center;">3</td>
  <td style="text-align: center;">3</td>
  <td style="text-align: center;">4</td>
</tr>
<tr>
  <td><code>col[k]</code></td>
  <td style="text-align: center;">0</td>
  <td style="text-align: center;">3</td>
  <td style="text-align: center;">0</td>
  <td style="text-align: center;">1</td>
  <td style="text-align: center;">3</td>
  <td style="text-align: center;">0</td>
  <td style="text-align: center;">2</td>
  <td style="text-align: center;">3</td>
  <td style="text-align: center;">4</td>
  <td style="text-align: center;">2</td>
  <td style="text-align: center;">3</td>
  <td style="text-align: center;">4</td>
</tr>
<tr>
  <td><code>val[k]</code></td>
  <td style="text-align: center;">1.0</td>
  <td style="text-align: center;">2.0</td>
  <td style="text-align: center;">3.0</td>
  <td style="text-align: center;">4.0</td>
  <td style="text-align: center;">5.0</td>
  <td style="text-align: center;">6.0</td>
  <td style="text-align: center;">7.0</td>
  <td style="text-align: center;">8.0</td>
  <td style="text-align: center;">9.0</td>
  <td style="text-align: center;">10.0</td>
  <td style="text-align: center;">11.0</td>
  <td style="text-align: center;">12.0</td>
</tr>
</tbody>
</table>

<p><strong>Option 2: Completely Unsorted COO (1-Indexed, Fortran/Course Lecture Convention):</strong></p>
<pre><code>Index:   1     2     3     4     5     6     7     8     9    10    11    12
─────────────────────────────────────────────────────────────────────────────
AA:   [ 12.0   9.0   7.0   5.0   1.0   2.0  11.0   3.0   6.0   4.0   8.0  10.0 ]
JR:   [    5     3     3     2     1     1     4     2     3     2     3     4 ]
JC:   [    5     5     3     4     1     4     4     1     1     2     4     3 ]</code></pre>
<p><em>Exam Tip:</em> Both representations encode the exact same mathematical matrix! Notice that entry $A(3, 4) = 8.0$ (row 3, col 4 in 1-based indexing) is located at index 11 in the unsorted lecture array ($JR[11]=3, JC[11]=4, AA[11]=8.0$) and at index 7 in the sorted 0-based array ($row[7]=2, col[7]=3, val[7]=8.0$).</p>

<h3 id="23-coo-sparse-matrix-vector-multiplication-spmv">2.3 COO Sparse Matrix-Vector Multiplication (SpMV)</h3>
<p>
The sequential SpMV algorithm for $y = \alpha A x + \beta y$ under COO is extraordinarily simple:
</p>

<pre><code class="language-c">// Sequential COO SpMV: y = A * x
void spmv_coo(int nnz, const int *row, const int *col, const double *val, 
              const double *x, double *y, int n) {
    // 1. Initialize output vector y to zero
    for (int i = 0; i &lt; n; i++) {
        y[i] = 0.0;
    }
    // 2. Stream through all non-zeros
    for (int k = 0; k &lt; nnz; k++) {
        y[row[k]] += val[k] * x[col[k]];
    }
}</code></pre>

<div class="cram-box trap-box" style="border-left: 4px solid #ef4444; background: var(--bg-alt); padding: 1.25rem; border-radius: 0 8px 8px 0; margin: 1.5rem 0;">
<h4 style="margin-top: 0; color: #ef4444;">⚠️ The Parallel Write-Collision Hazard in COO SpMV</h4>
<p>
Why don't high-performance linear solvers use COO for Krylov iterations (CG, GMRES)?
Look closely at the inner update:
<code>y[row[k]] += val[k] * x[col[k]];</code>
<br/>
If we parallelize this loop across multiple CPU threads or GPU cores using <code>#pragma omp parallel for</code>:
Multiple threads processing different non-zeros $k_1$ and $k_2$ may encounter the <strong>same row index</strong> ($\text{row}[k_1] == \text{row}[k_2]$). Both threads will attempt to read, add, and write back to the exact same memory location <code>y[i]</code> simultaneously!
<br/>
This causes a <strong>critical race condition</strong> (data hazard). To make it correct, you must use <code>#pragma omp atomic</code> or atomic hardware instructions (<code>atomicAdd</code> on CUDA).
Atomic instructions serialize memory writes, destroying parallel throughput! Therefore, COO is virtually <strong>never</strong> used for iterative solver execution.
</p>
</div>

<h3 id="24-summary-of-coo">2.4 Summary of COO</h3>
<ul>
  <li><strong>Memory Footprint:</strong> Exactly $3 N_{nz}$ words $= N_{nz} \times (\text{sizeof(double)} + 2 \times \text{sizeof(int)}) = \mathbf{16 N_{nz}\text{ bytes}}$.</li>
  <li><strong>Lookup Complexity $A(i, j)$:</strong> $O(N_{nz})$ for unsorted COO; $O(\log N_{nz})$ for sorted COO (binary search).</li>
  <li><strong>Ideal Use Cases:</strong> Finite element stiffness matrix assembly, mesh generation, reading/writing disk files in the Matrix Market (<code>.mtx</code>) standard.</li>
</ul>

</section>

<!-- SECTION 3 -->
<section id="ch20-sec3">
<h2 id="3-compressed-sparse-row-csr--crs--yale-format">3. Compressed Sparse Row (CSR / CRS / Yale Format)</h2>

<h3 id="31-csr-concept--compression-mechanism">3.1 Concept &amp; Compression Mechanism</h3>
<p>
The <strong>Compressed Sparse Row (CSR)</strong> format—also historically known as the <strong>Yale Sparse Matrix Format</strong> or <strong>Compressed Row Storage (CRS)</strong>—is the universal industry and academic standard for general unstructured sparse matrix computations.
</p>
<p>
<strong>The Fundamental Insight of CSR:</strong> In the COO format, row indices in the <code>row</code> array are repeated over and over for every non-zero in that row (e.g. if row 2 has 100 non-zeros, the number <code>2</code> is written 100 times!).
CSR eliminates this massive redundancy by <strong>compressing the row indices into an array of pointers</strong>.
</p>

<p>CSR represents an $M \times N$ matrix using <strong>three arrays</strong>:</p>
<ol>
  <li><code>val</code> (or <code>AA</code> / <code>A</code>): Floating-point array containing the numerical values of the non-zeros, stored <strong>row by row</strong> from top to bottom, and left-to-right within each row. (Size $= N_{nz}$).</li>
  <li><code>col_ind</code> (or <code>JA</code>): Integer array containing the column indices corresponding to each element in <code>val</code>. (Size $= N_{nz}$).</li>
  <li><code>row_ptr</code> (or <code>IA</code>): Integer array containing the starting offset/index in <code>val</code> and <code>col_ind</code> for each row. (Size $= M + 1$).</li>
</ol>

<div class="cram-box notation-box" style="border-left: 4px solid var(--primary); background: var(--bg-alt); padding: 1.25rem; border-radius: 0 8px 8px 0; margin: 1.5rem 0;">
<h4 style="margin-top: 0; color: var(--primary);">📐 The Mathematical Invariants of CSR</h4>
<p>For any valid CSR representation (0-indexed):</p>
<ul>
  <li><code>row_ptr[0] = 0</code> (the first row starts at index 0).</li>
  <li><code>row_ptr[M] = N_nz</code> (the sentinel element points to the total number of non-zeros).</li>
  <li><strong>Row Non-Zero Count:</strong> The number of non-zero entries in row $i$ is given instantaneously by:
  $$\text{nnz}_i = \text{row\_ptr}[i+1] - \text{row\_ptr}[i]$$
  </li>
  <li><strong>Range of Row $i$:</strong> The elements of row $i$ are stored in the index slice:
  $$k \in [\text{row\_ptr}[i], \; \text{row\_ptr}[i+1] - 1]$$
  </li>
  <li><strong>Empty Rows:</strong> If row $i$ contains only zeroes, then:
  $$\text{row\_ptr}[i] = \text{row\_ptr}[i+1] \iff \text{nnz}_i = 0$$
  </li>
</ul>
</div>

<h3 id="32-worked-example-converting-canonical-matrix-a-to-csr">3.2 Worked Example: Converting Canonical Matrix A to CSR</h3>
<p>Recall Canonical Matrix $A$ ($N = 5, N_{nz} = 12$):</p>
$$A = \begin{bmatrix}
1.0 & 0.0 & 0.0 & 2.0 & 0.0 \\
3.0 & 4.0 & 0.0 & 5.0 & 0.0 \\
6.0 & 0.0 & 7.0 & 8.0 & 9.0 \\
0.0 & 0.0 & 10.0 & 11.0 & 0.0 \\
0.0 & 0.0 & 0.0 & 0.0 & 12.0
\end{bmatrix}$$

<p><strong>Step-by-Step Derivation of CSR Arrays:</strong></p>
<ol>
  <li><strong>Row 0:</strong> Non-zeros are $A(0, 0) = 1.0$ and $A(0, 3) = 2.0$.
  <br/>Values: <code>[1.0, 2.0]</code>. Columns: <code>[0, 3]</code>. Starts at index $0$. Length $= 2$.
  </li>
  <li><strong>Row 1:</strong> Non-zeros are $A(1, 0) = 3.0, A(1, 1) = 4.0, A(1, 3) = 5.0$.
  <br/>Values: <code>[3.0, 4.0, 5.0]</code>. Columns: <code>[0, 1, 3]</code>. Starts at index $0 + 2 = 2$. Length $= 3$.
  </li>
  <li><strong>Row 2:</strong> Non-zeros are $A(2, 0) = 6.0, A(2, 2) = 7.0, A(2, 3) = 8.0, A(2, 4) = 9.0$.
  <br/>Values: <code>[6.0, 7.0, 8.0, 9.0]</code>. Columns: <code>[0, 2, 3, 4]</code>. Starts at index $2 + 3 = 5$. Length $= 4$.
  </li>
  <li><strong>Row 3:</strong> Non-zeros are $A(3, 2) = 10.0, A(3, 3) = 11.0$.
  <br/>Values: <code>[10.0, 11.0]</code>. Columns: <code>[2, 3]</code>. Starts at index $5 + 4 = 9$. Length $= 2$.
  </li>
  <li><strong>Row 4:</strong> Non-zero is $A(4, 4) = 12.0$.
  <br/>Value: <code>[12.0]</code>. Column: <code>[4]</code>. Starts at index $9 + 2 = 11$. Length $= 1$.
  </li>
  <li><strong>Sentinel:</strong> Total non-zeros $= 11 + 1 = 12$. Index $= 12$.
  </li>
</ol>

<p><strong>Complete CSR Arrays (0-Indexed, C/C++ Standard):</strong></p>
<pre><code>val     = [ 1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0, 11.0, 12.0 ]  (length 12)
col_ind = [   0,    3,    0,    1,    3,    0,    2,    3,    4,    2,    3,    4 ]  (length 12)
row_ptr = [   0,    2,    5,    9,   11,   12 ]                                        (length 6 = N+1)</code></pre>

<p><strong>Complete CSR Arrays (1-Indexed, Fortran/Saad Textbook Standard):</strong></p>
<pre><code>AA = [ 1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0, 11.0, 12.0 ]
JA = [   1,    4,    1,    2,    4,    1,    3,    4,    5,    3,    4,    5 ]
IA = [   1,    3,    6,   10,   12,   13 ]</code></pre>

<div class="cram-box walkthrough-box" style="border-left: 4px solid #10b981; background: var(--bg-alt); padding: 1.25rem; border-radius: 0 8px 8px 0; margin: 1.5rem 0;">
<h4 style="margin-top: 0; color: #10b981;">✍️ Exam Walkthrough: How to Query Element $A(i, j)$ in CSR</h4>
<p><strong>Problem:</strong> Using the 0-indexed CSR arrays above, determine whether element $A(2, 3)$ is non-zero, and if so, find its value.</p>
<ol>
  <li><strong>Step 1: Locate Row 2 Bounds:</strong>
  Look up <code>row_ptr[2]</code> and <code>row_ptr[3]</code>:
  $$\text{start} = \text{row\_ptr}[2] = 5, \quad \text{end} = \text{row\_ptr}[3] = 9$$
  Row 2 entries occupy indices $k \in [5, 8]$.
  </li>
  <li><strong>Step 2: Search Column Indices in <code>col_ind[5 : 8]</code>:</strong>
  Examine the column entries at indices $5, 6, 7, 8$:
  $$k = 5 \implies \text{col\_ind}[5] = 0$$
  $$k = 6 \implies \text{col\_ind}[6] = 2$$
  $$k = 7 \implies \text{col\_ind}[7] = 3 \quad \mathbf{\leftarrow \text{Match found for column } j = 3!}$$
  </li>
  <li><strong>Step 3: Retrieve Value:</strong>
  Fetch the value from <code>val[7]</code>:
  $$A(2, 3) = \text{val}[7] = \mathbf{8.0}$$
  <em>Complexity:</em> If columns are sorted within each row, search takes $O(\log(\text{nnz}_i))$ via binary search!
</ol>
</div>

<h3 id="33-high-performance-spmv-in-csr">3.3 High-Performance SpMV in CSR</h3>
<p>
The CSR Sparse Matrix-Vector Multiply ($y = A x$) is the central computational kernel of all modern numerical software packages (Intel MKL, cuSPARSE, PETSc, SciPy):
</p>

<pre><code class="language-c">// Production-Grade C99 CSR SpMV with OpenMP Multicore Parallelization
void spmv_csr(int m, const int *row_ptr, const int *col_ind, const double *val,
              const double *x, double *y) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i &lt; m; i++) {
        double sum = 0.0;
        int row_start = row_ptr[i];
        int row_end   = row_ptr[i + 1];
        
        #pragma GCC unroll 4
        for (int k = row_start; k &lt; row_end; k++) {
            sum += val[k] * x[col_ind[k]];
        }
        y[i] = sum;
    }
}</code></pre>

<div class="exam-highlight-box" style="border: 1px solid var(--primary-border); background: var(--primary-subtle); padding: 1rem 1.25rem; border-radius: 8px; margin: 1rem 0;">
  <strong>Why CSR SpMV is Embarrassingly Parallel on Multicore CPUs:</strong>
  Notice that thread $i$ computes only $y[i]$. <strong>Every single row writes to a distinct entry of $y$!</strong>
  There are strictly <strong>ZERO write collisions</strong>. No locks, no atomic directives, and no reduction trees are needed. This allows linear scaling up to dozens of CPU cores on structured grids.
</div>

<h3 id="34-gpu-pathology-why-csr-struggles-on-cuda">3.4 GPU Pathology: Why CSR Struggles on CUDA</h3>
<p>
On NVIDIA GPUs, assigning 1 thread per row in CSR causes devastating performance loss due to:
</p>
<ol>
  <li><strong>Warp Divergence:</strong> A GPU executes 32 threads in lockstep (a warp). If thread 0 processes a row with 2 entries and thread 1 processes a row with 32 entries, thread 0 must sit idle for 30 cycles waiting for thread 1!</li>
  <li><strong>Uncoalesced Memory Access:</strong> While <code>val[k]</code> is read consecutively by an individual thread across time, across the 32 threads in a warp at any single cycle, the threads are reading from widely separated memory addresses, shattering memory coalescing.</li>
</ol>

</section>

<!-- SECTION 4 -->
<section id="ch20-sec4">
<h2 id="4-compressed-sparse-column-csc--ccs--harwell-boeing-format">4. Compressed Sparse Column (CSC / CCS / Harwell-Boeing Format)</h2>

<h3 id="41-csc-concept--column-compression">4.1 Concept &amp; Column Compression</h3>
<p>
The <strong>Compressed Sparse Column (CSC)</strong> format—also called <strong>Compressed Column Storage (CCS)</strong> or the <strong>Harwell-Boeing format</strong>—is the exact mathematical dual of CSR.
Instead of compressing rows, CSC slices the matrix <strong>column by column</strong> and compresses the column indices into column pointers.
</p>

<p>CSC consists of <strong>three arrays</strong> for an $M \times N$ matrix:</p>
<ol>
  <li><code>val</code> (or <code>AA</code>): Non-zero numerical values stored <strong>column by column</strong> from left to right, top to bottom. (Size $= N_{nz}$).</li>
  <li><code>row_ind</code> (or <code>JA</code>): Integer row indices of each entry in <code>val</code>. (Size $= N_{nz}$).</li>
  <li><code>col_ptr</code> (or <code>IA</code>): Integer column pointers marking where each column begins in <code>val</code> and <code>row_ind</code>. (Size $= N + 1$).</li>
</ol>

<h3 id="42-the-fundamental-transpose-duality">4.2 The Fundamental Transpose Duality</h3>
<div class="cram-box intuition-box" style="border-left: 4px solid #f59e0b; background: var(--bg-alt); padding: 1.25rem; border-radius: 0 8px 8px 0; margin: 1.5rem 0;">
<h4 style="margin-top: 0; color: #f59e0b;">⭐ The Master Transpose Identity: CSR(Aᵀ) ≡ CSC(A)</h4>
<p>
A profound mathematical duality connects CSR and CSC:
<br/>
<strong>The CSC representation of any matrix $A$ is precisely identical to the CSR representation of its transpose $A^T$!</strong>
$$\mathbf{\text{CSC}(A) \equiv \text{CSR}(A^T)}$$
<br/>
<em>Operational Miracle:</em> If you have a matrix stored in CSR format and you need to compute with its transpose $A^T$, you do <strong>not</strong> need to move a single byte of data! You simply pass the CSR arrays into a CSC algorithm, and it acts as $A^T$ with zero computational cost.
</p>
</div>

<h3 id="43-worked-example-converting-canonical-matrix-a-to-csc">4.3 Worked Example: Converting Canonical Matrix A to CSC</h3>
<p>Recall Canonical Matrix $A$ ($N = 5, N_{nz} = 12$):</p>
$$A = \begin{bmatrix}
1.0 & 0.0 & 0.0 & 2.0 & 0.0 \\
3.0 & 4.0 & 0.0 & 5.0 & 0.0 \\
6.0 & 0.0 & 7.0 & 8.0 & 9.0 \\
0.0 & 0.0 & 10.0 & 11.0 & 0.0 \\
0.0 & 0.0 & 0.0 & 0.0 & 12.0
\end{bmatrix}$$

<p><strong>Step-by-Step Column Scan:</strong></p>
<ul>
  <li><strong>Col 0:</strong> $A(0,0)=1.0, A(1,0)=3.0, A(2,0)=6.0$. Values: <code>[1.0, 3.0, 6.0]</code>. Rows: <code>[0, 1, 2]</code>. Starts at $0$.</li>
  <li><strong>Col 1:</strong> $A(1,1)=4.0$. Values: <code>[4.0]</code>. Rows: <code>[1]</code>. Starts at $0 + 3 = 3$.</li>
  <li><strong>Col 2:</strong> $A(2,2)=7.0, A(3,2)=10.0$. Values: <code>[7.0, 10.0]</code>. Rows: <code>[2, 3]</code>. Starts at $3 + 1 = 4$.</li>
  <li><strong>Col 3:</strong> $A(0,3)=2.0, A(1,3)=5.0, A(2,3)=8.0, A(3,3)=11.0$. Values: <code>[2.0, 5.0, 8.0, 11.0]</code>. Rows: <code>[0, 1, 2, 3]</code>. Starts at $4 + 2 = 6$.</li>
  <li><strong>Col 4:</strong> $A(2,4)=9.0, A(4,4)=12.0$. Values: <code>[9.0, 12.0]</code>. Rows: <code>[2, 4]</code>. Starts at $6 + 4 = 10$.</li>
  <li><strong>Sentinel:</strong> Total non-zeros $= 10 + 2 = 12$.</li>
</ul>

<p><strong>Complete CSC Arrays (0-Indexed):</strong></p>
<pre><code>val     = [ 1.0, 3.0, 6.0, 4.0, 7.0, 10.0, 2.0, 5.0, 8.0, 11.0, 9.0, 12.0 ]  (length 12)
row_ind = [   0,   1,   2,   1,   2,    3,   0,   1,   2,    3,   2,    4 ]  (length 12)
col_ptr = [   0,   3,   4,   6,  10,   12 ]                                  (length 6 = N+1)</code></pre>

<h3 id="44-side-by-side-comparison-csr-vs-csc-arrays">4.4 Side-by-Side Comparison: CSR vs. CSC Arrays</h3>
<table>
<thead>
<tr>
  <th style="text-align: left;">Format</th>
  <th style="text-align: left;">Data Values Array (<code>val</code>)</th>
  <th style="text-align: left;">Coordinate Array (<code>col_ind</code> / <code>row_ind</code>)</th>
  <th style="text-align: left;">Pointer Array (<code>row_ptr</code> / <code>col_ptr</code>)</th>
</tr>
</thead>
<tbody>
<tr>
  <td><strong>CSR</strong></td>
  <td><code>[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]</code></td>
  <td><code>[0, 3, 0, 1, 3, 0, 2, 3, 4,  2,  3,  4]</code> (col)</td>
  <td><code>[0, 2, 5, 9, 11, 12]</code> (row start)</td>
</tr>
<tr>
  <td><strong>CSC</strong></td>
  <td><code>[1, 3, 6, 4, 7, 10, 2, 5, 8, 11, 9, 12]</code></td>
  <td><code>[0, 1, 2, 1, 2,  3, 0, 1, 2,  3,  2,  4]</code> (row)</td>
  <td><code>[0, 3, 4, 6, 10, 12]</code> (col start)</td>
</tr>
</tbody>
</table>

<h3 id="45-spmv-algorithm-in-csc-the-axpy-column-sweep">4.5 SpMV Algorithm in CSC (The AXPY / Column-Sweep)</h3>
<pre><code class="language-c">// CSC SpMV: y = A * x (Column-Oriented AXPY formulation)
void spmv_csc(int n, const int *col_ptr, const int *row_ind, const double *val,
              const double *x, double *y, int m) {
    for (int i = 0; i &lt; m; i++) y[i] = 0.0;
    
    for (int j = 0; j &lt; n; j++) {
        double xj = x[j];
        int col_start = col_ptr[j];
        int col_end   = col_ptr[j + 1];
        
        for (int k = col_start; k &lt; col_end; k++) {
            y[row_ind[k]] += val[k] * xj;  // &lt;-- Scatter write!
        }
    }
}</code></pre>

<div class="cram-box trap-box" style="border-left: 4px solid #ef4444; background: var(--bg-alt); padding: 1.25rem; border-radius: 0 8px 8px 0; margin: 1.5rem 0;">
<h4 style="margin-top: 0; color: #ef4444;">⚠️ Why CSC SpMV is Terrible for Parallel CPU Execution</h4>
<p>
In the CSC SpMV loop, the outer loop iterates over columns $j$. In the inner loop, non-zero entries in column $j$ update <code>y[row_ind[k]] += val[k] * xj</code>.
<br/>
If multiple threads process different columns $j_1$ and $j_2$ in parallel, both columns can have non-zeros in the <strong>same row</strong>, causing <strong>write contention on vector $y$</strong>!
Consequently, parallel CSC SpMV requires thread-private output vectors followed by a reduction, or heavy atomic locking.
<br/>
<strong>Where CSC Shines:</strong> CSC is the native storage format of <strong>MATLAB</strong> and direct sparse factorization libraries (<strong>SuperLU</strong>, <strong>UMFPACK</strong>, <strong>CHOLMOD</strong>) because Gaussian elimination on sparse matrices involves column-slicing, column scaling, and column rank-1 updates.
</p>
</div>

</section>

<!-- SECTION 5 -->
<section id="ch20-sec5">
<h2 id="5-diagonal-storage-format-dia--diag">5. Diagonal Storage Format (DIA / DIAG)</h2>

<h3 id="51-concept--mathematical-formulation">5.1 Concept &amp; Mathematical Formulation</h3>
<p>
The <strong>Diagonal Storage Format (DIA)</strong> is engineered specifically for banded matrices originating from finite difference (FDM) and finite volume (FVM) discretizations on structured orthogonal grids.
On regular grids, non-zeros do not appear at random coordinates; they lie strictly along continuous diagonal lines parallel to the main diagonal.
</p>

<p><strong>The Diagonal Offset Formulation:</strong></p>
<p>For any matrix entry $A_{ij}$, its <strong>diagonal offset</strong> $k$ is defined by the column-row difference:
$$\mathbf{k = j - i}$$
</p>
<ul>
  <li><strong>Main Diagonal ($k = 0$):</strong> All entries $A_{ii}$ where $j = i$.</li>
  <li><strong>Super-Diagonals ($k > 0$):</strong> Diagonals sitting strictly above the main diagonal (e.g. $k = +1$ is the first super-diagonal, $k = +N_x$ is the far block-diagonal).</li>
  <li><strong>Sub-Diagonals ($k < 0$):</strong> Diagonals sitting strictly below the main diagonal (e.g. $k = -1$ is the first sub-diagonal, $k = -N_x$ is the far lower block-diagonal).</li>
</ul>

<p>DIA compresses an $N \times N$ matrix with $N_{\text{diag}}$ active diagonals into <strong>two arrays</strong>:</p>
<ol>
  <li><code>offset</code> (or <code>IOFF</code>): 1D integer array of size $N_{\text{diag}}$ containing the active diagonal offsets sorted in increasing order.</li>
  <li><code>data</code> (or <code>DIAG</code>): 2D floating-point array of dimension $N \times N_{\text{diag}}$ (or $N_{\text{diag}} \times N$). Entry $\text{data}[i][d]$ stores element $A_{i, \; i + \text{offset}[d]}$.</li>
</ol>

<h3 id="52-the-boundary-padding-rule">5.2 The Boundary Padding Rule</h3>
<p>
Because sub-diagonals ($k < 0$) start below row 0, their top $|k|$ positions do not correspond to any valid matrix coordinate.
Similarly, super-diagonals ($k > 0$) terminate before row $N-1$, so their bottom $k$ positions have no valid matrix entries.
In DIA, these out-of-bounds positions are <strong>padded with dummy values</strong> (conventionally denoted by an asterisk <code>*</code> or zero <code>0.0</code>).
</p>

<h3 id="53-worked-example-converting-canonical-matrix-b-to-dia">5.3 Worked Example: Converting Canonical Matrix B to DIA</h3>
<p>Recall Canonical Matrix $B$ ($N = 5$ banded matrix from course lecture slides):</p>
$$B = \begin{bmatrix}
1.0 & 0.0 & 2.0 & 0.0 & 0.0 \\
3.0 & 4.0 & 0.0 & 5.0 & 0.0 \\
0.0 & 6.0 & 7.0 & 0.0 & 8.0 \\
0.0 & 0.0 & 9.0 & 10.0 & 0.0 \\
0.0 & 0.0 & 0.0 & 11.0 & 12.0
\end{bmatrix}$$

<p><strong>Step-by-Step Diagonal Identification:</strong></p>
<ol>
  <li><strong>Sub-diagonal ($k = -1$):</strong> Entries are $B(1,0)=3.0, B(2,1)=6.0, B(3,2)=9.0, B(4,3)=11.0$.
  <br/>Row 0 has no subdiagonal entry $\implies$ padded with <code>*</code>.
  </li>
  <li><strong>Main diagonal ($k = 0$):</strong> Entries are $B(0,0)=1.0, B(1,1)=4.0, B(2,2)=7.0, B(3,3)=10.0, B(4,4)=12.0$.
  </li>
  <li><strong>Second super-diagonal ($k = +2$):</strong> Entries are $B(0,2)=2.0, B(1,3)=5.0, B(2,4)=8.0$.
  <br/>Rows 3 and 4 have no entries on offset $+2 \implies$ padded with <code>*</code>.
  </li>
</ol>

<p><strong>Resulting DIA Arrays:</strong></p>
$$\text{offset} = \begin{bmatrix} -1 & 0 & +2 \end{bmatrix} \quad (N_{\text{diag}} = 3)$$

$$\text{data} = \begin{bmatrix}
* & 1.0 & 2.0 \\
3.0 & 4.0 & 5.0 \\
6.0 & 7.0 & 8.0 \\
9.0 & 10.0 & * \\
11.0 & 12.0 & *
\end{bmatrix}$$

<h3 id="54-the-zero-indirection-streaming-spmv-kernel">5.4 The Zero-Indirection Streaming SpMV Kernel</h3>
<p>
The greatest virtue of the DIA format is that it <strong>completely eliminates column index arrays</strong>!
The column index is computed via simple register arithmetic: $j = i + \text{offset}[d]$.
</p>

<pre><code class="language-c">// High-Performance DIA SpMV: y = B * x
void spmv_dia(int n, int ndiag, const int *offset, const double *data,
              const double *x, double *y) {
    // 1. Clear output vector
    for (int i = 0; i &lt; n; i++) y[i] = 0.0;
    
    // 2. Stream through active diagonals
    for (int d = 0; d &lt; ndiag; d++) {
        int k = offset[d];
        // Calculate valid row bounds for diagonal k
        int i_start = (k &lt; 0) ? -k : 0;
        int i_end   = (k &gt; 0) ? n - k : n;
        
        #pragma omp parallel for schedule(static)
        for (int i = i_start; i &lt; i_end; i++) {
            // Memory access data[i * ndiag + d] is stride-1 or prefetchable!
            // x[i + k] is regular stride-1 access!
            y[i] += data[i * ndiag + d] * x[i + k];
        }
    }
}</code></pre>

<div class="cram-box intuition-box" style="border-left: 4px solid #10b981; background: var(--bg-alt); padding: 1.25rem; border-radius: 0 8px 8px 0; margin: 1.5rem 0;">
<h4 style="margin-top: 0; color: #10b981;">🚀 Why DIA Delivers Up to 40% Speedup over CSR</h4>
<p>
In the lecture slides benchmark (2D Laplace Jacobi iterative solve on a structured grid):
<br/>
• <strong>CSR Jacobi Runtime:</strong> $11.57\text{ seconds}$
<br/>
• <strong>DIA Jacobi Runtime:</strong> $\mathbf{7.12\text{ seconds}}$ ($\mathbf{38.5\%}$ runtime reduction!)
<br/><br/>
<em>Why does this happen?</em>
1. <strong>Zero Indirection:</strong> No integer <code>col_ind</code> array is loaded from DRAM. This saves $33\%$ of the memory bandwidth!
<br/>
2. <strong>Hardware Prefetching:</strong> Both <code>data</code> and vector <code>x</code> are accessed with constant unit strides, allowing the CPU L1/L2 stream prefetchers to operate at $100\%$ efficiency with near-zero pipeline stalls.
</p>
</div>

<h3 id="55-the-pathological-failure-mode-of-dia">5.5 The Pathological Failure Mode of DIA</h3>
<div class="cram-box trap-box" style="border-left: 4px solid #ef4444; background: var(--bg-alt); padding: 1.25rem; border-radius: 0 8px 8px 0; margin: 1.5rem 0;">
<h4 style="margin-top: 0; color: #ef4444;">⚠️ The Outlier Diagonal Padding Catastrophe</h4>
<p>
DIA is strictly for structured banded systems. If a matrix has a single non-zero outlier in the top-right corner ($A(0, N-1) \ne 0 \implies k = N-1$), DIA must allocate an entire column of length $N$ in <code>data</code> just to store that <strong>one</strong> entry!
<br/>
Total storage becomes $N \times N_{\text{diag}}$. If an unstructured mesh has 1,000 distinct diagonal offsets, DIA requires $1000 \times N$ words—vastly larger than CSR, wasting gigabytes of memory on padded zeroes.
</p>
</div>

</section>

<!-- SECTION 6 -->
<section id="ch20-sec6">
<h2 id="6-ellpack-ell-format">6. ELLPACK (ELL Format)</h2>

<h3 id="61-concept--motivation-vector-supercomputers--gpus">6.1 Concept &amp; Motivation: Vector Supercomputers &amp; GPUs</h3>
<p>
The <strong>ELLPACK format</strong> (named after the ELLPACK numerical PDE package developed at Purdue University by John Rice et al.) was designed to solve the fatal flaw of CSR on SIMD (Single Instruction, Multiple Data) vector processors and modern GPU accelerators.
</p>
<p>
<strong>The Central Thesis of ELLPACK:</strong> In physical PDE grids, almost every node has at most $K_{\max}$ connections.
For example:
<ul>
  <li>2D 5-Point Stencil: $K_{\max} = 5$</li>
  <li>3D 7-Point Stencil: $K_{\max} = 7$</li>
  <li>2D 9-Point Stencil: $K_{\max} = 9$</li>
</ul>
Instead of using variable-length lists of pointers like CSR, <strong>ELLPACK compresses all rows into a fixed rectangular array of width $K_{\max}$</strong>.
Any row having fewer than $K_{\max}$ non-zeros is padded with zeroes.
</p>

<p>ELLPACK stores an $N \times N$ matrix using <strong>two 2D rectangular arrays</strong> of dimension $N \times K_{\max}$:</p>
<ol>
  <li><code>COEF</code> (or <code>AS</code> / <code>val</code>): Floating-point array of size $N \times K_{\max}$. Row $i$ contains the non-zero values of row $i$, padded with $0.0$ to length $K_{\max}$.</li>
  <li><code>JCOEF</code> (or <code>JA</code> / <code>col</code>): Integer array of size $N \times K_{\max}$. Stores the corresponding column indices for each entry in <code>COEF</code>. For padded positions, a dummy/sentinel column index (such as the row index $i$, 0, or $-1$) is stored.</li>
</ol>

<h3 id="62-worked-example-converting-canonical-matrix-a-to-ellpack">6.2 Worked Example: Converting Canonical Matrix A to ELLPACK</h3>
<p>Recall Canonical Matrix $A$ ($N = 5$):</p>
$$A = \begin{bmatrix}
1.0 & 0.0 & 0.0 & 2.0 & 0.0 \\
3.0 & 4.0 & 0.0 & 5.0 & 0.0 \\
6.0 & 0.0 & 7.0 & 8.0 & 9.0 \\
0.0 & 0.0 & 10.0 & 11.0 & 0.0 \\
0.0 & 0.0 & 0.0 & 0.0 & 12.0
\end{bmatrix}$$

<p><strong>Step 1: Determine Maximum Non-Zeros per Row ($K_{\max}$):</strong></p>
<ul>
  <li>Row 0: 2 entries $\implies [1.0, 2.0]$ at cols $[0, 3]$</li>
  <li>Row 1: 3 entries $\implies [3.0, 4.0, 5.0]$ at cols $[0, 1, 3]$</li>
  <li>Row 2: 4 entries $\implies [6.0, 7.0, 8.0, 9.0]$ at cols $[0, 2, 3, 4]$</li>
  <li>Row 3: 2 entries $\implies [10.0, 11.0]$ at cols $[2, 3]$</li>
  <li>Row 4: 1 entry $\implies [12.0]$ at col $[4]$</li>
</ul>
$$\mathbf{K_{\max} = \max(2, 3, 4, 2, 1) = 4}$$

<p><strong>Step 2: Pack into Rectangular $5 \times 4$ Arrays (0-Indexed):</strong></p>
$$\text{COEF} = \begin{bmatrix}
1.0 & 2.0 & 0.0 & 0.0 \\
3.0 & 4.0 & 5.0 & 0.0 \\
6.0 & 7.0 & 8.0 & 9.0 \\
10.0 & 11.0 & 0.0 & 0.0 \\
12.0 & 0.0 & 0.0 & 0.0
\end{bmatrix}, \qquad
\text{JCOEF} = \begin{bmatrix}
0 & 3 & 0 & 0 \\
0 & 1 & 3 & 1 \\
0 & 2 & 3 & 4 \\
2 & 3 & 2 & 2 \\
4 & 4 & 4 & 4
\end{bmatrix}$$
<p><em>Padding Convention Note:</em> For zero-padded slots in <code>COEF</code>, <code>JCOEF</code> is filled with a valid column index (often repeating the row index or the last valid column) to prevent out-of-bounds array reads in $x[\text{JCOEF}[i][c]] \times 0.0$.</p>

<p><strong>Course Lecture Representation (1-Indexed, Banded Example $Nd = 3$):</strong></p>
$$\text{COEF} = \begin{bmatrix}
1.0 & 2.0 & 0.0 \\
3.0 & 4.0 & 5.0 \\
6.0 & 7.0 & 8.0 \\
9.0 & 10.0 & 0.0 \\
11.0 & 12.0 & 0.0
\end{bmatrix}, \qquad
\text{JCOEF} = \begin{bmatrix}
1 & 3 & 1 \\
1 & 2 & 4 \\
2 & 3 & 5 \\
3 & 4 & 4 \\
4 & 5 & 5
\end{bmatrix}$$

<h3 id="63-gpu-memory-coalescing-row-major-vs-column-major-layout">6.3 GPU Memory Coalescing: Row-Major vs. Column-Major Layout</h3>
<p>
The true magic of ELLPACK on GPU hardware depends entirely on its <strong>in-memory array storage ordering</strong>:
</p>

<table>
<thead>
<tr>
  <th style="text-align: left;">Architecture / Target</th>
  <th style="text-align: left;">Array Memory Layout</th>
  <th style="text-align: left;">Indexing Formula</th>
  <th style="text-align: left;">Hardware Access Pattern</th>
</tr>
</thead>
<tbody>
<tr>
  <td><strong>CPU (Multithreaded)</strong></td>
  <td>Row-Major Layout</td>
  <td><code>COEF[i * K_max + c]</code></td>
  <td>Each thread iterates through contiguous memory for row $i$. High L1 spatial cache locality.</td>
</tr>
<tr>
  <td><strong>GPU (CUDA / SIMT Warps)</strong></td>
  <td><strong>Column-Major Layout</strong></td>
  <td><code>COEF[c * N + i]</code></td>
  <td><strong>Perfect 100% Memory Coalescing!</strong> 32 threads in a warp access 32 consecutive addresses in DRAM.</td>
</tr>
</tbody>
</table>

<div class="cram-box intuition-box" style="border-left: 4px solid #8b5cf6; background: var(--bg-alt); padding: 1.25rem; border-radius: 0 8px 8px 0; margin: 1.5rem 0;">
<h4 style="margin-top: 0; color: #8b5cf6;">💡 The GPU Memory Coalescing Mechanism Explained</h4>
<p>
In CUDA, assign Thread $i$ to process Row $i$.
<br/>
At iteration $c = 0$, all 32 threads in Warp 0 read their first element simultaneously:
<ul>
  <li>Thread 0 reads <code>COEF[0 * N + 0]</code></li>
  <li>Thread 1 reads <code>COEF[0 * N + 1]</code></li>
  <li>...</li>
  <li>Thread 31 reads <code>COEF[0 * N + 31]</code></li>
</ul>
Because the array is stored in <strong>column-major order</strong>, these 32 double-precision floating-point numbers sit in a <strong>single continuous 256-byte chunk of global memory</strong>!
The GPU memory controller issues a single 256-byte coalesced memory transaction to DRAM, fully saturating memory bus bandwidth.
<br/>
Furthermore, since every thread runs a loop from $c = 0$ to $K_{\max}-1$, <strong>there is zero warp divergence</strong>!
</p>
</div>

</section>

<!-- SECTION 7 -->
<section id="ch20-sec7">
<h2 id="7-ellpack-itpack-format-and-ellpack-ltpack-disambiguation">7. ELLPACK-ITPACK Format (and ELLPACK-LTPACK Disambiguation)</h2>

<h3 id="71-nomenclature-disambiguation-itpack-vs-ltpack">7.1 Nomenclature Disambiguation: ITPACK vs. "ltpack"</h3>
<div class="exam-highlight-box" style="border: 1px solid var(--primary-border); background: var(--primary-subtle); padding: 1rem 1.25rem; border-radius: 8px; margin: 1rem 0;">
  <strong>Pedagogical &amp; Typographical Disambiguation:</strong>
  In academic courses, examination papers, and lecture notes, students frequently encounter references to both <strong>"ELLPACK-ITPACK"</strong> and <strong>"Ellpack-ltpack"</strong>.
  <br/><br/>
  <strong>They are the EXACT SAME format!</strong>
  The format was developed as the storage standard for the famous <strong>ITPACK</strong> (Iterative Package) numerical linear algebra software library developed at the Center for Numerical Analysis at the University of Texas at Austin (David Kincaid, Thomas Grimes, David Young, and John Rice).
  <br/>
  In standard sans-serif computer typography (such as Arial, Helvetica, or Calibri), the uppercase letter <strong>'I'</strong> and the lowercase letter <strong>'l'</strong> are visually indistinguishable. Optical character recognition (OCR) software and handwritten lecture notes routinely transcribe "ITPACK" as "ltpack".
</div>

<h3 id="72-the-itpack-innovation-active-row-lengths">7.2 The ITPACK Innovation: Active Row Lengths</h3>
<p>
While standard ELLPACK achieves glorious memory coalescing on GPUs, it has a glaring inefficiency:
<strong>Every thread must perform $K_{\max}$ multiply-adds, even if its row only has 1 or 2 non-zeros!</strong>
In our Matrix $A$ example, Row 4 has only 1 non-zero ($A_{44}=12.0$), but standard ELLPACK forces it to execute 4 multiplications:
$$y_4 = 12.0 \times x_4 + 0.0 \times x_4 + 0.0 \times x_4 + 0.0 \times x_4$$
$75\%$ of the arithmetic operations on Row 4 are useless waste!
</p>

<p>
<strong>The ELLPACK-ITPACK Solution:</strong>
Introduce a third array, <code>len</code> (or <code>row_len</code> / <code>JLEN</code>), of size $N$:
$$\text{len}[i] = \text{number of actual non-zero elements in row } i \quad (0 \le \text{len}[i] \le K_{\max})$$
</p>

<p>The SpMV inner loop for row $i$ now terminates early at $\text{len}[i]$:</p>
<pre><code class="language-c">// ELLPACK-ITPACK SpMV with Early Exit
for (int c = 0; c &lt; len[i]; c++) {
    sum += COEF[c * N + i] * x[JCOEF[c * N + i]];
}</code></pre>

<h3 id="73-worked-example-converting-canonical-matrix-a-to-ellpack-itpack">7.3 Worked Example: Converting Canonical Matrix A to ELLPACK-ITPACK</h3>
<p>For Canonical Matrix $A$ ($N = 5, K_{\max} = 4$):</p>
<ul>
  <li><code>COEF</code> and <code>JCOEF</code> are identical to standard ELLPACK (dimension $5 \times 4$).</li>
  <li>The active length array is:
  $$\mathbf{\text{len} = [2, \; 3, \; 4, \; 2, \; 1]}$$
  </li>
</ul>

<p><strong>Quantitative Arithmetic Comparison (SpMV on Matrix A):</strong></p>
<ul>
  <li><strong>Dense Matrix-Vector Multiply:</strong> $N^2 = 5 \times 5 = \mathbf{25\text{ multiplications}}$</li>
  <li><strong>Standard ELLPACK SpMV:</strong> $N \times K_{\max} = 5 \times 4 = \mathbf{20\text{ multiplications}}$</li>
  <li><strong>ELLPACK-ITPACK SpMV:</strong> $\sum_{i=0}^4 \text{len}[i] = 2 + 3 + 4 + 2 + 1 = \mathbf{12\text{ multiplications}}$</li>
  <li><strong>FLOP Reduction:</strong> ELLPACK-ITPACK eliminates $\frac{20 - 12}{20} = \mathbf{40\%}$ of all arithmetic operations compared to standard ELLPACK!</li>
</ul>

<h3 id="74-the-permutation-enhancement-sorted-ellpack-itpack">7.4 The Permutation Enhancement: Sorted ELLPACK-ITPACK</h3>
<p>
On GPUs, early-exit loops can reintroduce warp divergence if adjacent threads have different <code>len</code> values.
To achieve the absolute pinnacle of GPU performance, ELLPACK-ITPACK introduces an optional <strong>permutation vector</strong> <code>perm[N]</code>:
</p>
<ol>
  <li>Sort all rows of matrix $A$ in <strong>descending order of non-zero count</strong>.</li>
  <li>Group rows of equal length into contiguous chunks of 32 rows (matching the GPU warp size).</li>
  <li>All 32 threads in each warp now execute the <strong>exact same loop bound</strong> $\text{len}[i]$, achieving $100\%$ warp occupancy, zero divergence, and zero wasted multiplications!</li>
</ol>

<h3 id="75-modern-gpu-descendants-sell-c-sigma-and-hyb">7.5 Modern GPU Descendants: SELL-C-σ and HYB</h3>
<ul>
  <li><strong>Sliced ELLPACK (SELL-C-$\sigma$):</strong> Breaks the matrix into horizontal slices of $C$ rows (e.g. $C = 32$ for a GPU warp). Each slice is padded only to the maximum row length <em>within that slice</em>, preventing an isolated outlier row from inflating the entire matrix.</li>
  <li><strong>Hybrid Format (HYB = ELL + COO):</strong> The default format in <strong>NVIDIA cuSPARSE</strong>. The regular portion of the matrix (up to a small threshold, e.g. $K=5$) is stored in ELLPACK for blazing GPU coalescing, while outlier entries beyond $K$ are spilled into an auxiliary COO buffer.</li>
</ul>

</section>

<!-- SECTION 8 -->
<section id="ch20-sec8">
<h2 id="8-modified-sparse-row-msr-format">8. Modified Sparse Row (MSR) Format</h2>

<h3 id="81-concept--the-diagonal-dominance-exploit">8.1 Concept &amp; The Diagonal Dominance Exploit</h3>
<p>
The <strong>Modified Sparse Row (MSR)</strong> format, developed by Yousef Saad, is a clever evolution of CSR designed specifically for stationary iterative solvers (Jacobi, Gauss-Seidel, SOR, SSOR).
</p>
<p>
<strong>The Physical Insight:</strong> In virtually every PDE discretization, the main diagonal entries $A_{ii}$ are <strong>strictly non-zero</strong> and are accessed on every single iteration step (e.g. in Jacobi relaxation: $x_i^{(k+1)} = \frac{1}{A_{ii}}(b_i - \dots)$).
MSR stores all $N$ main diagonal elements contiguously in the front of a single array, reducing memory overhead to <strong>exactly two arrays</strong> (eliminating the third array entirely!).
</p>

<p>Both arrays have length exactly equal to $N_{nz} + 1$:</p>
<ol>
  <li><code>AA</code> (Floating-point values, length $N_{nz} + 1$):
    <ul>
      <li>Positions $0 \dots N-1$: Main diagonal entries in order ($A_{00}, A_{11}, \dots, A_{N-1, N-1}$).</li>
      <li>Position $N$: Unused dummy placeholder (conventionally marked with <code>*</code>).</li>
      <li>Positions $N+1 \dots N_{nz}$: Strictly <strong>off-diagonal</strong> non-zeros stored row-by-row.</li>
    </ul>
  </li>
  <li><code>JA</code> (Integer indices &amp; pointers, length $N_{nz} + 1$):
    <ul>
      <li>Positions $0 \dots N$: Pointers marking where each row's off-diagonals start in <code>AA</code>.</li>
      <li>Positions $N+1 \dots N_{nz}$: Column indices of the off-diagonal elements in <code>AA</code>.</li>
    </ul>
  </li>
</ol>

<h3 id="82-worked-example-converting-canonical-matrix-a-to-msr">8.2 Worked Example: Converting Canonical Matrix A to MSR</h3>
<p>Recall Canonical Matrix $A$ ($N = 5, N_{nz} = 12$):</p>
$$A = \begin{bmatrix}
1.0 & 0.0 & 0.0 & 2.0 & 0.0 \\
3.0 & 4.0 & 0.0 & 5.0 & 0.0 \\
6.0 & 0.0 & 7.0 & 8.0 & 9.0 \\
0.0 & 0.0 & 10.0 & 11.0 & 0.0 \\
0.0 & 0.0 & 0.0 & 0.0 & 12.0
\end{bmatrix}$$

<p><strong>Array Length:</strong> $N_{nz} + 1 = 12 + 1 = 13$ elements.</p>
<ul>
  <li>Main diagonal entries ($A_{00} \dots A_{44}$): $[1.0, 4.0, 7.0, 11.0, 12.0]$.</li>
  <li>Off-diagonal entries row-by-row:
    <ul>
      <li>Row 0 off-diagonals: $A(0, 3) = 2.0$ (col 3)</li>
      <li>Row 1 off-diagonals: $A(1, 0) = 3.0$ (col 0), $A(1, 3) = 5.0$ (col 3)</li>
      <li>Row 2 off-diagonals: $A(2, 0) = 6.0$ (col 0), $A(2, 3) = 8.0$ (col 3), $A(2, 4) = 9.0$ (col 4)</li>
      <li>Row 3 off-diagonals: $A(3, 2) = 10.0$ (col 2)</li>
      <li>Row 4 off-diagonals: <strong>NONE!</strong> (Row 4 has only diagonal entry $12.0$)</li>
    </ul>
  </li>
</ul>

<p><strong>Complete MSR Arrays (0-Indexed, Size 13):</strong></p>
<pre><code>Index:   0     1     2     3     4     5   │   6     7     8     9    10    11    12
───────────────────────────────────────────┼─────────────────────────────────────────
AA:   [ 1.0   4.0   7.0  11.0  12.0    *   │  2.0   3.0   5.0   6.0   8.0   9.0  10.0 ]
        └──────────────┬───────────────┘   │  └──────────────────┬──────────────────┘
             Main Diagonal Entries         │             Off-Diagonal Entries
                                           │
JA:   [   6     7     9    12    13    13  │    3     0     3     0     3     4     2  ]
        └──────────────┬───────────────┘   │  └──────────────────┬──────────────────┘
           Row Off-Diagonal Pointers       │          Off-Diagonal Column Indices</code></pre>

<div class="cram-box walkthrough-box" style="border-left: 4px solid #10b981; background: var(--bg-alt); padding: 1.25rem; border-radius: 0 8px 8px 0; margin: 1.5rem 0;">
<h4 style="margin-top: 0; color: #10b981;">✍️ High-Yield Exam Trap: Why JA[4] == JA[5] == 13?</h4>
<p>
Look closely at indices 4 and 5 of array <code>JA</code>:
$$\text{JA}[4] = 13, \quad \text{JA}[5] = 13$$
<br/>
The number of off-diagonal entries in Row 4 is:
$$\text{count} = \text{JA}[5] - \text{JA}[4] = 13 - 13 = \mathbf{0}$$
<strong>Crucial Exam Takeaway:</strong> Whenever two consecutive pointer entries in <code>JA</code> are equal ($\text{JA}[i] == \text{JA}[i+1]$), it signifies that <strong>Row $i$ possesses ZERO off-diagonal elements</strong>!
</p>
</div>

</section>

<!-- SECTION 9 -->
<section id="ch20-sec9">
<h2 id="9-comprehensive-comparative-analysis--memory-engineering">9. Comprehensive Comparative Analysis &amp; Memory Engineering</h2>

<h3 id="91-analytical-memory-footprint-formulas">9.1 Analytical Memory Footprint Formulas</h3>
<p>
Let $N$ be the number of rows, $M$ the number of columns, and $N_{nz}$ the non-zero count.
Assuming standard IEEE 754 double-precision floating-point (<code>sizeof(double) = 8 bytes</code>) and 32-bit signed integers (<code>sizeof(int) = 4 bytes</code>):
</p>

<table>
<thead>
<tr>
  <th style="text-align: left;">Format</th>
  <th style="text-align: left;">Array Breakdown &amp; Sizes</th>
  <th style="text-align: left;">Total Memory Formula (Bytes)</th>
  <th style="text-align: left;">Memory Overhead per Non-Zero</th>
</tr>
</thead>
<tbody>
<tr>
  <td><strong>Dense</strong></td>
  <td>$1 \times [N \times M]$ double</td>
  <td>$\mathbf{8 \cdot N \cdot M}$</td>
  <td>$\infty$ (stores all zeroes)</td>
</tr>
<tr>
  <td><strong>COO</strong></td>
  <td>$1 \times [N_{nz}]$ double, $2 \times [N_{nz}]$ int</td>
  <td>$\mathbf{16 \cdot N_{nz}}$</td>
  <td>$16\text{ bytes / non-zero}$</td>
</tr>
<tr>
  <td><strong>CSR</strong></td>
  <td>$1 \times [N_{nz}]$ double, $1 \times [N_{nz}]$ int, $1 \times [N+1]$ int</td>
  <td>$\mathbf{12 \cdot N_{nz} + 4N + 4}$</td>
  <td>$\approx 12\text{ bytes} + \frac{4N}{N_{nz}}$</td>
</tr>
<tr>
  <td><strong>CSC</strong></td>
  <td>$1 \times [N_{nz}]$ double, $1 \times [N_{nz}]$ int, $1 \times [M+1]$ int</td>
  <td>$\mathbf{12 \cdot N_{nz} + 4M + 4}$</td>
  <td>$\approx 12\text{ bytes} + \frac{4M}{N_{nz}}$</td>
</tr>
<tr>
  <td><strong>MSR</strong></td>
  <td>$1 \times [N_{nz}+1]$ double, $1 \times [N_{nz}+1]$ int</td>
  <td>$\mathbf{12 \cdot N_{nz} + 12}$</td>
  <td>$\approx 12\text{ bytes / non-zero}$</td>
</tr>
<tr>
  <td><strong>DIA</strong></td>
  <td>$1 \times [N \times N_{\text{diag}}]$ double, $1 \times [N_{\text{diag}}]$ int</td>
  <td>$\mathbf{8 \cdot N \cdot N_{\text{diag}} + 4 \cdot N_{\text{diag}}}$</td>
  <td>$8 \cdot \frac{N \cdot N_{\text{diag}}}{N_{nz}}$ (zero index overhead!)</td>
</tr>
<tr>
  <td><strong>ELLPACK</strong></td>
  <td>$1 \times [N \times K_{\max}]$ double, $1 \times [N \times K_{\max}]$ int</td>
  <td>$\mathbf{12 \cdot N \cdot K_{\max}}$</td>
  <td>$12 \cdot \frac{N \cdot K_{\max}}{N_{nz}}$</td>
</tr>
<tr>
  <td><strong>ELLPACK-ITPACK</strong></td>
  <td>$1 \times [N \times K_{\max}]$ double, $1 \times [N \times K_{\max}]$ int, $1 \times [N]$ int</td>
  <td>$\mathbf{12 \cdot N \cdot K_{\max} + 4N}$</td>
  <td>$12 \cdot \frac{N \cdot K_{\max}}{N_{nz}} + \frac{4N}{N_{nz}}$</td>
</tr>
</tbody>
</table>

<h3 id="92-large-scale-hpc-case-studies">9.2 Large-Scale HPC Case Studies ($N = 10^6$ Equations)</h3>
<p>Let us compute the exact memory footprint in Megabytes (MB) across three distinct scientific workloads:</p>

<table>
<thead>
<tr>
  <th style="text-align: left;">Matrix Storage Format</th>
  <th style="text-align: center;">Case A: 2D Poisson (5-pt Stencil)<br/>$N=10^6, N_{nz}=5\times 10^6$<br/>$K_{\max}=5, N_{\text{diag}}=5$</th>
  <th style="text-align: center;">Case B: 3D Poisson (7-pt Stencil)<br/>$N=10^6, N_{nz}=7\times 10^6$<br/>$K_{\max}=7, N_{\text{diag}}=7$</th>
  <th style="text-align: center;">Case C: Unstructured Aerodynamics<br/>$N=10^6, N_{nz}=7\times 10^6$<br/>$K_{\max}=120 \text{ (1 hub node!)}, N_{\text{diag}}=800$</th>
</tr>
</thead>
<tbody>
<tr>
  <td><strong>Dense ($N \times N$)</strong></td>
  <td style="text-align: center; color: #ef4444;">8,000,000 MB (8 TB)</td>
  <td style="text-align: center; color: #ef4444;">8,000,000 MB (8 TB)</td>
  <td style="text-align: center; color: #ef4444;">8,000,000 MB (8 TB)</td>
</tr>
<tr>
  <td><strong>COO</strong></td>
  <td style="text-align: center;">80.0 MB</td>
  <td style="text-align: center;">112.0 MB</td>
  <td style="text-align: center;">112.0 MB</td>
</tr>
<tr>
  <td><strong>CSR</strong></td>
  <td style="text-align: center;"><strong>64.0 MB</strong></td>
  <td style="text-align: center;"><strong>88.0 MB</strong></td>
  <td style="text-align: center;"><strong style="color: #10b981;">88.0 MB (Optimal!)</strong></td>
</tr>
<tr>
  <td><strong>CSC</strong></td>
  <td style="text-align: center;">64.0 MB</td>
  <td style="text-align: center;">88.0 MB</td>
  <td style="text-align: center;">88.0 MB</td>
</tr>
<tr>
  <td><strong>DIA</strong></td>
  <td style="text-align: center;"><strong style="color: #10b981;">40.0 MB (Minimal RAM!)</strong></td>
  <td style="text-align: center;"><strong style="color: #10b981;">56.0 MB (Minimal RAM!)</strong></td>
  <td style="text-align: center; color: #ef4444;">6,400.0 MB (Disaster!)</td>
</tr>
<tr>
  <td><strong>ELLPACK</strong></td>
  <td style="text-align: center;">60.0 MB</td>
  <td style="text-align: center;">84.0 MB</td>
  <td style="text-align: center; color: #ef4444;">1,440.0 MB (95% zero padding!)</td>
</tr>
<tr>
  <td><strong>ELLPACK-ITPACK</strong></td>
  <td style="text-align: center;">64.0 MB</td>
  <td style="text-align: center;">88.0 MB</td>
  <td style="text-align: center;">1,444.0 MB</td>
</tr>
</tbody>
</table>

<div class="cram-box intuition-box" style="border-left: 4px solid var(--primary); background: var(--bg-alt); padding: 1.25rem; border-radius: 0 8px 8px 0; margin: 1.5rem 0;">
<h4 style="margin-top: 0; color: var(--primary);">📊 Crucial Engineering Insights from the Benchmark:</h4>
<ol>
  <li><strong>Structured Grids (Cases A &amp; B):</strong> DIA requires the least memory (only $56\text{ MB}$ for 7 million non-zeros) because it stores <em>zero</em> column index arrays! ELLPACK is a close second ($84\text{ MB}$) and delivers the fastest GPU SpMV.</li>
  <li><strong>Unstructured Grids with Hub Nodes (Case C):</strong> A single node connected to 120 neighbors inflates ELLPACK's $K_{\max}$ to 120, exploding memory by $16\times$ ($1.44\text{ GB}$). DIA suffers an even worse catastrophe ($6.4\text{ GB}$).
  In contrast, <strong>CSR remains rock-solid at 88 MB</strong>, completely immune to variations in row lengths! This is why CSR is the king of general unstructured finite-element packages.</li>
</ol>
</div>

<h3 id="93-the-master-sparse-matrix-rosetta-stone">9.3 The Master Sparse Matrix Rosetta Stone</h3>
<table>
<thead>
<tr>
  <th style="text-align: left;">Format</th>
  <th style="text-align: left;">Primary Data Arrays</th>
  <th style="text-align: center;">Total Storage Footprint</th>
  <th style="text-align: center;">Random Element Lookup $A_{ij}$</th>
  <th style="text-align: left;">SpMV Parallelization Strategy</th>
  <th style="text-align: center;">GPU Suitability</th>
  <th style="text-align: left;">Primary Practical Domain</th>
</tr>
</thead>
<tbody>
<tr>
  <td><strong>COO</strong></td>
  <td><code>val</code>, <code>row</code>, <code>col</code></td>
  <td style="text-align: center;">$16 N_{nz}$ B</td>
  <td style="text-align: center;">$O(N_{nz})$</td>
  <td>Poor (atomic write races on $y$)</td>
  <td style="text-align: center;">★☆☆☆☆</td>
  <td>FEM stiffness assembly, Matrix Market I/O</td>
</tr>
<tr>
  <td><strong>CSR</strong></td>
  <td><code>val</code>, <code>col_ind</code>, <code>row_ptr</code></td>
  <td style="text-align: center;">$12 N_{nz} + 4N$ B</td>
  <td style="text-align: center;">$O(\log(\text{nnz}_i))$</td>
  <td>Excellent (1 thread per row, no races)</td>
  <td style="text-align: center;">★★★☆☆</td>
  <td>General sparse linear solvers, Multicore CPUs</td>
</tr>
<tr>
  <td><strong>CSC</strong></td>
  <td><code>val</code>, <code>row_ind</code>, <code>col_ptr</code></td>
  <td style="text-align: center;">$12 N_{nz} + 4M$ B</td>
  <td style="text-align: center;">$O(\log(\text{nnz}_j))$</td>
  <td>Poor for $y=Ax$; Great for $y=A^Tx$</td>
  <td style="text-align: center;">★★☆☆☆</td>
  <td>Direct factorization (SuperLU, MATLAB)</td>
</tr>
<tr>
  <td><strong>DIA</strong></td>
  <td><code>data</code>, <code>offset</code></td>
  <td style="text-align: center;">$8 N \cdot N_{\text{diag}}$ B</td>
  <td style="text-align: center;">$O(1)$</td>
  <td>Optimal (stride-1 streaming, no indirection)</td>
  <td style="text-align: center;">★★★★☆</td>
  <td>Structured finite-difference stencils</td>
</tr>
<tr>
  <td><strong>ELLPACK</strong></td>
  <td><code>COEF</code>, <code>JCOEF</code></td>
  <td style="text-align: center;">$12 N \cdot K_{\max}$ B</td>
  <td style="text-align: center;">$O(K_{\max})$</td>
  <td>Optimal (100% coalesced column-major GPU loads)</td>
  <td style="text-align: center;">★★★★★</td>
  <td>GPUs, CUDA SIMT, bounded degree meshes</td>
</tr>
<tr>
  <td><strong>ELLPACK-ITPACK</strong></td>
  <td><code>COEF</code>, <code>JCOEF</code>, <code>len</code></td>
  <td style="text-align: center;">$12 N \cdot K_{\max} + 4N$ B</td>
  <td style="text-align: center;">$O(\text{len}_i)$</td>
  <td>Optimal + eliminates zero-padding flops</td>
  <td style="text-align: center;">★★★★★</td>
  <td>GPU vector computing, sorted row kernels</td>
</tr>
<tr>
  <td><strong>MSR</strong></td>
  <td><code>AA</code>, <code>JA</code> (2 arrays)</td>
  <td style="text-align: center;">$12 N_{nz} + 12$ B</td>
  <td style="text-align: center;">$O(1)$ diag, $O(\log k)$ off</td>
  <td>High (direct diagonal inversion)</td>
  <td style="text-align: center;">★★★☆☆</td>
  <td>Stationary relaxation (Jacobi, GS, SOR)</td>
</tr>
</tbody>
</table>

<h3 id="94-the-master-format-selection-decision-tree">9.4 The Master Format Selection Decision Tree</h3>
<pre><code>                     [What is your primary computational task?]
                                      │
           ┌──────────────────────────┴──────────────────────────┐
           ▼                                                     ▼
   [Matrix Assembly / I/O]                            [Iterative Linear Solver]
           │                                                     │
     Use ➔ COO                                  [What is your matrix structure?]
                                                                 │
                                ┌────────────────────────────────┴────────────────────────────────┐
                                ▼                                                                 ▼
                      [Structured Grid PDE]                                            [Unstructured Mesh / Graph]
                                │                                                                 │
                   [Are diagonals constant?]                                           [Target Hardware?]
                                │                                                                 │
                     ┌──────────┴──────────┐                                      ┌───────────────┴───────────────┐
                     ▼                     ▼                                      ▼                               ▼
                 {YES}                    {NO}                             [Multicore CPU]                  [NVIDIA GPU]
                   │                       │                                      │                               │
               Use ➔ DIA          [Are row counts bounded?]                   Use ➔ CSR                 [Is max degree bounded?]
                                           │                                                                      │
                                ┌──────────┴──────────┐                                                ┌──────────┴──────────┐
                                ▼                     ▼                                                ▼                     ▼
                              {YES}                  {NO}                                            {YES}                  {NO}
                                │                     │                                                │                     │
                        Use ➔ ELLPACK /          Use ➔ CSR                                      Use ➔ ELLPACK-        Use ➔ HYB
                              ELLPACK-ITPACK                                                          ITPACK          (cuSPARSE)</code></pre>

</section>

<!-- SECTION 10 -->
<section id="ch20-sec10">
<h2 id="10-production-grade-conversion-algorithms">10. Production-Grade Conversion Algorithms</h2>

<h3 id="101-coo-to-csr-conversion-the-prefix-sum-algorithm">10.1 COO to CSR Conversion: The Prefix-Sum Algorithm</h3>
<p>
Converting COO to CSR requires transforming unordered $(i, j, v)$ triplets into sorted, compressed row slices.
The standard high-performance algorithm executes in <strong>three linear-time passes ($O(N_{nz} + N)$)</strong> without full sorting:
</p>

<pre><code class="language-c">// Production-Grade C99 COO to CSR Converter
#include &lt;stdlib.h&gt;
#include &lt;string.h&gt;

void coo_to_csr(int n_rows, int nnz, 
                const int *coo_row, const int *coo_col, const double *coo_val,
                int *csr_row_ptr, int *csr_col_ind, double *csr_val) {
    // Step 1: Count non-zeros in each row (Histogram pass)
    memset(csr_row_ptr, 0, (n_rows + 1) * sizeof(int));
    for (int k = 0; k &lt; nnz; k++) {
        csr_row_ptr[coo_row[k] + 1]++;
    }
    
    // Step 2: Cumulative prefix sum to establish row pointers
    for (int i = 0; i &lt; n_rows; i++) {
        csr_row_ptr[i + 1] += csr_row_ptr[i];
    }
    
    // Step 3: Temporary tracking array for current row insert positions
    int *head = (int *)malloc(n_rows * sizeof(int));
    memcpy(head, csr_row_ptr, n_rows * sizeof(int));
    
    // Step 4: Scatter entries into their compressed slots
    for (int k = 0; k &lt; nnz; k++) {
        int r = coo_row[k];
        int dest = head[r]++;
        csr_val[dest] = coo_val[k];
        csr_col_ind[dest] = coo_col[k];
    }
    
    free(head);
}</code></pre>

<h3 id="102-csr-to-ellpack-and-ellpack-itpack-conversion">10.2 CSR to ELLPACK and ELLPACK-ITPACK Conversion</h3>
<pre><code class="language-c">// CSR to ELLPACK &amp; ELLPACK-ITPACK Converter (Column-Major GPU Layout)
void csr_to_ellpack(int n_rows, const int *row_ptr, const int *col_ind, const double *val,
                    int *K_max, double **COEF, int **JCOEF, int **row_len) {
    // 1. Determine maximum row length
    *K_max = 0;
    *row_len = (int *)malloc(n_rows * sizeof(int));
    for (int i = 0; i &lt; n_rows; i++) {
        int len = row_ptr[i + 1] - row_ptr[i];
        (*row_len)[i] = len;
        if (len &gt; *K_max) *K_max = len;
    }
    
    // 2. Allocate column-major 2D arrays: dimension [K_max * n_rows]
    int total_slots = (*K_max) * n_rows;
    *COEF  = (double *)calloc(total_slots, sizeof(double));
    *JCOEF = (int *)malloc(total_slots * sizeof(int));
    
    // 3. Populate arrays with padding
    for (int i = 0; i &lt; n_rows; i++) {
        int r_start = row_ptr[i];
        int r_len   = (*row_len)[i];
        
        for (int c = 0; c &lt; *K_max; c++) {
            int slot = c * n_rows + i; // Column-major index!
            if (c &lt; r_len) {
                (*COEF)[slot]  = val[r_start + c];
                (*JCOEF)[slot] = col_ind[r_start + c];
            } else {
                (*COEF)[slot]  = 0.0;
                // Pad with self-row index to guarantee safe in-bounds read
                (*JCOEF)[slot] = i; 
            }
        }
    }
}</code></pre>

</section>

<!-- SECTION 11 -->
<section id="ch20-sec11">
<h2 id="11-high-yield-examination-problems--model-solutions">11. High-Yield Examination Problems &amp; Model Solutions</h2>

<!-- Question 1 -->
<div class="cram-box walkthrough-box" style="border-left: 4px solid var(--primary); background: var(--bg-alt); padding: 1.5rem; border-radius: 0 8px 8px 0; margin: 2rem 0;">
<h3 id="exam-q1" style="margin-top: 0; color: var(--primary);">Exam Problem 1: Full Multi-Format Conversion of a $6 \times 6$ Matrix [8 Marks]</h3>
<p><strong>Problem Statement:</strong><br/>
A sparse linear system arising from an electrical network simulation has the following coefficient matrix $A \in \mathbb{R}^{6 \times 6}$ with $N_{nz} = 11$:</p>
$$A = \begin{bmatrix}
5 & 0 & 0 & 1 & 0 & 0 \\
0 & 8 & 2 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
3 & 0 & 0 & 9 & 0 & 4 \\
0 & 6 & 0 & 0 & 7 & 0 \\
0 & 0 & 0 & 0 & 0 & 2
\end{bmatrix}$$
<ol>
  <li>Write down the 0-indexed <strong>CSR</strong> arrays: <code>val</code>, <code>col_ind</code>, and <code>row_ptr</code>.</li>
  <li>Write down the 0-indexed <strong>CSC</strong> arrays: <code>val</code>, <code>row_ind</code>, and <code>col_ptr</code>.</li>
  <li>State the value of $K_{\max}$ and construct the <strong>ELLPACK</strong> arrays <code>COEF</code> and <code>JCOEF</code>.</li>
  <li>Construct the <strong>ELLPACK-ITPACK</strong> array <code>len</code>.</li>
  <li>Can this matrix be efficiently represented in the <strong>DIA</strong> format? Justify mathematically using diagonal offsets.</li>
</ol>

<h4 style="color: #10b981;">Complete Model Solution:</h4>
<ol>
  <li><strong>CSR Format (0-Indexed):</strong>
  <ul>
    <li>Row 0: $A(0,0)=5, A(0,3)=1$ (2 entries)</li>
    <li>Row 1: $A(1,1)=8, A(1,2)=2$ (2 entries)</li>
    <li>Row 2: <strong>EMPTY ROW!</strong> (0 entries)</li>
    <li>Row 3: $A(3,0)=3, A(3,3)=9, A(3,5)=4$ (3 entries)</li>
    <li>Row 4: $A(4,1)=6, A(4,4)=7$ (2 entries)</li>
    <li>Row 5: $A(5,5)=2$ (1 entry)</li>
  </ul>
  <pre><code>val     = [ 5, 1,  8, 2,  3, 9, 4,  6, 7,  2 ]           (length 10 ... wait, let's recount entries:
            (0,0)=5, (0,3)=1 [2]; (1,1)=8, (1,2)=2 [2]; (3,0)=3, (3,3)=9, (3,5)=4 [3]; 
            (4,1)=6, (4,4)=7 [2]; (5,5)=2 [1] -> Total N_nz = 10 entries)
val     = [ 5, 1, 8, 2, 3, 9, 4, 6, 7, 2 ]
col_ind = [ 0, 3, 1, 2, 0, 3, 5, 1, 4, 5 ]
row_ptr = [ 0, 2, 4, 4, 7, 9, 10 ]  &lt;-- Notice row_ptr[2] == row_ptr[3] == 4!</code></pre>
  </li>

  <li><strong>CSC Format (0-Indexed):</strong>
  Scan column by column:
  <ul>
    <li>Col 0: $A(0,0)=5, A(3,0)=3$</li>
    <li>Col 1: $A(1,1)=8, A(4,1)=6$</li>
    <li>Col 2: $A(1,2)=2$</li>
    <li>Col 3: $A(0,3)=1, A(3,3)=9$</li>
    <li>Col 4: $A(4,4)=7$</li>
    <li>Col 5: $A(3,5)=4, A(5,5)=2$</li>
  </ul>
  <pre><code>val     = [ 5, 3,  8, 6,  2,  1, 9,  7,  4, 2 ]
row_ind = [ 0, 3,  1, 4,  2,  0, 3,  4,  3, 5 ]
col_ptr = [ 0, 2,  4,     5,  7,     8, 10 ]</code></pre>
  </li>

  <li><strong>ELLPACK Format ($K_{\max} = \max(2, 2, 0, 3, 2, 1) = 3$):</strong>
  Dimension is $6 \times 3$:
  $$\text{COEF} = \begin{bmatrix}
  5 & 1 & 0 \\
  8 & 2 & 0 \\
  0 & 0 & 0 \\
  3 & 9 & 4 \\
  6 & 7 & 0 \\
  2 & 0 & 0
  \end{bmatrix}, \qquad
  \text{JCOEF} = \begin{bmatrix}
  0 & 3 & 0 \\
  1 & 2 & 1 \\
  0 & 0 & 0 \\
  0 & 3 & 5 \\
  1 & 4 & 1 \\
  5 & 5 & 5
  \end{bmatrix}$$
  </li>

  <li><strong>ELLPACK-ITPACK Active Length Array:</strong>
  $$\mathbf{\text{len} = [2, \; 2, \; 0, \; 3, \; 2, \; 1]}$$
  Notice that Row 2 has $\text{len}[2] = 0$, completely skipping all computation during SpMV!
  </li>

  <li><strong>DIA Evaluation:</strong>
  Let us compute the diagonal offsets $k = j - i$ for all entries:
  <ul>
    <li>$A(0,0) \implies k = 0$</li>
    <li>$A(0,3) \implies k = +3$</li>
    <li>$A(1,1) \implies k = 0$</li>
    <li>$A(1,2) \implies k = +1$</li>
    <li>$A(3,0) \implies k = -3$</li>
    <li>$A(3,3) \implies k = 0$</li>
    <li>$A(3,5) \implies k = +2$</li>
    <li>$A(4,1) \implies k = -3$</li>
    <li>$A(4,4) \implies k = 0$</li>
    <li>$A(5,5) \implies k = 0$</li>
  </ul>
  The active offsets are: $\{-3, 0, +1, +2, +3\} \implies N_{\text{diag}} = 5$.
  <br/>
  <strong>Conclusion:</strong> For an $N = 6$ matrix with only 10 non-zeros, storing 5 full diagonals requires $6 \times 5 = 30$ words! The density in DIA is $\frac{10}{30} \approx 33\%$, meaning <strong>$67\%$ of the storage is wasted on padded zeroes</strong>. DIA is highly inefficient for this matrix; CSR is vastly superior.
  </li>
</ol>
</div>

<!-- Question 2 -->
<div class="cram-box walkthrough-box" style="border-left: 4px solid #8b5cf6; background: var(--bg-alt); padding: 1.5rem; border-radius: 0 8px 8px 0; margin: 2rem 0;">
<h3 id="exam-q2" style="margin-top: 0; color: #8b5cf6;">Exam Problem 2: GPU Warp Divergence &amp; Bandwidth Analysis [6 Marks]</h3>
<p><strong>Problem Statement:</strong><br/>
An aerospace CFD code executes SpMV on an NVIDIA GPU (warp size $= 32$, DRAM bandwidth $= 900\text{ GB/s}$).
Matrix $A$ has $N = 10^6$ rows.
Consider Warp 0, which processes Rows $0 \dots 31$.
Suppose 31 of the rows have exactly 4 non-zeros, but Row 0 has 64 non-zeros.</p>
<ol>
  <li>Explain the execution behavior of Warp 0 under a 1-thread-per-row <strong>CSR</strong> kernel. What is the warp execution efficiency?</li>
  <li>Explain the execution behavior under a column-major <strong>ELLPACK</strong> kernel.</li>
  <li>How does <strong>ELLPACK-ITPACK with row sorting</strong> resolve this inefficiency?</li>
</ol>

<h4 style="color: #10b981;">Complete Model Solution:</h4>
<ol>
  <li><strong>CSR Kernel Behavior:</strong>
  In CUDA, all 32 threads in a warp share a single instruction decoder and execute in lockstep SIMT.
  Thread 0 requires 64 loop iterations. Threads $1 \dots 31$ require only 4 loop iterations.
  After 4 iterations, threads $1 \dots 31$ become inactive (diverged) and must <strong>idle for the remaining 60 iterations</strong>!
  $$\text{Warp Execution Efficiency} = \frac{\text{Useful Work}}{\text{Total Allocated Thread Cycles}} = \frac{64 + 31 \times 4}{32 \times 64} = \frac{64 + 124}{2048} = \frac{188}{2048} \approx \mathbf{9.18\%}$$
  More than $\mathbf{90.8\%}$ of the compute capability of Warp 0 is utterly wasted!
  </li>
  <li><strong>ELLPACK Kernel Behavior:</strong>
  In standard ELLPACK, $K_{\max} = 64$ across the warp. All 32 threads execute 64 iterations unconditionally. While there is no branch divergence at the hardware level, threads $1 \dots 31$ execute 60 redundant multiplications by $0.0$, burning arithmetic energy.
  </li>
  <li><strong>Remediation via ELLPACK-ITPACK + Permutation:</strong>
  By applying a permutation vector <code>perm</code> that sorts rows by non-zero count:
  Row 0 (length 64) is grouped with other dense rows into a high-density warp.
  The 31 rows of length 4 are grouped with other rows of length 4.
  In the length-4 warp, $K_{\max} = 4$, and the warp terminates after exactly 4 cycles with <strong>100% warp efficiency and zero idle threads</strong>!
  </li>
</ol>
</div>

<!-- Question 3 -->
<div class="cram-box walkthrough-box" style="border-left: 4px solid #f59e0b; background: var(--bg-alt); padding: 1.5rem; border-radius: 0 8px 8px 0; margin: 2rem 0;">
<h3 id="exam-q3" style="margin-top: 0; color: #f59e0b;">Exam Problem 3: Decoding MSR Arrays and Hand SpMV [5 Marks]</h3>
<p><strong>Problem Statement:</strong><br/>
A $4 \times 4$ linear system is stored in 0-indexed Modified Sparse Row (MSR) format:</p>
<pre><code>AA = [ 20.0,  30.0,  40.0,  50.0,    *,  -2.0,   1.0,  -4.0,   5.0 ]
JA = [    5,     6,     8,     8,    9,     1,     0,     3,     2 ]</code></pre>
<ol>
  <li>Reconstruct the full $4 \times 4$ matrix $A$.</li>
  <li>Given vector $x = [1, 2, 3, 4]^T$, compute the matrix-vector product $y = A x$ step-by-step using only the MSR arrays.</li>
</ol>

<h4 style="color: #10b981;">Complete Model Solution:</h4>
<ol>
  <li><strong>Matrix Reconstruction:</strong>
  $N = 4$.
  <ul>
    <li>Diagonal elements: $A_{00}=20.0, A_{11}=30.0, A_{22}=40.0, A_{33}=50.0$.</li>
    <li>Row 0: Off-diagonals start at $\text{JA}[0]=5$ up to $\text{JA}[1]-1=5$ (1 element).
    <br/>Index 5: $\text{val}=-2.0, \text{col}=\text{JA}[5]=1 \implies A_{01} = -2.0$.
    </li>
    <li>Row 1: Off-diagonals start at $\text{JA}[1]=6$ up to $\text{JA}[2]-1=7$ (2 elements).
    <br/>Index 6: $\text{val}=1.0, \text{col}=\text{JA}[6]=0 \implies A_{10} = 1.0$.
    <br/>Index 7: $\text{val}=-4.0, \text{col}=\text{JA}[7]=3 \implies A_{13} = -4.0$.
    </li>
    <li>Row 2: Off-diagonals start at $\text{JA}[2]=8$ up to $\text{JA}[3]-1=7 \implies \mathbf{JA}[2] == \text{JA}[3] == 8$.
    <br/><strong>Row 2 has zero off-diagonals!</strong>
    </li>
    <li>Row 3: Off-diagonals start at $\text{JA}[3]=8$ up to $\text{JA}[4]-1=8$ (1 element).
    <br/>Index 8: $\text{val}=5.0, \text{col}=\text{JA}[8]=2 \implies A_{32} = 5.0$.
    </li>
  </ul>
  $$A = \begin{bmatrix}
  20.0 & -2.0 & 0.0 & 0.0 \\
  1.0 & 30.0 & 0.0 & -4.0 \\
  0.0 & 0.0 & 40.0 & 0.0 \\
  0.0 & 0.0 & 5.0 & 50.0
  \end{bmatrix}$$
  </li>

  <li><strong>Hand SpMV Calculation ($y = Ax$):</strong>
  <ul>
    <li>$y_0 = \text{AA}[0] \cdot x_0 + \text{AA}[5] \cdot x_{\text{JA}[5]} = 20(1) + (-2)(2) = 20 - 4 = \mathbf{16.0}$</li>
    <li>$y_1 = \text{AA}[1] \cdot x_1 + \text{AA}[6] \cdot x_{\text{JA}[6]} + \text{AA}[7] \cdot x_{\text{JA}[7]} = 30(2) + 1(1) + (-4)(4) = 60 + 1 - 16 = \mathbf{45.0}$</li>
    <li>$y_2 = \text{AA}[2] \cdot x_2 + (\text{no off-diagonals}) = 40(3) = \mathbf{120.0}$</li>
    <li>$y_3 = \text{AA}[3] \cdot x_3 + \text{AA}[8] \cdot x_{\text{JA}[8]} = 50(4) + 5(3) = 200 + 15 = \mathbf{215.0}$</li>
  </ul>
  $$\mathbf{y = \begin{bmatrix} 16.0 \\ 45.0 \\ 120.0 \\ 215.0 \end{bmatrix}}$$
  </li>
</ol>
</div>

</section>
"""

def generate_standalone_html(body_content):
    return f"""<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=yes" />
  <title>Chapter 20: Sparse Matrix Storage Formats & High-Performance Kernels | HPSC Master Notes</title>
  <link rel="stylesheet" href="../css/style.css" />
  <style>
    .standalone-nav-bar {{
      position: fixed;
      top: 0;
      left: 0;
      right: 0;
      height: var(--topbar-height);
      background: rgba(255, 255, 255, 0.92);
      backdrop-filter: blur(12px);
      -webkit-backdrop-filter: blur(12px);
      border-bottom: 1px solid var(--border-color);
      display: flex;
      align-items: center;
      justify-content: space-between;
      padding: 0 1.25rem;
      z-index: 900;
      box-shadow: var(--shadow-sm);
    }}
    [data-theme="dark"] .standalone-nav-bar {{
      background: rgba(17, 24, 39, 0.92);
    }}
    .chapter-nav-cards {{
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 1rem;
      margin: 3rem 0 2rem 0;
      padding-top: 1.5rem;
      border-top: 1px solid var(--border-color);
    }}
    @media (max-width: 640px) {{
      .chapter-nav-cards {{ grid-template-columns: 1fr; }}
    }}
    .nav-card {{
      display: flex;
      flex-direction: column;
      padding: 1rem 1.25rem;
      border: 1px solid var(--border-color);
      border-radius: 8px;
      background: var(--bg-surface);
      text-decoration: none;
      transition: all 0.2s ease;
    }}
    .nav-card:hover {{
      border-color: var(--primary);
      transform: translateY(-2px);
      box-shadow: var(--shadow-md);
    }}
    .nav-card-label {{
      font-size: 0.78rem;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.5px;
      color: var(--primary);
      margin-bottom: 0.25rem;
    }}
    .nav-card-title {{
      font-size: 0.95rem;
      font-weight: 600;
      color: var(--text-main);
      line-height: 1.35;
    }}
    .nav-next {{
      text-align: right;
    }}
    .btn-hub-link {{
      display: inline-flex;
      align-items: center;
      gap: 0.4rem;
      padding: 0.4rem 0.85rem;
      background: var(--primary);
      color: #ffffff;
      border-radius: 6px;
      font-size: 0.85rem;
      font-weight: 600;
      text-decoration: none;
      transition: all 0.15s ease;
    }}
    .btn-hub-link:hover {{
      background: var(--primary-hover);
      color: #ffffff;
    }}
  </style>
  <script>
    window.MathJax = {{
      tex: {{
        inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
        displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
        processEscapes: true,
        processEnvironments: true
      }},
      options: {{
        skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code']
      }}
    }};
  </script>
  <script src="../js/mathjax-bundle.js" defer></script>
</head>
<body id="top">
  <div id="reading-progress"></div>

  <nav class="standalone-nav-bar" aria-label="Navigation">
    <div class="nav-left">
      <a href="../INDEX.html#chapter-20" class="btn-hub-link">
        <span>🏠</span>
        <span>Master Hub</span>
      </a>
      <a href="chapter-19.html" class="btn-nav-control" title="Previous: The Ultimate Last-Night Cramming Guide & Complete Notation Decryptor">← Prev</a>
      <a href="likely-exam-questions.html" class="btn-nav-control" title="Next: HPSC Examination Question Bank & Model Solutions Repository">Next →</a>
    </div>
    <div class="nav-right">
      <button id="btn-theme-toggle" class="btn-nav-control" onclick="toggleTheme()" title="Toggle dark/light theme" type="button">
        <span id="theme-icon">🌓</span>
      </button>
    </div>
  </nav>

  <main class="master-wrapper" style="margin-top: 1rem;">
    <article class="chapter-card" style="border-top: 4px solid var(--primary);">
      <header class="chapter-card-header">
        <div class="chapter-header-left">
          <div class="chapter-badge-row">
            <span class="hub-badge badge-primary">Chapter 20</span><span class="hub-badge badge-primary">Part I / Hardware Kernels</span>
          </div>
          <h1 class="chapter-title-heading" style="margin: 0; font-size: 1.6rem;">Sparse Matrix Storage Formats & High-Performance Kernels: CSR, CSC, DIA, ELLPACK, ELLPACK-ITPACK, and COO</h1>
        </div>
        <div class="chapter-actions">
          <a href="../INDEX.html#chapter-20" class="btn-action">🏠 Master Hub</a>
          <a href="#top" class="btn-action">⬆ Top</a>
        </div>
      </header>

      <div class="chapter-card-body">
        {body_content}
      </div>
    </article>

    <div class="chapter-nav-cards">
      <a href="chapter-19.html" class="nav-card nav-prev">
        <span class="nav-card-label">← Previous Chapter</span>
        <span class="nav-card-title">Chapter 19: The Ultimate Last-Night Cramming Guide &amp; Complete Notation Decryptor</span>
      </a>
      <a href="likely-exam-questions.html" class="nav-card nav-next">
        <span class="nav-card-label">Next Resource →</span>
        <span class="nav-card-title">HPSC Examination Question Bank &amp; Model Solutions Repository</span>
      </a>
    </div>
  </main>

  <script>
    // Theme Toggle
    function toggleTheme() {{
      var currentTheme = document.documentElement.getAttribute('data-theme') || 'light';
      var newTheme = currentTheme === 'light' ? 'dark' : 'light';
      document.documentElement.setAttribute('data-theme', newTheme);
      try {{
        localStorage.setItem('hpsc_theme', newTheme);
      }} catch(e) {{}}
    }}

    (function() {{
      try {{
        var savedTheme = localStorage.getItem('hpsc_theme');
        if (savedTheme) {{
          document.documentElement.setAttribute('data-theme', savedTheme);
        }} else if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {{
          document.documentElement.setAttribute('data-theme', 'dark');
        }}
      }} catch(e) {{}}
    }})();

    // Progress Bar
    window.addEventListener('scroll', function() {{
      var docEl = document.documentElement;
      var scrollTotal = docEl.scrollHeight - docEl.clientHeight;
      var scrollCurrent = docEl.scrollTop || document.body.scrollTop;
      var progress = (scrollTotal > 0) ? (scrollCurrent / scrollTotal) * 100 : 0;
      var progressBar = document.getElementById('reading-progress');
      if (progressBar) {{
        progressBar.style.width = Math.min(100, Math.max(0, progress)) + '%';
      }}
    }}, {{ passive: true }});
  </script>
</body>
</html>
"""

def generate_js_data(body_content):
    # Escape backticks and ${} for template literal, or output JSON string
    json_escaped = json.dumps(body_content)
    return f"""window.__HPSC_CHAPTERS = window.__HPSC_CHAPTERS || {{}};
window.__HPSC_CHAPTERS["chapter-20"] = {json_escaped};
"""

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    chapters_dir = os.path.join(base_dir, "chapters")
    data_dir = os.path.join(chapters_dir, "data")
    
    body = get_chapter20_body_html()
    
    # 1. Write chapters/chapter-20.html
    html_file = os.path.join(chapters_dir, "chapter-20.html")
    with open(html_file, "w", encoding="utf-8") as f:
        f.write(generate_standalone_html(body))
    print(f"Wrote {html_file} ({os.path.getsize(html_file)} bytes)")
    
    # 2. Write chapters/data/chapter-20.js
    js_file = os.path.join(data_dir, "chapter-20.js")
    with open(js_file, "w", encoding="utf-8") as f:
        f.write(generate_js_data(body))
    print(f"Wrote {js_file} ({os.path.getsize(js_file)} bytes)")

if __name__ == "__main__":
    main()
