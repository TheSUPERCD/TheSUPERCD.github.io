# Chapter 01: Foundations of Scientific Computing & High-Performance Computing (HPC)

---

## 1. Introduction and Course Overview

High Performance Scientific Computing (HPSC) lies at the intersection of three fundamental pillars:
1. **Domain Science & Engineering**: Physics, Chemistry, Biology, Mechanical Engineering, Fluid Dynamics, Materials Science, Astrophysics.
2. **Mathematics & Numerical Analysis**: Partial Differential Equations (PDEs), Numerical Linear Algebra, Matrix Equations ($Ax = b$), Optimization, Discretization Theory.
3. **Computer Science & Architecture**: High Performance Computing (HPC) infrastructure, Multi-core Processors, Graphics Processing Units (GPUs), Memory Hierarchies, Caching Strategies, Parallel Programming Models (MPI, OpenMP, CUDA).

The course is partitioned into two major intertwined modules:
- **Module A: Scientific Computing**: Formulation, discretization, storage, conditioning, and algorithmic solution of large-scale linear systems of equations ($Ax = b$) originating from physical continuum models.
- **Module B: High Performance & Parallel Computing**: Principles of parallel architecture, memory systems, interconnection networks, task decomposition, parallel algorithm design, scalability analysis, and parallel numerical linear algebra.

---

## 2. Serial vs. Parallel Computing & Sequential Architecture

### 2.1 The Sequential (Von Neumann) Computer Architecture
A standard sequential computing node consists of three tightly coupled components communicating across a system bus:

```
┌───────────────────────────────────────────────┐
│                     CPU                       │
│  ┌─────────────────────────┐  ┌─────────────┐ │      Instruction Stream
│  │        Registers        │  │   Control   │─┼──────────────────────┐
│  ├─────────────────────────┤  │    Unit     │ │                      ▼
│  │ Arithmetic & Logic Unit │  │             │ │◄───────────────►┌───────────────┐
│  │          (ALU)          │  │             │ │   Data Stream   │ Dynamic RAM   │
│  └─────────────────────────┘  └─────────────┘ │                 │    (DRAM)     │
└───────────────────────┬───────────────────────┘                 └───────────────┘
                        │
                        ▼
         ┌───────────────────────────────┐
         │ Input/Output (I/O) Controller │
         └──────────────┬────────────────┘
                        │
    ┌───────────┬───────┴───┬──────────────┬────────────┐
    ▼           ▼           ▼              ▼            ▼
Hard Disk   Keyboard/   Network/       Display/    Experimental
 / SSD        Mouse       Modem        Monitor      Instruments
```

#### Performance Driving Rules:
1. **Processor Speed**: For example, a $1.0\text{ GHz}$ processor cycle takes $10^{-9}\text{ s}$. Executing 4 floating-point operations per cycle yields an ideal peak speed of $4\text{ GigaFLOPS}$.
2. **RAM Bandwidth**: Rate at which data is transferred from DRAM into processor registers (bytes or words per cycle).
3. **RAM Latency**: Time elapsed between the CPU issuing a read/write memory request and the data arriving in registers.
4. **Governing Law of Real Performance**:
   > [!IMPORTANT]
   > *Computer performance is strictly driven by both processor speed and RAM latency, whichever is smaller!* A fast CPU stalled on slow RAM operates only at memory-bus speeds (the "Memory Wall").

### 2.2 Serial vs. Parallel Computing Analogy
- **Serial Paradigm**: A computational problem is broken into a linear sequence of discrete instructions executed one after another on a single core.
  - *Analogy*: One person solves a 2000-piece jigsaw puzzle alone $\implies$ takes $10\text{ hours}$.
- **Parallel Paradigm**: The problem is partitioned into concurrent tasks executed across multiple cooperative processing elements.
  - *Analogy*: Two people collaborate to solve the 2000-piece puzzle $\implies$ takes $5\text{ hours } 30\text{ minutes}$.
  - *Why not in 5 hours?* The extra $30\text{ minutes}$ is lost to **parallel overhead**: communication (discussing piece locations), coordination, synchronization, and contention for shared table space.

---

## 3. Physical Motivations for Large-Scale Computing

Scientific computing is indispensable because experimental investigations often reach fundamental physical, economic, or geometric limits:

### 3.1 Limitations of Experimental Methods vs. Theoretical/Numerical Studies
| Metric / Aspect | Experimental Methods | Theoretical & Numerical Simulation |
| :--- | :--- | :--- |
| **Intrusiveness** | Obtrusive (probes, Pitot tubes, thermocouples disrupt the flow field). | Completely non-intrusive (field variables computed everywhere). |
| **Scale & Scope** | Limited to prototypes, wind tunnels, or scaled-down test benches. | Capable of modeling full-scale, real-time astrophysical or nano-scale systems. |
| **Cost** | Extremely expensive; often destructive (e.g., crash tests, explosive testing). | Highly cost-effective once software and compute infrastructure are established. |
| **Measurement Errors** | Instrumentation drift, precision errors, calibration limits, human error. | Governed strictly by discretization, truncation, and round-off errors. |
| **Data Granularity** | Gross parametric results at a sparse set of probe locations. | High-fidelity, full-field spatial and temporal resolutions. |
| **Visualization** | Difficult to visualize internal 3D flow patterns and stresses. | Comprehensive 3D volume rendering, streamlines, vorticity contours. |

### 3.2 Multiscale and Multiphysics Cascades

#### 1. Multiscale Modeling of Brain Vasculature (Perdikaris et al., 2016)
Spans a profound cascade of spatial and temporal scales:
- **Protein Folding**: Length $\sim 10^{-8}\text{ m}$, Time $\sim 10^{-8}\text{ s}$ (Yotta/Zetta scale operations).
- **Cell Biomechanics**: Molecular Dynamics (MD), Length $\sim 10^{-6}\text{ m}$, Time $\sim 10^{-6}\text{ s}$.
- **Cell Interactions**: Dissipative Particle Dynamics (DPD), Length $\sim 10^{-4}\text{ m}$, Time $\sim 10^{-4}\text{ s}$.
- **Flow-Structure Interactions (FSI)**: 3D/1D Continuum Mechanics, Spectral/$hp$ Element Methods, Length $\sim 10^{-2}\text{ m}$, Time $\sim 10^{-2}\text{ s}$.

#### 2. Multiscale Nanochannel Flow Simulation (O'Connel & Thomson, 1995)
Couples three spatial subdomains seamlessly:
- $\Omega_C$: Continuum fluid mechanics subdomain (Navier-Stokes equations).
- $\Omega_O$: Overlap handshake domain where velocity and flux boundary conditions match.
- $\Omega_A$: Atomistic molecular dynamics subdomain resolving microscopic wall interactions.

---

## 4. Quantitative Demand Estimation: Turbulent Combustion Example

Consider the Direct Numerical Simulation (DNS) of high-speed turbulent jet flame premixed combustion (*Gauding et al., Physics of Fluids, 2007*):

### Mathematical Derivation of FLOP & Memory Requirements
1. **Spatial Resolution (Kolmogorov Scale $\eta$)**:
   - Smallest turbulent eddy scale scales as $\eta / L \sim Re^{-3/4}$.
   - To resolve all eddies in a 3D volume without empirical turbulence modeling:
     $$N_{grid} \sim \left(\frac{L}{\eta}\right)^3 \sim \left(Re^{3/4}\right)^3 = Re^{9/4}$$
   - For realistic turbulent combustion at $Re \sim 10^6$:
     $$N_{grid} \approx (10^6)^{9/4} = 10^{13.5} \approx 10^{13}\text{ to } 10^{14}\text{ grid points!}$$

2. **Per-Time-Step Work**:
   - Solving coupled Navier-Stokes, energy, and multi-species chemical kinetics at each node requires:
     $$\text{Work per step} \approx O(10^{13})\text{ FLOP}$$

3. **Temporal Resolution (CFL Limit)**:
   - For numerical stability, $\Delta t \sim 10^{-3}\text{ s}$. Simulating a modest $10\text{ seconds}$ physical event:
     $$N_{steps} = \frac{10\text{ s}}{10^{-3}\text{ s}} = 10^4\text{ time steps}$$

4. **Total Work & Serial Wall-Clock Time**:
   $$\text{Total FLOPs} = 10^4\text{ steps} \times 10^{13}\text{ FLOP/step} = 10^{17}\text{ FLOP}$$
   - On a top-tier uniprocessor delivering sustained $10\text{ GFLOPS}$ ($10^{10}\text{ FLOP/s}$):
     $$T_{serial} = \frac{10^{17}\text{ FLOP}}{10^{10}\text{ FLOP/s}} = 10^7\text{ seconds} \approx 317\text{ years!}$$

5. **Memory Capacity Wall**:
   - Storing $10^{13}$ double-precision complex or coupled state vectors:
     $$\text{Memory} = 10^{13} \times 16\text{ bytes} = 1.6 \times 10^{14}\text{ bytes} = 160\text{ Terabytes!}$$
   - No single workstation RAM can hold this matrix, necessitating large-scale distributed-memory clusters.

---

## 5. Overview of Discretization Methods

The journey from a continuum physical law to a digital solution proceeds by converting continuous differential equations into algebraic linear equations ($Ax = b$).

```
Physical Governing Laws (PDEs)
  (Mass, Momentum, Energy)
          │
          ▼
Discretization Scheme
  (FDM / FVM / FEM)
          │
          ▼
System of Algebraic Equations
         Ax = b
          │
          ▼
Sparse Matrix Linear Solver
  (Direct or Iterative)
          │
          ▼
Physical Solution Field T(x, y, z, t)
```

### 5.1 Finite Difference Method (FDM)
- **Basis**: Truncated Taylor series expansions approximating derivatives at discrete grid points.
- **Formulation**:
  $$\left.\frac{\partial^2 T}{\partial x^2}\right|_i \approx \frac{T_{i+1} - 2T_i + T_{i-1}}{\Delta x^2} + O(\Delta x^2)$$
- **Pros**: Conceptually straightforward; easy to implement high-order schemes on structured grids.
- **Cons**: Difficult to apply to complex, curved, or irregular geometries.

### 5.2 Finite Volume Method (FVM)
- **Basis**: Integral conservation laws applied across discrete control volumes.
- **Formulation**:
  $$\int_{\partial V} \mathbf{F} \cdot \mathbf{n} \, dA = \int_V S \, dV$$
- **Pros**: Strictly conservative locally and globally; well-suited for fluid dynamics and convective transport.
- **Cons**: Higher-order reconstructions on unstructured 3D meshes become mathematically complex.

### 5.3 Finite Element Method (FEM)
- **Basis**: Weak (variational) formulation of PDEs over unstructured subdivisions (elements).
- **Formulation**: Residual error weighted with test functions $\psi_i$ and integrated over elements:
  $$\int_\Omega \left( \nabla \psi_i \cdot \nabla T - \psi_i f \right) d\Omega = 0$$
- **Pros**: Exceptional geometric flexibility; mathematically rigorous; handles anisotropic boundary conditions naturally.
- **Cons**: High memory and computational cost for stiffness matrix assembly and solution.

---

## 6. Numerical Errors, Stability, and Optimal Discretization

Every numerical solution inevitably incurs errors. In an exam, distinguishing the exact categories of error is essential:

### 6.1 Classification of Numerical Errors
1. **Modeling Error**: Difference between the true physical phenomenon and the formulated mathematical equations (e.g., assuming laminar flow when turbulence exists, or assuming ideal gas behavior).
2. **Discretization / Truncation Error**: Error introduced by replacing continuous differential operators with discrete algebraic approximations (e.g., dropping higher-order terms in a Taylor series).
   $$\text{Error}_{trunc} = O(\Delta x^p)$$
   - *Behavior*: Decreases monotonically as grid spacing $\Delta x \to 0$ (or number of grid points $N \to \infty$).
3. **Round-Off Error**: Error caused by the finite precision representation of real numbers in computer hardware (IEEE 754 floating point standard).
   - *Behavior*: Accumulates with every floating-point operation. As $\Delta x$ decreases, the total number of operations $N$ surges, causing round-off error to **increase**!

### 6.2 Total Error and the Optimal Grid Spacing
The total computational error is the superposition of truncation error and round-off error:
$$\text{Error}_{total}(\Delta x) = \text{Error}_{discretization}(\Delta x) + \text{Error}_{round-off}(\Delta x) \approx C_1 (\Delta x)^p + \frac{C_2}{\Delta x^q}$$

```
Error
  ▲
  │ \                                           /  Round-off Error
  │  \                                         /   (accumulates with ops)
  │   \                                       /
  │    \             Total Error             /
  │     \           \___________/           /
  │      \                │                /
  │       \               ▼               /
  │        \_        Optimal Point       /
  │          \______       •       _____/
  │                 \____     ____/
  │                      \___/
  │                       │
  │ Discretization Error  │
  │ (vanishes as Δx -> 0) │
  └───────────────────────┼────────────────────────► Grid Spacing (Δx)
                      Δx_optimal
```

- **Optimal Solution ($\Delta x_{optimal}$)**: The grid spacing at which the total error achieves its global minimum. Refining beyond this threshold leads to numerical degradation due to round-off error domination.
- **Grid-Independent Solution**: A numerical solution where further refinement of $\Delta x$ changes the solution by less than a prescribed tolerance $\epsilon_{tol}$, confirming that discretization error is sufficiently minimized.

---

## 7. Verification vs. Validation (V&V)

```
               Real World Physical Phenomenon
                            │
               Modeling     ▼
               Conceptual Mathematical Model (PDEs)
                            │
            Discretization  ▼
               Computer Code / Numerical Solver
                            │
               Simulation   ▼
                    Numerical Result
```

### 7.1 Verification ("Solving the equations right")
- **Definition**: Confirming that the numerical code accurately reflects the conceptual mathematical model.
- **Grid Convergence Study**: Plotting error versus grid spacing on log-log coordinates:
  $$\log(\text{Error}) = p \cdot \log(\Delta x) + \text{constant}$$
- **Slide 14 Benchmark**: For 2D flow past a cylinder at $Re = 100$:
  - Measuring drag coefficient $C_d$ and recirculation length $x_s$ across varying $\Delta x$.
  - On a plot of $\log(|\text{error}|)$ vs. $\log(dx)$, the computed slopes verify the code's spatial accuracy:
    - **Slope = 1**: First-order scheme verification.
    - **Slope = 2**: Second-order central-difference scheme verification.

### 7.2 Validation ("Solving the right equations")
- **Definition**: Determining that the mathematical model adequately represents the true physical reality for its intended purpose.
- **Quantitative Comparison**:
  - Plotting normalized axial and radial velocity profiles $V/V_T$ vs. normalized radius $r/R$.
  - Direct overlay of laser Doppler velocimetry (LDV) experimental data points against computational curves across multiple axial stations.

---

## 8. High Performance Computing Infrastructure & Supercomputing Systems

### 8.1 Performance Metrics: FLOPS & Energy Efficiency
- **Ideal Peak Processor Speed**:
  $$\text{Peak FLOPS} = \text{Cores} \times \text{Clock Frequency (Hz)} \times \frac{\text{FLOPs}}{\text{Cycle}}$$
- **Power Usage Effectiveness (PUE)**:
  Supercomputers consume immense electrical power, generating massive thermal loads:
  $$\text{PUE} = \frac{\text{Total Facility Energy Requirement}}{\text{Energy Requirement of IT Components (CPUs, GPUs, Network, Storage)}}$$
  - For leading high-efficiency green supercomputers: $\text{PUE} \approx 1.1\text{ to } 1.2$.
  - For typical enterprise supercomputers: $\text{PUE} \approx 1.5$.

### 8.2 Top500 Global Supercomputers (July 2026 Standing)
| Rank | System | Specifications & Interconnect | Cores | $R_{max}$ (PFLOPS) | $R_{peak}$ (PFLOPS) | Power (kW) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | **LineShine** (China) | LingKun LX2 304C 1.55 GHz, LingQi, Kylin OS | 13,789,440 | **2,198.40** | 2,735.82 | 42,220 |
| **2** | **El Capitan** (LLNL, US) | HPE Cray EX255a, AMD EPYC 24C + MI300A, Slingshot-11 | 11,340,000 | **1,809.00** | 2,821.10 | 29,685 |
| **3** | **Frontier** (ORNL, US) | HPE Cray EX235a, AMD EPYC 64C + MI250X, Slingshot-11 | 9,066,176 | **1,353.00** | 2,055.72 | 24,607 |
| **4** | **Aurora** (ANL, US) | HPE Cray EX, Intel Xeon Max 9470 + Data Center GPU Max | 9,264,128 | **1,012.00** | 1,980.01 | 38,698 |
| **5** | **JUPITER Booster** (FZJ, EU) | BullSequana XH3000, NVIDIA GH200 Grace Hopper, NDR200 | 4,801,344 | **1,000.00** | 1,226.28 | 15,794 |

### 8.3 Indian Supercomputing Ecosystem (National Supercomputing Mission - NSM)
Funded at ₹4,500 Crores under NSM to deploy indigenous petascale infrastructure:
1. **AIRAWAT - PSAI** (C-DAC Pune): Rank 188, $R_{max} = 8.50\text{ PFLOPS}$, $R_{peak} = 13.17\text{ PFLOPS}$.
2. **Arka** (IITM Pune, Earth Sciences): Rank 196, $R_{max} = 5.94\text{ PFLOPS}$, $R_{peak} = 7.40\text{ PFLOPS}$.
3. **Arunika** (NCMRWF): Rank 251, $R_{max} = 5.94\text{ PFLOPS}$.
4. **Pratyush** (Cray XC40, IITM): Rank 338, $R_{max} = 3.76\text{ PFLOPS}$.
5. **Mihir** (Cray XC40, NCMRWF): $R_{max} = 2.57\text{ PFLOPS}$.
6. **PARAM Shivay** (IIT BHU): First NSM supercomputer inaugurated Feb 2019.
7. **PARAM Shakti (IIT Kharagpur)**:
   - Peak Performance: **$1.66\text{ PFLOPS}$**.
   - **CPU-Only Compute Nodes (384 nodes)**:
     - $2\times$ Intel Xeon SKL G-6148 (40 cores/node, 2.4 GHz).
     - Memory: 192 GB DDR4 2666 MHz per node. Total Cores = 15,360; Total RAM = 73,728 GB. Local SSD scratch = 480 GB.
   - **High-Memory Compute Nodes (36 nodes)**:
     - $2\times$ Intel Xeon SKL G-6148. Memory: 768 GB DDR4 per node. Total Cores = 1,440; Total RAM = 27,648 GB.
   - **GPU Compute Nodes (22 nodes)**:
     - $2\times$ Intel Xeon SKL G-6148 (880 total CPU cores), 192 GB RAM.
     - $2\times$ NVIDIA Tesla V100 (16 GB HBM2) per node (44 GPUs total).
     - GPU Cores per node = $2 \times 5,120 = 10,240$ CUDA cores; $2 \times 640 = 1,280$ Tensor cores.

---

## 9. Real-World HPC Scientific & Biomedical Case Studies

### 9.1 Vesicle Fusion in Binary Lipid Mixtures (Kohlmayer, 2011)
- **Problem**: Simulating vesicle fusion and lipid phase separation in computational biochemistry.
- **Workload**: 4 million particles per vesicle + solvent system using the **LAMMPS** molecular dynamics engine.
- **Compute Resources**: 30 million CPU-hours on the **Cray XT5 Jaguar** supercomputer (224,256 cores, 1.75 PFLOPS).

### 9.2 Black Hole Grazing Collisions & Gravitational Waves (Alcubierre et al., 2001; Ott et al., 2008)
- **Problem**: 3D numerical relativity solving Einstein field equations for binary black hole mergers.
- **Code**: $100^3$ grid WaveToy benchmark on the **Abe** supercomputer.
- **Scaling Comparison**: Pure MPI (**PUGH**) vs. Hybrid MPI+OpenMP (**Carpet** adaptive mesh refinement) demonstrating sustained parallel scaling up to 4,096 processor cores.

### 9.3 GPU Interactive Cloud Dynamics (Harris et al., 2003, ACM SIGGRAPH)
Demonstrated early GPU acceleration for solving the 2D/3D pressure Poisson equation ($\nabla^2 p = S$):
| Poisson Solver Algorithm | Convergence / msec | Relative Convergence (100 iters) | Time for 100 iters (msec) |
| :--- | :--- | :--- | :--- |
| **Jacobi 2D** | 0.50 | 0.078 | 45.9 ms |
| **Red-Black Gauss-Seidel 2D** | 0.85 | 0.124 | 45.3 ms |
| **Vectorized Jacobi 2D** | **1.00** | **0.079** | **17.3 ms** ($2.65\times$ speedup) |
| **Jacobi 3D** | N/C | N/C | 110.0 ms |
| **Vectorized Jacobi 3D** | N/C | N/C | **49.0 ms** ($2.24\times$ speedup) |

### 9.4 Impeller Mixing in Chemical Stirred Tanks
- **Problem**: 3D turbulent scalar mixing across axial planes ($z/R = -0.5, 0.0, 1.0, 2.0$) with rotating blades.
- **Discretization**: 3.8 million grid points executed across 1,000 parallel processors over 60,000 iterations.
- **Runtime**: 2 months continuous execution. Illustrates sublinear speedup scaling due to all-to-all communication at blade-baffle interfaces.

### 9.5 Biomedical Pulsatile Blood Flow in Stenosed Arteries
- **Clinical Motivation**: Severe arterial narrowing ($58\%$ stenosis) induces disturbed post-stenotic vortex shedding and abnormal Time-Averaged Wall Shear Stress (TAWSS: $0.002 - 0.022\text{ Pa}$), triggering thrombus formation.
- **Computational Benchmark**: Immersed Boundary (IB) solver comparing $2\times 16$-core Intel Xeon E5-2698 v3 (2.3 GHz) vs. NVIDIA Tesla P100 (16 GB):
  - At 4,000 time steps: Sequential CPU ($\sim 30,000\text{ s}$) $\to$ Multicore OpenACC ($\sim 5,000\text{ s}$) $\to$ GPU OpenACC ($\sim 500\text{ s}$, **$60\times$ overall acceleration**).

---

## 10. Course Structure & Reference Textbooks

### Textbooks for Module 1 (Scientific Computing & Linear Algebra)
1. **Gilbert Strang**, *Introduction to Linear Algebra*, 4th Edition, Cengage / Wellesley-Cambridge Press (2006).
2. **Yousef Saad**, *Iterative Methods for Sparse Linear Systems*, 2nd Edition, SIAM (2003).

### Textbooks for Module 2 (Parallel Computing & GPU Programming)
1. **Ananth Grama, Anshul Gupta, George Karypis, Vipin Kumar**, *Introduction to Parallel Computing*, 2nd Edition, Addison-Wesley / Pearson (2003).
2. **Michael J. Quinn**, *Parallel Programming in C with MPI and OpenMP*, McGraw-Hill (2004).
3. **William Gropp, Ewing Lusk, Anthony Skjellum**, *Using MPI: Portable Parallel Programming with the Message-Passing Interface*, 3rd Edition, MIT Press.
4. **Miguel Hermanns**, *Parallel Programming in Fortran 95 using OpenMP*, Universidad Politécnica de Madrid.
5. **Blaise Barney**, *Introduction to Parallel Computing Tutorial*, Lawrence Livermore National Laboratory (LLNL).

