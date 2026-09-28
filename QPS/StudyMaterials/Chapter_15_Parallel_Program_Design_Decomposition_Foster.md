# Chapter 15: Parallel Program Design, Decomposition & Foster's Methodology

---

## 1. Principles of Parallel Algorithm Design

Transforming a sequential scientific computation into an efficient parallel program requires orchestrating four fundamental activities:
1. **Decomposition**: Dividing total computational work into concurrent primitive tasks.
2. **Communication**: Coordinating data exchanges and dependencies between tasks.
3. **Agglomeration**: Clustering primitive tasks into larger grains to optimize locality.
4. **Mapping**: Assigning agglomerated composite tasks to physical hardware cores.

---

## 2. Task Dependency Graphs & Directed Acyclic Graphs (DAGs)

The structural workflow of a parallel algorithm is formally modeled as a **Task Dependency Graph** (a Directed Acyclic Graph, DAG).

```
       Start: Task 1 (Initialize Grid)
                     │
         ┌───────────┴───────────┐
         ▼                       ▼
      Task 2 (North Half)     Task 3 (South Half)    ◄── Concurrency = 2
         │                       │
         └───────────┬───────────┘
                     ▼
             Task 4 (Exchange Boundaries)
                     │
                     ▼
             Task 5 (Compute Global Residual)
```

### 2.1 Formal Components of a DAG
- **Nodes ($V$)**: Individual, indivisible units of computation called **tasks**.
- **Directed Edges ($E$)**: Data or control dependencies. A directed edge from Task $A$ to Task $B$ ($A \to B$) indicates that Task $B$ cannot begin execution until Task $A$ completes and its output data becomes available.

### 2.2 Crucial DAG Metrics
1. **Total Work ($W$ or $T_1$)**:
   The total computational time spent executing all tasks sequentially on a single processor:
   $$W = \sum_{v \in V} t(v)$$
2. **Critical Path**:
   The longest directed path from an entry node to an exit node in the DAG.
3. **Critical Path Length ($T_\infty$)**:
   The execution time of the tasks along the critical path.
   - **Theoretical Foundation**: Even if an infinite number of processors are available ($p = \infty$), the parallel execution time can **never be less than $T_\infty$**:
     $$\mathbf{T_p \ge T_\infty}$$
4. **Maximum Degree of Concurrency**:
   The maximum number of independent tasks that can execute simultaneously at any point in time.
5. **Average Degree of Concurrency**:
   The ratio of total sequential work to the critical path length:
   $$\mathbf{\text{Average Concurrency} = \frac{W}{T_\infty}}$$
   - **Exam Significance**: The average degree of concurrency represents an upper bound on the number of processors that can be productively utilized before diminishing returns dominate!

---

### 2.3 Worked Case Study: Database Query DAG Optimization (Slides 316–321)

Consider processing a complex database query across a student registry:
$$\text{Query: } (\text{Fee} \le 12000) \;\mathbf{AND}\; (\text{Board} = \text{State}) \;\mathbf{AND}\; \Big(\text{Medium} = \text{'Hindi'} \;\mathbf{OR}\; \text{Medium} = \text{'English'}\Big)$$

The same algorithm can be structured into two distinct task dependency graphs with identical operational outputs:

```
        Abstraction-1 (Balanced Tree)                    Abstraction-2 (Cascade)

  Task 4     Task 3     Task 2     Task 1          Task 4     Task 3     Task 2     Task 1
 (Fee≤12k)  (State)    (Hindi)    (English)       (Fee≤12k)  (State)    (Hindi)    (English)
  [t = 10]   [t = 10]   [t = 10]   [t = 10]        [t = 10]   [t = 10]   [t = 10]   [t = 10]
     │          │          │          │               │          │          │          │
     └────┬─────┘          └────┬─────┘               │          │          └────┬─────┘
          ▼                     ▼                     │          │               ▼
        Task 6                Task 5                  │          │             Task 5
    (Intersection)           (Union)                  │          │         (Union, t = 6)
       [t = 9]               [t = 6]                  │          │               │
          │                     │                     │          └───────┬───────┘
          └──────────┬──────────┘                     │                  ▼
                     ▼                                │                Task 6
                   Task 7                             │            (Intersection)
             (Final Filter)                           │               [t = 11]
                [t = 8]                               │                  │
                                                      └──────────┬───────┘
                                                                 ▼
                                                               Task 7
                                                           (Final Filter)
                                                              [t = 7]
```

#### Quantitative Performance Comparison:
1. **Abstraction-1**:
   - **Total Work ($W$)**: $10 + 10 + 10 + 10 + 6 + 9 + 8 = \mathbf{63}$.
   - **Critical Path**: $\text{Task 4} \to \text{Task 6} \to \text{Task 7}$ (or $\text{Task 3} \to 6 \to 7$).
   - **Critical Path Length ($T_\infty$)**: $10 + 9 + 8 = \mathbf{27}$.
   - **Average Degree of Concurrency**:
     $$\frac{W}{T_\infty} = \frac{63}{27} = \mathbf{2.33}$$
   - **Maximum Concurrency**: 4 simultaneous tasks (Tasks 1, 2, 3, 4 at start).

2. **Abstraction-2**:
   - **Total Work ($W$)**: $10 + 10 + 10 + 10 + 6 + 11 + 7 = \mathbf{64}$.
   - **Critical Path**: $\text{Task 1} \to \text{Task 5} \to \text{Task 6} \to \text{Task 7}$.
   - **Critical Path Length ($T_\infty$)**: $10 + 6 + 11 + 7 = \mathbf{34}$.
   - **Average Degree of Concurrency**:
     $$\frac{W}{T_\infty} = \frac{64}{34} = \mathbf{1.88}$$
   - **Maximum Concurrency**: 2 simultaneous tasks (Tasks 1 and 2; Task 3 must idle waiting for Task 5, Task 4 idles waiting for Task 6).

> [!IMPORTANT]
> **Exam Conclusion**: Abstraction-1 is significantly superior for parallel execution because its **shorter critical path ($27 < 34$)** provides a **higher average degree of concurrency ($2.33 > 1.88$)**, maximizing parallel efficiency.


---

## 3. Decomposition Techniques

Parallel algorithms employ four foundational decomposition strategies:

### 3.1 Recursive Decomposition (Divide-and-Conquer)
- **Concept**: A problem is recursively decomposed into smaller, independent subproblems of the same type until a primitive base case is reached. Subproblem solutions are subsequently combined.
- **Examples**:
  - Parallel Merge Sort and Quick Sort.
  - Tree-based global reductions (e.g., computing the sum or maximum of an array across $p$ processors).

### 3.2 Data Decomposition (Domain Decomposition)
The most ubiquitous strategy in High Performance Scientific Computing. The problem is partitioned by slicing the underlying data structures.

```
       Input Data Decomposition                Output Data Decomposition
       (e.g., Numerical Quadrature)            (e.g., Matrix Multiplication C = A x B)
       ┌──────┬──────┬──────┬──────┐           ┌───────────┬───────────┐
       │P0    │P1    │P2    │P3    │           │ Task 1:   │ Task 2:   │
       │x0..x1│x1..x2│x2..x3│x3..x4│           │ C11 = ... │ C12 = ... │
       └──────┴──────┴──────┴──────┘           ├───────────┼───────────┤
         ▲      ▲      ▲      ▲                │ Task 3:   │ Task 4:   │
         └──────┴──────┼──────┘                │ C21 = ... │ C22 = ... │
                       ▼                       └───────────┴───────────┘
                Global Reduction (Sum)
```

1. **Partitioning the Output Data**:
   - Each task is responsible for computing a designated subset of the output data.
   - *Example: 2D Block Matrix Multiplication*:
     $$\begin{bmatrix} C_{11} & C_{12} \\ C_{21} & C_{22} \end{bmatrix} = \begin{bmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix} \begin{bmatrix} B_{11} & B_{12} \\ B_{21} & B_{22} \end{bmatrix}$$
     - Task 1 computes $C_{11} = A_{11} B_{11} + A_{12} B_{21}$.
     - Task 2 computes $C_{12} = A_{11} B_{12} + A_{12} B_{22}$.
     - Task 3 computes $C_{21} = A_{21} B_{11} + A_{22} B_{21}$.
     - Task 4 computes $C_{22} = A_{21} B_{12} + A_{22} B_{22}$.
2. **Partitioning the Input Data**:
   - Input data is divided among processes, each performing local calculations on its partition, followed by an aggregation step.
   - *Example: Numerical Integration (Trapezoidal Rule)*:
     $$I = \int_a^b f(x) \, dx = \sum_{i=1}^p \text{LocalSum}_i \quad \text{where } \text{LocalSum}_i = \int_{x_{i-1}}^{x_i} f(x) \, dx$$
3. **The Owner Computes Rule**:
   A cardinal rule in distributed-memory programming: **A process executes only the computations that determine data elements it officially owns in local memory**.

### 3.3 Exploratory Decomposition
- Used when searching through large, irregular discrete search spaces (e.g., game-tree search, graph backtracking, integer programming).
- Different processors search distinct branches of the state space concurrently. As soon as one processor finds the optimal solution, it broadcasts an abort signal to prune remaining branches.

### 3.4 Speculative Decomposition
- Used when dependencies between tasks are uncertain. A processor begins executing a tentative task before knowing whether its input data will be valid.
- If the speculation proves correct, latency is saved; if incorrect, the work is rolled back.

---

## 4. Mapping & Load Balancing

**Load Imbalance** is a primary source of parallel performance degradation. If $p-1$ processors finish in $1\text{ second}$ but 1 processor requires $10\text{ seconds}$, the overall parallel execution time is $10\text{ seconds}$, and the $p-1$ processors sit idle ($90\%$ wasted capacity).

```
       Load Imbalance:                       Balanced Load:
       P0: [████████████████████] (10s)      P0: [██████████] (5s)
       P1: [████] idle............ (2s)      P1: [██████████] (5s)
       P2: [██████] idle.......... (3s)      P2: [██████████] (5s)
       P3: [██] idle.............. (1s)      P3: [██████████] (5s)
```

### 4.1 Static Mapping Strategies
Decisions are made prior to execution and remain fixed throughout runtime.

1. **Block Distribution**:
   - Continuous arrays of size $N$ are partitioned into contiguous blocks of size $N/p$ assigned to processor $k$:
     $$\text{Elements on Processor } k: \quad \left[ k \cdot \frac{N}{p}, \; (k+1) \frac{N}{p} - 1 \right]$$
   - *Pros*: Minimizes communication surface area in grid stencils; excellent spatial cache locality.
2. **Cyclic (Round-Robin) Distribution**:
   - Elements are assigned round-robin: element $i$ maps to processor $i \pmod p$.
   - *Best Use Case*: Used when computational effort per index increases or decreases monotonically (e.g., LU decomposition or Gaussian elimination, where rows become progressively shorter as elimination proceeds).
3. **Block-Cyclic Distribution**:
   - Divides the array into blocks of size $b$, and distributes these blocks cyclically across processors.
   - Balances the trade-off between cache locality (block size $b$) and load balancing (cyclic distribution). Foundation of ScaLAPACK.
4. **Graph Partitioning (e.g., METIS, ParMETIS)**:
   - For unstructured FEM/FVM meshes, the dual connectivity graph is partitioned such that:
     1. Each processor receives an equal number of mesh elements (load balance).
     2. The number of cut edges between partitions (inter-processor communication volume) is strictly minimized.

### 4.2 Dynamic Mapping Strategies
Used when task runtimes are unpredictable, highly irregular, or data-dependent.

1. **Centralized Master-Slave / Work-Pool Model**:
   - A master process maintains a queue of available tasks.
   - Idle worker processes send a request to the master, receive a task, compute the result, and return for the next task.
   - *Limitation*: The master becomes a communication bottleneck if tasks are too fine-grained.
2. **Distributed Work Stealing**:
   - Each processor maintains a local double-ended queue (deque) of tasks.
   - When a processor exhausts its local work, it becomes a "thief" and steals tasks from the deque of a randomly chosen neighboring processor.

---

## 5. Parallel Algorithm Models

1. **Data-Parallel Model**:
   - Processors execute uniform operations simultaneously on distinct partitions of a large dataset.
   - High data locality, minimal synchronization.
2. **Task-Graph Model**:
   - Tasks with arbitrary dependencies are represented as a DAG and executed as their predecessor tasks complete.
3. **Work-Pool Model**:
   - Unordered collection of tasks dynamically assigned to available workers.
4. **Master-Slave Model**:
   - One master coordinates execution, manages I/O, and assigns work to worker processes.
5. **Pipeline / Producer-Consumer Model**:
   - Computation is partitioned into sequential stages. Data streams continuously through the stages in an assembly-line fashion.
6. **Hybrid Model**:
   - Hierarchical composition of multiple models (e.g., pipelined stages, where each individual stage is internally data-parallel).

---

## 6. Foster's Design Methodology (PCAM)

In 1995, Prof. Ian Foster formulated the **PCAM Methodology**, a systematic 4-stage engineering pipeline for parallel algorithm design:

```
            Problem Specification
                      │
        1. PARTITIONING │ (Divide into maximum fine-grained tasks)
                      ▼
        2. COMMUNICATION│ (Define local vs global data exchanges)
                      ▼
        3. AGGLOMERATION│ (Combine tasks to increase granularity)
                      ▼
        4. MAPPING      │ (Assign to physical cores for load balance)
                      ▼
             Parallel Algorithm
```

### Stage 1: Partitioning
- **Objective**: Decompose total computation and data into the **maximum possible number of primitive tasks**.
- Ignore machine architecture, number of physical cores, and practical constraints.
- Maximize potential concurrency.

### Stage 2: Communication
- **Objective**: Analyze the data dependencies between primitive tasks.
- Classify communications:
  - **Local Communication**: Tasks exchange data with immediate neighbors (e.g., stencil neighbor exchanges).
  - **Global Communication**: All tasks participate in collective data operations (e.g., global sum reduction, broadcasts).
- Identify communication bottlenecks and eliminate unnecessary data transfers.

### Stage 3: Agglomeration
- **Objective**: Group primitive tasks into larger composite tasks to **increase granularity** and **reduce communication overhead**.
- *Design Trade-off*:
  - If tasks are too fine-grained: Communication overhead overwhelms arithmetic execution.
  - If tasks are too coarse: Available concurrency drops, leaving physical cores idle.
- Agglomeration combines tasks that share high communication volumes, converting inter-processor network traffic into zero-cost internal memory accesses!

### Stage 4: Mapping
- **Objective**: Assign agglomerated composite tasks to physical processing cores.
- Goals:
  1. Equalize computational load across all cores (Load Balancing).
  2. Map tasks that communicate frequently to the same processor socket or neighboring physical nodes to minimize network hop distance and latency.

---

## 7. Concrete Examples of Foster's Methodology (Slides 338–339)

### Example 1: Domain Decomposition for 2D Heat Conduction
- **Step 1 (Partitioning)**:
  Every discrete grid point $(i, j)$ in the 2D domain is treated as a single primitive task (updating temperature $T_{i, j}$ via the 5-point stencil).
- **Step 2 (Communication)**:
  Each grid point $(i, j)$ requires values from its 4 immediate spatial neighbors $(i \pm 1, j)$ and $(i, j \pm 1)$ at every iteration (Local Communication).
- **Step 3 (Agglomeration)**:
  Primitive grid nodes are grouped into contiguous 2D rectangular sub-domains of size $m \times m$ (agglomerating $m^2$ grid points per sub-domain).
  - *Surface-to-Volume Ratio*:
    - Internal computations: $O(m^2)$ operations.
    - Boundary communication: Only $4m$ edge boundary values need to be communicated across the network!
    - As sub-domain size $m$ grows, the communication-to-computation ratio drops as $\frac{4m}{m^2} = \frac{4}{m} \to 0$.
- **Step 4 (Mapping)**:
  Each agglomerated rectangular block is mapped to one processor core. Neighboring sub-domains are mapped to neighboring processors on the physical network mesh/torus.

---

### Example 2: Finding the Maximum of $n$ Numbers (Reduction Tree)
- **Step 1 (Partitioning)**:
  Computing the pairwise comparison $\max(x_a, x_b)$ is a primitive task. A total of $n-1$ comparison tasks are required.
- **Step 2 (Communication)**:
  Each comparison task requires inputs from two previous tasks, forming a directed binary tree.
- **Step 3 (Agglomeration)**:
  If $p \ll n$, assigning 1 comparison per task is wildly inefficient. Agglomerate $n/p$ numbers onto each of the $p$ processors. Each processor sequentially finds its local maximum in $(n/p) - 1$ steps with **zero communication**.
- **Step 4 (Mapping)**:
  The $p$ local maxima are reduced to the global maximum using a binary reduction tree across the $p$ processors in $\log_2 p$ steps, terminating at root processor $P_0$.
