# Chapter 16: Parallel Performance Analysis, Scalability, Amdahl's Law & Isoefficiency

---

## 1. Quantitative Parallel Performance Metrics

To evaluate whether a parallel implementation delivers real performance gains, parallel computer scientists utilize a formal framework of metrics:

### 1.1 Sequential Execution Time ($T_s$ or $W$)
- **Definition**: The wall-clock execution time taken by the **fastest known sequential algorithm** executing on a single core of the target architecture.
- **Exam Rule**: $T_s$ must **never** be measured by simply running the parallel code on one processor ($p = 1$), because the parallel code contains communication setup and indexing overheads that artificially inflate the baseline.

### 1.2 Parallel Execution Time ($T_p$)
- **Definition**: The elapsed wall-clock time from the moment the parallel algorithm initiates to the moment the very last processor completes its execution.

### 1.3 Speedup ($S$ or $S_p$)
The ratio of sequential execution time to parallel execution time:
$$\mathbf{S(p) = \frac{T_s}{T_p}}$$
- **Linear Speedup**: $S(p) = p$ (ideal benchmark: $p$ processors execute $p$ times faster).
- **Sublinear Speedup**: $S(p) < p$ (the normal real-world regime, caused by communication, idling, and synchronization).
- **Superlinear Speedup**: $S(p) > p$ (occurs under special architectural or algorithmic conditions).

#### Why Does Superlinear Speedup Occur?
1. **Aggregate Cache Memory Effects (Most Common)**:
   - When a massive problem of size $N$ runs on a single core, its data cannot fit into L1/L2/L3 caches, forcing continuous high-latency DRAM page thrashing.
   - When distributed across $p$ cores, each core handles a sub-domain of size $N/p$.
   - The aggregated cache capacity ($p \times \text{Cache}$) is large enough to fit the entire dataset in fast on-chip SRAM cache! Memory stall cycles vanish, yielding $S(p) > p$.
2. **Exploratory Search Pruning**:
   - In combinatorial or state-space searches, parallel threads explore different sub-trees. One thread may encounter a pruning solution early, eliminating massive sub-trees that the serial code would have systematically explored.

---

### 1.4 Efficiency ($E$ or $E_p$)
Efficiency measures the fraction of time processors spend engaged in productive computation:
$$\mathbf{E(p) = \frac{S(p)}{p} = \frac{T_s}{p T_p}}$$
- Ideal linear speedup corresponds to **$100\%$ efficiency ($E = 1.0$)**.
- As the number of processors $p$ increases for a fixed problem size, efficiency monotonically declines toward zero ($E \to 0$).

---

### 1.5 Total Parallel Overhead ($T_o$)
The total overhead $T_o$ represents the collective wasted time spent by all $p$ processors on non-serial activities (communication latency, idling, synchronization barriers, and redundant calculations):
$$\mathbf{T_o(W, p) = p T_p - T_s}$$
Rearranging for parallel execution time:
$$\mathbf{T_p = \frac{T_s + T_o(W, p)}{p}}$$
Substituting into the efficiency equation:
$$\mathbf{E = \frac{T_s}{p T_p} = \frac{T_s}{T_s + T_o(W, p)} = \frac{1}{1 + \frac{T_o(W, p)}{T_s}} = \frac{1}{1 + \frac{T_o(W, p)}{W}}}$$
And speedup:
$$\mathbf{S = \frac{p}{1 + \frac{T_o(W, p)}{W}}}$$

---

## 2. Benchmark Case Study: Parallel Reduction (Adding $n$ Numbers on $p$ Processors)

The lecture slides detail a fundamental benchmark: summing $n$ numbers distributed across $p$ processors (Slides 376–377).

```
    Processor Work:                  Binary Tree Reduction Across Processors:
    P0: Adds n/p numbers locally      P0 ◄── P1   P2 ◄── P3   P4 ◄── P5   P6 ◄── P7  (Step 1)
    P1: Adds n/p numbers locally       \    /      \    /      \    /      \    /
    ...                                 P0   ◄────── P2         P4   ◄────── P6      (Step 2)
    Pp-1: Adds n/p numbers locally        \            /          \            /
                                           P0         ◄──────────── P4                   (Step 3)
```

### 2.1 Derivation of Runtimes and Overhead
1. **Local Sequential Phase**:
   Each processor adds $n/p$ local numbers in $O(n/p)$ time.
2. **Global Reduction Phase**:
   The $p$ partial sums are aggregated across a binary tree network of depth $\log_2 p$. At each tree level, a communication and addition occur, taking $2 \log_2 p$ time.
3. **Parallel Execution Time**:
   $$T_p = \frac{n}{p} + 2 \log_2 p$$
4. **Total Overhead**:
   $$T_o = p T_p - T_s = p \left( \frac{n}{p} + 2 \log_2 p \right) - n = \mathbf{2 p \log_2 p}$$

### 2.2 Analytical Efficiency and Speedup Expressions
$$E = \frac{n}{p T_p} = \frac{n}{p \left( \frac{n}{p} + 2 \log_2 p \right)} = \mathbf{\frac{1}{1 + \frac{2 p \log_2 p}{n}}}$$
$$S = p E = \mathbf{\frac{p}{1 + \frac{2 p \log_2 p}{n}} = \frac{1}{\frac{1}{p} + \frac{2 \log_2 p}{n}}}$$

### 2.3 The Two Fundamental Observations (Slide 377)
1. **Processor Scaling at Fixed Problem Size ($n = \text{const}, p \to \infty$)**:
   $$\lim_{p \to \infty} \frac{2 p \log_2 p}{n} = \infty \implies E \to 0, \quad S \to \frac{n}{2 \log_2 p} \quad (\text{Speedup flattens out!})$$
   Adding processors indefinitely to a fixed-size problem causes communication overhead to dominate completely.
2. **Problem Size Scaling at Fixed Processor Count ($p = \text{const}, n \to \infty$)**:
   $$\lim_{n \to \infty} \frac{2 p \log_2 p}{n} = 0 \implies \mathbf{E \to 1, \quad S \to p}$$
   As problem size increases, the computational volume ($O(n/p)$) overwhelms the communication overhead ($O(\log p)$), yielding near-perfect linear speedup!

---

## 3. Strong Scaling vs. Weak Scaling

| Feature | Strong Scaling (Amdahl's Metric) | Weak Scaling (Gustafson's Metric) |
| :--- | :--- | :--- |
| **Problem Size ($W$)** | **Fixed** as $p$ increases | **Grows proportionally** with $p$ ($W/p = \text{const}$) |
| **Primary Goal** | Minimize time to solution for a fixed problem | Maximize simulation resolution in fixed time |
| **Typical Behavior** | Efficiency declines rapidly beyond a critical $p$ | Efficiency remains near constant up to massive $p$ |
| **Scaling Limit** | Dictated by Amdahl's serial fraction $f$ | Dictated by memory bandwidth and network latency |

```
    Speedup
       ▲
       │                                     Weak Scaling (S ~ p)
       │                                   .´
       │                                 .´
       │                               .´
       │                             .´
       │                           .´  Strong Scaling (Saturation curve)
       │                      ...─´
       │               .─-─´´
       │        .─-─´´
       │  .─-─´´
       └────────────────────────────────────────► Number of Processors (p)
```

---

## 4. The Isoefficiency Metric & Function

### 4.1 Concept and Purpose
Because efficiency drops as $p$ increases for a fixed $W$, and efficiency rises as $W$ increases for a fixed $p$, a parallel system can maintain a **constant target efficiency $E$** if the problem size $W$ is scaled upwards simultaneously with $p$.
The **Isoefficiency Function** $W = f(p)$ defines the exact rate at which problem size $W$ must grow as a function of $p$ to keep efficiency $E$ strictly constant.

### 4.2 Complete Mathematical Derivation (Slides 378–380)
From the fundamental relation:
$$E = \frac{1}{1 + \frac{T_o(W, p)}{W}}$$
Rearranging:
$$1 + \frac{T_o(W, p)}{W} = \frac{1}{E} \implies \frac{T_o(W, p)}{W} = \frac{1 - E}{E}$$
Inverting:
$$\frac{W}{T_o(W, p)} = \frac{E}{1 - E}$$
Define the constant $K = \frac{E}{1 - E}$:
$$\mathbf{W = K \cdot T_o(W, p)}$$

### 4.3 Isoefficiency Example: Adding $n$ Numbers
For our reduction example, $T_o = 2 p \log_2 p$:
$$W = K (2 p \log_2 p) \implies \mathbf{W = O(p \log p)}$$
- **Interpretation**: To maintain an efficiency of $E = 0.80$ while doubling the number of processors from $p$ to $2p$, the problem size $W$ only needs to increase by slightly more than double ($O(p \log p)$).
- **Scalability Ranking**:
  - $W = O(p)$: Perfectly scalable.
  - $W = O(p \log p)$: **Highly scalable** (e.g., FFT, Reduction).
  - $W = O(p^{1.5})$: Moderately scalable (e.g., Cannon’s Matrix Multiply).
  - $W = O(p^2)$ or higher: **Poorly scalable** (requires massive data explosion to maintain efficiency).

### 4.4 Worked Isoefficiency Problem from Lecture (Slide 380)
> **Exam Problem (Slide 380)**: Consider a hypothetical parallel program running on a given architecture with $p$ processors with an overhead function:
> $$T_o = p^{3/4} W^{3/4}$$
> 1. Determine the isoefficiency function relating $W$ and $p$.
> 2. If the problem size increases by a factor of $m$, by what factor must the number of processors be increased to maintain constant efficiency?

#### Analytical Solution
1. Using the isoefficiency relationship $W = k T_o(W, p)$ where $k = \frac{E}{1 - E}$:
   $$W = k \left( p^{3/4} W^{3/4} \right)$$
   Dividing both sides by $W^{3/4}$:
   $$W^{1/4} = k p^{3/4}$$
   Raising both sides to the 4th power:
   $$\mathbf{W = k^4 p^3 \implies W = O(p^3)}$$
2. Inverting to express processors as a function of workload:
   $$p = \left( \frac{W}{k^4} \right)^{1/3} \propto W^{1/3}$$
   Therefore, if the problem workload increases by a factor of $m$ ($W_{\text{new}} = m W$):
   $$p_{\text{new}} \propto (m W)^{1/3} = m^{1/3} W^{1/3} = \mathbf{m^{1/3} p}$$
   **Conclusion**: To maintain constant efficiency, the number of processors must be increased by a factor of **$m^{1/3}$**.

---

## 5. Amdahl's Law

Gene Amdahl (1967) established the theoretical upper bound on speedup for fixed-size workloads.

### 5.1 Step-by-Step Mathematical Derivation (Slides 381–382)
Let a serial program on data size $n$ run in time $T_s$, decomposed into two parts:
1. **$\Phi_s(n)$**: The inherently sequential component that cannot be parallelized (thread joining, barrier synchronization, sequential I/O).
2. **$\Psi(n)$**: The perfectly parallelizable component.
$$T_s = \Phi_s(n) + \Psi(n)$$

When executed on $p$ processors with overhead $T_o(n, p)$:
$$T_p = \Phi_s(n) + \frac{\Psi(n)}{p} + T_o(n, p)$$
The actual parallel speedup is:
$$S = \frac{T_s}{T_p} = \frac{\Phi_s(n) + \Psi(n)}{\Phi_s(n) + \frac{\Psi(n)}{p} + T_o(n, p)}$$

Since overhead $T_o(n, p) \ge 0$, the upper bound on speedup is:
$$S \le \frac{\Phi_s(n) + \Psi(n)}{\Phi_s(n) + \frac{\Psi(n)}{p}}$$

Define the fraction of the serial algorithm that cannot be parallelized as $f$ (Slide 382):
$$f = \frac{\Phi_s(n)}{\Phi_s(n) + \Psi(n)} \implies \Phi_s(n) + \Psi(n) = \frac{\Phi_s(n)}{f}$$
From which:
$$\frac{\Psi(n)}{\Phi_s(n)} = \frac{1}{f} - 1 \implies \Psi(n) = \Phi_s(n) \left(\frac{1}{f} - 1\right)$$

Substituting into the speedup bound:
$$S \le \frac{\frac{\Phi_s(n)}{f}}{\Phi_s(n) + \frac{\Phi_s(n) \left( \frac{1}{f} - 1 \right)}{p}} = \frac{\frac{1}{f}}{1 + \frac{1 - f}{f p}} = \frac{1}{f \left( 1 + \frac{1 - f}{f p} \right)} = \mathbf{\frac{1}{f + \frac{1 - f}{p}}}$$

$$\mathbf{S(p) \le \frac{1}{f + \frac{1 - f}{p}}}$$

### 5.2 The Asymptotic Speedup Limit
Taking the mathematical limit as the number of processors approaches infinity ($p \to \infty$):
$$\mathbf{S_{\max} = \lim_{p \to \infty} S(p) = \frac{1}{f + 0} = \frac{1}{f}}$$

#### Concrete Numerical Examples
- If **$5\%$ of a code is sequential** ($f = 0.05$):
  $$S_{\max} = \frac{1}{0.05} = \mathbf{20\times}$$
  Even on a supercomputer with $1,000,000$ cores, the maximum theoretical speedup can **never exceed 20**!
- If **$1\%$ of a code is sequential** ($f = 0.01$):
  $$S_{\max} = \frac{1}{0.01} = \mathbf{100\times}$$

### 5.3 Limiting Cases Analyzed in Lecture (Slide 383)
1. **Strictly Sequential Algorithm ($f \to 1$)**:
   $$S \le 1$$
   Because real systems incur overhead ($T_o > 0$), the parallel execution time is strictly worse than running on a single sequential processor ($S < 1$)!
2. **Negligible Sequential Component ($f \to 0$)**:
   $$S \le p$$
   Speedup is bounded strictly by communication overhead and hardware limits.

---

## 6. Gustafson-Barsis's Law (Scaled Speedup)

John Gustafson (1988) pointed out the fundamental flaw in Amdahl’s premise:
- Amdahl assumed that problem size $W$ remains fixed as we acquire larger supercomputers.
- In reality, scientists use more powerful supercomputers to solve **larger, higher-resolution problems**, not to run tiny problems in milliseconds!

### 6.1 Derivation of Scaled Speedup
Let $s$ be the fraction of time spent executing the sequential part on the **parallel machine** on a scaled problem:
$$T_p = s + (1 - s) = 1 \quad (\text{normalized parallel time})$$
To run this scaled problem sequentially on a single core, the sequential fraction requires time $s$, while the parallel fraction requires $p \times (1 - s)$:
$$T_s = s + p (1 - s)$$
The **Scaled Speedup** is:
$$S_{\text{scaled}}(p) = \frac{T_s}{T_p} = \frac{s + p (1 - s)}{1}$$
$$\mathbf{S_{\text{scaled}}(p) = p - (p - 1) s}$$

#### Comparison Example
Suppose $s = 0.05$ ($5\%$ serial) on $p = 1024$ processors:
- **Amdahl’s Law (Fixed Workload)**:
  $$S \le \frac{1}{0.05 + \frac{0.95}{1024}} \approx \frac{1}{0.05 + 0.00093} \approx \mathbf{19.64\times}$$
- **Gustafson’s Law (Scaled Workload)**:
  $$S_{\text{scaled}} = 1024 - (1023)(0.05) = 1024 - 51.15 = \mathbf{972.85\times}$$
**Core Lesson**: Parallelism is viable and massively scalable when problem sizes are scaled to match the compute capacity of the machine!

---

> [!IMPORTANT]
> **Past Exam Focus on Parallel Performance & Communication Overhead (2025 Midsem Q4)**:
> In the 2025 examination, a sequential program ran in $T_s = 80\text{ s}$ with $80\%$ of execution time parallelizable. On $p = 8$ processors, parallel runtime was $T_p(8) = 32\text{ s}$. Given that communication time scales linearly with processor count, find the parallel efficiency on $p = 16$ processors:
> 1. Serial component: $\Phi_s = 0.20 \times 80 = 16\text{ s}$; Parallel component: $\Psi = 64\text{ s}$.
> 2. Computation time on $8$ processors: $\frac{64}{8} = 8\text{ s} \implies T_{\text{comm}}(8) = 32 - (16 + 8) = 8\text{ s}$.
> 3. Linear scaling of communication: $T_{\text{comm}}(16) = \frac{16}{8} \times 8 = 16\text{ s}$.
> 4. Parallel runtime on $16$ processors: $T_p(16) = 16 + \frac{64}{16} + 16 = 36\text{ s}$.
> 5. Parallel efficiency: $E(16) = \frac{T_s}{16 \times T_p(16)} = \frac{80}{16 \times 36} = \frac{5}{36} \approx \mathbf{13.89\%}$.
>
> See the complete step-by-step derivation:  
> 👉 [**2025 Exam Q4 Detailed Solution**](file:///home/thesupercd/Documents/HPSC_TutorialsNAssignments/ExamNotes/Previous_Years_Question_Papers_Solutions.md#question-4-4-marks).
