# Chapter 14: Parallel Computer Architectures, Flynn's Taxonomy, Interconnects & Communication

---

## 1. Flynn's Taxonomy of Computer Architectures

Michael J. Flynn (1966) classified all computer architectures according to the multiplicity of hardware instruction streams and data streams:

```
                              Instruction Streams
                             Single         Multiple
                         ┌─────────────┬─────────────┐
                  Single │    SISD     │    MISD     │
     Data Streams        ├─────────────┼─────────────┤
                Multiple │    SIMD     │    MIMD     │
                         └─────────────┴─────────────┘
```

### 1.1 Single Instruction, Single Data stream (SISD)
- **Description**: The classic uniprocessor sequential architecture (Von Neumann).
- **Operation**: A single control unit decodes a single instruction stream from memory per unit time, operating on a single scalar data element in registers/memory.
- **Example**: Legacy single-core CPUs (e.g., Intel 8086, Pentium uniprocessors).

### 1.2 Single Instruction, Multiple Data stream (SIMD)
- **Description**: Vector processors and data-parallel accelerators.
- **Operation**: A single control unit fetches and decodes an instruction, and then **broadcasts that exact same instruction to multiple ALUs**, which execute it synchronously in strict lockstep on distinct, independent data elements.
- **Applications**:
  - Vector CPU extensions: Intel SSE, AVX-512, ARM Neon (processing 8 or 16 floats per cycle in vector registers).
  - Graphics Processing Units (GPUs): Streaming Multiprocessors (SMs) executing warp/wavefront instructions across thousands of threads.
  - Matrix-vector multiplication: Each processing element computes one row or one element of the vector dot product simultaneously.

### 1.3 Multiple Instruction, Single Data stream (MISD)
- **Description**: Multiple autonomous processing elements execute distinct instructions on the exact same shared stream of data.
- **Applications**: Rare in scientific computing.
  - Fault-tolerant aerospace systems: Redundant flight computers executing distinct algorithms on the same sensor data stream, using voting logic to detect component failure.
  - Systolic arrays and pipelined streaming data filters.

### 1.4 Multiple Instruction, Multiple Data stream (MIMD)
- **Description**: The dominant architecture for general high-performance parallel computing.
- **Operation**: Multiple independent processors operate asynchronously, each possessing its own control unit and ALU, decoding different instruction streams on different data streams.
- **Examples**: Multi-core microprocessors, distributed HPC clusters, Top500 supercomputers.

---

## 2. MIMD Architectural Classifications: Memory Organization

MIMD systems are categorized by how physical and logical memory are organized across processing elements.

```
       Shared Memory (Multiprocessor)        Distributed Memory (Multicomputer)
       ┌─────┐ ┌─────┐ ┌─────┐ ┌─────┐       ┌─────┐┌─────┐   ┌─────┐┌─────┐
       │CPU 0│ │CPU 1│ │CPU 2│ │CPU 3│       │CPU 0││Mem 0│   │CPU 1││Mem 1│
       └──┬──┘ └──┬──┘ └──┬──┘ └──┬──┘       └──┬──┘└─────┘   └──┬──┘└─────┘
          └───────┼───────┼───────┘             │                │
                  ▼                             ▼ Interconnect   ▼
       ┌─────────────────────────────┐       ──────────────────────────────
       │    Global Shared Memory     │         Network (Message Passing)
       │           (RAM)             │       ──────────────────────────────
       └─────────────────────────────┘          ▲                ▲
                                             ┌──┴──┐┌─────┐   ┌──┴──┐┌─────┐
                                             │CPU 2││Mem 2│   │CPU 3││Mem 3│
                                             └─────┘└─────┘   └─────┘└─────┘
```

### 2.1 Shared Memory Systems (Multiprocessors)
In a shared memory system, all CPUs share a single, unified physical address space:
- Any processor can directly read or write any memory location via conventional load/store instructions.
- Communication is **implicit**: Processor 0 writes data to address $X$, and Processor 1 subsequently reads from address $X$.

#### Sub-classifications
1. **Uniform Memory Access (UMA / Symmetric Multiprocessing - SMP)**:
   - All processors connect to central RAM via a shared bus or crossbar switch.
   - Every processor experiences identical access latency and memory bandwidth to all addresses.
   - Scalability limit: Typically limited to $16 - 64$ sockets due to physical bus contention.
2. **Non-Uniform Memory Access (NUMA)**:
   - Memory is physically distributed across individual processor sockets, but interconnected hardware provides a logically unified global address space.
   - Access to **local memory** (on the same socket) is extremely fast. Access to **remote memory** (on another socket across QPI/UPI links) has significantly higher latency and lower bandwidth.
   - Programmers must enforce NUMA data affinity to prevent performance degradation.

### 2.2 Distributed Memory Systems (Multicomputers)
In a distributed memory system, each node consists of an autonomous processor and private local memory:
- No shared global address space exists. Processor 0 cannot directly address or load memory from Processor 1.
- Communication is **explicit**: Data exchange occurs exclusively by sending and receiving discrete packets across an interconnection network via **Message Passing** (e.g., MPI - Message Passing Interface).

### 2.3 Hybrid Architecture: Distributed-Shared Memory
Modern supercomputing clusters (including PARAM Shakti at IIT Kharagpur, Frontier, and Summit) combine both paradigms:
- **Node-Level**: Multi-socket, multi-core CPUs and GPUs sharing local node RAM (Shared Memory / NUMA).
- **Cluster-Level**: Thousands of nodes interconnected by high-speed, low-latency fabric like InfiniBand or Slingshot (Distributed Memory).
- **Hybrid Programming Model**: MPI across nodes + OpenMP / CUDA within each node.

> [!IMPORTANT]
> **Past Exam Focus (2025 Midsem Q5: Shared vs Distributed Memory Comparison)**:
> In the 2025 exam, students were asked to contrast Shared-Memory and Distributed-Memory architectures. For an exam-ready comparative table across memory space, interconnect, communication, scalability, bottlenecks (cache coherency vs latency), and APIs (OpenMP vs MPI):  
> 👉 [**2025 Exam Q5 Detailed Solution**](file:///home/thesupercd/Documents/HPSC_TutorialsNAssignments/ExamNotes/Previous_Years_Question_Papers_Solutions.md#question-5-2-marks).

### 2.4 Comparative Summary: Shared vs. Distributed Memory
| Architectural Feature | Shared Memory (Multiprocessor) | Distributed Memory (Multicomputer) |
| :--- | :--- | :--- |
| **Address Space** | Single, globally unified address space | Multiple, disjoint private address spaces |
| **Communication Mechanism** | Implicit (Loads, Stores, Locks, Mutexes) | Explicit Message Passing (`MPI_Send`, `MPI_Recv`) |
| **Hardware Complexity** | High (Complex cache-coherency logic) | Low (Commodity compute nodes + network) |
| **Scalability Limit** | Moderate ($\le 64$ cores before bus saturates) | Massive ($> 1,000,000$ cores) |
| **Programming Difficulty** | Easier for basic tasks (OpenMP directives) | Higher (Requires explicit domain decomposition) |
| **Major Bottlenecks** | Cache coherency, false sharing, bus contention | Network latency, message startup overhead ($t_s$) |

---

## 3. Cache Coherency and False Sharing in Shared Memory

### 3.1 The Cache Coherency Problem (Slides 458–459)
Because each processor in an SMP/NUMA system possesses private L1 and L2 caches, multiple cached copies of the exact same memory address can exist simultaneously.

#### The Hazard: Concurrent Inconsistent Updates (Slide 458)
Consider an initial shared variable $x = 5$ stored in main memory:

```
        Processor-1                     Cache-1
        [ x = x + 3 ] ───────────────► [ x = 8  ] ────┐
                                                      │ (Write: x = 8)
                                                      ▼
                                                ┌─────────────┐
                                                │ Main Memory │
                                                │   x = 5     │
                                                └─────────────┘
                                                      ▲
                                                      │ (Write: x = 10?)
        Processor-2                     Cache-2       │
        [ x = x + 5 ] ───────────────► [ x = 10 ] ────┘
```
1. Processor-1 reads $x = 5$ into Cache-1, performs $x \leftarrow x + 3 = 8$, and writes $x = 8$ to memory.
2. Simultaneously, Processor-2 reads $x = 5$ into Cache-2, performs $x \leftarrow x + 5 = 10$, and attempts to write $x = 10$ to memory.
3. **The Catastrophic Failure**: Depending on bus arbitration timing, either $x = 8$ or $x = 10$ overwrites memory, completely losing the other arithmetic update (the mathematically correct sequential result is $5 + 3 + 5 = 13$)! Furthermore, each processor's local cache now contradicts main memory.

#### The Solution: Hardware Serialization via Coherency (Slide 459)
To prevent inconsistent states, hardware controllers serialize memory operations:

```
    Sequential Serialization:
    1. Processor-1 executes: x = x + 3 = 8.
       Cache-1 updates to x = 8; broadcasts update/invalidation to memory (x = 8).
    2. Processor-2 stalls / waits until cache update is completed across the bus.
    3. Processor-2 resumes, fetches latest value x = 8 into Cache-2.
    4. Processor-2 executes: x = 8 + 5 = 13.
       Cache-2 and Main Memory are updated to x = 13!
```
Processors must operate sequentially on shared memory addresses to guarantee coherency, enforced via **Snoopy Cache** (broadcast bus) or **Directory-Based** protocols!

### 3.2 Coherency Protocols
Hardware controllers enforce coherency through two primary mechanisms:
1. **Snooping Protocols (Bus-Based Systems)**:
   All cache controllers constantly monitor ("snoop") the common broadcast bus.
   - **Write-Invalidate (Most Common)**: When Processor A writes to $X$, it broadcasts an invalidation signal. All other caches holding $X$ change its status to **Invalid**. If Processor B subsequently accesses $X$, it suffers a cache miss and re-fetches the updated value from A.
   - **Write-Update**: Processor A broadcasts the newly written data, updating all caches holding $X$.
   - **MESI Protocol**: State machine tracking cache lines: **M**odified (dirty, exclusive), **E**xclusive (clean, only in this cache), **S**hared (clean, present in multiple caches), **I**nvalid.
2. **Directory-Based Protocols (NUMA / Scalable Systems)**:
   Because broadcasting on a bus does not scale to hundreds of nodes, a centralized or distributed directory tracks which nodes hold copies of each memory block. Point-to-point invalidation messages are sent only to the specific nodes sharing the line.

### 3.3 The False Sharing Pathology
> **Definition**: **False sharing** occurs when two distinct processors independently modify completely different logical variables that happen to reside within the **same 64-byte physical cache line**.

```
       Physical 64-Byte Cache Line in Memory:
       ┌────────────────────────────────────────────────────────┐
       │  A[0]  │  A[1]  │  A[2]  │  A[3]  │ ... │  A[7]        │
       └────────────────────────────────────────────────────────┘
           ▲        ▲
           │        │
         CPU 0    CPU 1
       (writes) (writes)
```

- **Mechanism of Performance Collapse**:
  - Thread 0 on CPU 0 is assigned to update scalar $A[0]$.
  - Thread 1 on CPU 1 is assigned to update scalar $A[1]$.
  - Logically, there is zero data sharing or race condition.
  - However, because $A[0]$ and $A[1]$ sit on the same 64-byte line:
    1. CPU 0 writes to $A[0]$ $\implies$ marks the entire 64-byte line **Invalid** in CPU 1's cache!
    2. CPU 1 tries to write to $A[1]$ $\implies$ encounters a cache miss, stalls, and invalidates CPU 0's cache line!
    3. The cache line bounces endlessly back and forth across the interconnect bus ("ping-ponging" / cache thrashing).
- **HPC Solution**: Pad variables or allocate thread-local variables with at least $64\text{ bytes}$ of padding (`alignas(64)`), ensuring distinct threads never write to the same cache line.

---

## 4. Interconnection Networks & Topological Metrics

In distributed multicomputers, nodes are connected via an interconnection network.

### 4.1 Key Topological Metrics
In an exam, evaluating network topologies requires computing 5 quantitative metrics:
1. **Node Degree ($d$)**: The number of physical communication links incident on a single node. Determines hardware complexity and router pin count.
2. **Diameter ($D$)**: The maximum shortest-path distance (number of hops) between any pair of nodes in the network. Represents the worst-case communication latency.
3. **Bisection Width ($B_w$)**: The minimum number of physical links that must be severed to partition the network into two equal halves with $\lfloor p/2 \rfloor$ nodes each.
4. **Bisection Bandwidth**: The minimum data transfer rate (in GB/s) across the bisection cut ($B_w \times \text{link bandwidth}$). Determines vulnerability to global communication bottlenecks.
5. **Network Cost (Link Count, $C$)**: The total number of communication links in the network.

---

### 4.2 Standard Network Topologies and Metric Comparisons

```
    Ring (p=8)             2D Mesh (4x4, p=16)             Hypercube (d=3, p=8)
      0 ─── 1                0 ─── 1 ─── 2 ─── 3                 011 ───── 111
    /       \                │     │     │     │                / │       / │
   7         2               4 ─── 5 ─── 6 ─── 7               001 ───── 101│
   │         │               │     │     │     │               │  010 ───│─ 110
   6         3               8 ─── 9 ─── 10─── 11              │ /       │ /
    \       /                │     │     │     │               000 ───── 100
      5 ─── 4                12─── 13─── 14─── 15
```

| Topology | Node Degree ($d$) | Diameter ($D$) | Bisection Width ($B_w$) | Total Link Count ($C$) |
| :--- | :---: | :---: | :---: | :---: |
| **Linear Array** | 2 | $p - 1$ | 1 | $p - 1$ |
| **Ring** | 2 | $\lfloor p/2 \rfloor$ | 2 | $p$ |
| **2D Mesh ($\sqrt{p} \times \sqrt{p}$)** | 4 | $2(\sqrt{p} - 1)$ | $\sqrt{p}$ | $2(p - \sqrt{p})$ |
| **2D Torus (Wrap-around)** | 4 | $2 \lfloor \sqrt{p}/2 \rfloor$ | $2\sqrt{p}$ | $2p$ |
| **Hypercube ($d$-cube, $p = 2^d$)**| $d = \log_2 p$ | $d = \log_2 p$ | $p/2$ | $\frac{p}{2} \log_2 p$ |
| **Complete Binary Tree** | 3 (root has 2) | $2 \log_2 p$ | 1 (severe bottleneck) | $p - 1$ |
| **Fat-Tree** | Variable | $2 \log_2 p$ | $O(p)$ (no bottleneck!) | $O(p \log p)$ |
| **Fully Connected (Crossbar)** | $p - 1$ | 1 | $p^2/4$ | $\frac{p(p-1)}{2}$ |

- **Hypercube Property**: In an $n$-dimensional hypercube, two nodes are connected if and only if their binary address representations differ by **exactly one bit**. (E.g., node $000$ connects to $001, 010, 100$).
- **Fat-Tree Principle**: In a regular binary tree, the bisection width is 1 (the single link at the root), causing massive traffic jams. A **Fat-Tree** increases link capacity and channel width progressively as one ascends toward the root, maintaining constant bisection bandwidth.

---

## 5. Analytical Communication Cost Models

When two processes communicate across a network, message delivery incurs latency.

### 5.1 Three Components of Communication Time (Slide 356)
1. **Startup Time (Latency, $t_s$)**:
   - The time required to prepare the message at the sending node (packetization, appending headers/checksums, routing algorithm execution, initializing network DMA) plus unpacking at the receiver.
   - Independent of message size.
2. **Per-Hop Time (Switch Latency, $t_h$)**:
   - The time required for the message header to travel between two directly connected router switches across a single network link.
3. **Per-Word Transfer Time ($t_w$)**:
   - The time required to transmit a single word of data across a link.
   - Governed by link bandwidth $r$ (words/second):
     $$t_w = \frac{1}{r}$$
   - Transmitting a message of $m$ words takes $m \cdot t_w = m / r$ seconds.

### 5.2 Store-and-Forward vs. Cut-Through (Wormhole) Routing
1. **Store-and-Forward Routing**:
   - Each intermediate switch along the path must receive and buffer the entire message of $m$ words before forwarding it to the next hop.
   - Total time over $l$ links (hops):
     $$\mathbf{T_{comm} = t_s + (l \times m \cdot t_w)}$$
   - *Severe Disadvantage*: Communication time scales multiplicatively with distance $l$.

2. **Cut-Through / Wormhole Routing (Modern Standard)**:
   - A message is divided into tiny flow control units called **flits** (a few bytes).
   - As soon as the header flit arrives at an intermediate switch and determines the outgoing port, it is immediately forwarded. The remaining data flits follow in a pipelined stream.
   - Total time over $l$ links for an $m$-word message:
     $$\mathbf{T_{comm} = t_s + l \cdot t_h + m \cdot t_w}$$
   - In modern optical interconnects (InfiniBand), the switch latency is in nanoseconds ($t_h \approx 10^{-8}\text{ s}$), while startup latency is in microseconds ($t_s \approx 10^{-6}\text{ s}$).
   - Because $l \cdot t_h \ll t_s$, the per-hop latency is completely negligible!
   - **Canonical Analytical Model for HPC Communication**:
     $$\mathbf{T_{comm} \approx t_s + m \cdot t_w}$$
     Communication time is effectively **distance-independent**!

---

## 6. Parallel Random Access Machine (PRAM) Models

The **PRAM** is the theoretical foundation for analyzing parallel algorithmic complexity. It assumes $p$ processors connected to an ideal shared memory where all memory locations can be accessed in **unit time ($O(1)$)**.

### 6.1 The Four PRAM Variants
1. **Exclusive Read, Exclusive Write (EREW)**:
   - No two processors can read or write to the same memory cell simultaneously.
   - The most restrictive model.
2. **Concurrent Read, Exclusive Write (CREW)**:
   - Multiple processors can read the same memory location simultaneously.
   - Only one processor can write to a given location at any time.
3. **Exclusive Read, Concurrent Write (ERCW)**:
   - Multiple writes allowed, but reads must be exclusive (rarely used).
4. **Concurrent Read, Concurrent Write (CRCW)**:
   - Multiple processors can concurrently read and write to the same memory cell.

### 6.2 Conflict Resolution Rules for CRCW PRAM
When multiple processors attempt to write to the same cell simultaneously, the conflict is resolved by one of three protocols:
1. **Common CRCW**: Concurrent writes are permitted only if all processors write the **exact same value**. If values differ, an error occurs.
2. **Arbitrary CRCW**: One arbitrary processor's write succeeds, and its value is stored; all other writes are ignored.
3. **Priority CRCW**: Processors are assigned unique priority IDs (e.g., processor rank). The processor with the **highest priority** successfully writes its value.
