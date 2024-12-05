# [Computing Approximate Centerpoints in Polynomial Time](https://ieeexplore.ieee.org/abstract/document/10756130)

Yeshwanth Cherapanamjeri, Massachusetts Institute of Technology.



### Problem Settings

Given a set of points, a centerpoint is a point such that any hyperplane passing through it divides the point set into roughly balanced pieces.

Formally, given a point set, $X = \{x_i \in \mathbb{R}^d\}^n_{i=1}$, the centerpoint is any point $x$ which satisfies the following: 
$\forall v \in \mathbb{R}^d : \text{Depth}(X,x,v) \geq \frac{1}{d+1}$ where 
$Depth(X,x,v) := \min \{\frac{1}{n}\sum_{i=1}^n\mathbb{1}\{⟨x_i,v⟩ \leq ⟨x,v⟩\}, \frac{1}{n}\sum_{i=1}^n\mathbb{1}\{⟨x_i,v⟩ ⩾ ⟨x,v⟩\}\}$.

You can understand the conception in the one-dimensional setting easily, because the centerpoint in 1-dimension just divides all points into fifty-fifty, just the median. 

A natural question is that can we get a point with depth $1/2$ for a higher dimension? The answer is NO. Helly proves that $\frac{1}{d+1}$ is the best achievable in the worst case in 1923.



As a consequence of Helly’s theorem in convex geometry, such a point actually exists for all point sets. Then the problem is that how can we efficiently compute the centerpoints.



### History

The best algorithms for this problem only guarantee a point with depth $Ω(1/d^2)$. (Approximating center points with iterated radon points, 1993)



### Preliminaries: approximate centerpoints

Define the directional quantiles, $\forall ||v||=1:q_v:=\min\{x:\#\{i:⟨x_i,v⟩ \leq x\}\geq \Omega(\frac{n}{d+1})\}$. Then a valid point with depth $\Omega(1/d)$ can be defined by $\forall ||v||=1:⟨x,v⟩ \geq q_v$. For error parameter $\varepsilon \gt 0$, an approximate centerpoint will satisfy $\forall ||v||=1:⟨x,v⟩ \geq q_v-\varepsilon$.



### Results

There exists an algorithm computing an $\varepsilon$-approximate centerpoint in time poly($n,d,1/\varepsilon$).



### Approach

1. Compute an approximate gradient
2. Use the approximate gradient to update the estimate

#### A simple example to show how an iteration works

Suppose $0$ is a point of low depth, and let $v$ be a direction with low depth such that $\#\{i:⟨X_i,v⟩<0\}\leq\frac{n}{4d}$. We should compute an approximate gradient $\hat{g}$ such that $⟨\hat{g},v⟩\geq\gamma>0$ for $\gamma$ bounded away from $0$ for all such $v$.

The writer uses **Radial Isotropic Transformation** to compute $\hat{g}$.

We say that a linear transformation $R : \mathbb{R}^d \rightarrow \mathbb{R}^d$ places a set of points $u_1,\dots,u_n \in \mathbb{R}^d$ in Radial Isotropic Position if $\frac{1}{n}\sum^n_{i=1}\frac{Ru_i}{||Ru_i||}⊗ \frac{Ru_i}{||Ru_i||} = \frac{I}{d}$.

They prove that the Radial Isotropic Transformation exists almost everywhere and can be (approximately) computable in polynomial time.

Then $\hat{g}=\frac{1}{n}\sum_{i=1}^n\frac{x_i}{||Rx_i||}$ is the required approximate gradient.

The process is like gradient descent.



### Open questions

1. Can the poly($1/\varepsilon$) be improved to logarithmic dependency?
2. Beyond the worst-case analysis, can we get better depth guarantees for non-worst case point sets?





# [Near-Optimal $(1+\varepsilon)$-Approximate Fully-Dynamic All-Pairs Shortest Paths in Planar Graphs](https://eprints.cs.univie.ac.at/8169/1/Planar_dynamic_APSP_FOCS_2024_Filtser_Goranci_Patel_Probst_Gutenberg.pdf)

Arnold Filtser, Bar Ilan University Ramat Gan, Israel
Neel Patel, University of Southern California, USA
Gramoz Goranci, University of Vienna, Austria
Maximilian Probst Gutenberg, ETH Zurich, Switzerland



### Problem Setting (fully dynamic all-pair shortest path problem, fully dynamic APSP)

Given an $n$-vertex undirected weighted planar graph $G = (V,E,w)$ undergoing edge insertions and deletions, find a structure and corresponding algorithm to support distance queries about the shortest path between any source-target vertex pair as quickly as possible.



### Related Work and History

#### Exact Algorithms

1997, Faster shortest-path algorithms for planar graphs, $SSSP$ algorithm, O(n) for answer queries.

SODA 2021, Planar distance oracles with better time-space tradeoffs, update time $n^{1+o(1)}$, query time $\tilde{O}(1)$.

ICALP 2018, Improved bounds for shortest paths in dense distance graphs, update and query time $O(n^{2/3}\frac{\log^{5/3}n}{\log^{4/3}\log n}$).

FOCS 2016, Popular conjectures as a barrier for dynamic planar graph algorithms, lower bound $O(n^{1/2-\varepsilon})$.

#### Approximation Algorithm

1998, A fully dynamic approximation scheme for shortest paths in planar graphs, $\tilde{O}(n^{2/3})$ update and query time.

STOC 2012, A deterministic near-linear time algorithm for finding minimun cuts in planar graphs, $\tilde{O}(n^{1/2})$ update and query time.



### Result

Given an $n$-vertex undirected planar graph $G = (V,E,w)$ with weights in $[1,W]$ undergoing edge insertions and deletions that preserve the embedding, and a precision parameter $1/\text{polylog}(n) < \varepsilon < 1$, there is a deterministic data structure that supports edge updates and queries for the $(1+\varepsilon)$-approximate distance between any pair $s,t$ of vertices in the current graph $G$. The data structure can be initialized in  $\tilde{O}(n\log W)$ time and achieves $n^{o(1)} \log W$ amortized update and query time.



#### Def 1

Given a planar graph $G = (V,E,w)$ and a subset of terminals $T \subseteq V$, we define an instance as a tuple $(G,T)$. We say that an instance $(H,T)$ is an $\varepsilon$-emulator of $(G,T)$ if $H = (V_H,E_H,w_H)$ is a planar graph over a set of vertices $V_H$ containing the terminals $T \subseteq V_H$, that preserves distances between terminals up to a small multiplicative factor: $\forall u,v \in T : d_G(u,v) \leq d_H(u,v) \leq (1+\varepsilon)\cdot d_G(u,v)$.

Several work developed algorithms to construct an $\varepsilon$-emulator of any given instance $(G,T)$ denoted as $(H,T)$ of size $|V (H)| ≤ O (\frac{|T|}{\varepsilon^{O(1)}})$ which is independent of the size of the underlying planar graph $G$ in $\tilde{O}(\frac{n}{\varepsilon^{O(1)}})$.




#### Key 1 Initialize the $\varepsilon$-emulator

Divide and conquer to initialize the $\varepsilon$-emulator for the instance $(G, T=\emptyset)$. Use a decomposition tree $\mathcal{T}$.

Each node in $\mathcal{T}$ is an instance. All nodes are divided into four categories. For an instance $(R,S)$,

##### Type-1: $|V(R)| \leq \frac{n}{\tau}$ and all terminals in $S$ are on a single face of $R$. (no hole, a leaf node)

##### Type-2: $|V(R)| \leq \frac{n}{\tau}$ and one hole.

##### Type-3: $|V(R)| \leq \frac{n}{\tau}$ and a constant number of holes.

##### Type-4: $|V(R)|>\frac{n}{\tau}$ and a constant number of holes. $(G, T=\emptyset)$ is a Type-4 instance.

Here, a hole means its a face of $R$ but not a face of the origin graph $G$.

##### Bound the depth of $\mathcal{T}$ and the number of copies: use the decomposition relationship between these four types. 

The result is that the depth of $\mathcal{T}$ is bounded by $O(\log n)$ and each edge can belong to at most $O(\log^2 n)$ different instance. 

#### Key 2 Maintain the $\varepsilon$-emulator



### Open problem

1. For the approximation algorithm, whether we can improve subpoly update to actual polylog update time.

2. Whether we can obtain $\sqrt{n}$ APAS exact distance oracles.



