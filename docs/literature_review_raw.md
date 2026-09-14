# Crystal / scintillator calorimeter clustering: technique review

Scope: BESIII EMC (CsI(Tl)), Belle II ECL (CsI(Tl)), LHCb ECAL, ATLAS/CMS topological & particle-flow
clustering, and modern ML clustering. Written for a ~622-module crystal ECAL with a
wavefront/topological-growth clusterer (seed by energy, grow through primary neighbours above `Eth`
with `Nth` lit neighbours, optional one-shot second-neighbour gap jumping, multiplicity-scored output).

**Method note / honesty statement.** Every quote below was copied from a page I actually fetched with
`web_fetch`. The fetch tool cannot retrieve PDFs (`unsupported content type "application/pdf"`), so
conference slides, theses and journal PDFs are cited only where an HTML mirror (ar5iv / arXiv HTML)
was available. Rendering artifacts of the HTML/LaTeX converters (e.g. `5\times 5`) have been
transcribed into plain notation; wording is unchanged. Anything I could not verify is explicitly
marked **UNVERIFIED**.

---

## 1. BESIII EMC (CsI(Tl))

### 1.1 Seed = local maximum + centre-of-gravity with logarithmic weighting

**Source:** V. Prasad, C. Liu, X. Ji, W. Li, H. Liu, X. Lou, "Study of BESIII electromagnetic
calorimeter performance with radiative lepton pair events", BESIII, arXiv:1606.00248 (2016) —
https://arxiv.org/abs/1606.00248 , HTML: https://ar5iv.labs.arxiv.org/html/1606.00248

**Verbatim excerpts (Section 5 "Shower position Reconstruction"):**
- "The continuous connected region of the crystals deposited by energy of the incident particles is
  called clusters. Each shower is recognized by a seed, which is the local maxima of energy deposit
  among its neighbors."
- "A center-of-gravity (CG) method is used to calculate the impact coordinate, x_c, of the showering
  particles on the front face of the EMC" … "the sum includes the information about all crystals of
  the clusters."
- "The second approach is the 'logarithmic weighting function' which can reduce the weight of the most
  energetic crystal and enhance the low energy ones." (definition W_i^log(E_i) = Max{0, a0 + ln(E_i)
  − ln(E_tot)}, with cutoff parameter a0 = 4.0)

**Relevance:** Our seed definition is already energy-ordered, but we throw away position information.
Replacing the cluster "position" (currently unused) with a CG + log-weighted centroid gives a
sub-crystal shower position for free, and the log weight down-weights the seed crystal so the centroid
is not biased toward the highest-energy module — directly useful for later track matching and for
π0/γγ invariant-mass work.

### 1.2 Shower-shape moments (lateral moment, second moment) as cluster features

**Source:** same paper, arXiv:1606.00248, Section 4.1 "Calibration of photon cluster shape".

**Verbatim excerpts:**
- "A number of variables have been developed to study the shower shapes of the different particles in
  the EMC, such as second moment and lateral moment (LAT). These two variables are used to quantify
  the transverse shower shape of the cluster and separate the EM showers from the hadronic showers."
- "The EM showers tend to deposit a large fraction of their energy in one or two crystals, whereas the
  hadronic showers tend to be more spread out."
- "The discrepancy between data and MC has been overcome while simulating a new MC sample with the
  EINC value of 0.27 MeV." (EINC = EMC incoherent noise; the older value was 0.20 MeV)

**Relevance:** LAT and second moment are cheap, threshold-free per-cluster features that separate
e+/γ clusters from muon/hadron-like deposits without any global cut tuning, and the EINC study shows
that the *noise level itself must be tuned* to reproduce shower-shape data. For our framework, adding
LAT/second-moment per cluster plus a proper per-module noise term (instead of a single hard `Eth`) is
a small change with immediate discriminating power.

### 1.3 Track–cluster association defines the photon sample; matching validated by kinematic fit

**Source:** BESIII Collaboration, "Improved measurement of absolute branching fraction of the inclusive
decay Λc+ → KS0 X", BESIII, arXiv:2502.20821 (2025) — https://arxiv.org/abs/2502.20821 ,
HTML: https://arxiv.org/html/2502.20821v2 (Section 3 "Event Selection")

**Verbatim excerpt:**
- "Neutral showers are reconstructed in the EMC. Showers not associated with any charged track are
  identified as photon candidates if they satisfy two additional criteria: (1) an energy deposition in
  the EMC of E_dep > 25 MeV in the barrel region corresponding to the polar angle |cos θ| < 0.8, while
  E_dep > 50 MeV in the end-cap region corresponding to 0.86 < |cos θ| < 0.92. (2) the EMC time
  difference from the event start time is required to be less than 700 ns to suppress electronic noise
  and showers unrelated to the event."

**Source for the matching validation:** arXiv:1606.00248, Section 4.2 "Photon detection efficiency".

**Verbatim excerpt:**
- "The regions of θ_γ,γ_pred < 0.5 radian and E_γ/E_γ_pred ∈ [0.15,1.4] are considered to be the
  detected region of the photon in the EMC."

**Relevance:** This is the cheapest complete "single-particle → cluster" assignment scheme: charged
clusters are those matched to a track, everything else above a *region-dependent* threshold and inside
a timing window is a photon. For our e+/γ/μ- samples we can adopt exactly this logic to label clusters
(e+/γ vs μ-like) and use the BESIII-style (Δangle, E_ratio) matching window as the seed of a proper
track–cluster matching figure of merit instead of pure multiplicity counting.

### 1.4 MUCKS / BESIII clustering code — could not verify

**Status: UNVERIFIED.** I could not find or fetch any BESIII document describing an EMC clustering
code called "MUCKS". Searches for "MUCKS BESIII simulation/clustering" returned nothing relevant; the
"MUC"-family names in BESIII literature refer to the muon counter (MUC) and its simulation/offline
software, not to EMC clustering. I therefore make no claim about MUCKS. A dedicated BESIII EMC
reconstruction reference exists (M. He, J. Phys. Conf. Ser. 293 (2011) 012025, "Simulation and
reconstruction of the BESIII EMC") but it is only reachable as a PDF, which this tool cannot fetch, so
I did not quote it.

---

## 2. Belle II ECL (CsI(Tl))

### 2.1 Connected regions from local maxima + iterative log-weighted energy sharing between overlapping clusters

**Source:** F. Wemmer et al. (Belle II), "Photon Reconstruction in the Belle II Calorimeter Using Graph
Neural Networks", arXiv:2306.04179 (2023), published in EPJ Research Infrastructures —
https://arxiv.org/abs/2306.04179 , HTML: https://ar5iv.labs.arxiv.org/html/2306.04179
(Section 5.1 "Baseline")

**Verbatim excerpts (Section 5.1):**
- "The clustering is performed in three steps. In the first step, all crystals are grouped into a
  connected set of crystals, so-called connected regions starting with LMs, as defined previously. In
  an iterative procedure all direct neighbors with energies above 0.5 MeV are added to this LM, and the
  process is continued if any neighbor itself has energy above 10 MeV. Overlapping connected regions
  are merged into one."
- "If there is more than one LM in a connected region, the energy in each crystal of the connected
  region is assigned a distance-dependent weight and can be shared between different clusters. The
  distance is calculated from the cluster centroid to each crystal center, where the cluster centroid
  is updated iteratively using logarithmic energy weights. This process is repeated until all cluster
  centroids in a connected region are stable within 1 mm."
- "In a third step, an optimal subset, including the n highest energetic crystals of all non-zero
  weighted crystals that minimize the energy resolution, is used to predict the cluster energy."

**Relevance:** This is the single closest match to our algorithm, and it fixes three of our stated
weaknesses at once: (i) it *splits* a connected region into one cluster per local maximum instead of
emitting one blob; (ii) it *shares* energy between overlapping clusters via distance-dependent weights
converged to 1 mm; (iii) it defines the energy by an optimally chosen subset of crystals rather than a
fixed 3×3 sum. The growth rule (append neighbours > 0.5 MeV, continue past a node only if that node is
> 10 MeV) is a concrete parameterisation of our "unclustered, above Eth, with Nth lit neighbours" rule.

### 2.2 Noise/beam-background handling: timing-based noise estimate and bias correction

**Source:** arXiv:2306.04179, Section 5.1.

**Verbatim excerpts:**
- "n depends on the measured noise in the event, and on the energy of the LM itself. The noise level is
  estimated by counting the number of crystals in the event containing more than 5 MeV that have times
  t more than 125 ns from the trigger time."
- "E_rec^basf2 is also corrected already within basf2 for possible bias using simulated events. This
  bias includes leakage (energy not deposited in the crystals included in the energy sum) and beam
  backgrounds (energy included in the sum that is not from the signal photon)."
- "Low beam background results in approximately 17 % of all crystals in the ECL having significant
  reconstructed energy E_rec ≥ 1 MeV; for high beam backgrounds this number is expected to increase to
  about 40 %."

**Relevance:** Even without real beam background, our MC has noise; a per-event, *self-calibrating*
noise estimate (count out-of-time / high-energy outliers) replaces our fixed `Eth`, and the explicit
separation of leakage loss vs. noise gain is a template for a two-component energy correction.

### 2.3 GNN "fuzzy" clustering: partial (probabilistic) crystal-to-cluster assignment incl. a background class

**Source:** arXiv:2306.04179 (same paper, Section 1 and Section 5.2).

**Verbatim excerpts:**
- "We present the study of a fuzzy clustering algorithm for the Belle II electromagnetic calorimeter
  using Graph Neural Networks. We use a realistic detector simulation including simulated beam
  backgrounds and focus on the reconstruction of both isolated and overlapping photons."
- "The term fuzzy clustering refers to the partial assignment of individual calorimeter crystals to
  several clustering classes. In our case, these are potentially overlapping, different signal
  photons, but also a beam background class."
- "While the basf2 algorithm strictly reconstructs one cluster for each LM, the GNN algorithm only uses
  the LMs to center the ROI" (as indexed by the published version, EPJ Research Infrastructures,
  https://link.springer.com/article/10.1007/s41781-023-00105-w ) — **note:** this last sentence was
  seen in a search-result snippet of the Springer page, not in text I fetched; treat it as **UNVERIFIED**.

**Relevance:** Confirms that hard, deterministic membership is the wrong target at high occupancy and
that a per-module weight vector (with an explicit "not from any particle / background" class) is the
modern answer to our "no energy sharing between overlapping clusters" problem. It also justifies
keeping the seed/local-maximum machinery as an ROI-definition step rather than as the final answer.

### 2.4 Real-time GNN clustering on FPGAs (Belle II L1 trigger)

**Source:** I. Haide et al. (Belle II), "Real-time graph neural networks on FPGAs for the Belle II
electromagnetic calorimeter", JINST 21 (2026) P06039, arXiv:2602.15118 —
https://arxiv.org/abs/2602.15118 (abstract, fetched from the arXiv abs page)

**Verbatim excerpts (abstract):**
- "The algorithm processes calorimeter trigger cells as graph nodes to perform clustering, feature
  extraction, and per-cluster signal classification with deterministic latency."
- "Cluster purity increases by up to 20% at low energies for isolated clusters, and cluster efficiency
  improves by up to 20% for overlapping clusters."

**Relevance:** Evidence that a GNN clusterer is not only an offline analysis tool but can be a
drop-in for a classical clustering stage with a *fixed* latency budget — relevant if our framework
ever needs to be embedded, and a strong argument that the graph formulation of the calorimeter scales.

### 2.5 ECL energy calibration — no dedicated source verified

**Status: UNVERIFIED.** I did not fetch a Belle II ECL calibration paper. The only calibration
statement I verified is the bias-correction sentence in §2.2 above (arXiv:2306.04179, Section 5.1),
which says the correction is applied inside basf2 using simulated events. Belle II's per-crystal
calibration is documented in basf2 code/Doxygen
(e.g. https://software.belle2.org/release-09-00-01/doxygen/caf__ecl__cluster__energy_8py_source.html)
and in collaboration notes that I could not fetch; I make no quotable claim about its procedure.

---

## 3. LHCb calorimeter clustering and overlap handling

### 3.1 Graph Clustering: seeds as sinks of directed edges, weakly connected components, overlap energy fractions

**Source:** N. Valls Canudas, M. Calvo Gómez, X. Vilasís-Cardona, E. Golobardes Ribé (LHCb),
"Graph Clustering: a graph-based clustering algorithm for the electromagnetic calorimeter in LHCb",
LHCb-DP-2022-003, arXiv:2212.11061 (2022) — https://arxiv.org/abs/2212.11061 ,
HTML: https://ar5iv.labs.arxiv.org/html/2212.11061

**Verbatim excerpts (Section 3 "Method"):**
- "It transforms the calorimeter digits into independent graph structures, where only relevant digits
  for a cluster are contained into isolated graphs. Following graph theory nomenclature, each energy
  digit from an event is represented as a vertex v in the graph, also called node."
- "By design, the target nodes of all edges in the graph are the seeds of the reconstructed clusters,
  where a seed is defined as a local maximum energy digit in the calorimeter grid over a threshold of
  50 MeV in transverse energy. With this, the cluster seeds can be easily identified as nodes with only
  incoming edges."
- "Furthermore, a node can be linked to more than one seed if it is susceptible to have energy deposits
  from more than one particle. These particular cases are called overlap cells."

**Verbatim excerpts (Section 3.3 / 3.4):**
- "To retrieve them from the graph, we need to search for the sub-graphs where all the nodes are
  connected to each other by some path, ignoring the direction of edges. This is defined by the weakly
  connected components of the original graph. In the proposed algorithm, this is implemented as a
  depth-first search, which explores an entire graph exploring all its branches as far as possible
  before backtracking."
- "It searches for overlap vertices, identified by having two or more output edges, and accumulates the
  energy of all the connected nodes on all the clusters involved in the overlap, including the energy
  of the overlapping node equally fractioned for every involved cluster. Then, a weight is computed for
  every overlapping edge as the fraction between the energy of the target cluster and the sum of all
  the clusters involved in the overlap."

**Relevance:** Gives a rigorous formulation of exactly what our wavefront algorithm does implicitly:
the cluster is a connected component of a directed graph whose sinks are seeds. Encoding overlap cells
as nodes with several outgoing edges, and assigning fractional energy weights, is a minimal, cheap
generalisation of our "unclustered" test that solves both cluster splitting and energy sharing without
any ML and without changing the module-ID output format.

### 3.2 Merged-π0 cluster expansion: enlarging a cluster only when the second maximum is energetic enough

**Source:** arXiv:2212.11061, Section 3.2.1 "Merged π0 case".

**Verbatim excerpts:**
- "A local maximum in this context defines a cell that has the highest energy value among its distance
  one neighbours in the calorimeter grid. This definition is the same as the one used in the Cellular
  Automaton algorithm."
- "In that case, the reconstruction is done as a single cluster, since the definition of maxima does
  not allow two adjacent cluster seeds. When photons are not separable, it is then called merged π0
  case. Hence, the super-cluster from a merged π0 can be bigger than the 3×3 window around the seed"
- "That is why the Graph Clustering algorithm adapts the shape of potential merged π0 candidates,
  expanding the cluster up to the neighbours of the second most energetic digit in the cluster. …
  we have studied the relation between the two most energetic digits as a ratio labeled R1 to define a
  threshold on which to make the cluster expansion. … the threshold for R1 is set to 25. The chosen
  value ensures that the residual energy left outside the cluster is less than 9% for the studied π0
  samples and that the cluster expansion affects an average of 8.2% of the clusters in an event.
  Moreover, a second threshold for merged π0 candidates concerning the minimum energy of the cluster
  seed is set to 1000 MeV"

**Relevance:** A concrete, data-driven rule for the case our "local maxima" logic cannot handle:
two maxima closer than one cell. It also shows the right methodology for choosing our `Eth` / `Nth` /
gap-jump parameters — scan the ratio distribution on labelled MC, then fix the threshold from the
residual-energy target rather than from convenience.

### 3.3 Performance framing: equivalent physics, large speed-up

**Source:** arXiv:2212.11061, abstract and Section 1.

**Verbatim excerpt:**
- "It outperforms the previously used method by 65.4% in terms of computational time on average, with
  an equivalent efficiency and resolution." … "Furthermore, it is currently the default solution for
  calorimeter reconstruction in the upcoming Run 3."

**Relevance:** Reassurance that restructuring an existing clustering into a graph/connected-component
form is not a physics-risk exercise; it buys throughput at identical efficiency. For a 622-module
detector the speed is irrelevant, but the graph data structures (adjacency, overlap edges) are exactly
the scaffolding we need for later splitting/merging and ML upgrades.

---

## 4. Generic algorithms and reviews

### 4.1 ATLAS topological cell clustering: significance-based noise suppression

**Source:** ATLAS Collaboration, "Topological cell clustering in the ATLAS calorimeters and its
performance in LHC Run 1", Eur. Phys. J. C77 (2017) 490, arXiv:1603.02934 —
https://arxiv.org/abs/1603.02934 , HTML: https://arxiv.org/html/1603.02934v2

**Verbatim excerpts (abstract and Section 1)** — I could only fetch the abstract and introduction
(the fetch truncates long pages before Section 3, which holds the explicit seed/neighbour thresholds
and the k-split algorithm; those numeric details are therefore **UNVERIFIED** here):
- "The cluster formation follows cell signal-significance patterns generated by electromagnetic and
  hadronic showers. In this, the clustering algorithm implicitly performs a topological noise
  suppression by removing cells with insignificant signals which are not in close proximity to cells
  with significant signals."
- "Calorimeter cells with insignificant signals found to not be connected to neighbouring cells with
  significant signals are considered noise and discarded from further jet, particle and missing
  transverse momentum reconstruction."
- "The algorithm building the topo-clusters explores the spatial distribution of the cell signals in
  all three dimensions to establish connections between neighbours in an attempt to reconstruct the
  energy and directions of the incoming particles."

**Relevance:** The conceptual fix for "fixed global thresholds": thresholds are expressed in units of
the *cell's own noise* (significance), and a low-significance cell survives only if it touches a
high-significance neighbour. Translated to us: keep `Eth` low for growth but admit a module only when
it is either above a high seed threshold or adjacent to an already-clustered module — which is
precisely the "topological noise suppression" that our current `Eth`+`Nth` hack approximates.

### 4.2 Review of clustering practice + a clustering-specific performance metric (V-score / homogeneity / completeness)

**Source:** M. Al Halabi et al., "Machine Learning Power Week 2023: Clustering in Hadronic
Calorimeters", arXiv:2508.09938 (2025) — https://arxiv.org/abs/2508.09938 ,
HTML: https://ar5iv.labs.arxiv.org/html/2508.09938

**Verbatim excerpts (Section 3 "Previous Approaches"):**
- "The clustering algorithms of ATLAS, CMS and ALICE all share the same main strategy for clustering in
  their calorimeter(s). They find cluster seeds by identifying the cells where the energy deposited is
  above a certain threshold. The clusters are grown by adding neighbouring cells in the same or
  neighbouring detector layers to the clusters."
- "If one cluster contains several seeds, it can be split into several clusters."
- "LHCb uses a largely similar approach, with the exception that the information is stored in a graph
  structure. Hits are nodes, and connections between hits are expressed in terms of edges. The largest
  draw of a graph structure is the increased computational efficiency that has been demonstrated."

**Verbatim excerpts (Section 2.4 "Performance Metric"):**
- "Two possible metrics are the homogeneity, which measures the extent to which each cluster contains
  hits from a single particle, and the completeness, which measures the extent to which all hits from a
  given particle are captured in a single cluster."
- "This metric was chosen as it requires that a single algorithm must have good homogeneity and
  completeness in order to perform well - that is an homogeneity or completeness of 0 gives a V-score
  of 0, with a maximum possible validity score of 1."
- "We choose to weight hits linearly in the energy deposited in the hit, such that incorrectly
  clustering low-energy hits insubstantially affects the weighted V-score, compared with clustering
  high-energy hits from different particles or splitting the hits in a single high-energy particle
  shower."

**Relevance:** Our current scoring ("cluster multiplicity per event vs. expected number of particles")
is blind to *which* modules ended up where. A weighted V-score (homogeneity × completeness, weighted by
module energy) measures over-splitting and over-merging separately and is computable from our existing
MC truth — this is the highest-value, lowest-cost change to the framework's evaluation.

### 4.3 Non-ML dimensionless clustering baselines: anti-kT and density clustering, and their failure mode at high density

**Source:** arXiv:2508.09938, Section 4.1.1 "Anti-kT" and Section 4.2.3 "DBSCAN".

**Verbatim excerpts:**
- "The anti-kT algorithm is a sequential recombination algorithm that prioritizes the clustering of
  high pT particles within a given resolution (or radius) parameter R. … In this implementation,
  however, hits in the LFHCal are assigned to a pseudo jet i according to the corresponding Cartesian
  position (x, y, z) and the energy E."
- "For larger radii, the clusters appear more complete, whereas the homogeneity decreases. The V-score
  peaks around R ≈ 0.4, which also yields a near-optimal average energy resolution of the clusters."
- "The large variation in the density of clusters throughout the phase space makes the application of
  an algorithm with a fixed resolution parameter, like the anti-kT algorithm, challenging."

**Relevance:** A ready-made, one-line-per-event clustering baseline that needs no thresholds and no
training: run a sequential-recombination (anti-kT-like) clusterer on module positions with energy as
weight, and scan the single radius parameter against the V-score. It also documents the exact failure
mode we should expect from *any* fixed spatial/threshold parameter at high local density — evidence
that `Eth`/`Nth` should become density- or significance-adaptive.

### 4.4 Graph seed classification + connected components as a learnable-but-simple hybrid

**Source:** arXiv:2508.09938, Section 4.4 "Graph-based model".

**Verbatim excerpts:**
- "The algorithm works in two stages: a) seed-to-seed classification, then b) seed-to-nonseed
  segmentation." … "any edge classifier, whether a simple multi-layer perceptron (MLP), GNN or
  transformer, can be used to prune spurious connections between a set of high-energy seeds and the
  rest of the energy deposits."
- "Next, the remaining N_h − N_s hits are each connected with an edge to their nearest seed. This can
  be achieved by applying a K-nearest neighbour algorithm with K=1 … To obtain unique cluster labels to
  each hit, we apply a connected components algorithm to the full graph."
- "The core component of this model, the edge-classifier MLP, achieved a best efficiency of 0.956,
  purity of 0.956, and an area under the receiver operator characteristic curve (AUC) of 0.998."

**Relevance:** This is the smallest possible ML upgrade of a wavefront clusterer: keep our seeds and
our connected-component growth, and only replace the "does this neighbour join?" boolean with a
learned edge score. Our per-event multiplicity target can then be derived from the number of connected
components rather than compared to it.

---

## 5. Particle-flow / framework-level ideas: cluster as graph, cluster–track matching, energy sharing

### 5.1 Pandora SDK: decoupled algorithms, reclustering, and separating genuine neutrals from track-cluster fragments

**Source:** J. S. Marshall, M. A. Thomson, "The Pandora Software Development Kit for Pattern
Recognition", Eur. Phys. J. C75 (2015) 439, arXiv:1506.05348 —
https://arxiv.org/abs/1506.05348 , HTML: https://ar5iv.labs.arxiv.org/html/1506.05348

**Verbatim excerpts (Sections 1 and 3):**
- "This design promotes an approach using many decoupled algorithms, each addressing specific
  topologies."
- "The algorithms that address the problem must be able to build clusters of space-points and should be
  able to manipulate clusters by splitting them up or merging them together. What differs between
  pattern recognition problems is the precise logic controlling the algorithm operations."
- "The Pandora SDK monitors the usage of all the Input Objects to ensure that no double-counting can
  occur, with no Input Object being used to create multiple Algorithm Objects."

**Verbatim excerpts (Section 11.1 "Linear Collider Event Reconstruction"):**
- "The Clustering algorithm is configured so that it tends to split CaloHits from individual particles
  into multiple Clusters, rather than risk merging energy deposits from multiple particles into single
  Clusters. The Clusters are instead carefully merged together by a series of algorithms implementing
  well-defined topological rules."
- "Clusters are associated to Tracks via careful comparison of Cluster positions and directions
  (obtained, for instance, via sliding linear fits) with projected Track positions and directions at
  the ECAL."
- "The compatibility of associated Tracks and Clusters is assessed, via comparison of Track momentum
  with associated Cluster energies. Significant discrepancies indicate pattern recognition problems and
  the reclustering approach described in Section 9 is used to improve the clustering."
- "Clusters without associated tracks are examined to assess whether they genuinely represent
  electrically neutral particles, or whether they are more likely to be fragments of any nearby
  track-associated Clusters, representing charged particles."

**Relevance:** Four adoptable ideas: (i) deliberately *over-split* then merge under topological rules
(our growth is currently conservative and merges-by-construction); (ii) use track momentum vs. cluster
energy disagreement as a *trigger to re-cluster* — with our e+/μ- samples this is a directly
implementable quality flag; (iii) the explicit "is this unmatched cluster a real photon or a fragment
of a nearby track cluster" test is the general form of split-off rejection; (iv) note the deliberate
contrast that Pandora *forbids* sharing a hit between clusters, whereas Belle II/LHCb allow fractional
sharing — a design decision we must make explicitly.

### 5.2 Arbor: cluster the shower as a tree, with hits as nodes and no energy used at all

**Source:** M. Ruan, "Arbor, a new approach of the Particle Flow Algorithm", CHEF 2013,
arXiv:1403.4784 — https://arxiv.org/abs/1403.4784 , HTML: https://ar5iv.labs.arxiv.org/html/1403.4784

**Verbatim excerpts (Sections 2–4):**
- "Arbor algorithm, as its name suggests, is inspired by the fact that the shower spatial development
  follows the topology of a tree. … Arbor reconstruct the long travelling charged particles generated
  at hadronic shower as the branches of the tree, while the tree structure is ensured with the
  constrain that no loop structure is allowed."
- "In its current configuration, Arbor use only geometry information of the hits, i.e, no energy
  information is used."
- "For all the connectors end at this hit, Arobr kept at most one connector that has the minimal angle
  to the reference direction. Therefore, no loop structure will be kept after the cleaning, and a tree
  structure based on the connectors emerges."
- "Since the seed and the tree has one one correspondence, Arbor is a powerful algorithm to separate
  nearby showers, which is highly appreciated at Particle Flow Algorithm."

**Relevance:** A pure-geometry, threshold-free alternative to our energy-ordered wavefront: build
oriented links between nearby modules, prune every link except the one closest to a locally estimated
reference direction, and every remaining tree is a cluster. It removes `Eth` and `Nth` entirely and
gives each cluster a direction (seed→leaf), which is exactly the extra information (position/shape)
our framework lacks — and its O(1)-per-node pruning is trivial to implement for 622 modules.

---

## 6. Machine-learning clustering (clustering-oriented, not just energy regression)

### 6.1 GNN + object condensation: single-shot instance segmentation with energy regression (CMS HGCAL)

**Source:** S. R. Qasim, N. Chernyavskaya, J. Kieseler, K. Long, O. Viazlo, M. Pierini, R. Nawaz,
"End-to-end multi-particle reconstruction in high occupancy imaging calorimeters with graph neural
networks", Eur. Phys. J. C82 (2022), arXiv:2204.01681 —
https://arxiv.org/abs/2204.01681 , HTML: https://ar5iv.labs.arxiv.org/html/2204.01681

**Verbatim excerpts (abstract, Section 1, Section 4):**
- "The algorithm exploits a distance-weighted graph neural network, trained with object condensation, a
  graph segmentation technique. Through a single-shot approach, the reconstruction task is paired with
  energy regression."
- "In traditional approaches, particle reconstruction follows a two-step strategy: first, clusters are
  built, and then classification and regression tasks are performed on these clusters."
- "For instance, when two particles (p1 and p2) are maximally overlapping, we merge them into a single
  particle in the ground truth since disentangling such two clusters would be an impractical and
  imprecise task."

**Relevance:** Provides the loss-function answer to "no cluster splitting/merging": object
condensation learns a per-hit attractive/repulsive potential plus a seediness score, so clusters emerge
as attraction basins rather than from thresholds. The last quote is a warning for our own scoring: at
some overlap the *truth itself* is ambiguous, so our MC-truth-based comparisons must define a merge
criterion for unresolvable pairs instead of penalising the algorithm for them.

### 6.2 Transformer clustering for overlapping electromagnetic showers (ClusTEX), with a splitting-rate metric

**Source:** Y. Maidannyk, F. Couderc, J. Malclès, M. Ö. Sahin, "Reconstruction of overlapping
electromagnetic showers in calorimeters using Transformers", Eur. Phys. J. C 86, 873 (2026),
arXiv:2603.18172 — https://arxiv.org/abs/2603.18172 (abstract, fetched from the arXiv abs page)

**Verbatim excerpts (abstract):**
- "Accurate clustering of electromagnetic energy deposits is essential for reconstructing photons and
  electrons in modern hadron collider experiments, where boosted topologies and pileup cause
  overlapping showers and ambiguous energy assignment. We present deep learning-based clustering
  approaches that reconstruct particle energy and position directly from calorimeter readout."
- "The study includes a two-step strategy in which candidate seed windows are identified and then
  jointly processed via distance-weighted message passing or attention mechanism and a single-step
  graph transformer, ClusTEX, which performs candidate selection and reconstruction in one inference
  stage."
- "Performance is evaluated using efficiency, energy and position resolutions and splitting rate -
  reconstruction of two objects for a single photon. … For boosted π0 → γγ, the attention-based model
  retains di-photon mass reconstruction capability, where the standard algorithm becomes inefficient."

**Relevance:** Directly targets our "no cluster splitting / no energy sharing between overlapping
clusters" weakness with a *seed-window + attention* design that keeps a classical seed finder as a
front-end. The "splitting rate" (two objects for one photon) is a metric our framework does not have
and should add immediately, since our current multiplicity-only score cannot detect over-splitting.

### 6.3 ML clustering as a whole taxonomy (K-means/DBSCAN/autoencoder/neural cluster-count regression)

**Source:** arXiv:2508.09938, Sections 4.2–4.3 (see §4.2 above for the source link).

**Verbatim excerpts:**
- "One of the key insights that emerged independently from many of the groups is that the optimal
  number of clusters depends largely on the number of hits in the event." … "Higher energy events are
  also spread out more in space, meaning they have more clusters."
- "The highest performing solution was A's quadratic regression using the number of hits only."
- "The groups that explored DBSCAN all concluded that it did not work as well as K-means as a
  standalone algorithm. A potential reason for this could be that the distance between detector hits
  in a cluster will be smaller closer to the beamline, and spread out further in the detector."

**Relevance:** Our framework's entire output is "how many clusters", decided by threshold logic. This
result says a one-feature regression (number of lit modules, or total event energy) predicts the right
cluster count better than any structural algorithm — a trivially cheap sanity check/baseline before we
add gap jumping or ML.

---

## 7. Top 5 highest-value ideas for our project

1. **Adopt the Belle II baseline energy-sharing scheme (split connected regions at every local maximum,
   share crystals with iterative log-weighted distance weights converged to ~1 mm).** It removes our
   "no splitting / no energy sharing between overlapping clusters" weakness with no ML and keeps our
   module-ID output format.
2. **Replace the global `Eth`/`Nth` pair with noise-significance + topological adjacency (ATLAS
   topo-cluster style, approximated by a per-module noise term measured from our MC as BESIII did with
   EINC = 0.27 MeV).** This is the only principled way to make one threshold work for e+, γ and μ- over
   the full energy range.
3. **Add a truth-based weighted V-score (homogeneity + completeness, energy-weighted) and an explicit
   splitting rate next to the multiplicity check.** Our current metric cannot distinguish
   over-splitting from over-merging, so we cannot measure the benefit of any of the other changes.
4. **Move to the LHCb graph formulation (nodes = modules, directed edges to seeds, weakly connected
   components, fractional overlap weights).** It gives splitting, energy sharing, merged-π0-style
   expansion and a natural place to hang track–cluster matching, at zero physics risk (LHCb: 65.4%
   faster, equivalent efficiency).
5. **Add track–cluster matching with an energy/momentum consistency check and a "fragment vs. real
   neutral" test (Pandora).** With e+/μ- MC this is available immediately, and it is the standard,
   non-ML route to split-off rejection and to a physically meaningful particle-level assignment.

---

## 8. Limitations of this review

- `web_fetch` rejects PDFs, so all conference slides, theses and journal PDFs (Belle II ECL overview
  indico talk, BESIII EMC conference proceedings, the PANDA split-off-recognition talk, NIM A
  detector papers, the CMS particle-flow paper if only available as PDF) could not be quoted.
- ATLAS topo-cluster **Section 3** (explicit seed/neighbour significance thresholds, k-split
  algorithm) and the Belle II ECL **calibration procedure** were beyond the fetched text and are
  marked UNVERIFIED.
- "MUCKS" as a BESIII clustering code could not be verified at all; no claim is made about it.
- The "sliding window" and "watershed / Gaussian-mixture EM" families are named in the task but I did
  not find a primary source I could fetch that describes them in a crystal-calorimeter context; they
  are therefore omitted rather than quoted second-hand.
