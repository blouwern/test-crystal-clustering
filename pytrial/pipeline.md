# Pipeline — Reconstruction & Scoring

ECAL crystal clustering pipeline: wavefront-propagation reconstruction followed by statistical or event-by-event scoring against ground truth.

---

## Reconstruction (`reconstruction.py`)

Wavefront-propagation clustering algorithm operating on per-event crystal energy depositions.

### Algorithm

1. **Load topology** — reads crystal neighbor graph from `ecal_neighbor_info.root`
2. **Seed selection** — sorts crystals by energy descending; the highest-energy unclustered crystal above `energy_threshold` becomes a seed
3. **Wavefront propagation** — from the seed, iteratively propagates outward through neighbors. A crystal joins the cluster if it is:
   - Unclustered
   - Above `energy_threshold`
   - Has at least `n_neighbor_threshold` lit neighbors
4. **Repeat** — continue finding seeds until no remaining crystals qualify
5. **Output** — writes `ECALClusters` tree (branches: `EvtID/I`, `ModIDList`) to a ROOT file

### API

```python
def reconstruction(filename_edep_of_each_evt,
                  energy_threshold, n_neighbor_threshold,
                  data_dir="$ECAL_CLUSTERING_DATA_DIR",
                  output_path=None) -> Path
```

| Parameter | Type | Default | Description |
|:----------|:-----|:--------|:------------|
| `filename_edep_of_each_evt` | `str` | required | Input file name under `data_dir/processed/` |
| `energy_threshold` | `float` | required | Minimum crystal energy (MeV) for seeding & cluster membership |
| `n_neighbor_threshold` | `int` | required | Minimum lit neighbors required for wavefront propagation |
| `data_dir` | `str` | `$ECAL_CLUSTERING_DATA_DIR` | Root data directory (env var) |
| `output_path` | `str\|None` | `None` | Output ROOT path; `None` → auto-name `ecal_clusters_Eth-X_Nth-Y.root` |

Returns the output `Path` for chaining into scoring.

### CLI

```
python reconstruction.py <filename> [Eth=10] [Nth=2] [output_name]
                         [expected_lambda] [--event-by-event] [--data-path path]
```

Runs reconstruction then automatically calls `scoring.run()` on the output.

### Output tree schema

```
ECALClusters
├── EvtID/I           event index
└── ModIDList         vector<int> crystal IDs in this cluster
```

---

## Scoring (`scoring.py`)

Evaluates clustering output against ground truth or expected distributions. Two assessment modes.

### Modes

#### 1. Statistical mode (default)

Compares the observed cluster-count distribution to a reference:

- **Chi-squared goodness-of-fit** — automated bin merging (expected < 5 merged into adjacent bins)
- **KS test** — one-sample vs Poisson, or two-sample vs ground-truth NTrack distribution

Reference selection (priority order):
1. **Empirical ground-truth** — if `data_path` has an `NTrack` branch → two-sample KS, chi-squared vs empirical distribution
2. **Explicit lambda** — `expected_lambda=X` → Poisson(λ=X), ddof=0
3. **Auto-estimation** — MLE from cluster data, ddof=1

#### 2. Event-by-event mode (`event_by_event=True`)

Requires `data_path` with `NTrack` branch. Per-event cluster-vs-truth comparison:

| Metric | Description |
|:-------|:------------|
| Exact match rate | Fraction where `n_clusters == n_tracks` |
| Over/under-clustering | Fraction where `n_clusters >` or `< n_tracks` |
| MAE | Mean absolute error |
| RMSE | Root mean squared error |
| Pearson r | Correlation between predicted and true counts |

### API

```python
def run(input_path, expected_lambda=None, output_plot=None,
        event_by_event=False, data_path=None)
```

| Parameter | Type | Default | Description |
|:----------|:-----|:--------|:------------|
| `input_path` | `str\|Path` | required | Reconstruction output ROOT file (`ECALClusters` tree) |
| `expected_lambda` | `float\|None` | `None` | Poisson λ for statistical comparison; `None` = MLE from data |
| `output_plot` | `str\|None` | `None` | Plot path; `None` = auto-name from input |
| `event_by_event` | `bool` | `False` | Enable per-event comparison |
| `data_path` | `str\|None` | `None` | Data ROOT with `NTrack` branch in `EdepOfEachEvt` tree |

### Internal helpers

```python
_read_cluster_counts(input_path) -> np.ndarray    # cluster count per event
_read_ntrack(data_path, n_expected) -> np.ndarray  # ground-truth track count per event
```

---

## Examples

```python
from reconstruction import reconstruction
import scoring

# Reconstruct with Eth=10, Nth=2
out = reconstruction("child_stitched_edeps.root", 10.0, 2.0)

# Statistical scoring, auto-estimate lambda
scoring.run(out)

# Statistical scoring, compare against ground-truth distribution
scoring.run(out, data_path="data/processed/child_stitched_edeps.root")

# Event-by-event scoring
scoring.run(out, data_path="data/processed/child_stitched_edeps.root",
            event_by_event=True)
```

```bash
# CLI: reconstruct + event-by-event scoring against ground truth
python reconstruction.py child_stitched_edeps.root 10 2 \
    --event-by-event --data-path $ECAL_CLUSTERING_DATA_DIR/processed/child_stitched_edeps.root
```

---

## Data flow

```
child_stitched_edeps.root          ecal_neighbor_info.root
(EdeOfEachEvt: Edeps, NTrack)      (ECALCrystalNeighbors: neighbors)
         │                                    │
         └────────────┬───────────────────────┘
                      │
              reconstruction()
                      │
                      ▼
           ecal_clusters_Eth-X_Nth-Y.root
           (ECALClusters: EvtID, ModIDList)
                      │
                      ▼
                scoring.run()
                      │
            ┌─────────┴──────────┐
            ▼                    ▼
     Statistical mode     Event-by-event mode
     (distribution)       (per-event accuracy)
```
