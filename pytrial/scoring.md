# Scoring

Evaluates ECAL clustering output against ground truth or expected distributions. Supports two assessment modes.

## Modes

### 1. Statistical mode (default)
Compares the observed cluster-count distribution to a reference using:
- **Chi-squared goodness-of-fit** with automated bin merging (expected < 5 merged into adjacent bins)
- **Kolmogorov-Smirnov test** (one-sample vs Poisson, or two-sample vs ground-truth)

The reference comes from, in priority order:
1. **Empirical ground-truth** — if `data_path` points to a file with an `NTrack` branch in its `EdepOfEachEvt` tree, uses that distribution directly (two-sample KS, chi-squared vs empirical)
2. **Explicit lambda** — `expected_lambda` parameter specifies a Poisson λ (ddof=0 in chi-squared)
3. **Auto-estimation** — MLE from the cluster data itself (ddof=1)

### 2. Event-by-event mode (`event_by_event=True`)
Requires `data_path` with `NTrack` branch. For each event, compares the number of reconstructed clusters to the known number of tracks:

| Metric | Description |
|:-------|:------------|
| Exact match rate | Fraction of events where `n_clusters == n_tracks` |
| Over/under-clustering rate | Fraction where `n_clusters >` or `< n_tracks` |
| MAE | Mean absolute error |
| RMSE | Root mean squared error |
| Pearson r | Linear correlation between predicted and true counts |
| Difference histogram | Per-event `n_clusters - n_tracks` distribution |

## API

```python
def run(input_path, expected_lambda=None, output_plot=None,
        event_by_event=False, data_path=None)
```

| Parameter | Type | Default | Description |
|:----------|:-----|:--------|:------------|
| `input_path` | `str\|Path` | required | Reconstruction output ROOT file (`ECALClusters` tree) |
| `expected_lambda` | `float\|None` | `None` | Poisson λ for statistical comparison; `None` = estimate from data |
| `output_plot` | `str\|None` | `None` | Path for distribution plot; `None` = auto-name from input |
| `event_by_event` | `bool` | `False` | Enable per-event ground-truth comparison |
| `data_path` | `str\|None` | `None` | Path to data ROOT file with `NTrack` branch in `EdepOfEachEvt` tree |

### Internal helpers (usable standalone)

```python
_read_cluster_counts(input_path) -> np.ndarray    # cluster count per event
_read_ntrack(data_path, n_expected) -> np.ndarray  # ground-truth track count per event
```

## Example

```python
import scoring

# Statistical, auto-estimate lambda from clusters
scoring.run("ecal_clusters_Eth-10_Nth-2.root")

# Statistical, compare against ground-truth NTrack distribution
scoring.run("ecal_clusters_Eth-10_Nth-2.root",
            data_path="data/processed/child_stitched_edeps.root")

# Event-by-event precision/recall scoring
scoring.run("ecal_clusters_Eth-10_Nth-2.root",
            data_path="data/processed/child_stitched_edeps.root",
            event_by_event=True)
```
