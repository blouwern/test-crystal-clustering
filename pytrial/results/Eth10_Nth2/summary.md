# Reconstruction Summary — Eth=10, Nth=2

Reconstruction of single-particle datasets with Eth=10.0, Nth=2.0, expected lambda=1.

| Dataset | Events | Mean | 0-cl | 1-cl | 2-cl | 3+ | chi²(Pois1) | KS D |
|:--------|-------:|-----:|-----:|-----:|-----:|---:|---------:|-----:|
| ep_iso | 1,000,000 | 1.244 | 11.9% | 54.4% | 31.1% | 2.6% | 378,088 | 0.616 |
| gamma_iso | 1,000,000 | 1.184 | 12.7% | 58.1% | 27.3% | 1.9% | 371,852 | 0.609 |
| mu_rest | 1,000,000 | 1.030 | 16.3% | 65.1% | 17.9% | 0.7% | 400,037 | 0.573 |

- All datasets are strongly 1-cluster dominant (54–65%), consistent with single-particle inputs
- mu_rest is cleanest: lowest mean (1.03), highest 1-cluster fraction (65.1%), lowest 2+ cluster rate (18.6%)
- All sharply non-Poisson (p ≈ 0) — distributions are tightly peaked, not Poisson-like

### Input

- `$ECAL_CLUSTERING_DATA_DIR/processed/edep_of_each_evt_*.root`
- Neighbor topology: `$ECAL_CLUSTERING_DATA_DIR/utilities/ecal_neighbor_info.root`

### Parameters

- Energy threshold: 10.0
- Neighbor threshold: 2
- Neighbor file: ecal_neighbor_info.root (primary neighbors only)
