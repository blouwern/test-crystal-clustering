# Reconstruction Summary — Eth=10, Nth=2 (v2, gap-jumping)

Reconstruction with **second-neighbor gap-jumping** mechanism enabled.

| Dataset | Events | Mean | 0-cl | 1-cl | 2-cl | 3+ | chi²(Pois1) | KS D |
|:--------|-------:|-----:|-----:|-----:|-----:|---:|---------:|-----:|
| ep_iso | 1M | 1.244 | 11.9% | 54.4% | 31.1% | 2.6% | 378,088 | 0.616 |
| gamma_iso | 1M | 1.184 | 12.7% | 58.1% | 27.3% | 1.9% | 371,852 | 0.609 |
| mu_rest | 1M | 1.030 | 16.3% | 65.1% | 17.9% | 0.7% | 400,037 | 0.573 |

### vs v1 (no gap-jumping)

Results are **identical** to the original Nth=2 run — single-particle events produce contiguous energy depositions, so the gap-jumping mechanism (secondary neighbors) has no effect on these datasets.

### Parameters

- Energy threshold: 10.0
- Neighbor threshold: 2
- Neighbor file: ecal_neighbor_info_added.root (primary + secondary neighbors)
- Gap-jumping: enabled (if wavefront dies, fallback to second-nearest neighbors)
