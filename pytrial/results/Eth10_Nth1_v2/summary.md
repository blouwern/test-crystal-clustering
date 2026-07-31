# Reconstruction Summary — Eth=10, Nth=1 (v2, gap-jumping)

Reconstruction with **second-neighbor gap-jumping** mechanism enabled.

| Dataset | Events | Mean | 0-cl | 1-cl | 2-cl | 3+ | chi²(Pois1) | KS D |
|:--------|-------:|-----:|-----:|-----:|-----:|---:|---------:|-----:|
| ep_iso | 1M | 1.217 | 11.9% | 56.6% | 29.2% | 2.2% | 381,623 | 0.616 |
| gamma_iso | 1M | 1.161 | 12.7% | 60.1% | 25.6% | 1.6% | 385,158 | 0.609 |
| mu_rest | 1M | 1.016 | 16.3% | 66.4% | 16.7% | 0.6% | 423,566 | 0.573 |

### vs v1 (no gap-jumping)

Results are **identical** to the original Nth=1 run — single-particle datasets have no gaps for the mechanism to bridge.

### Parameters

- Energy threshold: 10.0
- Neighbor threshold: 1
- Neighbor file: ecal_neighbor_info_added.root (primary + secondary neighbors)
- Gap-jumping: enabled (if wavefront dies, fallback to second-nearest neighbors)
