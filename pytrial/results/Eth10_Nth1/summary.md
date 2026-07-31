# Reconstruction Summary — Eth=10, Nth=1

Reconstruction of single-particle datasets with Eth=10.0, Nth=1.0, expected lambda=1.

| Dataset | Events | Mean | 0-cl | 1-cl | 2-cl | 3+ | chi²(Pois1) | KS D |
|:--------|-------:|-----:|-----:|-----:|-----:|---:|---------:|-----:|
| ep_iso | 1,000,000 | 1.217 | 11.9% | 56.6% | 29.2% | 2.2% | 381,625 | 0.616 |
| gamma_iso | 1,000,000 | 1.161 | 12.7% | 60.1% | 25.6% | 1.6% | 385,147 | 0.609 |
| mu_rest | 1,000,000 | 1.016 | 16.3% | 66.4% | 16.7% | 0.6% | 423,552 | 0.573 |

### Comparison vs Nth=2

| Dataset | Nth=2 mean | Nth=1 mean | Δmean | Nth=2 1-cl | Nth=1 1-cl |
|:--------|----------:|----------:|------:|----------:|----------:|
| ep_iso | 1.244 | 1.217 | −0.027 | 54.4% | 56.6% |
| gamma_iso | 1.184 | 1.161 | −0.023 | 58.1% | 60.1% |
| mu_rest | 1.030 | 1.016 | −0.014 | 65.1% | 66.4% |

Nth=1 produces fewer clusters than Nth=2 across all datasets — looser propagation requirement reduces cluster splitting, pushing mean closer to 1.0 and increasing 1-cluster fraction.

### Parameters

- Energy threshold: 10.0
- Neighbor threshold: 1
- Neighbor file: ecal_neighbor_info.root (primary neighbors only)
