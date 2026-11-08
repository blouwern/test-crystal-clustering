# Background scan — best points

## Best point per neighbour threshold (mean over the three signal samples)

| Nth | seed (MeV) | wave (MeV) | exact-1 | efficiency | over-split | fake/event | purity | signal capture | bkg capture | bias (E>20 MeV) |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 0 | 10 | 0.02 | 0.8601 | 0.8656 | 0.0054 | 0.003 | 0.9909 | 0.9684 | 0.0319 | -0.0175 |
| 1 | 8 | 0.02 | 0.8614 | 0.8668 | 0.0054 | 0.001 | 0.9913 | 0.9687 | 0.0309 | -0.0174 |
| 2 | 6 | 0.02 | 0.8379 | 0.8421 | 0.0042 | 0.001 | 0.9914 | 0.9477 | 0.0298 | -0.0164 |

## Best point per dataset (any Nth)

| dataset | Nth | seed | wave | exact-1 | efficiency | purity | fake/event | signal capture |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| ep_iso | 0 | 13 | 0.02 | 0.8783 | 0.8802 | 0.9925 | 0.000 | 0.9746 |
| gamma_iso | 0 | 14 | 0.02 | 0.8739 | 0.8747 | 0.9928 | 0.000 | 0.9799 |
| mu_rest | 1 | 7 | 0.02 | 0.8423 | 0.8550 | 0.9876 | 0.007 | 0.9557 |

## Same-code no-background comparison at the best bkg points

| dataset | bkg exact | no-bkg exact | Δ exact | bkg purity | no-bkg purity |
|:--|--:|--:|--:|--:|--:|
| ep_iso | 0.8783 | 0.8790 | -0.0007 | 0.9925 | nan |
| gamma_iso | 0.8739 | 0.8744 | -0.0005 | 0.9928 | nan |
| mu_rest | 0.8423 | 0.8483 | -0.0061 | 0.9876 | nan |
