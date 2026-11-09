"""Multi-particle evaluation modules.

The package is split by concern so the source files stay small:

* :mod:`evaluation.counts`     - count-based metrics
* :mod:`evaluation.matching`   - overlap matrix, Hungarian assignment, merge/split
* :mod:`evaluation.background` - signal purity, fake clusters, background capture
* :mod:`evaluation.energy`     - energy capture, containment, per-primary response
* :mod:`evaluation.io`         - stitched-sample and cluster-file readers
"""

from . import background, counts, energy, matching, separation

__all__ = ["background", "counts", "energy", "matching", "separation"]
