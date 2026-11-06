"""Multi-particle evaluation modules.

The package is split by concern so the source files stay small:

* :mod:`evaluation.counts`   - count-based metrics (implemented in Stage 1)
* :mod:`evaluation.matching` - truth matching, efficiency/purity, merge/split
  (Stage 2/3)
* :mod:`evaluation.energy`   - per-primary energy metrics (Stage 2)
"""

from . import counts

__all__ = ["counts"]
