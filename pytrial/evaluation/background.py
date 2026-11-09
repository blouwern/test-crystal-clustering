"""Background-aware metrics for multi-particle clustering.

Definitions follow ``perf_report_bkg_20261106.md``:

* fake cluster: a cluster with less than ``fake_threshold`` MeV of signal
  energy (i.e. essentially a background-only cluster);
* signal purity: cluster signal energy / cluster total energy;
* background capture: background energy inside all clusters / total
  background energy;
* leading-cluster background fraction: background energy in the leading
  cluster (highest total energy) / that cluster's total energy.
"""

from __future__ import annotations

import numpy as np


def event_background(cluster_sig, cluster_bkg, e_sig_total, e_bkg_total,
                     fake_threshold=1.0):
    """Background metrics for one event.

    ``cluster_sig`` / ``cluster_bkg`` are arrays with one entry per
    reconstructed cluster; ``e_sig_total`` / ``e_bkg_total`` are the event
    totals.  ``bkg_capture`` is NaN when the event has no background.
    """
    cluster_sig = np.asarray(cluster_sig, dtype=np.float64)
    cluster_bkg = np.asarray(cluster_bkg, dtype=np.float64)
    cluster_tot = cluster_sig + cluster_bkg

    sig_in = float(cluster_sig.sum())
    bkg_in = float(cluster_bkg.sum())
    fake = cluster_sig < float(fake_threshold)

    lead = int(np.argmax(cluster_tot)) if cluster_tot.size else -1
    lead_tot = float(cluster_tot[lead]) if lead >= 0 else 0.0
    lead_bkg = float(cluster_bkg[lead] / lead_tot) if lead_tot > 0 else 0.0

    return {
        "n_fake": int(fake.sum()),
        "has_fake": bool(fake.any()),
        "signal_purity": (sig_in / (sig_in + bkg_in)
                          if sig_in + bkg_in > 0 else 0.0),
        "bkg_capture": (bkg_in / float(e_bkg_total)
                        if e_bkg_total > 0 else float("nan")),
        "lead_bkg_fraction": lead_bkg,
    }
