"""Energy metrics for multi-particle samples.

For a multi-particle event the leading cluster is only one shower, so a
"leading cluster vs total signal energy" ratio is not a per-shower energy
response.  The response is therefore computed per primary from the matched
overlap entry, while the containment/capture metrics are event-level.
"""

from __future__ import annotations

import numpy as np


def event_energy(cluster_sig, cluster_bkg, e_sig_total, overlap, assignment,
                 primary_totals, min_primary_energy=20.0):
    """Energy metrics for one event.

    ``overlap`` is the ``n_clusters x n_primary`` signal-overlap matrix and
    ``assignment`` maps clusters to primaries (Hungarian).  The response is
    only accumulated for primaries with at least ``min_primary_energy`` MeV of
    deposited signal energy.
    """
    cluster_sig = np.asarray(cluster_sig, dtype=np.float64)
    cluster_bkg = np.asarray(cluster_bkg, dtype=np.float64)
    cluster_tot = cluster_sig + cluster_bkg
    overlap = np.asarray(overlap, dtype=np.float64)

    lead = int(np.argmax(cluster_tot)) if cluster_tot.size else -1
    signal_capture = (float(cluster_sig.sum()) / float(e_sig_total)
                      if e_sig_total > 0 else 0.0)
    lead_containment = (float(cluster_sig[lead]) / float(e_sig_total)
                        if lead >= 0 and e_sig_total > 0 else 0.0)

    biases = []
    for i, j in enumerate(assignment):
        if j >= 0 and primary_totals[j] >= float(min_primary_energy):
            biases.append(float((overlap[i, j] - primary_totals[j])
                                / primary_totals[j]))

    return {
        "signal_capture": signal_capture,
        "lead_signal_containment": lead_containment,
        "energy_biases": biases,
        "n_energy": len(biases),
    }
