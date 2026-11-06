"""Multi-particle event stitching.

A stitched event is the per-crystal sum of N single-particle events.  The
detector response is assumed linear and the particles are assumed not to
interact, which is adequate for energy-map clustering studies (not for physics
yields).

Stage 1 implements the count-only variant (Scheme A): only the summed ``Edeps``
and ``NTrack`` (plus the source indices for reproducibility) are stored.
Stage 2 will add per-primary truth bookkeeping (Scheme B).
"""

from __future__ import annotations

import array
from pathlib import Path

import numpy as np
import ROOT
from ROOT import vector

N_MOD = 568
TREE_NAME = "EdepOfEachEvt"


def load_edeps(path, n_events=None, tree_name=TREE_NAME, n_mod=N_MOD):
    """Load the first ``n_events`` ``Edeps`` vectors as a float32 matrix."""
    with ROOT.TFile.Open(str(path), "READ") as f:
        tree = f.Get(tree_name)
        total = tree.GetEntries()
        n = total if n_events is None else min(int(n_events), total)
        edeps = vector("float")()
        tree.SetBranchAddress("Edeps", edeps)
        out = np.zeros((n, n_mod), dtype=np.float32)
        for i in range(n):
            tree.GetEntry(i)
            m = min(edeps.size(), n_mod)
            out[i, :m] = np.asarray(edeps, dtype=np.float32)[:m]
    return out


def positive_energy_mask(edeps, eps=0.0):
    """Mask of events with total deposited energy above ``eps``."""
    return np.asarray(edeps, dtype=np.float64).sum(axis=1) > eps


def stitch_count_only(edeps, n_primary, n_events, rng):
    """Sum ``n_primary`` single-particle events into ``n_events`` events.

    Parameters
    ----------
    edeps : ndarray, shape (pool, n_mod)
        Source single-particle events.
    n_primary : int
        Number of particles per stitched event (N).
    n_events : int
        Number of stitched events to generate.
    rng : numpy.random.Generator
        Seeded generator, so the sample is reproducible.

    Returns
    -------
    stitched : ndarray, shape (n_events, n_mod), float32
    ntrack : ndarray, shape (n_events,), int32
    source_idx : ndarray, shape (n_events, n_primary), int32
        Indices into the source pool, useful for debugging and re-running.
    """
    edeps = np.asarray(edeps, dtype=np.float32)
    pool = len(edeps)
    if n_primary < 1:
        raise ValueError("n_primary must be >= 1")
    if n_primary > pool:
        raise ValueError(f"n_primary={n_primary} exceeds pool size {pool}")

    stitched = np.zeros((int(n_events), edeps.shape[1]), dtype=np.float32)
    source_idx = np.zeros((int(n_events), int(n_primary)), dtype=np.int32)
    for i in range(int(n_events)):
        idx = rng.choice(pool, size=int(n_primary), replace=False)
        source_idx[i] = idx
        stitched[i] = edeps[idx].sum(axis=0)
    ntrack = np.full(int(n_events), int(n_primary), dtype=np.int32)
    return stitched, ntrack, source_idx


def write_stitched(path, edeps, ntrack, source_idx=None, tree_name=TREE_NAME):
    """Write a stitched sample as ``EdepOfEachEvt`` with ``Edeps``/``NTrack``.

    ``SourceIdx`` (one entry per primary) is written when ``source_idx`` is
    given.  The reconstruction only needs ``Edeps``; the other branches are
    truth/bookkeeping.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    edeps = np.asarray(edeps, dtype=np.float32)
    ntrack = np.asarray(ntrack, dtype=np.int32)
    source_idx = None if source_idx is None else np.asarray(source_idx, dtype=np.int32)

    f = ROOT.TFile.Open(str(path), "RECREATE")
    tree = ROOT.TTree(tree_name, "stitched multi-particle events")
    vec = vector("float")()
    n = array.array("i", [0])
    src = vector("int")()
    tree.Branch("Edeps", vec)
    tree.Branch("NTrack", n, "NTrack/I")
    if source_idx is not None:
        tree.Branch("SourceIdx", src)

    for i in range(len(edeps)):
        vec.clear()
        for value in edeps[i]:
            vec.push_back(float(value))
        n[0] = int(ntrack[i])
        if source_idx is not None:
            src.clear()
            for value in source_idx[i]:
                src.push_back(int(value))
        tree.Fill()

    tree.Write()
    f.Close()
    return path


def stitch_with_truth(*args, **kwargs):
    """Per-primary truth stitching (Scheme B, implemented in Stage 2)."""
    raise NotImplementedError(
        "Per-primary truth bookkeeping is the Stage 2 deliverable; "
        "Stage 1 only provides stitch_count_only()."
    )
