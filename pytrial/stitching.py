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


def event_energies(edeps):
    """Total deposited energy per event (MeV)."""
    return np.asarray(edeps, dtype=np.float64).sum(axis=1)


def background_vectors(counts, mods, energies, n_mod=N_MOD):
    """Split flat background draws into per-event vectors and sparse lists.

    ``counts``/``mods``/``energies`` are the flat arrays returned by
    ``bkg_est_add.add_background.draw_injections``.  Returns
    ``(vectors, mods_list, energies_list)`` where ``vectors[i]`` is a float32
    array of length ``n_mod`` and the lists hold the raw per-event hits.
    """
    counts = np.asarray(counts, dtype=np.int64)
    mods = np.asarray(mods, dtype=np.int64)
    energies = np.asarray(energies, dtype=np.float64)
    vectors = np.zeros((len(counts), n_mod), dtype=np.float32)
    mods_list, energies_list = [], []
    start = 0
    for i, k in enumerate(counts):
        k = int(k)
        stop = start + k
        m = mods[start:stop].astype(np.int64)
        e = energies[start:stop].astype(np.float32)
        mods_list.append(m)
        energies_list.append(e)
        if k:
            np.add.at(vectors[i], m, e)
        start = stop
    return vectors, mods_list, energies_list


def write_multiparticle_sample(
    path,
    edeps,
    ntrack,
    species,
    source_idx,
    primary_energy,
    source_energy,
    signal_energy,
    bkg_energy_total,
    nbkg,
    truth_crystal,
    truth_primary,
    truth_energy,
    bkg_crystal,
    bkg_energy,
    separation=None,
    sep_bin=None,
    tree_name=TREE_NAME,
):
    """Write a multi-particle sample with per-primary and background truth.

    ``edeps`` is the summed signal+background map used by the reconstruction.
    The sparse ``Truth*`` branches hold the signal contribution of every
    primary (``TruthPrimary`` is the primary index), and ``Bkg*`` holds the
    background hits.  Set ``nbkg=0`` and empty background lists for a
    signal-only control sample.

    For separation studies, pass ``separation`` (mm) and ``sep_bin``; the
    writer then adds ``Separation`` and ``SepBin`` branches.
    """
    if (separation is None) != (sep_bin is None):
        raise ValueError("separation and sep_bin must be provided together")
    has_sep = separation is not None

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    edeps = np.asarray(edeps, dtype=np.float32)
    ntrack = np.asarray(ntrack, dtype=np.int32)
    species = np.asarray(species, dtype=np.int32)
    source_idx = np.asarray(source_idx, dtype=np.int32)
    primary_energy = np.asarray(primary_energy, dtype=np.float32)
    source_energy = np.asarray(source_energy, dtype=np.float32)
    signal_energy = np.asarray(signal_energy, dtype=np.float32)
    bkg_energy_total = np.asarray(bkg_energy_total, dtype=np.float32)
    nbkg = np.asarray(nbkg, dtype=np.int32)
    if has_sep:
        separation = np.asarray(separation, dtype=np.float32)
        sep_bin = np.asarray(sep_bin, dtype=np.int32)

    f = ROOT.TFile.Open(str(path), "RECREATE")
    tree = ROOT.TTree(tree_name, "stitched multi-particle events")
    vec = vector("float")()
    nt = array.array("i", [0])
    spec = vector("int")()
    src = vector("int")()
    pe = vector("float")()
    se = vector("float")()
    sig_e = array.array("f", [0.0])
    bkg_e = array.array("f", [0.0])
    nb = array.array("i", [0])
    tc = vector("int")()
    tp = vector("int")()
    te = vector("float")()
    bc = vector("int")()
    be = vector("float")()
    sep_arr = array.array("f", [0.0]) if has_sep else None
    bin_arr = array.array("i", [0]) if has_sep else None

    tree.Branch("Edeps", vec)
    tree.Branch("NTrack", nt, "NTrack/I")
    tree.Branch("Species", spec)
    tree.Branch("SourceIdx", src)
    tree.Branch("PrimaryEnergy", pe)
    tree.Branch("SourceEnergy", se)
    tree.Branch("SignalEnergy", sig_e, "SignalEnergy/F")
    tree.Branch("BkgEnergyTotal", bkg_e, "BkgEnergyTotal/F")
    tree.Branch("NBkg", nb, "NBkg/I")
    tree.Branch("TruthCrystal", tc)
    tree.Branch("TruthPrimary", tp)
    tree.Branch("TruthEnergy", te)
    tree.Branch("BkgCrystal", bc)
    tree.Branch("BkgEnergy", be)
    if has_sep:
        tree.Branch("Separation", sep_arr, "Separation/F")
        tree.Branch("SepBin", bin_arr, "SepBin/I")

    for i in range(len(edeps)):
        vec.clear()
        for value in edeps[i]:
            vec.push_back(float(value))
        nt[0] = int(ntrack[i])

        for name, buf, row in (("species", spec, species[i]),
                               ("source", src, source_idx[i])):
            buf.clear()
            for value in row:
                buf.push_back(int(value))
        for buf, row in ((pe, primary_energy[i]), (se, source_energy[i])):
            buf.clear()
            for value in row:
                buf.push_back(float(value))

        sig_e[0] = float(signal_energy[i])
        bkg_e[0] = float(bkg_energy_total[i])
        nb[0] = int(nbkg[i])
        if has_sep:
            sep_arr[0] = float(separation[i])
            bin_arr[0] = int(sep_bin[i])

        tc.clear()
        for value in truth_crystal[i]:
            tc.push_back(int(value))
        tp.clear()
        for value in truth_primary[i]:
            tp.push_back(int(value))
        te.clear()
        for value in truth_energy[i]:
            te.push_back(float(value))
        bc.clear()
        for value in bkg_crystal[i]:
            bc.push_back(int(value))
        be.clear()
        for value in bkg_energy[i]:
            be.push_back(float(value))
        tree.Fill()

    tree.Write()
    f.Close()
    return path

