"""I/O helpers for multi-particle samples and cluster files.

The stitched sample schema is written by :func:`pytrial.stitching.write_multiparticle_sample`:

* ``Edeps``, ``NTrack``
* ``Species``, ``SourceIdx``, ``PrimaryEnergy``, ``SourceEnergy``
* ``SignalEnergy``, ``BkgEnergyTotal``, ``NBkg``
* ``TruthCrystal`` / ``TruthPrimary`` / ``TruthEnergy`` (sparse signal truth)
* ``BkgCrystal`` / ``BkgEnergy`` (sparse background truth)

``load_multiparticle_sample`` turns the sparse truth into one ``(N, 568)``
signal array and one ``(568,)`` background array per event, which is what the
evaluation modules consume.
"""

from __future__ import annotations

import array
from pathlib import Path

import numpy as np
import ROOT
from ROOT import vector

N_MOD = 568


def load_multiparticle_sample(path, n_mod=N_MOD):
    """Load a stitched sample with per-primary and background truth."""
    n_events = 0
    with ROOT.TFile.Open(str(path), "READ") as f:
        tree = f.Get("EdepOfEachEvt")
        n_events = tree.GetEntries()

        ntrack = array.array("i", [0])
        nbkg = array.array("i", [0])
        sig_e = array.array("f", [0.0])
        bkg_e = array.array("f", [0.0])
        species = vector("int")()
        source = vector("int")()
        primary_e = vector("float")()
        source_e = vector("float")()
        tc = vector("int")()
        tp = vector("int")()
        te = vector("float")()
        bc = vector("int")()
        be = vector("float")()

        tree.SetBranchAddress("NTrack", ntrack)
        tree.SetBranchAddress("NBkg", nbkg)
        tree.SetBranchAddress("SignalEnergy", sig_e)
        tree.SetBranchAddress("BkgEnergyTotal", bkg_e)
        tree.SetBranchAddress("Species", species)
        tree.SetBranchAddress("SourceIdx", source)
        tree.SetBranchAddress("PrimaryEnergy", primary_e)
        tree.SetBranchAddress("SourceEnergy", source_e)
        tree.SetBranchAddress("TruthCrystal", tc)
        tree.SetBranchAddress("TruthPrimary", tp)
        tree.SetBranchAddress("TruthEnergy", te)
        tree.SetBranchAddress("BkgCrystal", bc)
        tree.SetBranchAddress("BkgEnergy", be)

        branch_names = {b.GetName() for b in tree.GetListOfBranches()}
        has_sep = ("Separation" in branch_names) and ("SepBin" in branch_names)
        sep_arr = array.array("f", [0.0])
        bin_arr = array.array("i", [0])
        if has_sep:
            tree.SetBranchAddress("Separation", sep_arr)
            tree.SetBranchAddress("SepBin", bin_arr)

        out = {
            "path": str(path),
            "n_events": int(n_events),
            "ntrack": np.zeros(n_events, dtype=np.int32),
            "species": [],
            "source_idx": [],
            "primary_energy": [],
            "source_energy": [],
            "signal_energy": np.zeros(n_events, dtype=np.float32),
            "bkg_energy_total": np.zeros(n_events, dtype=np.float32),
            "nbkg": np.zeros(n_events, dtype=np.int32),
            "signal_arrays": [],
            "bkg_arrays": [],
            "separation": np.full(n_events, np.nan, dtype=np.float32),
            "sep_bin": np.full(n_events, -1, dtype=np.int32),
        }

        for i in range(n_events):
            tree.GetEntry(i)
            nt = int(ntrack[0])
            out["ntrack"][i] = nt
            out["nbkg"][i] = int(nbkg[0])
            out["signal_energy"][i] = float(sig_e[0])
            out["bkg_energy_total"][i] = float(bkg_e[0])
            out["species"].append(np.array(species, dtype=np.int32))
            out["source_idx"].append(np.array(source, dtype=np.int32))
            out["primary_energy"].append(np.array(primary_e, dtype=np.float32))
            out["source_energy"].append(np.array(source_e, dtype=np.float32))
            if has_sep:
                out["separation"][i] = float(sep_arr[0])
                out["sep_bin"][i] = int(bin_arr[0])

            sig = np.zeros((nt, n_mod), dtype=np.float32)
            for c, p, e in zip(tc, tp, te):
                sig[int(p), int(c)] += float(e)
            out["signal_arrays"].append(sig)

            bkg = np.zeros(n_mod, dtype=np.float32)
            for c, e in zip(bc, be):
                bkg[int(c)] += float(e)
            out["bkg_arrays"].append(bkg)

    return out


def load_clusters(cluster_path, n_events):
    """Return ``clusters[event] = [[mod ids], ...]`` for ``n_events`` events."""
    clusters = [[] for _ in range(int(n_events))]
    with ROOT.TFile.Open(str(cluster_path), "READ") as f:
        tree = f.Get("ECALClusters")
        for entry in tree:
            evt = int(entry.EvtID)
            if 0 <= evt < len(clusters):
                clusters[evt].append([int(m) for m in entry.ModIDList])
    return clusters
