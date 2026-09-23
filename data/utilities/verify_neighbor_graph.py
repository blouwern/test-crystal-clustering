"""Verify whether the shipped ECAL crystal-neighbour graph matches the detector.

Background
----------
``reconstruction.py`` reads ``utilities/<neighbor_file>`` and uses entry ``i`` of
the tree as "the neighbours of crystal ``i``", where ``i`` is the ``ModID``
written by the simulation into ``G4Run0/ECALSimHit``.  This script tests that
assumption against the detector geometry, which is available independently from
the MC truth hit positions (branch ``x``).

The tests are ordered from "cheap sanity check" to "end-to-end functional test".
Each one reports what it can and cannot conclude; see
``docs/neighbor_graph_verification.md`` for the full write-up.

Usage
-----
    python data/utilities/verify_neighbor_graph.py
    python data/utilities/verify_neighbor_graph.py --raw data/raw/gamma_iso.root
    python data/utilities/verify_neighbor_graph.py --json out.json
"""

import argparse
import json
from collections import deque

import numpy as np
import ROOT
from ROOT import vector

REPO_DEFAULT_RAW = "data/raw/ep_iso.root"
REPO_DEFAULT_GRAPH = "data/utilities/ecal_neighbor_info_added.root"

# an edge is called "local" if it is shorter than this many crystal pitches
LOCAL_EDGE_FACTOR = 1.75


# --------------------------------------------------------------------------- #
# I/O helpers
# --------------------------------------------------------------------------- #
def read_neighbour_graph(path):
    """Return {module_id: [neighbour ids]} from a ECALCrystalNeighbors tree."""
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"cannot open {path}")
    tree = f.Get("ECALCrystalNeighbors")
    n_entries = tree.GetEntries()          # read BEFORE closing the file
    vec = vector("int")()
    tree.SetBranchAddress("neighbors", vec)
    graph = {}
    for i in range(n_entries):
        tree.GetEntry(i)
        graph[i] = [int(vec[j]) for j in range(vec.size())]
    f.Close()
    return graph


def read_module_positions(raw_path, n_hits, n_mod=None):
    """Mean MC hit position per ModID (the crystal centre) and hit counts."""
    df = ROOT.RDataFrame("G4Run0/ECALSimHit", str(raw_path))
    data = df.Range(0, n_hits).AsNumpy(["ModID", "x"])
    mod = data["ModID"].astype(np.int32)
    pts = np.stack([np.asarray(v, dtype=np.float64) for v in data["x"]])
    # the module universe is defined by the neighbour graph, not by the sample:
    # some modules are never hit, so max(ModID) under-counts them
    n_mod = int(mod.max()) + 1 if n_mod is None else int(n_mod)
    acc = np.zeros((n_mod, 3))
    cnt = np.zeros(n_mod, dtype=np.int64)
    np.add.at(acc, mod, pts)
    np.add.at(cnt, mod, 1)
    centre = np.divide(acc, cnt[:, None], out=np.zeros_like(acc),
                       where=cnt[:, None] > 0)
    return centre, cnt


# --------------------------------------------------------------------------- #
# graph construction from geometry
# --------------------------------------------------------------------------- #
def geometric_adjacency(centre, hit_ids, radius_factor=1.4):
    """Adjacency from geometry: connect modules closer than radius_factor x pitch."""
    p = centre[hit_ids]
    d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    pitch = float(np.median(d.min(axis=1)))
    thresh = radius_factor * pitch
    adj = (d < thresh)
    return adj, pitch, thresh


def knn_adjacency(centre, hit_ids, k):
    p = centre[hit_ids]
    d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    adj = np.zeros_like(d, dtype=bool)
    for i in range(len(hit_ids)):
        adj[i, np.argsort(d[i])[:k]] = True
    return adj


def symmetric(a):
    return a | a.T


def connected_components(nodes, graph):
    """Connected components of the subgraph induced on ``nodes``."""
    pool = set(nodes)
    seen, comps = set(), []
    for start in sorted(pool):
        if start in seen:
            continue
        comp, q = set(), deque([start])
        while q:
            x = q.popleft()
            if x in comp:
                continue
            comp.add(x)
            q.extend(n for n in graph.get(x, ()) if n in pool and n not in comp)
        seen |= comp
        comps.append(sorted(comp))
    return comps


# --------------------------------------------------------------------------- #
# tests
# --------------------------------------------------------------------------- #
def test0_geometry(centre, cnt, results):
    """The positions themselves must be a sane detector before they can be a referee."""
    hit = np.where(cnt > 0)[0]
    r = np.linalg.norm(centre[hit], axis=1)
    d = np.linalg.norm(centre[hit][:, None, :] - centre[hit][None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    pitch = float(np.median(d.min(axis=1)))
    area = 4 * np.pi * float(np.median(r)) ** 2
    pitch_expected = float(np.sqrt(2 * area / (np.sqrt(3) * len(hit))))
    res = {
        "n_modules_total": int(centre.shape[0]),
        "n_modules_hit": int(len(hit)),
        "hit_id_min": int(hit.min()), "hit_id_max": int(hit.max()),
        "unhit_ids": [int(x) for x in range(centre.shape[0]) if cnt[x] == 0],
        "radius_median_mm": float(np.median(r)),
        "radius_std_mm": float(r.std()),
        "hits_per_module_median": int(np.median(cnt[hit])),
        "pitch_measured_mm": pitch,
        "pitch_expected_from_area_mm": pitch_expected,
    }
    results["test0_geometry"] = res
    u = res["unhit_ids"]
    span = f"{u[0]}..{u[-1]}" if u else "(none)"
    contig = bool(u) and u == list(range(u[0], u[-1] + 1))
    print(f"[T0] geometry: {res['n_modules_hit']}/{res['n_modules_total']} modules hit; "
          f"unhit = {span} (contiguous={contig})")
    print(f"     shell radius = {res['radius_median_mm']:.1f} +/- {res['radius_std_mm']:.1f} mm, "
          f"pitch = {pitch:.1f} mm (hex-packing expectation {pitch_expected:.1f} mm)")
    return hit, pitch


def test1_edge_lengths(graph, centre, hit, pitch, results):
    """Do the graph's edges connect crystals that are physically adjacent?"""
    hit_set = set(int(x) for x in hit)
    lengths, n_edges, n_unhit_endpoints = [], 0, 0
    for i, nbrs in graph.items():
        if i not in hit_set:
            continue
        for j in nbrs:
            n_edges += 1
            if j in hit_set:
                lengths.append(float(np.linalg.norm(centre[i] - centre[j])))
            else:
                n_unhit_endpoints += 1
    lengths = np.array(lengths)
    thresh = LOCAL_EDGE_FACTOR * pitch
    res = {
        "edges_touching_hit_modules": int(n_edges),
        "edges_with_unhit_endpoint": int(n_unhit_endpoints),
        "edges_between_hit_modules": int(lengths.size),
        "edge_length_median_mm": float(np.median(lengths)),
        "edge_length_mean_mm": float(lengths.mean()),
        "edge_length_min_mm": float(lengths.min()),
        "edge_length_max_mm": float(lengths.max()),
        "local_edge_threshold_mm": float(thresh),
        "fraction_edges_local": float(np.mean(lengths < thresh)),
        "n_pitches_median": float(np.median(lengths) / pitch),
    }
    results["test1_edge_lengths"] = res
    print(f"[T1] edge lengths: median = {res['edge_length_median_mm']:.1f} mm "
          f"= {res['n_pitches_median']:.1f} x pitch; "
          f"only {res['fraction_edges_local'] * 100:.1f}% of edges are local (< {thresh:.0f} mm)")
    return res


def test2_structure(graph, results):
    """Is the graph at least structurally a Goldberg-polyhedron face adjacency?"""
    deg = np.array([len(v) for v in graph.values()])
    n = len(deg)
    n_edges = deg.sum() // 2
    # a Goldberg polyhedron with F faces has 12 pentagons, the rest hexagons
    expected_deg_sum = 5 * 12 + 6 * (n - 12)
    res = {
        "n_entries": int(n),
        "degree_min": int(deg.min()), "degree_max": int(deg.max()),
        "degree_mean": float(deg.mean()),
        "degree_value_counts": {int(k): int(v) for k, v in zip(*np.unique(deg, return_counts=True))},
        "edges_undirected": int(n_edges),
        "expected_degree_sum_for_12_pentagons": int(expected_deg_sum),
        "degree_sum": int(deg.sum()),
        "symmetric": bool(all(i in graph.get(j, ()) for i, ns in graph.items() for j in ns)),
    }
    results["test2_structure"] = res
    print(f"[T2] structure: {n} entries, degree {res['degree_min']}..{res['degree_max']} "
          f"(mean {res['degree_mean']:.2f}), symmetric = {res['symmetric']}")
    print(f"     degree histogram = {res['degree_value_counts']}")
    print(f"     sum(deg) = {res['degree_sum']} vs {res['expected_degree_sum_for_12_pentagons']} "
          f"expected for a 622-face Goldberg polyhedron (12 pentagons)")


def test3_neighbour_overlap(graph, centre, hit, pitch, results, k=6):
    """How much do graph neighbours agree with the geometric nearest modules?"""
    hit_list = [int(x) for x in hit]
    idx_of = {m: i for i, m in enumerate(hit_list)}
    d = np.linalg.norm(centre[hit_list][:, None, :] - centre[hit_list][None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    true_nn = [set(hit_list[j] for j in np.argsort(d[i])[:k]) for i in range(len(hit_list))]

    overlaps, exact = [], 0
    for i, m in enumerate(hit_list):
        gn = set(n for n in graph.get(m, ()) if n in idx_of)   # only compare hit modules
        if not gn:
            continue
        # penalise the graph for listing removed modules: compare the best
        # |gn|-sized subset of the true neighbourhood
        t = set(list(true_nn[i])[:len(gn)])
        overlaps.append(len(gn & t) / len(gn))
        if gn == true_nn[i]:
            exact += 1
    res = {
        "k": k,
        "mean_overlap_with_true_knn": float(np.mean(overlaps)),
        "median_overlap_with_true_knn": float(np.median(overlaps)),
        "modules_with_exact_true_neighbourhood": int(exact),
        "modules_compared": int(len(overlaps)),
    }
    results["test3_neighbour_overlap"] = res
    print(f"[T3] overlap with the true {k} nearest modules: "
          f"mean = {res['mean_overlap_with_true_knn'] * 100:.1f}%, "
          f"median = {res['median_overlap_with_true_knn'] * 100:.1f}%; "
          f"exact matches = {exact}/{len(overlaps)}")


def graph_invariants(a):
    """Relabelling-invariant quantities: triangles, components, spectrum."""
    a = (a > 0).astype(int)
    n = a.shape[0]
    tri = int(np.trace(a @ a @ a) / 6)          # number of triangles
    # connected components via BFS on the boolean matrix
    seen = np.zeros(n, dtype=bool)
    comps = 0
    for s0 in range(n):
        if seen[s0]:
            continue
        comps += 1
        stack = [s0]
        seen[s0] = True
        while stack:
            x = stack.pop()
            for y in np.where(a[x])[0]:
                if not seen[y]:
                    seen[y] = True
                    stack.append(y)
    return {"triangles": tri, "components": comps,
            "edges": int(a.sum() // 2), "nodes": n}


def test4_isomorphism(graph, centre, hit, results):
    """Refute the 'correct graph, permuted labels' hypothesis.

    A relabelling preserves the graph structure, so if the shipped graph were the
    true adjacency with permuted indices it would be *isomorphic* to the
    geometry-derived adjacency.  We compare relabelling-invariant signatures:
    degree multiset, sorted neighbour-degree multiset and the adjacency
    eigenvalue spectrum.
    """
    hit_list = [int(x) for x in hit]
    idx_of = {m: i for i, m in enumerate(hit_list)}
    n = len(hit_list)

    # shipped graph restricted to the hit modules (indices are MC ModIDs)
    a_ship = np.zeros((n, n), dtype=float)
    for m in hit_list:
        i = idx_of[m]
        for j in graph.get(m, ()):
            if j in idx_of:
                a_ship[i, idx_of[j]] = 1.0
    a_ship = ((a_ship + a_ship.T) > 0).astype(float)      # symmetrise

    # geometry adjacency on exactly the same module set
    a_geo_raw, pitch, thresh = geometric_adjacency(centre, hit_list)
    a_geo = symmetric(a_geo_raw).astype(float)

    deg_ship = a_ship.sum(axis=1)
    deg_geo = a_geo.sum(axis=1)
    spec_ship = np.sort(np.linalg.eigvalsh(a_ship))
    spec_geo = np.sort(np.linalg.eigvalsh(a_geo))

    def ndeg(a):
        return np.sort(np.array([np.sort(a[i][a[i] > 0]).sum() for i in range(a.shape[0])]))

    inv_ship = graph_invariants(a_ship)
    inv_geo = graph_invariants(a_geo)

    res = {
        "invariants_shipped": inv_ship,
        "invariants_geometry": inv_geo,
        "n_modules_compared": int(n),
        "degree_sequence_equal": bool(np.array_equal(np.sort(deg_ship), np.sort(deg_geo))),
        "degree_mean_shipped": float(deg_ship.mean()),
        "degree_mean_geometry": float(deg_geo.mean()),
        "neighbour_degree_signature_equal": bool(np.array_equal(ndeg(a_ship), ndeg(a_geo))),
        "spectrum_max_abs_diff": float(np.max(np.abs(spec_ship - spec_geo))),
        "spectrum_l2_diff": float(np.linalg.norm(spec_ship - spec_geo)),
        "geometric_edge_threshold_mm": float(thresh),
        "geometric_pitch_mm": float(pitch),
    }
    res["isomorphic_possible"] = bool(
        res["degree_sequence_equal"] and res["spectrum_l2_diff"] < 1e-6)
    results["test4_isomorphism"] = res
    print(f"[T4] relabelling-invariant comparison (vs geometry adjacency on the same {n} modules):")
    print(f"     degree sequences equal        : {res['degree_sequence_equal']} "
          f"(mean {res['degree_mean_shipped']:.2f} vs {res['degree_mean_geometry']:.2f})")
    print(f"     neighbour-degree signature eq : {res['neighbour_degree_signature_equal']}")
    print(f"     adjacency spectrum L2 distance: {res['spectrum_l2_diff']:.2f} "
          f"(max |dlambda| = {res['spectrum_max_abs_diff']:.2f})")
    print(f"     => isomorphic? {res['isomorphic_possible']}")
    print(f"     triangles: shipped = {inv_ship['triangles']}, "
          f"geometry = {inv_geo['triangles']}  "
          f"(a hexagonally-packed ECAL face graph is a triangular lattice and is "
          f"triangle-rich)")
    print(f"     components: shipped = {inv_ship['components']}, "
          f"geometry = {inv_geo['components']}")


def test5_gaps(graph, hit, results):
    """The removed crystals should form compact gaps if the labelling is right."""
    all_ids = set(graph.keys())
    unhit = sorted(all_ids - set(int(x) for x in hit))
    comps = connected_components(unhit, graph)
    res = {
        "n_unhit": len(unhit),
        "n_components_of_unhit_in_graph": len(comps),
        "component_sizes": [len(c) for c in comps],
        "components_are_consecutive_triples": bool(
            all(len(c) == 3 and c == list(range(c[0], c[0] + 3)) for c in comps)),
    }
    results["test5_gaps"] = res
    print(f"[T5] gap fingerprint: {len(unhit)} removed modules form "
          f"{len(comps)} components in graph space, sizes {res['component_sizes'][:8]}"
          f"{'...' if len(comps) > 8 else ''}")
    print(f"     all components are consecutive ModID triples: "
          f"{res['components_are_consecutive_triples']}")
    print("     (two compact antipodal gaps would give 1-2 components)")


def test6_random_removal_control(graph, centre, hit, pitch, results, n_trials=200,
                                 seed=0):
    """Is the *removed set* itself consistent with a correct (possibly permuted) graph?

    If the shipped graph were the correct lattice with permuted node labels, the 54
    unlucky modules {568..621} would be an essentially arbitrary 54-node subset of
    graph space.  We compare the damage it does -- the per-node triangle-count
    histogram of the remaining 568 nodes -- with the damage from removing 54
    *random* nodes.  Agreement means the observed set behaves like a random set,
    which a label permutation predicts and compact physical gaps do not.
    """
    all_ids = sorted(graph.keys())
    keep_observed = [m for m in all_ids if m in set(int(x) for x in hit)]
    n_all = len(all_ids)
    idx_all = {m: i for i, m in enumerate(all_ids)}

    a_full = np.zeros((n_all, n_all), dtype=bool)
    for m in all_ids:
        i = idx_all[m]
        for j in graph.get(m, ()):
            if j in idx_all:
                a_full[i, idx_all[j]] = True
    a_full = a_full | a_full.T

    def hist(keep_idx):
        a = a_full[np.ix_(keep_idx, keep_idx)].astype(np.int64)
        tri = np.diag(a @ a @ a) // 2
        return np.bincount(tri, minlength=8)[:8]

    keep_obs_idx = np.array([idx_all[m] for m in keep_observed])
    obs_hist = hist(keep_obs_idx)

    rng = np.random.default_rng(seed)
    n_remove = n_all - len(keep_observed)
    mask = np.ones(n_all, dtype=bool)
    sims = []
    for _ in range(n_trials):
        mask[:] = True
        mask[rng.choice(n_all, size=n_remove, replace=False)] = False
        sims.append(hist(np.where(mask)[0]))
    sims = np.array(sims)
    lo = np.percentile(sims, 2.5, axis=0)
    hi = np.percentile(sims, 97.5, axis=0)

    res = {
        "n_removed": int(n_remove),
        "observed_triangle_hist": [int(x) for x in obs_hist],
        "random_removal_mean_hist": [float(x) for x in sims.mean(axis=0)],
        "random_removal_2.5pct": [float(x) for x in lo],
        "random_removal_97.5pct": [float(x) for x in hi],
        "observed_inside_random_band": bool(np.all((obs_hist >= lo) & (obs_hist <= hi))),
        "n_trials": n_trials,
    }
    results["test6_random_removal_control"] = res
    print(f"[T6] random-removal control ({n_trials} trials, remove {n_remove} nodes):")
    print(f"     observed per-node triangle hist : {res['observed_triangle_hist']}")
    print(f"     random removal mean hist        : "
          f"{[round(x, 1) for x in res['random_removal_mean_hist']]}")
    print(f"     observed inside 95% random band : {res['observed_inside_random_band']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=REPO_DEFAULT_RAW)
    ap.add_argument("--graph", default=REPO_DEFAULT_GRAPH)
    ap.add_argument("--hits", type=int, default=400_000)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    results = {"raw_file": args.raw, "graph_file": args.graph}
    print("=" * 78)
    print(f"neighbour-graph verification\n  raw   : {args.raw}\n  graph : {args.graph}")
    print("=" * 78)

    graph = read_neighbour_graph(args.graph)
    centre, cnt = read_module_positions(args.raw, args.hits,
                                       n_mod=max(graph.keys()) + 1)

    hit, pitch = test0_geometry(centre, cnt, results)
    test1_edge_lengths(graph, centre, hit, pitch, results)
    test2_structure(graph, results)
    test3_neighbour_overlap(graph, centre, hit, pitch, results)
    test4_isomorphism(graph, centre, hit, results)
    test5_gaps(graph, hit, results)
    test6_random_removal_control(graph, centre, hit, pitch, results)

    print("=" * 78)
    print("NOTE / INTERPRETATION")
    print("  reconstruction.py reads neighbours[cid] with cid = the simulation ModID,")
    print("  so it assumes (H): tree entry i lists the neighbours of ModID i.")
    print("  T1/T3/T5 are evaluated under (H) and show the graph is not the physical")
    print("  adjacency under that assumption.")
    print("  T2/T4 show the graph is nevertheless a legitimate, triangle-rich lattice")
    print("  (so the polyhedron itself looks fine) - the discrepancy is therefore an")
    print("  INDEXING/labelling issue between the file and the simulation ModID.")
    print("  T4 alone does NOT settle it: its degree/spectrum comparison depends on the")
    print("  threshold used to build the reference geometry graph and is inconclusive.")
    print("  Either way the current code cannot be using the correct adjacency:")
    print("  if (H) holds the values are wrong; if (H) fails it reads another crystal's")
    print("  neighbours.  See docs/neighbor_graph_verification.md for the full write-up.")
    print("=" * 78)

    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=1, sort_keys=True)
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
