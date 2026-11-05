#!/usr/bin/env python3
"""Independent geometric cross-check of the GP(8,0) crystal volume.

The estimate uses the average face area of an ideal GP(8,0) polyhedron because
that is what the task prescribes.  The actual detector, however, is described
by the GDML geometry embedded in the raw simulations
(``SimMACEPhaseI_gdml`` in ``data/raw/*.root``).  This script parses that GDML,
computes the inner/outer face area and the volume of every crystal from its
tessellated solid, and compares the averages with the model numbers.

Run from the repository root or from this directory:

    python check_geometry_gdml.py
    python check_geometry_gdml.py --raw ../../data/raw/gamma_iso.root
"""

from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path

import numpy as np
import ROOT

import bkg_model

ROOT.gROOT.SetBatch(True)
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def load_gdml(raw_path: Path) -> str:
    f = ROOT.TFile.Open(str(raw_path), "READ")
    macro = f.Get("SimMACEPhaseI_gdml")
    if macro is None:
        raise RuntimeError(f"{raw_path} has no SimMACEPhaseI_gdml object")
    lines = macro.GetListOfLines()
    text = "\n".join(lines[i].GetName() for i in range(lines.GetEntries()))
    f.Close()
    return text


def _polygon_area(points: list[np.ndarray]) -> float:
    s = np.zeros(3)
    for i, p in enumerate(points):
        s += np.cross(p, points[(i + 1) % len(points)])
    return 0.5 * float(np.linalg.norm(s))


def parse_crystals(text: str):
    """Return per-crystal dicts with inner/outer face areas and volume."""
    positions = {}
    for m in re.finditer(
        r'<position name="(ECALCrystal_\d+)_v(\d+)" unit="mm" '
        r'x="([^"]+)" y="([^"]+)" z="([^"]+)"', text
    ):
        positions[(m.group(1), int(m.group(2)))] = np.array(
            [float(m.group(3)), float(m.group(4)), float(m.group(5))]
        )

    solids = re.findall(
        r'<tessellated[^>]*name="(ECALCrystal_\d+)0x[0-9a-f]+">(.*?)</tessellated>',
        text, re.S,
    )

    crystals = []
    for name, body in solids:
        faces = []
        for _, verts in re.findall(r'<(triangular|quadrangular)([^/]*)/>', body):
            idx = [int(x) for x in re.findall(r'_v(\d+)"', verts)]
            faces.append([positions[(name, i)] for i in idx])

        volume = 0.0
        area_in = area_out = 0.0
        for pts in faces:
            p0 = pts[0]
            for i in range(1, len(pts) - 1):
                volume += float(np.dot(p0, np.cross(pts[i], pts[i + 1]))) / 6.0
            radii = [np.linalg.norm(p) for p in pts]
            if all(r < 400.0 for r in radii):        # inner end cap
                area_in += _polygon_area(pts)
            elif all(r >= 400.0 for r in radii):     # outer end cap
                area_out += _polygon_area(pts)
        crystals.append({
            "name": name,
            "n_vertices": len({k for k in positions if k[0] == name}),
            "area_in": area_in,
            "area_out": area_out,
            "volume": abs(volume),
        })
    return crystals


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--raw", type=Path, default=REPO / "data/raw/ep_iso.root")
    args = parser.parse_args()

    crystals = parse_crystals(load_gdml(args.raw))
    vol = np.array([c["volume"] for c in crystals])
    a_in = np.array([c["area_in"] for c in crystals])
    a_out = np.array([c["area_out"] for c in crystals])

    geo = bkg_model.geometry()
    print(f"GDML source                      {args.raw}")
    print(f"crystals in GDML                 {len(crystals)}")
    print(f"  vertices per crystal           {dict(sorted(Counter(c['n_vertices'] for c in crystals).items()))}")
    print()
    print(f"{'quantity':<32s}{'GDML mean':>14s}{'GP(8,0) model':>18s}{'ratio':>10s}")
    print(f"{'inner face area [mm^2]':<32s}{a_in.mean():>14.2f}{geo['s1_mm2']:>18.2f}"
          f"{a_in.mean() / geo['s1_mm2']:>10.3f}")
    print(f"{'outer face area [mm^2]':<32s}{a_out.mean():>14.2f}{geo['s2_mm2']:>18.2f}"
          f"{a_out.mean() / geo['s2_mm2']:>10.3f}")
    print(f"{'crystal volume [mm^3]':<32s}{vol.mean():>14.1f}{geo['volume_mm3']:>18.1f}"
          f"{vol.mean() / geo['volume_mm3']:>10.3f}")
    print()
    print(f"GDML volume spread               min {vol.min():.0f}  max {vol.max():.0f} mm^3")
    print(f"GDML inner area spread           min {a_in.min():.1f}  max {a_in.max():.1f} mm^2")
    print()
    print("The GDML detector is a 568-crystal subset of the ideal 642-cell")
    print("polyhedron (74 cells absent) and every cell is cut by a plane midway")
    print("to its neighbours, so its mean volume is smaller than the ideal average.")


if __name__ == "__main__":
    main()
