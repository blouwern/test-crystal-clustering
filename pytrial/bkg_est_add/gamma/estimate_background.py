#!/usr/bin/env python3
"""Steps 1-4 of the gamma-background estimate: run the model and record it.

Usage
-----
    python estimate_background.py                 # print report + write JSON
    python estimate_background.py --plot bkg_spectrum.png
    python estimate_background.py --json out.json

The default JSON output is ``bkg_estimation.json`` next to this script; it is
the authoritative record of the numbers used by ``add_background.py``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import bkg_model

HERE = Path(__file__).resolve().parent


def print_report(result: dict) -> None:
    geo = result["geometry"]
    spec = result["spectrum"]
    dose = result["dose"]
    rates = result["rates"]

    line = "=" * 78
    print(line)
    print(" Gamma-background estimate for the MACE ECAL  (steps 1-4)")
    print(line)

    print(" 1. Detector geometry - GP(8,0) first-kind Goldberg polyhedron")
    print(f"      faces                                 {geo['n_faces']:>10d}"
          f"  ({geo['n_pentagons']} pentagons + {geo['n_hexagons']} hexagons)")
    print(f"      vertices / edges                      {geo['n_vertices']:>10d}"
          f" / {geo['n_edges']}   (Euler = {geo['euler_characteristic']})")
    print(f"      inner radius                          {geo['inner_radius_mm']:>10.1f} mm")
    print(f"      crystal length                        {geo['crystal_length_mm']:>10.1f} mm")
    print(f"      mean inner face area S1               {geo['s1_mm2']:>10.3f} mm^2")
    print(f"      mean outer face area S2 (R=540 mm)    {geo['s2_mm2']:>10.3f} mm^2")
    print(f"      crystal volume (exact frustum)        {geo['volume_mm3']:>10.1f} mm^3"
          f" = {geo['volume_cm3']:.3f} cm^3")
    print(f"      CsI mass per crystal                  {geo['mass_kg']:>10.4f} kg")
    print(f"      instrumented crystals (2 gaps)        {geo['n_crystals']:>10d}")
    print(f"      total crystal volume                  {geo['volume_total_cm3']:>10.1f} cm^3")
    print(f"      total CsI mass                        {geo['mass_total_kg']:>10.1f} kg")

    print(" 2. Energy-deposition spectrum - normal, truncated at E > 0")
    print(f"      nominal mean / variance               {spec['nominal_mean_mev']:>10.3f} MeV"
          f" / {spec['variance_mev2']:.3f} MeV^2")
    print(f"      sigma                                 {spec['sigma_mev']:>10.4f} MeV")
    print(f"      mean deposition E_bkg                 {spec['mean_deposition_mev']:>10.4f} MeV")

    print(" 3. Dose rate -> background interaction rate")
    print(f"      P_bkg                                 {dose['p_bkg_usv_per_h']:>10.1f} uSv/h"
          f"   (1 Sv = 1 Gy)")
    print(f"      dose rate                             {dose['dose_rate_gy_per_s']:>10.4e} Gy/s")
    print(f"      power, one crystal                    {dose['power_per_crystal_w']:>10.4e} W")
    print(f"      power, all crystals                   {dose['power_total_w']:>10.4e} W")
    print(f"      r_bkg = P_total / E_bkg               {rates['r_bkg_per_s']:>10.4e} s^-1")

    print(" 4. Event window -> expected hits per event")
    print(f"      window                                {rates['event_window_ns']:>10.1f} ns")
    print(f"      n_bkg = r_bkg * window                {rates['n_bkg_per_event']:>10.4f}")
    print(f"      P(0), P(1), P(2)                      "
          f"{rates['poisson_probability']['0']:.6f}, "
          f"{rates['poisson_probability']['1']:.6f}, "
          f"{rates['poisson_probability']['2']:.6f}")
    print(f"      expected hits in 1e6 events           {rates['expected_hits_per_1e6_events']:>10.1f}")
    print(line)
    print(" NOTE: the rate is computed for the summed volume of all 568")
    print("       instrumented crystals (the ideal GP(8,0) has 642 cells, but")
    print("       the detector has two upstream/downstream gaps).  Step 5 draws")
    print("       n_bkg hits per event and assigns each to a random crystal.")
    print(line)


def save_plot(result: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    spec = bkg_model.DepositionSpectrum(
        result["spectrum"]["nominal_mean_mev"], result["spectrum"]["variance_mev2"]
    )
    rng = np.random.default_rng(20261105)
    samples = spec.rvs(200_000, rng)
    x = np.linspace(0.0, 14.0, 800)

    fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=140)
    ax.hist(samples, bins=140, range=(0, 14), density=True, alpha=0.35,
            color="tab:blue", label="sampled spectrum")
    ax.plot(x, spec.pdf(x), color="tab:blue", lw=1.6, label="normal (5 MeV, var 2 MeV$^2$)")
    ax.axvline(spec.mean_deposition_mev, color="tab:red", ls="--", lw=1.2,
               label=f"$E_{{bkg}}$ = {spec.mean_deposition_mev:.3f} MeV")
    ax.set_xlabel("energy deposition per background gamma [MeV]")
    ax.set_ylabel("probability density [MeV$^{-1}$]")
    ax.set_title("Assumed background gamma deposition spectrum in one CsI crystal")
    ax.set_xlim(0, 14)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json", type=Path, default=HERE / "bkg_estimation.json",
                        help="output JSON (default: bkg_estimation.json next to this script)")
    parser.add_argument("--plot", type=Path, default=None,
                        help="optional output PNG of the assumed deposition spectrum")
    parser.add_argument("--dose", type=float, default=bkg_model.P_BKG_USV_PER_H,
                        help="ambient dose rate [uSv/h] (default: 30)")
    parser.add_argument("--spectrum-mean", type=float, default=bkg_model.SPECTRUM_MEAN_MEV,
                        help="mean deposition energy [MeV] (default: 5)")
    parser.add_argument("--spectrum-variance", type=float, default=bkg_model.SPECTRUM_VAR_MEV2,
                        help="deposition energy variance [MeV^2] (default: 2)")
    parser.add_argument("--window", type=float, default=bkg_model.WINDOW_NS,
                        help="event window [ns] (default: 100)")
    args = parser.parse_args()

    result = bkg_model.estimate(
        dose_usv_per_h=args.dose,
        spectrum_mean_mev=args.spectrum_mean,
        spectrum_var_mev2=args.spectrum_variance,
        window_ns=args.window,
    )
    print_report(result)
    out = bkg_model.write_json(result, args.json)
    print(f"wrote {out}")
    if args.plot is not None:
        save_plot(result, args.plot)
        print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()
