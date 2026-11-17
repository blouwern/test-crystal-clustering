#!/usr/bin/env python3
"""Steps 1-6 of the neutron-background estimate: run the model and record it.

Usage
-----
    python estimate_neutron.py                  # print report + write JSON
    python estimate_neutron.py --plot neutron_depth.png
    python estimate_neutron.py --json out.json

The default JSON output is ``neutron_estimation.json`` next to this script; it
is the authoritative record of the numbers used by
``add_neutron_background.py``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import neutron_model

HERE = Path(__file__).resolve().parent


def print_report(result: dict) -> None:
    geo = result["geometry"]
    rate = result["rate"]
    cap = result["capture"]
    gam = result["gamma"]
    dep = result["deposition"]
    en = result["energy"]

    line = "=" * 78
    print(line)
    print(" Neutron-background estimate for the MACE ECAL  (steps 1-6)")
    print(" CONSERVATIVE UPPER BOUND: all neutrons thermal, captured with P = 1")
    print(line)

    print(" 1. Detector geometry")
    print(f"      inner / outer radius                  {geo['inner_radius_mm']:>10.1f}"
          f" / {geo['outer_radius_mm']:.1f} mm")
    print(f"      crystal length                        {geo['crystal_length_mm']:>10.1f} mm")
    print(f"      shell surface S_out                   {geo['shell_surface_cm2']:>10.1f} cm^2")
    print(f"      effective area S/4                    {geo['effective_area_cm2']:>10.1f} cm^2")
    print(f"      crystals / ideal faces                {geo['n_crystals']:>10d}"
          f" / {geo['n_faces_ideal']}")
    print(f"      geometric coverage eps                {geo['geometric_coverage']:>10.4f}")
    print(f"      molecular density                     {geo['molecular_density_cm3']:>10.4e} cm^-3")
    print(f"      total CsI mass                        {geo['total_mass_kg']:>10.1f} kg")

    print(" 2. Incident rate -> captures per event  (P_cap = 1)")
    print(f"      flux Phi                              {rate['flux_n_cm2_s']:>10.1f} n/cm^2/s")
    print(f"      r_incident = Phi * S/4                {rate['r_incident_per_s']:>10.4e} s^-1")
    print(f"      r_capture  = r_incident * eps         {rate['r_capture_per_s']:>10.4e} s^-1")
    print(f"      window                                {rate['event_window_ns']:>10.1f} ns")
    print(f"      n_capture = r_capture * window        {rate['n_capture_per_event']:>10.4f}")

    print(" 3. Capture mean free path -> depth distribution")
    print(f"      sigma_cap (Cs + I)                    {cap['sigma_cap_b_per_molecule']:>10.1f} b"
          f"   ({cap['sigma_cap_components_b']['Cs-133']} + "
          f"{cap['sigma_cap_components_b']['I-127']})")
    print(f"      Sigma_cap                             {cap['macro_cross_section_cm_inv']:>10.4f} cm^-1")
    print(f"      lambda_capture                        {cap['lambda_capture_cm']:>10.3f} cm")
    print(f"      fraction within 1 / 3 lambda          "
          f"{cap['fraction_within_1_lambda']:.3f} / {cap['fraction_within_3_lambda']:.3f}")

    print(" 4. Capture-gamma mean free path")
    print(f"      mu/rho at 6.86 MeV                    {gam['mu_over_rho_cm2_g']:>10.4f} cm^2/g")
    print(f"      mu = rho * mu/rho                     {gam['mu_cm_inv']:>10.4f} cm^-1")
    print(f"      lambda_gamma                          {gam['lambda_gamma_cm']:>10.3f} cm")

    print(" 5. Deposition fraction  f_dep = int g(d) p(d) dd")
    print(f"      inward hemisphere                     {dep['inward_hemisphere']:>10.4f}")
    print(f"      outward hemisphere                    {dep['outward_hemisphere']:>10.4f}")
    print(f"      f_dep (closed form)                   {dep['f_dep']:>10.4f}")
    print(f"      f_dep (numeric cross-check)           {dep['f_dep_numeric']:>10.4f}")

    print(" 6. Energies")
    print(f"      Q_gamma                               {en['q_gamma_mev']:>10.3f} MeV")
    print(f"      E_gamma = Q_gamma * f_dep             {en['energy_per_capture_mev']:>10.4f} MeV"
          f"   (per capture)")
    print(f"      E_evt = n_capture * E_gamma           {en['energy_per_event_mev']:>10.4f} MeV"
          f"   (per event)")
    print(f"      total gamma energy per event          {en['total_gamma_energy_per_event_mev']:>10.4f} MeV")
    print(f"      escaped energy per event              {en['escaped_energy_per_event_mev']:>10.4f} MeV")
    print(f"      hard ceiling per capture              {en['hard_ceiling_per_capture_mev']:>10.3f} MeV")
    print(line)
    print(" NOTE: the result is an upper bound.  Every incident neutron is assumed")
    print("       thermal and captured once; no hadronic shower, no elastic recoil")
    print("       and no (n,n'gamma) channel is included.")
    print(line)


def save_plot(result: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    cap = result["capture"]
    dep = result["deposition"]
    geo = result["geometry"]

    lam_c = cap["lambda_capture_cm"]
    mu = result["gamma"]["mu_cm_inv"]
    L = geo["crystal_length_cm"]
    d = np.linspace(0.0, L, 400)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8), dpi=140)

    ax = axes[0]
    ax.plot(d, np.exp(-d / lam_c) / lam_c, color="tab:blue", lw=1.6)
    for k in (1, 3):
        ax.axvline(k * lam_c, color="tab:red", ls="--", lw=1.0)
        ax.text(k * lam_c + 0.4, 0.9 / lam_c, f"{k}$\\lambda_c$",
                color="tab:red", fontsize=8, va="top")
    ax.set_xlabel("capture depth d [cm]")
    ax.set_ylabel("p(d) [cm$^{-1}$]")
    ax.set_title(f"Capture depth: $\\lambda_c$ = {lam_c:.2f} cm")
    ax.set_xlim(0, L)

    ax = axes[1]
    g = 0.5 * (1.0 - np.exp(-mu * (L - d))) + 0.5 * (1.0 - np.exp(-mu * d))
    ax.plot(d, g, color="tab:green", lw=1.6, label="$g(d)$")
    ax.axhline(dep["f_dep"], color="tab:red", ls="--", lw=1.2,
               label=f"$f_{{dep}}$ = {dep['f_dep']:.4f}")
    ax.set_xlabel("capture depth d [cm]")
    ax.set_ylabel("deposited fraction")
    ax.set_title("Deposition fraction vs capture depth")
    ax.set_xlim(0, L)
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json", type=Path, default=HERE / "neutron_estimation.json",
                        help="output JSON (default: neutron_estimation.json next to this script)")
    parser.add_argument("--plot", type=Path, default=None,
                        help="optional output PNG (capture depth and deposition fraction)")
    parser.add_argument("--flux", type=float, default=neutron_model.FLUX_N_CM2_S,
                        help="neutron flux [n/cm^2/s] (default: 200)")
    parser.add_argument("--window", type=float, default=neutron_model.WINDOW_NS,
                        help="event window [ns] (default: 100)")
    args = parser.parse_args()

    result = neutron_model.estimate(flux=args.flux, window_ns=args.window)
    print_report(result)
    out = neutron_model.write_json(result, args.json)
    print(f"wrote {out}")
    if args.plot is not None:
        save_plot(result, args.plot)
        print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()
