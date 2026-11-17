"""Conservative estimate of the neutron background of the MACE ECAL (steps 1-6).

The estimate deliberately avoids any simulation of the hall neutron field.  It
uses the *most conservative* set of assumptions -- every incident neutron is
thermal and is captured with unit probability -- so that the result is a hard
upper bound on the neutron background rather than a best estimate.

Physics chain
-------------
1. **Geometry** -- the crystal array is a spherical shell of inner radius
   260 mm and outer radius 540 mm.  For an isotropic flux the rate at which
   neutrons enter a convex body is ``Phi * S/4`` (for a sphere this equals the
   unidirectional ``Phi * pi R^2``).  The ideal GP(8, 0) polyhedron has 642
   cells but only 568 crystals are instrumented, so ``eps = 568/642``.
2. **Capture rate** -- with ``P_cap = 1`` every entering neutron is captured
   once, giving ``n = Phi * (S/4) * eps * dt`` captures per event.
3. **Capture depth** -- the thermal capture cross section of CsI is
   ``sigma_cap = 29.0 b (Cs-133) + 6.2 b (I-127) = 35.2 b`` per molecule, hence
   ``Sigma_cap = N_mol * sigma_cap`` and ``lambda_c = 1/Sigma_cap = 2.72 cm``.
   The capture depth therefore follows ``p(d) = exp(-d/lambda_c)/lambda_c``.
4. **Gamma attenuation** -- the 6.86 MeV capture cascade has a mass attenuation
   coefficient ``mu/rho = 0.0353 cm^2/g`` in CsI, so ``mu = 0.159 cm^-1`` and
   ``lambda_gamma = 6.3 cm``.
5. **Deposition fraction** -- a capture at depth ``d`` emits gammas
   isotropically: half travel inwards (through ``L - d``) and half outwards
   (through ``d`` and then out of the array).  Averaging over ``p(d)`` gives
   ``f_dep = 0.6405``.
6. **Energies** -- one capture deposits ``E_gamma = Q_gamma * f_dep = 4.39 MeV``
   and one event deposits ``E_evt = n * E_gamma = 0.71 MeV``.

Why the other channels are negligible
-------------------------------------
Elastic scattering deposits at most ``0.031 E_n`` per collision (about 31 keV
for a 1 MeV neutron) and pure CsI quenches Cs/I recoils by a factor 0.02-0.1,
so the visible energy is at the keV level.  Inelastic ``(n, n'gamma)`` produces
0.1-0.6 MeV gammas a few times per neutron (tens of keV per event).  Both are
1-2 orders of magnitude below the capture channel, and neither can approach the
7 MeV seed threshold, so they are not modelled.

Units
-----
Energies in MeV (same as ``Edeps``), lengths in cm unless a name says ``_mm``,
cross sections in barn (1 b = 1e-24 cm^2).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- #
# Inputs of the estimate
# --------------------------------------------------------------------------- #

#: Hall neutron flux limit, n cm^-2 s^-1 (the only source parameter given).
FLUX_N_CM2_S = 200.0

#: Inner / outer radius of the crystal shell.
R_IN_MM = 260.0
R_OUT_MM = 540.0
CRYSTAL_LENGTH_MM = R_OUT_MM - R_IN_MM

#: Ideal GP(8,0) cells and instrumented crystals (two gaps, upstream/downstream).
N_FACES_IDEAL = 642
N_CRYSTALS = 568

#: Event time window.
WINDOW_NS = 100.0

#: CsI properties.
CSI_DENSITY_G_CM3 = 4.51
CSI_MOLAR_MASS_G_MOL = 259.809          # Cs (132.905) + I (126.904)

#: Thermal capture cross section per CsI molecule, barn (ENDF/B-VIII.0, 0.0253 eV).
SIGMA_CAP_CS_B = 29.0                   # Cs-133
SIGMA_CAP_I_B = 6.2                     # I-127
SIGMA_CAP_TH_B = SIGMA_CAP_CS_B + SIGMA_CAP_I_B

#: Total gamma energy released per capture (Cs 6.89, I 6.82, abundance weighted).
Q_GAMMA_MEV = 0.5115 * 6.89 + 0.4885 * 6.82

#: Mass attenuation coefficient of CsI at 6.86 MeV (NIST XCOM).
MU_RHO_CM2_G = 0.0353

#: Physical constants.
N_A = 6.02214076e23


# --------------------------------------------------------------------------- #
# Step 1: geometry and material constants
# --------------------------------------------------------------------------- #

def geometry() -> dict:
    """Return the geometric and material constants used by the estimate."""
    r_out = R_OUT_MM / 10.0                                  # cm
    r_in = R_IN_MM / 10.0                                    # cm
    length = CRYSTAL_LENGTH_MM / 10.0                        # cm

    s_out = 4.0 * math.pi * r_out ** 2
    s_eff = s_out / 4.0
    eps = N_CRYSTALS / N_FACES_IDEAL
    n_mol = CSI_DENSITY_G_CM3 / CSI_MOLAR_MASS_G_MOL * N_A

    # Crystal volume: exact frustum between the mean inner/outer face areas of
    # GP(8,0).  Re-derived here only for the mass bookkeeping; the detailed
    # polyhedron construction lives in ../gamma/bkg_model.py.
    s1_mm2 = 1319.3796636114332
    s2_mm2 = s1_mm2 * (R_OUT_MM / R_IN_MM) ** 2
    volume_mm3 = CRYSTAL_LENGTH_MM / 3.0 * (s1_mm2 + s2_mm2 + math.sqrt(s1_mm2 * s2_mm2))
    volume_cm3 = volume_mm3 / 1000.0
    mass_kg = volume_cm3 * CSI_DENSITY_G_CM3 / 1000.0

    return {
        "inner_radius_mm": R_IN_MM,
        "outer_radius_mm": R_OUT_MM,
        "crystal_length_mm": CRYSTAL_LENGTH_MM,
        "crystal_length_cm": length,
        "shell_surface_cm2": s_out,
        "effective_area_cm2": s_eff,
        "n_crystals": N_CRYSTALS,
        "n_faces_ideal": N_FACES_IDEAL,
        "geometric_coverage": eps,
        "molecular_density_cm3": n_mol,
        "crystal_volume_cm3": volume_cm3,
        "crystal_mass_kg": mass_kg,
        "total_volume_cm3": volume_cm3 * N_CRYSTALS,
        "total_mass_kg": mass_kg * N_CRYSTALS,
        "density_g_cm3": CSI_DENSITY_G_CM3,
    }


# --------------------------------------------------------------------------- #
# Steps 2-5: rate, capture depth, gamma attenuation, deposition fraction
# --------------------------------------------------------------------------- #

def rates(flux: float = FLUX_N_CM2_S, window_ns: float = WINDOW_NS,
          capture_probability: float = 1.0) -> dict:
    """Step 2: incident rate and captures per event (conservative: P_cap = 1)."""
    geo = geometry()
    r_in = flux * geo["effective_area_cm2"]
    r_cap = r_in * geo["geometric_coverage"] * capture_probability
    window_s = window_ns * 1e-9
    n_capture = r_cap * window_s

    poisson = {
        str(k): float(math.exp(-n_capture) * n_capture ** k / math.factorial(k))
        for k in range(4)
    }
    return {
        "flux_n_cm2_s": float(flux),
        "event_window_ns": float(window_ns),
        "capture_probability": float(capture_probability),
        "r_incident_per_s": float(r_in),
        "r_capture_per_s": float(r_cap),
        "n_capture_per_event": float(n_capture),
        "poisson_probability": poisson,
        "expected_captures_per_1e6_events": float(n_capture * 1e6),
    }


def capture_mean_free_path(sigma_cap_b: float = SIGMA_CAP_TH_B) -> dict:
    """Step 3: macroscopic capture cross section and capture depth distribution."""
    geo = geometry()
    sigma_cm2 = sigma_cap_b * 1e-24
    sigma_macro = geo["molecular_density_cm3"] * sigma_cm2
    lam_c = 1.0 / sigma_macro
    return {
        "sigma_cap_b_per_molecule": float(sigma_cap_b),
        "sigma_cap_components_b": {"Cs-133": SIGMA_CAP_CS_B, "I-127": SIGMA_CAP_I_B},
        "sigma_cap_cm2": float(sigma_cm2),
        "macro_cross_section_cm_inv": float(sigma_macro),
        "lambda_capture_cm": float(lam_c),
        "depth_pdf": "p(d) = exp(-d/lambda_c)/lambda_c,  0 <= d <= L",
        "fraction_within_1_lambda": float(1.0 - math.exp(-1.0)),
        "fraction_within_3_lambda": float(1.0 - math.exp(-3.0)),
    }


def gamma_mean_free_path(mu_rho: float = MU_RHO_CM2_G) -> dict:
    """Step 4: linear attenuation coefficient and gamma mean free path."""
    mu = CSI_DENSITY_G_CM3 * mu_rho
    return {
        "mu_over_rho_cm2_g": float(mu_rho),
        "density_g_cm3": CSI_DENSITY_G_CM3,
        "mu_cm_inv": float(mu),
        "lambda_gamma_cm": float(1.0 / mu),
    }


def deposit_fraction(lam_c: float | None = None, mu: float | None = None,
                     thickness_cm: float | None = None) -> dict:
    """Step 5: average fraction of the capture gammas deposited in the array.

    ``f_dep = integral_0^L g(d) p(d) dd`` with
    ``g(d) = 1/2 [1 - exp(-b(L-d))] + 1/2 [1 - exp(-b d)]``.
    The integral is done in closed form (see the README) and cross-checked
    numerically.
    """
    geo = geometry()
    if lam_c is None:
        lam_c = capture_mean_free_path()["lambda_capture_cm"]
    if mu is None:
        mu = gamma_mean_free_path()["mu_cm_inv"]
    if thickness_cm is None:
        thickness_cm = geo["crystal_length_cm"]

    a, b, L = 1.0 / lam_c, mu, thickness_cm

    # inward hemisphere
    inward = 0.5 * (1.0 - (a / (a - b)) * math.exp(-b * L)
                    * (1.0 - math.exp(-(a - b) * L)))
    # outward hemisphere
    outward = 0.5 * (1.0 - (a / (a + b)) * (1.0 - math.exp(-(a + b) * L)))

    # independent numerical cross-check
    d = np.linspace(0.0, L, 200001)
    g = 0.5 * (1.0 - np.exp(-b * (L - d))) + 0.5 * (1.0 - np.exp(-b * d))
    p = a * np.exp(-a * d)
    numeric = float(np.trapezoid(g * p, d))

    return {
        "a_cm_inv": float(a),
        "b_cm_inv": float(b),
        "thickness_cm": float(L),
        "inward_hemisphere": float(inward),
        "outward_hemisphere": float(outward),
        "f_dep": float(inward + outward),
        "f_dep_numeric": float(numeric),
        "closed_form": "1/2[1 - a/(a-b) e^{-bL}(1-e^{-(a-b)L})] "
                       "+ 1/2[1 - a/(a+b)(1-e^{-(a+b)L})]",
    }


# --------------------------------------------------------------------------- #
# Step 6: energies
# --------------------------------------------------------------------------- #

def energies(n_capture: float | None = None, f_dep: float | None = None,
             q_mev: float = Q_GAMMA_MEV) -> dict:
    """Step 6: deposited energy per capture and per event."""
    if n_capture is None:
        n_capture = rates()["n_capture_per_event"]
    if f_dep is None:
        f_dep = deposit_fraction()["f_dep"]
    e_gamma = q_mev * f_dep
    return {
        "q_gamma_mev": float(q_mev),
        "f_dep": float(f_dep),
        "energy_per_capture_mev": float(e_gamma),
        "n_capture_per_event": float(n_capture),
        "energy_per_event_mev": float(n_capture * e_gamma),
        "total_gamma_energy_per_event_mev": float(n_capture * q_mev),
        "escaped_energy_per_event_mev": float(n_capture * q_mev * (1.0 - f_dep)),
        "hard_ceiling_per_capture_mev": float(q_mev),
        "hard_ceiling_per_event_mev": float(n_capture * q_mev),
    }


# --------------------------------------------------------------------------- #
# Capture-depth sampler used by the injection
# --------------------------------------------------------------------------- #

class CaptureEnergy:
    """Per-capture deposited energy, sampled from the capture depth.

    A capture at depth ``d`` deposits the fraction
    ``g(d) = 1/2[1-exp(-b(L-d))] + 1/2[1-exp(-bd)]`` of ``Q_gamma``.  Drawing
    ``d ~ Exp(lambda_c)`` reproduces the depth-integrated mean ``Q*f_dep`` and
    gives the event-by-event spread for free.
    """

    def __init__(self, q_mev: float = Q_GAMMA_MEV,
                 lam_c_cm: float | None = None,
                 mu_cm_inv: float | None = None,
                 thickness_cm: float | None = None):
        geo = geometry()
        self.q_mev = float(q_mev)
        self.lam_c_cm = float(lam_c_cm if lam_c_cm is not None
                              else capture_mean_free_path()["lambda_capture_cm"])
        self.mu_cm_inv = float(mu_cm_inv if mu_cm_inv is not None
                               else gamma_mean_free_path()["mu_cm_inv"])
        self.thickness_cm = float(thickness_cm if thickness_cm is not None
                                  else geo["crystal_length_cm"])
        self.f_dep = float(deposit_fraction(self.lam_c_cm, self.mu_cm_inv,
                                            self.thickness_cm)["f_dep"])

    def fraction(self, depth_cm) -> np.ndarray:
        d = np.asarray(depth_cm, dtype=np.float64)
        b, L = self.mu_cm_inv, self.thickness_cm
        return 0.5 * (1.0 - np.exp(-b * (L - d))) + 0.5 * (1.0 - np.exp(-b * d))

    def rvs(self, size: int, rng: np.random.Generator) -> np.ndarray:
        d = rng.exponential(self.lam_c_cm, size=int(size))
        d = np.minimum(d, self.thickness_cm)
        return (self.q_mev * self.fraction(d)).astype(np.float32)

    @property
    def mean_mev(self) -> float:
        """Mean deposited energy per capture (equals ``Q_gamma * f_dep``)."""
        return self.q_mev * self.f_dep

    def summary(self) -> dict:
        return {
            "distribution": "Q_gamma * g(d),  d ~ Exp(lambda_c)",
            "q_gamma_mev": self.q_mev,
            "lambda_capture_cm": self.lam_c_cm,
            "mu_cm_inv": self.mu_cm_inv,
            "thickness_cm": self.thickness_cm,
            "mean_deposition_mev": self.mean_mev,
            "energy_range_mev": [float(self.q_mev * self.fraction(0.0)),
                                 float(self.q_mev * self.fraction(self.thickness_cm / 2))],
        }


# --------------------------------------------------------------------------- #
# Full estimate
# --------------------------------------------------------------------------- #

def estimate(flux: float = FLUX_N_CM2_S, window_ns: float = WINDOW_NS,
             sigma_cap_b: float = SIGMA_CAP_TH_B,
             mu_rho: float = MU_RHO_CM2_G,
             q_mev: float = Q_GAMMA_MEV) -> dict:
    """Run steps 1-6 and return every intermediate number."""
    geo = geometry()
    r = rates(flux, window_ns, 1.0)
    cap = capture_mean_free_path(sigma_cap_b)
    gam = gamma_mean_free_path(mu_rho)
    dep = deposit_fraction(cap["lambda_capture_cm"], gam["mu_cm_inv"],
                           geo["crystal_length_cm"])
    en = energies(r["n_capture_per_event"], dep["f_dep"], q_mev)

    return {
        "title": "Conservative upper bound on the neutron background in the MACE ECAL",
        "geometry": geo,
        "rate": r,
        "capture": cap,
        "gamma": gam,
        "deposition": dep,
        "energy": en,
        "assumptions": {
            "spectrum": "all incident neutrons are thermal and are captured; "
                        "this is an upper bound because each neutron can be "
                        "captured at most once",
            "capture_probability": "P_cap = 1 (real thermal value ~0.9)",
            "hadronic_showers": "not modelled (excluded by decision); they are "
                                "the only channel that could create fake clusters",
            "elastic_inelastic": "not modelled; quenched recoils are at the "
                                 "keV level and (n,n'gamma) gives tens of keV "
                                 "per event, 1-2 orders below the capture channel",
            "geometry": "GP(8,0) ideal 642 cells; 568 instrumented crystals",
            "window": "100 ns, continuous beam (steady state)",
        },
        "interpretation": {
            "n_capture": "expected number of captures per event over the whole "
                         "detector, upper bound",
            "f_dep": "fraction of the 6.86 MeV cascade deposited in the CsI "
                     "array; the remaining %d%% escapes through the outer surface"
                     % round(100 * (1 - dep["f_dep"])),
            "energy_per_event": "n_capture * Q_gamma * f_dep = n_capture * "
                                "energy_per_capture",
        },
    }


def write_json(result: dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")
    return path


if __name__ == "__main__":  # pragma: no cover - convenience only
    print(json.dumps(estimate(), indent=2, ensure_ascii=False))
