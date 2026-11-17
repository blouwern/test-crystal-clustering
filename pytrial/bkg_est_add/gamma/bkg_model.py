"""Conservative estimate of the gamma background of the MACE ECAL (steps 1-4).

The estimate deliberately avoids a full Geant4 simulation of the laboratory
background.  Instead it combines four ingredients:

1. the ambient gamma dose rate in the hall, ``P_bkg = 30 uSv/h``;
2. the assumed energy-deposition spectrum of one background gamma in a CsI
   crystal -- a normal distribution with mean 5 MeV and variance 2 MeV^2,
   truncated at zero so that only physical (positive) deposits are drawn;
3. the average crystal volume obtained from a first-kind ``GP(8, 0)`` Goldberg
   polyhedron with inner radius 260 mm and crystal length 280 mm;
4. a 100 ns event window, which turns the interaction rate into the expected
   number of background hits per event, ``n_bkg``.

The output of :func:`estimate` is a plain dictionary (JSON friendly).  The
injection script ``add_background.py`` imports the same functions, so the
number used to generate the data can never drift away from the documented
calculation.

Geometry
--------
``GP(8, 0)`` is the dual of the class-I (first-kind) geodesic sphere of
frequency 8 built on the regular icosahedron.  It has

    12 pentagons + 630 hexagons = 642 faces,  1280 vertices, 1920 edges,

so Euler's formula gives 1280 - 1920 + 642 = 2.  The average face area at
radius R is ``S(R)``.  Because the outer polyhedron is an exact radial scaling
of the inner one,

    S_2 = S_1 * (540 / 260)^2,

and the crystal volume is the exact frustum volume

    V_cry = L / 3 * (S_1 + S_2 + sqrt(S_1 * S_2)),   L = 280 mm.

The ideal polyhedron has 642 cells, but the detector has two gaps (upstream and
downstream), so only ``N_CRYSTALS = 568`` crystals are instrumented.  The rate
is computed for their summed volume, ``V_total = N_CRYSTALS * V_cry``.

Physics conversion
------------------
The dose rate is converted to an energy-deposition rate with the standard
photon approximation ``1 Sv ~ 1 Gy`` (radiation weighting factor 1) and
``1 Gy = 1 J/kg``.  The whole crystal array of mass ``M = N_CRYSTALS * m``
therefore absorbs

    P_total = D * M        [W],     D = 30e-6 Gy/h / 3600 s/h

and the interaction rate is ``r_bkg = P_total / E_bkg``.  Note that this treats
the ambient dose equivalent rate as the absorbed dose rate in CsI.  The true
dose conversion depends on the (unknown) background spectrum through the
CsI/air mass-energy-absorption ratio, so this is an O(1) approximation (see
the README for the caveat).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import truncnorm

# --------------------------------------------------------------------------- #
# Inputs of the estimate (steps 1-4 of the task)
# --------------------------------------------------------------------------- #

#: Ambient gamma dose rate in the laboratory hall.
P_BKG_USV_PER_H = 30.0

#: First-kind Goldberg polyhedron GP(FREQUENCY, 0), inner spherical radius.
R_IN_MM = 260.0

#: Single crystal length along the radial direction.
CRYSTAL_LENGTH_MM = 280.0

#: Outer radius of the crystal array.
R_OUT_MM = R_IN_MM + CRYSTAL_LENGTH_MM

#: Subdivision frequency of the GP(8, 0) polyhedron.
FREQUENCY = 8

#: Instrumented crystals in the real detector.  The ideal GP(8, 0) has 642
#: faces, but the detector has two gaps (upstream and downstream), so only 568
#: crystals are instrumented.  Confirmed by the GDML geometry embedded in the
#: raw files and by the length of the ``Edeps`` vector in
#: ``data/processed/*.root``.
N_CRYSTALS = 568

#: Assumed background gamma deposition spectrum: normal(mu, sigma^2), E > 0.
SPECTRUM_MEAN_MEV = 5.0
SPECTRUM_VAR_MEV2 = 2.0

#: Event time window.
WINDOW_NS = 100.0

#: Physical constants.
CSI_DENSITY_G_CM3 = 4.51          # CsI(Tl) / CsI, g cm^-3
MEV_TO_J = 1.602176634e-13        # J per MeV
SV_TO_GY = 1.0                    # photon radiation weighting factor


# --------------------------------------------------------------------------- #
# Step 1: GP(8, 0) Goldberg polyhedron
# --------------------------------------------------------------------------- #

_PHI = (1.0 + 5.0**0.5) / 2.0

_ICO_VERTICES = np.array(
    [
        [-1, _PHI, 0], [1, _PHI, 0], [-1, -_PHI, 0], [1, -_PHI, 0],
        [0, -1, _PHI], [0, 1, _PHI], [0, -1, -_PHI], [0, 1, -_PHI],
        [_PHI, 0, -1], [_PHI, 0, 1], [-_PHI, 0, -1], [-_PHI, 0, 1],
    ],
    dtype=np.float64,
)
_ICO_FACES = np.array(
    [
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
    ],
    dtype=np.int64,
)


def _normalise(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def geodesic_sphere(frequency: int = FREQUENCY) -> tuple[np.ndarray, np.ndarray]:
    """Class-I geodesic sphere of the given frequency.

    Every icosahedron face is subdivided into ``frequency**2`` congruent
    triangles with a barycentric grid, and the grid points are projected onto
    the unit sphere.  Returns ``(vertices, triangles)`` with unit-length
    vertices.
    """
    if frequency < 1:
        raise ValueError("frequency must be >= 1")

    ico = _normalise(_ICO_VERTICES.astype(np.float64))
    coords: list[np.ndarray] = []
    index: dict[tuple[float, float, float], int] = {}

    def vertex(p: np.ndarray) -> int:
        p = p / np.linalg.norm(p)
        key = (round(float(p[0]), 12), round(float(p[1]), 12), round(float(p[2]), 12))
        idx = index.get(key)
        if idx is None:
            idx = len(coords)
            index[key] = idx
            coords.append(p)
        return idx

    f = int(frequency)
    triangles: list[tuple[int, int, int]] = []
    for a, b, c in _ICO_FACES:
        A, B, C = ico[a], ico[b], ico[c]
        grid: dict[tuple[int, int, int], int] = {}
        for i in range(f + 1):
            for j in range(f + 1 - i):
                k = f - i - j
                grid[(i, j, k)] = vertex((i * A + j * B + k * C) / f)
        for i in range(f):
            for j in range(f - i):
                k = f - i - j
                triangles.append((grid[(i, j, k)], grid[(i + 1, j, k - 1)], grid[(i, j + 1, k - 1)]))
                if j <= f - i - 2:
                    triangles.append(
                        (grid[(i + 1, j, k - 1)], grid[(i + 1, j + 1, k - 2)], grid[(i, j + 1, k - 1)])
                    )

    return np.asarray(coords), np.asarray(triangles, dtype=np.int64)


def dual_faces(vertices: np.ndarray, triangles: np.ndarray) -> list[np.ndarray]:
    """Dual of a triangular mesh.

    Each mesh vertex becomes one polygon whose corners are the (normalised)
    centroids of the incident triangles, ordered by azimuth in the tangent
    plane of the vertex.  For a geodesic sphere this is exactly the Goldberg
    polyhedron: pentagons at the 12 icosahedron vertices, hexagons elsewhere.
    """
    incident: list[list[np.ndarray]] = [[] for _ in range(len(vertices))]
    for tri in triangles:
        for v in tri:
            incident[int(v)].append(tri)

    faces: list[np.ndarray] = []
    for v, tris in enumerate(incident):
        normal = vertices[v]
        helper = np.array([0.0, 0.0, 1.0]) if abs(normal[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        e1 = np.cross(normal, helper)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(normal, e1)

        corners = []
        for tri in tris:
            centre = vertices[tri].mean(axis=0)
            centre /= np.linalg.norm(centre)
            corners.append((math.atan2(float(centre @ e2), float(centre @ e1)), centre))
        corners.sort(key=lambda item: item[0])
        faces.append(np.asarray([c for _, c in corners], dtype=np.float64))

    return faces


def polygon_area(points: np.ndarray) -> float:
    """Area of a (possibly slightly non-planar) polygon: Newell's method."""
    cross = np.cross(points, np.roll(points, -1, axis=0))
    return 0.5 * float(np.linalg.norm(cross.sum(axis=0)))


def goldberg_faces(radius_mm: float, frequency: int = FREQUENCY) -> list[np.ndarray]:
    """Faces of GP(frequency, 0) with all vertices on a sphere of ``radius_mm``."""
    vertices, triangles = geodesic_sphere(frequency)
    return [face * radius_mm for face in dual_faces(vertices, triangles)]


def geometry(radius_in_mm: float = R_IN_MM, length_mm: float = CRYSTAL_LENGTH_MM) -> dict:
    """Build GP(8, 0) and return S_1, S_2 and the average crystal volume."""
    faces_in = goldberg_faces(radius_in_mm)
    faces_out = goldberg_faces(radius_in_mm + length_mm)

    areas_in = np.array([polygon_area(f) for f in faces_in])
    areas_out = np.array([polygon_area(f) for f in faces_out])

    sizes = np.array([len(f) for f in faces_in])
    _, triangles = geodesic_sphere(FREQUENCY)
    n_edges = len({tuple(sorted((int(tri[0]), int(tri[1])))) for tri in triangles}
                  | {tuple(sorted((int(tri[1]), int(tri[2])))) for tri in triangles}
                  | {tuple(sorted((int(tri[2]), int(tri[0])))) for tri in triangles})
    n_dual_vertices = len(triangles)   # one dual vertex per geodesic triangle

    s1 = float(areas_in.mean())
    s2 = float(areas_out.mean())

    # Exact volume of the frustum between the two end faces.  (The trapezoidal
    # (S1+S2)/2*L form is kept only as a reference: it overestimates a tapered
    # crystal by ~8%.)
    volume_mm3 = length_mm / 3.0 * (s1 + s2 + math.sqrt(s1 * s2))
    volume_trapezoid_mm3 = (s1 + s2) / 2.0 * length_mm

    euler = n_dual_vertices - n_edges + len(faces_in)
    if euler != 2:
        raise AssertionError(f"polyhedron is not closed: Euler characteristic {euler}")
    if not np.isclose(s2, s1 * ((radius_in_mm + length_mm) / radius_in_mm) ** 2, rtol=1e-12):
        raise AssertionError("outer face area does not follow the radial scaling law")

    volume_cm3 = volume_mm3 / 1000.0
    mass_kg = volume_cm3 * CSI_DENSITY_G_CM3 / 1000.0
    return {
        "polyhedron": "GP(8,0) first-kind Goldberg polyhedron (dual of class-I geodesic sphere)",
        "frequency": FREQUENCY,
        "inner_radius_mm": radius_in_mm,
        "outer_radius_mm": radius_in_mm + length_mm,
        "crystal_length_mm": length_mm,
        "n_faces": int(len(faces_in)),
        "n_pentagons": int(np.count_nonzero(sizes == 5)),
        "n_hexagons": int(np.count_nonzero(sizes == 6)),
        "n_vertices": int(n_dual_vertices),
        "n_edges": int(n_edges),
        "euler_characteristic": int(euler),
        "total_surface_mm2": float(areas_in.sum()),
        "sphere_surface_mm2": float(4.0 * math.pi * radius_in_mm**2),
        "s1_mm2": s1,
        "s2_mm2": s2,
        "s1_hexagon_mm2": float(areas_in[sizes == 6].mean()),
        "s1_pentagon_mm2": float(areas_in[sizes == 5].mean()),
        "volume_formula": "exact frustum: L/3 * (S1 + S2 + sqrt(S1*S2))",
        "volume_mm3": float(volume_mm3),
        "volume_cm3": float(volume_cm3),
        "volume_trapezoid_mm3": float(volume_trapezoid_mm3),
        "mass_kg": float(mass_kg),
        "n_crystals": N_CRYSTALS,
        "volume_total_cm3": float(volume_cm3 * N_CRYSTALS),
        "mass_total_kg": float(mass_kg * N_CRYSTALS),
        "density_g_cm3": CSI_DENSITY_G_CM3,
    }


# --------------------------------------------------------------------------- #
# Step 2: energy-deposition spectrum
# --------------------------------------------------------------------------- #

class DepositionSpectrum:
    """Normal spectrum truncated at zero, describing one background gamma hit.

    The nominal distribution is ``N(mean, variance)``.  Deposits below zero are
    unphysical, so the distribution is truncated (and renormalised) at E = 0.
    The truncation changes the mean by 0.02 % for the values used here.
    """

    def __init__(self, mean_mev: float = SPECTRUM_MEAN_MEV, var_mev2: float = SPECTRUM_VAR_MEV2):
        if mean_mev <= 0 or var_mev2 <= 0:
            raise ValueError("spectrum mean and variance must be positive")
        self.mean_mev = float(mean_mev)
        self.var_mev2 = float(var_mev2)
        self.sigma_mev = math.sqrt(self.var_mev2)
        self._dist = truncnorm(
            a=(0.0 - self.mean_mev) / self.sigma_mev,
            b=np.inf,
            loc=self.mean_mev,
            scale=self.sigma_mev,
        )

    @property
    def mean_deposition_mev(self) -> float:
        """Mean energy deposited per background gamma (truncated at zero)."""
        return float(self._dist.mean())

    def rvs(self, size: int, rng: np.random.Generator) -> np.ndarray:
        return self._dist.rvs(size=int(size), random_state=rng)

    def pdf(self, energy_mev: np.ndarray) -> np.ndarray:
        return self._dist.pdf(energy_mev)

    def summary(self) -> dict:
        return {
            "distribution": "normal truncated at E>0",
            "nominal_mean_mev": self.mean_mev,
            "variance_mev2": self.var_mev2,
            "sigma_mev": self.sigma_mev,
            "mean_deposition_mev": self.mean_deposition_mev,
        }


# --------------------------------------------------------------------------- #
# Steps 3-4: rate and expected hits per event
# --------------------------------------------------------------------------- #

def estimate(
    dose_usv_per_h: float = P_BKG_USV_PER_H,
    spectrum_mean_mev: float = SPECTRUM_MEAN_MEV,
    spectrum_var_mev2: float = SPECTRUM_VAR_MEV2,
    window_ns: float = WINDOW_NS,
    sv_to_gy: float = SV_TO_GY,
    density_g_cm3: float = CSI_DENSITY_G_CM3,
) -> dict:
    """Run steps 1-4 and return every intermediate number."""
    geo = geometry()
    geo["density_g_cm3"] = float(density_g_cm3)
    geo["mass_kg"] = float(geo["volume_cm3"] * density_g_cm3 / 1000.0)
    geo["mass_total_kg"] = float(geo["volume_total_cm3"] * density_g_cm3 / 1000.0)

    spectrum = DepositionSpectrum(spectrum_mean_mev, spectrum_var_mev2)
    e_bkg = spectrum.mean_deposition_mev

    dose_gy_per_h = dose_usv_per_h * 1e-6 * sv_to_gy
    dose_gy_per_s = dose_gy_per_h / 3600.0

    # The dose rate is absorbed over the whole instrumented crystal volume
    # (N_CRYSTALS crystals), not over a single crystal.
    power_per_crystal_w = dose_gy_per_s * geo["mass_kg"]
    power_total_w = dose_gy_per_s * geo["mass_total_kg"]
    energy_per_hit_j = e_bkg * MEV_TO_J

    r_bkg = power_total_w / energy_per_hit_j
    window_s = window_ns * 1e-9
    n_bkg = r_bkg * window_s

    # Poisson probabilities for the per-event hit multiplicity.
    poisson = {
        str(k): float(math.exp(-n_bkg) * n_bkg**k / math.factorial(k))
        for k in range(4)
    }

    return {
        "title": "Estimated gamma background in the MACE ECAL crystal array",
        "geometry": geo,
        "spectrum": spectrum.summary(),
        "dose": {
            "p_bkg_usv_per_h": float(dose_usv_per_h),
            "sv_to_gy": float(sv_to_gy),
            "dose_rate_gy_per_h": float(dose_gy_per_h),
            "dose_rate_gy_per_s": float(dose_gy_per_s),
            "power_per_crystal_w": float(power_per_crystal_w),
            "power_total_w": float(power_total_w),
            "energy_per_hit_j": float(energy_per_hit_j),
        },
        "rates": {
            "r_bkg_per_s": float(r_bkg),
            "event_window_ns": float(window_ns),
            "n_bkg_per_event": float(n_bkg),
            "poisson_probability": poisson,
            "expected_hits_per_1e6_events": float(n_bkg * 1e6),
        },
        "interpretation": {
            "volume": "the rate uses the summed volume of all N_CRYSTALS = 568 "
                      "instrumented crystals (the ideal GP(8,0) has 642 faces; the "
                      "detector has two gaps, upstream and downstream, so 74 cells "
                      "are absent), i.e. V_total = N_CRYSTALS * V_cry.",
            "n_bkg": "n_bkg = r_bkg * 100 ns is the expected number of background "
                     "hits per event over the whole detector; step 5 draws that "
                     "number and assigns each hit to a random crystal.",
        },
        "assumptions": {
            "dose": "1 Sv = 1 Gy for photons; the ambient dose rate is applied as "
                    "the absorbed dose rate in CsI (no CsI/air mass-energy-absorption correction)",
            "spectrum": "normal distribution, mean 5 MeV, variance 2 MeV^2, truncated at 0",
            "geometry": "crystal = frustum between the mean inner and outer face areas "
                        "of GP(8,0), V = L/3 * (S1 + S2 + sqrt(S1*S2)); total volume "
                        "= 568 * V",
            "statistics": "background hits are Poisson distributed in a 100 ns window",
        },
    }


def write_json(result: dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


if __name__ == "__main__":  # pragma: no cover - convenience only
    print(json.dumps(estimate(), indent=2, ensure_ascii=False))
