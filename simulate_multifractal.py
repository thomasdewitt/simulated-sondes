"""Synthetic multifractal sonde profiles.

Generate FIF, coarsen, fractionally integrate (uniform or broken), smooth to
mimic sonde inertia, trim to break periodicity, and save one NetCDF per case.
"""

import argparse
import os
import time as timer

import numpy as np
import xarray as xr
from scipy.interpolate import make_lsq_spline

import scaleinvariance
from scaleinvariance.simulation.FIF import extremal_levy


# Backend: float64 on CUDA via torch
scaleinvariance.set_backend("torch")
scaleinvariance.set_device("cuda")
scaleinvariance.set_numerical_precision("float64")


# Geometry
FIF_FINE_SIZE = 4_000_000       # 40,000 * 100, periodic
COARSEN = 100                    # → 40,000 flux samples at 1 m
N_FLUX = FIF_FINE_SIZE // COARSEN          # 40,000 (= 40 km at 1 m)
TOP_DELETE = 20_000              # delete top 20 km after smoothing → keep 20 km
N_KEEP = N_FLUX - TOP_DELETE     # 20,000

# Smoothing: cubic B-spline low-pass (Ooyama / NCAR ASPEN style), 50 m cutoff
# wavelength. Knot spacing = wavelength/2 so the spline basis cannot represent
# features shorter than the cutoff. See https://github.com/NCAR/bspline.
SMOOTH_WAVELENGTH_M = 50.0

# Post-processing
SCALE = 10.0

# FIF parameters
ALPHA = 1.8
C1 = 0.02
H_FIF = 0.0

# Three fractional-integral regimes
CASES = [
    {"name": "uniform_H06",          "kind": "uniform", "H": 0.6, "smooth": True},
    {"name": "uniform_H06_nosmooth", "kind": "uniform", "H": 0.6, "smooth": False},
    {"name": "broken_10m",           "kind": "broken",  "H_small": 0.3, "H_large": 1.0, "transition_m": 10.0,   "smooth": True},
    {"name": "broken_1km",           "kind": "broken",  "H_small": 0.3, "H_large": 1.0, "transition_m": 1000.0, "smooth": True},
]

DEFAULT_N_SONDES = 1000
DEFAULT_SEED = 42
DEFAULT_OUTPUT_DIR = "output"


def bspline_smooth_periodic(y, cutoff_wavelength):
    """Cubic B-spline low-pass (Ooyama/ASPEN style), treating y as periodic.

    Knots placed every cutoff_wavelength/2 units, so the spline basis cannot
    represent features shorter than the cutoff. Periodicity is handled by
    padding with wrapped data and trimming after evaluation.
    """
    n = len(y)
    k = 3
    dx = cutoff_wavelength / 2.0
    pad = int(round(5 * cutoff_wavelength))  # 5 wavelengths of padding
    y_pad = np.concatenate([y[-pad:], y, y[:pad]])
    x = np.arange(n + 2 * pad, dtype=np.float64)
    n_interior = max(1, int((x[-1] - x[0]) / dx) - 1)
    t_int = np.linspace(x[0], x[-1], n_interior + 2)[1:-1]
    t = np.concatenate((np.full(k + 1, x[0]), t_int, np.full(k + 1, x[-1])))
    spl = make_lsq_spline(x, y_pad, t, k=k)
    return spl(np.arange(pad, pad + n, dtype=np.float64))


def simulate_one(case, seed):
    """Run the full per-sonde pipeline and return a 1D profile of length N_KEEP."""
    noise = extremal_levy(alpha=ALPHA, size=FIF_FINE_SIZE, seed=seed)

    flux_fine = scaleinvariance.FIF_1D(
        FIF_FINE_SIZE,
        alpha=ALPHA, C1=C1, H=H_FIF,
        levy_noise=noise,
        periodic=True, causal=False,
    )

    flux = flux_fine.reshape(N_FLUX, COARSEN).mean(axis=1)

    if case["kind"] == "uniform":
        profile = scaleinvariance.fractional_integral_spectral(flux, H=case["H"])
    else:
        profile = scaleinvariance.broken_fractional_integral_spectral(
            flux,
            H_small_scale=case["H_small"],
            H_large_scale=case["H_large"],
            transition_wavelength=case["transition_m"],
        )

    profile = scaleinvariance.to_numpy(profile)
    if case["smooth"]:
        profile = bspline_smooth_periodic(profile, cutoff_wavelength=SMOOTH_WAVELENGTH_M)

    profile = profile[:N_KEEP]
    profile = profile / profile.mean() * SCALE

    return profile.astype(np.float32)


def run_case(case, n_sondes, base_seed, output_dir):
    t0 = timer.time()
    print(f"\n[{case['name']}] simulating {n_sondes} profiles...")

    u = np.empty((n_sondes, N_KEEP), dtype=np.float32)
    for i in range(n_sondes):
        u[i] = simulate_one(case, seed=base_seed + i)
        if (i + 1) % 50 == 0 or i + 1 == n_sondes:
            dt = timer.time() - t0
            print(f"  {i+1}/{n_sondes}  ({dt:.1f}s, {dt/(i+1):.2f}s/sonde)")

    altitude = np.arange(N_KEEP, dtype=np.float64) + 0.5  # cell centers, m

    attrs = {
        "source": "simulate_multifractal.py",
        "case_name": case["name"],
        "case_kind": case["kind"],
        "alpha": ALPHA,
        "C1": C1,
        "H_fif": H_FIF,
        "fif_fine_size": FIF_FINE_SIZE,
        "coarsen": COARSEN,
        "smooth_kind": "cubic_bspline_lowpass_ooyama" if case["smooth"] else "none",
        "smooth_cutoff_wavelength_m": SMOOTH_WAVELENGTH_M if case["smooth"] else 0.0,
        "top_delete_m": TOP_DELETE,
        "scale": SCALE,
        "seed_base": base_seed,
    }
    if case["kind"] == "uniform":
        attrs["case_H"] = case["H"]
    else:
        attrs["case_H_small"] = case["H_small"]
        attrs["case_H_large"] = case["H_large"]
        attrs["case_transition_m"] = case["transition_m"]

    ds = xr.Dataset(
        data_vars={
            "u": (("sonde", "altitude"), u, {"units": "m/s", "long_name": "multifractal wind-like field"}),
        },
        coords={
            "altitude": ("altitude", altitude, {"units": "m", "long_name": "altitude above surface"}),
            "sonde": ("sonde", np.arange(n_sondes, dtype=np.int32)),
        },
        attrs=attrs,
    )

    filepath = os.path.join(output_dir, f"simulated_multifractal_{case['name']}.nc")
    ds.to_netcdf(filepath, encoding={"u": {"zlib": True, "complevel": 4}})
    print(f"  Saved {filepath}  ({timer.time() - t0:.1f}s total)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-sondes", type=int, default=DEFAULT_N_SONDES)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--only", choices=[c["name"] for c in CASES],
                        help="Run only one case")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    for ci, case in enumerate(CASES):
        if args.only and case["name"] != args.only:
            continue
        # Distinct seed range per case so files don't share noise.
        case_seed_base = args.seed + ci * 1_000_000
        run_case(case, args.n_sondes, case_seed_base, args.output_dir)

    print("\nAll done.")


if __name__ == "__main__":
    main()
