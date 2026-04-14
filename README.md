# Simulated Sonde Measurements in TWPICE LES

Simulates dropsonde, radiosonde, and instantaneous vertical profiles through a TWPICE large-eddy simulation to test the effect of non-instantaneous measurements and horizontal drift on vertical structure functions.

## Method

We create simulated dropsonde and radiosonde profiles in a large-eddy simulation of tropical convection (SAM TWPICE). The simulation has a horizontal grid spacing of 100 m and a maximum vertical grid spacing of 100 m, which is similar to the effective resolution of real sondes due to their vertical inertia. After spin-up, 5-minute snapshots spanning a total of two hours are used.

Sondes are initialized at random times and horizontal locations within the first hour of the period considered, and at initial altitudes of 10 km (dropsondes) or the surface (radiosondes). Profiles are constructed by numerically integrating along the sonde path. Radiosondes ascend at 5 m/s, while dropsondes fall at a rate that linearly decreases from 20 m/s at 10 km to 11 m/s at the surface, approximating the observed fall rates during the ACTIVATE field campaign. As the sondes pass through the LES volume, they are advected by the local 3D wind field, defined as the nearest wind vector in time and space at each measurement location. The timestep for the numerical integration is 0.1 s, and measurements are saved every other timestep. The profiles are then regridded to a common 10 m height grid by bin averaging. A total of 1000 dropsonde and 1000 radiosonde profiles are simulated.

As a control, 1000 instantaneous vertical columns are also extracted at random horizontal locations, distributed evenly across all available timesteps. These columns are interpolated onto the same 10 m grid and serve as a reference unaffected by sonde drift or temporal inconsistency.

As a synthetic control with known scaling properties, 1000 multifractal profiles are also generated for each of four fractional-integral regimes: uniform `H=0.6`, uniform `H=0.6` without the inertia-smoothing step (for reference), a broken regime with `H=0.3` below 10 m and `H=1` above, and a broken regime with `H=0.3` below 1 km and `H=1` above. Each profile starts as an FIF flux (`alpha=1.8`, `C1=0.05`, `H=0`) of length 4,000,000 that is coarsened by averaging into 40,000 samples at 1 m spacing, fractionally integrated, smoothed with a cubic B-spline low-pass filter (Ooyama / NCAR ASPEN style, 50 m cutoff wavelength) to mimic sonde inertia and processing, truncated to the bottom 20 km (to break periodicity), normalised to unit mean, and multiplied by 10.

## Usage

Set up the environment with [uv](https://docs.astral.sh/uv/):

```bash
uv sync
```

Full LES-sonde run (requires machine with enough RAM for 2048x2048x255 arrays):

```bash
uv run python simulate.py
```

Test run on limited-memory machine:

```bash
uv run python simulate.py --subset-xy 128 --n-sondes 10
```

Multifractal synthetic sondes (GPU-accelerated via `scaleinvariance` torch/cuda backend):

```bash
uv run python simulate_multifractal.py
```

## Data

LES data must be locally available.

## Output

Seven NetCDF files in `output/`: `simulated_dropsondes.nc`, `simulated_radiosondes.nc`, `instantaneous_columns.nc` (each dims `(sonde, altitude)` with U, V, QV, P on a 10 m vertical grid), plus `simulated_multifractal_uniform_H06.nc`, `simulated_multifractal_uniform_H06_nosmooth.nc`, `simulated_multifractal_broken_10m.nc`, `simulated_multifractal_broken_1km.nc` (each dims `(sonde, altitude)` with variable `u` on a 1 m vertical grid from 0–20 km).
