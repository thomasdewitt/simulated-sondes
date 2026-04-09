"""Main simulation script for simulated sonde measurements in TWPICE LES."""

import os
import time as timer
import numpy as np
import netCDF4 as nc
from scipy.interpolate import interp1d

from config import (
    N_DROPSONDES, N_RADIOSONDES, N_INSTANTANEOUS,
    DT, SAVE_INTERVAL, DROP_START_ALT, RADIO_START_ALT, RADIO_END_ALT,
    LAUNCH_WINDOW, SAMPLE_VARIABLES, OUTPUT_DIR,
)
from les_data import LESDataManager
from sonde_physics import (
    dropsonde_fall_speed, radiosonde_rise_speed, compute_grid_indices,
)
from regrid import make_grid, regrid_profiles


def simulate_sondes(les, sonde_type, n_sondes, rng):
    """Simulate sonde trajectories through the LES domain.

    Parameters
    ----------
    les : LESDataManager
    sonde_type : "dropsonde" or "radiosonde"
    n_sondes : number of sondes
    rng : numpy random generator

    Returns
    -------
    z_meas, data_meas, meas_count, launch_times, launch_x, launch_y
    """
    is_drop = sonde_type == "dropsonde"

    # Initialize positions
    launch_times = rng.uniform(0, LAUNCH_WINDOW, size=n_sondes)
    x = rng.uniform(0, les.Lx, size=n_sondes)
    y = rng.uniform(0, les.Ly, size=n_sondes)
    launch_x = x.copy()
    launch_y = y.copy()

    if is_drop:
        z = np.full(n_sondes, DROP_START_ALT)
        max_meas = int(DROP_START_ALT / (11.0 * DT * SAVE_INTERVAL)) + 500
    else:
        z = np.full(n_sondes, RADIO_START_ALT)
        max_meas = int((RADIO_END_ALT - RADIO_START_ALT)
                       / (radiosonde_rise_speed() * DT * SAVE_INTERVAL)) + 500

    # State
    launched = np.zeros(n_sondes, dtype=bool)
    active = np.zeros(n_sondes, dtype=bool)

    # Pre-allocate measurement storage
    z_meas = np.full((n_sondes, max_meas), np.nan, dtype=np.float32)
    data_meas = {var: np.full((n_sondes, max_meas), np.nan, dtype=np.float32)
                 for var in SAMPLE_VARIABLES}
    meas_count = np.zeros(n_sondes, dtype=np.int32)

    # Time bounds
    if is_drop:
        max_traverse_time = DROP_START_ALT / 11.0  # slowest fall speed
    else:
        max_traverse_time = (RADIO_END_ALT - RADIO_START_ALT) / radiosonde_rise_speed()
    sim_time_max = LAUNCH_WINDOW + max_traverse_time + 60.0  # small buffer

    sim_time = 0.0
    step = 0
    t_start = timer.time()
    last_print = 0.0

    print(f"\nSimulating {n_sondes} {sonde_type}s...")
    print(f"  Max sim time: {sim_time_max:.0f}s, max measurements/sonde: {max_meas}")

    while sim_time <= sim_time_max:
        # Check if any sondes are still active or unlaunched
        if not np.any(~launched) and not np.any(active):
            break

        # Load nearest LES snapshot
        les.get_nearest(sim_time)

        # Activate newly launched sondes
        newly_launched = (~launched) & (launch_times <= sim_time)
        if np.any(newly_launched):
            launched[newly_launched] = True
            active[newly_launched] = True

        # Get active sonde mask
        active_idx = np.where(active)[0]
        if len(active_idx) == 0:
            sim_time += DT
            step += 1
            continue

        # Compute grid indices for active sondes
        iy, ix, iz = compute_grid_indices(
            x[active_idx], y[active_idx], z[active_idx], les
        )

        # Sample wind for advection
        u_local = les.sample("U", iy, ix, iz)
        v_local = les.sample("V", iy, ix, iz)
        w_local = les.sample("W", iy, ix, iz)

        # Compute sonde vertical velocity
        if is_drop:
            w_sonde = -dropsonde_fall_speed(z[active_idx])
        else:
            w_sonde = radiosonde_rise_speed()

        # Update positions
        x[active_idx] += u_local * DT
        y[active_idx] += v_local * DT
        z[active_idx] += (w_local + w_sonde) * DT

        # Save measurements
        if step % SAVE_INTERVAL == 0:
            # Recompute grid indices at new positions
            iy, ix, iz = compute_grid_indices(
                x[active_idx], y[active_idx], z[active_idx], les
            )
            mc = meas_count[active_idx]
            can_write = mc < max_meas
            write_idx = active_idx[can_write]
            write_mc = mc[can_write]
            if len(write_idx) > 0:
                iy_w, ix_w, iz_w = iy[can_write], ix[can_write], iz[can_write]
                for var in SAMPLE_VARIABLES:
                    if var == "P":
                        # Total pressure = base state (mb→Pa) + perturbation (Pa)
                        # Use same iz level for both to stay consistent
                        pp = les.sample("PP", iy_w, ix_w, iz_w)
                        pbase = les.pres[iz_w] * 100.0
                        vals = pbase + pp
                    else:
                        vals = les.sample(var, iy_w, ix_w, iz_w)
                    data_meas[var][write_idx, write_mc] = vals
                z_meas[write_idx, write_mc] = z[write_idx]
                meas_count[write_idx] = write_mc + 1

        # Terminate sondes
        if is_drop:
            terminate = z[active_idx] <= 0.0
        else:
            terminate = z[active_idx] >= RADIO_END_ALT

        if np.any(terminate):
            active[active_idx[terminate]] = False

        # Progress
        if sim_time - last_print >= 300.0:
            n_active = np.sum(active)
            n_done = np.sum(launched & ~active)
            elapsed = timer.time() - t_start
            print(f"  t={sim_time:.0f}s: {n_active} active, "
                  f"{n_done} done, elapsed={elapsed:.1f}s")
            last_print = sim_time

        sim_time += DT
        step += 1

    elapsed = timer.time() - t_start
    print(f"  Done. {step} steps, {elapsed:.1f}s elapsed.")

    return z_meas, data_meas, meas_count, launch_times, launch_x, launch_y


def simulate_instantaneous(les, n_columns, rng):
    """Extract instantaneous vertical columns from LES, distributed across timesteps.

    Returns regridded profiles on a 10m grid.
    """
    print(f"\nSimulating {n_columns} instantaneous columns...")

    centers, _ = make_grid()
    n_bins = len(centers)
    z_les = les.z
    z_mask = z_les <= 10000.0  # only below 10 km

    # Base-state pressure on the masked LES z-grid (Pa)
    pres_base_col = les.pres[z_mask] * 100.0

    # Distribute columns across all available timesteps
    n_steps = len(les.timesteps)
    per_step = n_columns // n_steps
    remainder = n_columns % n_steps
    counts_per_step = [per_step + (1 if i < remainder else 0)
                       for i in range(n_steps)]

    # Output arrays
    regridded = {var: np.full((n_columns, n_bins), np.nan, dtype=np.float32)
                 for var in SAMPLE_VARIABLES}
    col_times = np.zeros(n_columns, dtype=np.float64)
    col_x = np.zeros(n_columns, dtype=np.float64)
    col_y = np.zeros(n_columns, dtype=np.float64)

    # LES variables to load (PP instead of P)
    les_vars = [v if v != "P" else "PP" for v in SAMPLE_VARIABLES]

    col_offset = 0
    for si, step in enumerate(les.timesteps):
        nc_this = counts_per_step[si]
        if nc_this == 0:
            continue

        print(f"  Timestep {step}: extracting {nc_this} columns")

        data = les.load_timestep(step, les_vars)

        # Random horizontal locations
        ix = rng.integers(0, les.nx, size=nc_this)
        iy = rng.integers(0, les.ny, size=nc_this)

        for j in range(nc_this):
            idx = col_offset + j
            col_times[idx] = les.times_sec[si]
            col_x[idx] = les.x[ix[j]]
            col_y[idx] = les.y[iy[j]]

            for var, les_var in zip(SAMPLE_VARIABLES, les_vars):
                profile = data[les_var][iy[j], ix[j], z_mask].astype(np.float64)
                if var == "P":
                    profile = pres_base_col + profile  # base state + PP
                f_interp = interp1d(
                    z_les[z_mask], profile,
                    kind="linear", bounds_error=False, fill_value=np.nan,
                )
                regridded[var][idx, :] = f_interp(centers).astype(np.float32)

        # Free memory
        del data
        col_offset += nc_this

    return centers, regridded, col_times, col_x, col_y


def save_netcdf(filepath, centers, regridded,
                launch_times, launch_x, launch_y, sonde_type):
    """Save regridded profiles to NetCDF."""
    n_sondes = regridded["U"].shape[0]

    with nc.Dataset(filepath, "w", format="NETCDF4") as ds:
        ds.createDimension("sonde", n_sondes)
        ds.createDimension("altitude", len(centers))

        # Coordinates
        alt_var = ds.createVariable("altitude", "f8", ("altitude",))
        alt_var[:] = centers
        alt_var.units = "m"
        alt_var.long_name = "altitude above surface"

        # Metadata
        lt = ds.createVariable("launch_time", "f8", ("sonde",))
        lt[:] = launch_times
        lt.units = "seconds since simulation start"

        lx = ds.createVariable("launch_x", "f8", ("sonde",))
        lx[:] = launch_x
        lx.units = "m"

        ly = ds.createVariable("launch_y", "f8", ("sonde",))
        ly[:] = launch_y
        ly.units = "m"

        # Data variables
        var_meta = {
            "U": ("m/s", "zonal wind"),
            "V": ("m/s", "meridional wind"),
            "QV": ("g/kg", "water vapor mixing ratio"),
            "P": ("Pa", "total pressure"),
        }

        for var in regridded:
            units, long_name = var_meta.get(var, ("", var))
            v = ds.createVariable(var, "f4", ("sonde", "altitude"),
                                  zlib=True, complevel=4)
            v[:] = regridded[var]
            v.units = units
            v.long_name = long_name

        ds.description = f"Simulated {sonde_type} profiles in TWPICE LES"

    print(f"  Saved {filepath}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Simulate sonde measurements in TWPICE LES")
    parser.add_argument("--subset-xy", type=int, default=None,
                        help="Load only first N×N horizontal points (for testing)")
    parser.add_argument("--n-sondes", type=int, default=None,
                        help="Override number of sondes per type (for testing)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--only", choices=["drop", "radio", "inst"],
                        help="Run only one sonde type")
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    rng = np.random.default_rng(seed=args.seed)

    n_drop = args.n_sondes or N_DROPSONDES
    n_radio = args.n_sondes or N_RADIOSONDES
    n_inst = args.n_sondes or N_INSTANTANEOUS

    les = LESDataManager(subset_xy=args.subset_xy)

    run = args.only  # None means run all

    if run in (None, "drop"):
        z_meas, data_meas, meas_count, lt, lx, ly = simulate_sondes(
            les, "dropsonde", n_drop, rng
        )
        centers, regridded = regrid_profiles(
            z_meas, data_meas, n_drop, z_meas.shape[1], meas_count
        )
        save_netcdf(
            os.path.join(OUTPUT_DIR, "simulated_dropsondes.nc"),
            centers, regridded, lt, lx, ly, "dropsonde"
        )
        del z_meas, data_meas, regridded

    if run in (None, "radio"):
        z_meas, data_meas, meas_count, lt, lx, ly = simulate_sondes(
            les, "radiosonde", n_radio, rng
        )
        centers, regridded = regrid_profiles(
            z_meas, data_meas, n_radio, z_meas.shape[1], meas_count
        )
        save_netcdf(
            os.path.join(OUTPUT_DIR, "simulated_radiosondes.nc"),
            centers, regridded, lt, lx, ly, "radiosonde"
        )
        del z_meas, data_meas, regridded

    # Clear cached LES data before instantaneous (which loads its own)
    les._data.clear()

    if run in (None, "inst"):
        centers, regridded, col_times, col_x, col_y = simulate_instantaneous(
            les, n_inst, rng
        )
        save_netcdf(
            os.path.join(OUTPUT_DIR, "instantaneous_columns.nc"),
            centers, regridded, col_times, col_x, col_y, "instantaneous column"
        )

    print("\nAll done.")


if __name__ == "__main__":
    main()
