"""Lazy-loading manager for TWPICE LES data."""

import glob
import re
import numpy as np
import netCDF4 as nc

from config import DATA_DIR, LOAD_VARIABLES


class LESDataManager:
    """Manages LES data with lazy loading of the nearest timestep."""

    def __init__(self, data_dir=DATA_DIR, subset_xy=None):
        """
        Parameters
        ----------
        data_dir : str
            Path to TWPICE data directory.
        subset_xy : int or None
            If given, only load the first subset_xy points in x and y.
            Useful for testing on machines with limited memory.
        """
        self.data_dir = data_dir
        self._loaded_step = None
        self._data = {}
        self._subset_xy = subset_xy

        # Read coordinates from any file (they're the same across all)
        sample_file = sorted(
            glob.glob(f"{data_dir}/OUT_3D.U/TWPICE_LPT_3D_U_*.nc")
        )[0]
        with nc.Dataset(sample_file) as ds:
            x_full = ds.variables["x"][:].astype(np.float64)
            y_full = ds.variables["y"][:].astype(np.float64)
            self.z = ds.variables["z"][:].astype(np.float64)
            self.pres = ds.variables["pres"][:].astype(np.float64)  # mb, f(z)

        if subset_xy is not None:
            self.x = x_full[:subset_xy]
            self.y = y_full[:subset_xy]
        else:
            self.x = x_full
            self.y = y_full

        self.nx = len(self.x)
        self.ny = len(self.y)
        self.nz = len(self.z)
        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.x0 = self.x[0]
        self.y0 = self.y[0]
        # Domain length (periodic): from first cell center to last + one spacing
        self.Lx = self.nx * self.dx
        self.Ly = self.ny * self.dy

        # Build time mapping: scan all U files for available timesteps
        u_files = sorted(glob.glob(f"{data_dir}/OUT_3D.U/TWPICE_LPT_3D_U_*.nc"))
        self._timesteps = []  # list of (step_number, time_in_seconds)
        times_days = []
        for f in u_files:
            step_num = int(re.search(r"_(\d{10})\.nc$", f).group(1))
            with nc.Dataset(f) as ds:
                t_days = float(ds.variables["time"][:].item())
            self._timesteps.append(step_num)
            times_days.append(t_days)

        # Convert to seconds relative to first timestep
        t0 = times_days[0]
        self._times_sec = np.array([(t - t0) * 86400.0 for t in times_days])
        self._timesteps = np.array(self._timesteps)

        print(f"LES data: {len(self._timesteps)} timesteps, "
              f"time range 0 to {self._times_sec[-1]:.0f} s")
        print(f"Grid: {self.nx}x{self.ny}x{self.nz}, "
              f"dx={self.dx:.0f}m, Lx={self.Lx:.0f}m")

    def get_nearest(self, sim_time):
        """Load the LES timestep nearest to sim_time (in seconds)."""
        idx = np.argmin(np.abs(self._times_sec - sim_time))
        step = self._timesteps[idx]

        if step == self._loaded_step:
            return  # already loaded

        # Unload previous
        self._data.clear()

        print(f"  Loading LES timestep {step} "
              f"(t={self._times_sec[idx]:.0f}s, sim_t={sim_time:.0f}s)")

        for var in LOAD_VARIABLES:
            fpath = (f"{self.data_dir}/OUT_3D.{var}/"
                     f"TWPICE_LPT_3D_{var}_{step:010d}.nc")
            with nc.Dataset(fpath) as ds:
                if self._subset_xy is not None:
                    n = self._subset_xy
                    self._data[var] = ds.variables[var][0, :n, :n, :].astype(np.float32)
                else:
                    self._data[var] = ds.variables[var][0].astype(np.float32)

        self._loaded_step = step

    def load_timestep(self, step, variables):
        """Load specific variables for a specific timestep. Returns dict."""
        data = {}
        for var in variables:
            fpath = (f"{self.data_dir}/OUT_3D.{var}/"
                     f"TWPICE_LPT_3D_{var}_{step:010d}.nc")
            with nc.Dataset(fpath) as ds:
                if self._subset_xy is not None:
                    n = self._subset_xy
                    data[var] = ds.variables[var][0, :n, :n, :].astype(np.float32)
                else:
                    data[var] = ds.variables[var][0].astype(np.float32)
        return data

    def sample(self, var, iy, ix, iz):
        """Sample a loaded variable at integer grid indices."""
        return self._data[var][iy, ix, iz]

    @property
    def timesteps(self):
        return self._timesteps

    @property
    def times_sec(self):
        return self._times_sec
