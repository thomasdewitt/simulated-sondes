"""Bin-average regridding to a uniform vertical grid."""

import numpy as np
from config import REGRID_DZ, REGRID_ZMIN, REGRID_ZMAX


def make_grid():
    """Create the uniform regridding altitude grid. Returns bin centers."""
    n_bins = int((REGRID_ZMAX - REGRID_ZMIN) / REGRID_DZ)
    edges = REGRID_ZMIN + np.arange(n_bins + 1) * REGRID_DZ
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, edges


def regrid_profiles(z_meas, data_meas, n_sondes, max_meas, meas_count):
    """Bin-average measurements onto a uniform vertical grid.

    Parameters
    ----------
    z_meas : (n_sondes, max_meas) array of measurement altitudes
    data_meas : dict of (n_sondes, max_meas) arrays for each variable
    n_sondes : number of sondes
    max_meas : max measurements per sonde
    meas_count : (n_sondes,) array of actual measurement counts

    Returns
    -------
    centers : 1D altitude grid
    regridded : dict of (n_sondes, n_bins) arrays, NaN where no data
    """
    centers, edges = make_grid()
    n_bins = len(centers)

    regridded = {var: np.full((n_sondes, n_bins), np.nan, dtype=np.float32)
                 for var in data_meas}

    for i in range(n_sondes):
        n = meas_count[i]
        if n == 0:
            continue
        z_i = z_meas[i, :n]
        bin_idx = np.floor((z_i - REGRID_ZMIN) / REGRID_DZ).astype(int)
        valid = (bin_idx >= 0) & (bin_idx < n_bins)
        bi = bin_idx[valid]

        for var, arr in data_meas.items():
            vals = arr[i, :n][valid].astype(np.float64)
            sums = np.bincount(bi, weights=vals, minlength=n_bins)
            counts = np.bincount(bi, minlength=n_bins)
            occupied = counts > 0
            regridded[var][i, occupied] = (sums[occupied] / counts[occupied]).astype(np.float32)

    return centers, regridded
