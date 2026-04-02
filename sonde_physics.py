"""Sonde physics: fall/rise speeds and grid index computation."""

import numpy as np


def dropsonde_fall_speed(z):
    """Dropsonde fall speed (m/s, positive downward).

    Linear from 20 m/s at 10 km to 11 m/s at surface.
    """
    speed = 11.0 + 9.0 * (np.clip(z, 0, 10000) / 10000.0)
    return speed


def radiosonde_rise_speed():
    """Radiosonde ascent speed (m/s, positive upward). Constant 5 m/s."""
    return 5.0


def compute_grid_indices(x, y, z, les):
    """Convert physical positions to nearest LES grid indices.

    Horizontal positions are wrapped for periodicity.
    Vertical index is clipped to valid range.
    """
    # Wrap horizontal positions into [0, Lx)
    x_wrapped = x % les.Lx
    y_wrapped = y % les.Ly

    # Nearest grid index
    ix = np.round((x_wrapped - les.x0) / les.dx).astype(np.int32) % les.nx
    iy = np.round((y_wrapped - les.y0) / les.dy).astype(np.int32) % les.ny

    # Vertical: nearest index, clipped
    iz = np.searchsorted(les.z, z, side="right") - 1
    # searchsorted - 1 gives the index of the z level just below; round to nearest
    iz_below = np.clip(iz, 0, les.nz - 1)
    iz_above = np.clip(iz + 1, 0, les.nz - 1)
    d_below = np.abs(z - les.z[iz_below])
    d_above = np.abs(z - les.z[iz_above])
    iz = np.where(d_below <= d_above, iz_below, iz_above)

    return iy, ix, iz
