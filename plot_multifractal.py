"""Quick verification plot: 2 sondes per case, thin black line, high DPI."""

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

CASES = ["uniform_H06", "uniform_H06_nosmooth", "broken_10m", "broken_1km"]
SONDE_INDICES = [0, 1]

fig, axes = plt.subplots(
    len(CASES), len(SONDE_INDICES),
    figsize=(16, 24), sharey=True,
)

for i, case in enumerate(CASES):
    path = f"output/simulated_multifractal_{case}.nc"
    with nc.Dataset(path) as ds:
        alt = ds.variables["altitude"][:]
        u = ds.variables["u"][SONDE_INDICES, :]
    for j, idx in enumerate(SONDE_INDICES):
        ax = axes[i, j]
        ax.plot(u[j], alt / 1000.0, color="black", linewidth=0.3)
        ax.set_title(f"{case}  |  sonde {idx}", fontsize=10)
        if j == 0:
            ax.set_ylabel("altitude (km)")
        if i == len(CASES) - 1:
            ax.set_xlabel("u (m/s)")
        ax.grid(True, linewidth=0.2, alpha=0.4)

fig.tight_layout()
out = "output/multifractal_verification.png"
fig.savefig(out, dpi=400, bbox_inches="tight")
print(f"Saved {out}")
