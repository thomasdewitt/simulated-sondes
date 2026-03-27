# Simulated Sonde Measurements in TWPICE LES

Simulates dropsonde, radiosonde, and instantaneous vertical profiles through a TWPICE large-eddy simulation to test the effect of non-instantaneous measurements and horizontal drift on vertical structure functions.

## Usage

Full run (requires machine with enough RAM for 2048x2048x255 arrays):
```bash
python simulate.py
```

Test run on limited-memory machine:
```bash
python simulate.py --subset-xy 128 --n-sondes 10
```

## Data

LES data on `/Volumes/BLUE/TWPICE/OUT_3D.{U,V,W,TABS,QV,PP}/`. Currently 3 timesteps; will work with all 23 once downloaded.

## Output

Three NetCDF files in `output/`: `simulated_dropsondes.nc`, `simulated_radiosondes.nc`, `instantaneous_columns.nc`. Each has dims `(sonde, altitude)` with variables U, V, TABS, QV, P on a 10m vertical grid.
