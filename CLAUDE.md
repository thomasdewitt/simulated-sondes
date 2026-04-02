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

LES data on `/media/thomasdewitt/Expansion/LES-data/SAM-TWPICE/OUT_3D.{U,V,W,TABS,QV,PP}/`. 23 timesteps available.

## Output

Three NetCDF files in `output/`: `simulated_dropsondes.nc`, `simulated_radiosondes.nc`, `instantaneous_columns.nc`. Each has dims `(sonde, altitude)` with variables U, V, QV, P on a 10m vertical grid.
