"""Configuration for simulated sonde measurements in TWPICE LES."""

DATA_DIR = "/media/thomasdewitt/Expansion/LES-data/SAM-TWPICE"

# Variables saved per sonde (P = PP + base state, computed at sample time)
SAMPLE_VARIABLES = ["U", "V", "QV", "P"]

# Variables needed for trajectory advection
ADVECT_VARIABLES = ["U", "V", "W"]

# All variables that must be loaded from LES
LOAD_VARIABLES = list(set(["U", "V", "W", "QV", "PP"]))

N_DROPSONDES = 1000
N_RADIOSONDES = 1000
N_INSTANTANEOUS = 1000

DT = 0.1            # integration timestep (s)
SAVE_INTERVAL = 2   # save every 2nd step (0.2s)

DROP_START_ALT = 10000.0    # m
RADIO_START_ALT = 25.0      # m (lowest LES level)
RADIO_END_ALT = 10000.0     # m

REGRID_DZ = 10.0            # m
REGRID_ZMIN = 0.0           # m
REGRID_ZMAX = 10000.0       # m

LAUNCH_WINDOW = 3600.0      # sondes launched within first hour (s)

OUTPUT_DIR = "output"
