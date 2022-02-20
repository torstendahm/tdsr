import sys
sys.path.insert(0, "../")
from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional, RSM
from tdsm.plotting import plot
from tdsm.loading import StepLoading, FourPointLoading, BackgroundLoading
import matplotlib.pyplot as plt

current_dir = Path(__file__).parent.parent
config_file = current_dir / "config.toml"

print("Initialise a linear coulomb failure model using simple of loop")
trad = Traditional(config=Config.open(config_file))

# ---- define paramerter -----R

strend = 7.E-5
sstep  = 1.E0

chi0   = 1.E4
depthS = -0.5E0

deltaS = -depthS/60.
sigma_max = 3000.*deltaS

print("Calculate rate for step loading")
loading = StepLoading(_config=trad.config, sstep=sstep, strend=strend)
config, t, chiz, cf, ratez, neqz = trad(loading=loading, chi0=chi0)
plot(config, t, cf, ratez, neqz)
