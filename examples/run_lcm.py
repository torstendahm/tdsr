import sys
sys.path.insert(0, "../")
from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional, RSM
from tdsm.plotting import plot
from tdsm.loading import StepLoading, FourPointLoading, BackgroundLoading
import matplotlib.pyplot as plt

current_dir = Path(__file__).parent.parent
config_file = current_dir / "config.toml"

print("Initialise a tdsm lcm run")
lcm = LCM(config=Config.open(config_file))

# ---- Define parameters ------
strend = 7.E-5
sstep  = 1.E0

chi0   = 1.E4
depthS = -0.5E0

deltaS = -depthS/100.
sigma_max = 6000.*deltaS
precision = 18

print("Calculate rate for step loading starting from equilibrium distribution")
loading = StepLoading(_config=lcm.config, sstep=sstep, strend=strend)
config, t, chiz, cf, ratez, neqz = lcm(loading=loading, chi0=chi0, deltaS=deltaS, sigma_max=sigma_max, precision=18)
plot(config, t, cf, ratez, neqz)
