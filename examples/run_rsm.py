import sys
sys.path.insert(0, "../")
from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional, RSM
from tdsm.plotting import plot
from tdsm.loading import StepLoading, FourPointLoading, BackgroundLoading, CyclicLoading
import matplotlib.pyplot as plt

current_dir = Path(__file__).parent.parent
config_file = current_dir / "config.toml"

print("Initialise a run using the rsm model")
rsm = RSM(config=Config.open(config_file))

# ---- define parameter ----
sstep  = 1.E0
tstep = 10_000

ampsin = 2.0
Tsin = 43_200

chi0   = 1.E4

strend = 7.E-5
depthS = -0.5E0

deltaS = -depthS/60.
sigma_max = 3000.*deltaS

print("Calculate rate for step excitation")
loading = StepLoading(_config=rsm.config, sstep=sstep, tstep=tstep, strend=strend)
config, t, chiz, cf, ratez, neqz = rsm(loading=loading, chi0=chi0, depthS=depthS)
plot(config, t, cf, ratez, neqz)
