import sys
sys.path.insert(0, "../")
from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional, RSM
from tdsm.plotting import plot
from tdsm.loading import StepLoading, FourPointLoading, BackgroundLoading
import matplotlib.pyplot as plt

current_dir = Path(__file__).parent
config_file = current_dir / "config.toml"

print("Initialise tdsm runs")
tdsm = TDSM(config=Config.open(config_file))

# ---- define parameter ----
strend = 7.E-5
sstep  = 1.E0

chi0   = 1.E4
depthS = -0.5E0

deltaS = -depthS/60.
sigma_max = 3000.*deltaS

print("First estimate distribution chiz assuming steady state background")
loading = BackgroundLoading(_config=tdsm.config, strend=strend)
config, t, chiz_background, cf, ratez, neqz = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max)
plot(config, t, cf, ratez, neqz)

print("Calculate rate for step function starting from undepleted model")
loading = StepLoading(_config=tdsm.config, sstep=sstep, strend=strend)
config, t, chiz, cf, ratez, neqz = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max)
plot(config, t, cf, ratez, neqz)

print("Calculate rate for step function starting from equilibrium distribution, and use step time short time after start")
loading = StepLoading(_config=tdsm.config, sstep=sstep, strend=strend, tstep=10_000)
config, t, chiz, cf, ratez, neqz = tdsm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max)
plot(config, t, cf, ratez, neqz)
