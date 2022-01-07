from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional, RSM
from tdsm.plotting import plot
from tdsm.loading import StepLoading, FourPointLoading, BackgroundLoading
import matplotlib.pyplot as plt

current_dir = Path(__file__).parent
config_file = current_dir / "config.toml"

print("Aufsetzen eines tdsm runs")
tdsm = TDSM(config=Config.open(config_file))

# ----Rechne t, cf, ratez, neqz aus mit Werten aus config, wobei Parameter aus config ueberschrieben werden koennen
# config, t, chiz, cf, ratez, neqz = tdsm()
# config, t, chiz, cf, ratez, neqz = tdsm(chi0=10)
# config, t, chiz, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.1)
# config, t, chiz, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.8)
# config, t, chiz, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.2, tend=86_400)
# config, t, chiz, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.2, tend=200_000)
# loading = FourPointLoading(_config=tdsm.config)


<<<<<<< HEAD
strend = 7.E-5
sstep  = 1.E0
#sstep  = 1.E-1
#sstep  = 1.E+1
#sstep  = 2.E0

chi0   = 1.E4
depthS = -0.5E0

#deltaS = -depthS/40.
#sigma_max = 2000.*deltaS
deltaS = -depthS/60.
sigma_max = 3000.*deltaS

print("Rechnung Background")
# ---- chiz mit BackgroundLoading ausrechnen (vom Ende der Iteration) 
loading = BackgroundLoading(_config=tdsm.config, strend=strend)
#config, t, chiz_background, cf, ratez, neqz = tdsm(loading=loading, chi0=10.0, depthS=-0.2)
config, t, chiz_background, cf, ratez, neqz = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max)
plot(config, t, cf, ratez, neqz)

print("Rechnung mit step function wie bisher")
#loading = StepLoading(_config=tdsm.config, tstep=10_000)
loading = StepLoading(_config=tdsm.config, sstep=sstep, strend=strend)
#config, t, chiz, cf, ratez, neqz = tdsm(loading=loading, chi0=10.0, depthS=-0.2)
config, t, chiz, cf, ratez, neqz = tdsm(loading=loading, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max)
plot(config, t, cf, ratez, neqz)

print("Rechnung mit step function und Background chiz")
# ---- chiz mit BackgroundLoading ausrechnen (vom Ende der Iteration) 
loading = StepLoading(_config=tdsm.config, sstep=sstep, strend=strend, tstep=10_000)
#config, t, chiz, cf, ratez, neqz = tdsm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=10.0, depthS=-0.2)
config, t, chiz, cf, ratez, neqz = tdsm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=chi0, depthS=depthS, deltaS=deltaS, sigma_max=sigma_max)
=======
print("Rechnung Background")
# ---- chiz mit BackgroundLoading ausrechnen (vom Ende der Iteration) 
loading = BackgroundLoading(_config=tdsm.config)
config, t, chiz_background, cf, ratez, neqz = tdsm(loading=loading, chi0=10.0, depthS=-0.2)
plot(config, t, cf, ratez, neqz)

print("Rechnung mit step function wie bisher")
#loading = StepLoading(_config=tdsm.config, tstep=10_000)
loading = StepLoading(_config=tdsm.config)
config, t, chiz, cf, ratez, neqz = tdsm(loading=loading, chi0=10.0, depthS=-0.2)
plot(config, t, cf, ratez, neqz)

print("Rechnung mit step function und Background chiz")
# ---- chiz mit BackgroundLoading ausrechnen (vom Ende der Iteration) 
loading = StepLoading(_config=tdsm.config, tstep=20_000)
config, t, chiz, cf, ratez, neqz = tdsm(loading=loading, equilibrium=True, chiz=chiz_background, chi0=10.0, depthS=-0.2)
>>>>>>> refs/remotes/origin/master
plot(config, t, cf, ratez, neqz)

