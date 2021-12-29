from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional
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
plot(config, t, cf, ratez, neqz)

