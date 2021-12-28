from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional
from tdsm.plotting import plot
from tdsm.loading import StepLoading
import matplotlib.pyplot as plt

current_dir = Path(__file__).parent
config_file = current_dir / "config.toml"

print("unklar was dieser Befehl macht")
tdsm = TDSM(config=Config.open(config_file))

# ----Rechne t, cf, ratez, neqz aus mit Werten aus config, wobei zwei Parameter ueberschrieben werden
print("Rechnung")
# config, t, cf, ratez, neqz = tdsm()
# config, t, cf, ratez, neqz = tdsm(chi0=10)
# config, t, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.1)
loading = StepLoading(_config=tdsm.config, tstep=10_000)
config, t, cf, ratez, neqz = tdsm(loading=loading, chi0=10.0, depthS=-0.2)
# config, t, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.8)
# config, t, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.2, tend=86_400)
# config, t, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.2, tend=200_000)

plot(config, t, cf, ratez, neqz)
