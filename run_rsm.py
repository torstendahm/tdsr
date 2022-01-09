from pathlib import Path
from tdsm import Config, TDSM, LCM, Traditional, RSM
from tdsm.plotting import plot
from tdsm.loading import StepLoading, FourPointLoading, BackgroundLoading, CyclicLoading
import matplotlib.pyplot as plt

current_dir = Path(__file__).parent
config_file = current_dir / "config.toml"

print("Aufsetzen eines tdsm runs")
rsm = RSM(config=Config.open(config_file))

# ----Rechne t, cf, ratez, neqz aus mit Werten aus config, wobei Parameter aus config ueberschrieben werden koennen
# config, t, chiz, cf, ratez, neqz = tdsm()
# config, t, chiz, cf, ratez, neqz = tdsm(chi0=10)
# config, t, chiz, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.1)
# config, t, chiz, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.8)
# config, t, chiz, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.2, tend=86_400)
# config, t, chiz, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.2, tend=200_000)
# loading = FourPointLoading(_config=tdsm.config)


sstep  = 1.E0
#sstep  = 1.E-5
#sstep  = 0.0
#sstep  = 1.E+3
#sstep  = 2.E0
tstep = 80_000
tstep = 10_000
#tstep = 42_000

ampsin = 2.0
Tsin = 43_200

chi0   = 1.E4

strend = 7.E-5
depthS = -0.5E0
#strend = 0.01
#depthS = -0.2E0

#deltaS = -depthS/40.
#sigma_max = 2000.*deltaS
deltaS = -depthS/60.
sigma_max = 3000.*deltaS

print("Rechnung Background")
# ---- chiz mit BackgroundLoading ausrechnen (vom Ende der Iteration) 
#loading = BackgroundLoading(_config=rsm.config, strend=strend)
#print('loading=',loading.strend)
#config, t, chiz_background, cf, ratez, neqz = rsm(loading=loading, chi0=chi0, depthS=depthS)
#plot(config, t, cf, ratez, neqz)

#loading = CyclicLoading(_config=rsm.config, strend=strend, ampsin=ampsin, Tsin=Tsin)
#print('loading=',loading.strend)
#config, t, chiz, cf, ratez, neqz = rsm(loading=loading, chi0=chi0, depthS=depthS)
#plot(config, t, cf, ratez, neqz)

print("Rechnung mit step function wie bisher")
loading = StepLoading(_config=rsm.config, sstep=sstep, tstep=tstep, strend=strend)
#loading = StepLoading(_config=rsm.config, sstep=sstep, strend=strend)
config, t, chiz, cf, ratez, neqz = rsm(loading=loading, chi0=chi0, depthS=depthS)
plot(config, t, cf, ratez, neqz)
