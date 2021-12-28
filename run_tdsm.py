from tdsm import TDSM, LCM, Traditional
from tdsm.plotting import plot
import matplotlib.pyplot as plt

print('unklar was dieser Befehl macht')
tdsm = TDSM()

# ----Rechne t, cf, ratez, neqz aus mit Werten aus config, wobei zwei Parameter ueberschrieben werden
print('Rechnung')
#config, t, cf, ratez, neqz = tdsm()
#config, t, cf, ratez, neqz = tdsm(chi0=10)
#config, t, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.1)
config, t, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.2)
#config, t, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.8)
#config, t, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.2, tend=86_400)
#config, t, cf, ratez, neqz = tdsm(chi0=10.0, depthS=-0.2, tend=200_000)


print('begin plotting')
plt.plot(t,cf)
plt.xlabel("$t$")
plt.ylabel('$\sigma_c$')
plt.show()
plt.close()

plt.plot(t,ratez)
plt.xlabel("$t$")
plt.ylabel('$r$')
plt.show()
plt.close()

plt.plot(t[0:-2],neqz[0:-1])
plt.xlabel("$t$")
plt.ylabel('$n$')
plt.show()
plt.close()
