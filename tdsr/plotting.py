################################
# Time Dependent Seismicity Model - Plotting results
# T. Dahm, R. Dahm 26.12.2021
################################

import numpy as np
import numpy.typing as npt

from tdsr.config import Config


def plot(
    config: Config,
    t: npt.NDArray[np.float64],
    cf: npt.NDArray[np.float64],
    ratez: npt.NDArray[np.float64],
    neqz: npt.NDArray[np.float64],
) -> None:
    from matplotlib import pyplot as plt

    plt.plot(t, cf)
    plt.xlabel("$t$")
    plt.ylabel("$\sigma_c$")
    plt.show()
    plt.close()

    plt.plot(t, ratez)
    plt.xlabel("$t$")
    plt.ylabel("$r$")
    plt.show()
    plt.close()

    plt.plot(t[0:-2], neqz[0:-1])
    plt.xlabel("$t$")
    plt.ylabel("$n$")
    plt.show()
    plt.close()
