################################
# Time Dependent Seismicity Model - Plotting results
# T. Dahm, R. Dahm 26.12.2021
################################

import numpy as np
import numpy.typing as npt
from tdsm.config import Config


def plot(
    config: Config,
    t: npt.NDArray[np.float64],
    cf: npt.NDArray[np.float64],
    ratez: npt.NDArray[np.float64],
    neqz: npt.NDArray[np.float64],
) -> None:
    from matplotlib import pyplot as plt

    plt.plot(cf)
    plt.show()
    plt.close()
    plt.plot(ratez)
    plt.show()
    plt.close()
    plt.plot(neqz)
    plt.show()
    plt.close()
