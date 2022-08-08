#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test custom loading"""

import numpy as np

from .utils import DATA_DIR
from tdsr.loading import CustomLoading
from tests.test_fig5bc import compute_groningen_stress_loading


def test_custom_loading_file_and_data(fs):
    fs.add_real_directory(DATA_DIR)

    dsig = 1.0  # Asig = dsig in the RS model
    ta = 3e5  # [years] relaxation time ta=dsig/dotsigc
    dotsigc = dsig / ta
    tstart = 1960  # start time for model simulation
    tend = 2022.0  # end time for model simulation
    _, data = compute_groningen_stress_loading(tstart=tstart, tend=tend)
    tgrid = data[:, 0]
    dt = tgrid[1] - tgrid[0]

    loading_params = dict(
        scal_t=1.0,
        scal_cf=1.0,
        strend=dotsigc,
        c_tstart=tstart,
        tstart=tstart,
        tend=tend,
        deltat=dt,
    )

    ascii_data_file = "/tmp/ascii_loading.dat"
    fs.create_file(ascii_data_file, contents="test")
    with open(ascii_data_file, "wb") as f:
        np.savetxt(
            f,
            data,
            header="\n".join(
                [
                    "Ascii loading file is expected to have two header lines",
                    "1st column: time, 2nd column: coulomb stress",
                ]
            ),
        )

    expected = CustomLoading(data=data, **loading_params)
    ascii_loading = CustomLoading(file=ascii_data_file, **loading_params)

    assert np.allclose(expected.values(0), ascii_loading.values(0))
    # assert np.allclose(expected.stress_rate, ascii_loading.stress_rate)
