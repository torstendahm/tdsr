#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test scipy pdf drop in replacement"""

import numpy as np
from scipy.stats import norm
from tdsr.utils import pdf


def test_pdf():
    for loc in np.linspace(-3.0, 3.0, 10):
        for scale in np.linspace(-1.0, 10.0, 10):
            for x in np.linspace(-10.0, 10.0, 100):
                print(x, loc, scale)
                expected = norm.pdf(x, loc=loc, scale=scale)
                given = pdf(x, loc=loc, scale=scale)
                assert np.allclose([expected], [given], equal_nan=True)
