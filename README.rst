
===============================
proto-compile
===============================

.. image:: https://github.com/romnn/proto-compile/workflows/test/badge.svg
        :target: https://github.com/romnn/proto-compile/actions
        :alt: Build Status

.. image:: https://img.shields.io/pypi/v/proto-compile.svg
        :target: https://pypi.python.org/pypi/proto-compile
        :alt: PyPI version

.. image:: https://img.shields.io/github/license/romnn/proto-compile
        :target: https://github.com/romnn/proto-compile
        :alt: License

.. image:: https://codecov.io/gh/romnn/proto-compile/branch/master/graph/badge.svg
        :target: https://codecov.io/gh/romnn/proto-compile
        :alt: Test Coverage

""""""""

Your short description here. 

.. code-block:: console

    $ pip install proto-compile

Usage
-----

.. code-block:: python

    import proto_compile


Files:
-----

directories example:

run_Fig3b.py
run_Fig5a.py
run_Fig_appendix_A1.py
run_rsm.py
run_traditional.py
run_Fig3a+c.py
run_Fig4.py
run_Fig7c.py
run_lcm.py
run_tdsm.py

config.toml
Model parameter defined on two levels:
(1) general parameters, defined for all simulatiions moduls, including tdsm, rs,. lcm and traditional
(2) parameter for specific stress loading scenarios, as e.g. step, background, cyclic, 4points, trendchange, ramp

Development
-----------

For detailed instructions see the invoke pre-commit
