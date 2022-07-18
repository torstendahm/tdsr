
===============================
tdsm
===============================

.. image:: https://github.com/torstendahm/tdsm/workflows/test/badge.svg
        :target: https://github.com/torstendahm/tdsm/actions
        :alt: Build Status

.. image:: https://img.shields.io/pypi/v/tdsm.svg
        :target: https://pypi.python.org/pypi/tdsm
        :alt: PyPI version

.. image:: https://img.shields.io/github/license/torstendahm/tdsm
        :target: https://github.com/torstendahm/tdsm
        :alt: License

""""""""

Your short description here. 

.. code-block:: console

    $ pip install tdsm

Usage
-----

Usage instructions t.b.a.

.. code-block:: python

    import tdsm


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

If you want to contribute to ``tdsm``, we highly recommend setting up a virtual environment for local development. You can easily do so with ``pipenv``, which can be installed with ``pip install --user pipenv`` (`instructions <https://pipenv.pypa.io/en/latest/install/>`_) by following the steps below. Note that all commands should be run in the top level directory of the cloned ``tdsm`` repository.

.. code-block:: bash

    $ git clone https://github.com/torstendahm/tdsm.git
    $ cd tdsm
    $ pipenv install --dev

With the virtual environment set up, activate it with ``pipenv shell``. To exit the virtual environment, run ``exit``.

invoke, linting, formatting, code checks t.b.a
