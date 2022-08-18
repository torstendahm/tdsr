===============================
Citations
===============================
Dahm, Torsten; Hainzl, Sebastian; Dahm, Roman Albert (2022): Time-dependent stress response seismicity models (TDSR). GFZ Data Services. https://doi.org/10.5880/GFZ.2.1.2022.002
Dahm, Torsten and hainzl, Sebastian (2022): A Coulomb Stress response model for time-dependent earthquake forecasts. Journal of Geophysical Research, Solid Earth,  http://doi.org/10.1029/2022JB024443

===============================
tdsr
===============================

.. image:: https://github.com/torstendahm/tdsr/workflows/test/badge.svg
        :target: https://github.com/torstendahm/tdsr/actions
        :alt: Build Status

.. image:: https://img.shields.io/pypi/v/tdsr.svg
        :target: https://pypi.python.org/pypi/tdsr
        :alt: PyPI version

.. image:: https://img.shields.io/github/license/torstendahm/tdsr
        :target: https://github.com/torstendahm/tdsr
        :alt: License

.. image:: https://img.shields.io/badge/docs-tdsr-green
        :target: https://torstendahm.github.io/tdsr
        :alt: Documentation

""""""""

Time dependent stress response seismicity model.

A modified effective Coulomb failure model, that calculates earthquake rate as a function of stress loading and model parameter

.. code-block:: console

    $ pip install tdsr

Usage
-----

Usage instructions t.b.a.

.. code-block:: python

    import tdsr


Development
-----------

If you want to contribute to ``tdsr``, we highly recommend setting up a virtual environment for local development. You can easily do so with ``pipenv``, which can be installed with ``pip install --user pipenv`` (`instructions <https://pipenv.pypa.io/en/latest/install/>`_) by following the steps below. Note that all commands should be run in the top level directory of the cloned ``tdsr`` repository.

.. code-block:: bash

    $ git clone https://github.com/torstendahm/tdsr.git
    $ cd tdsr
    $ pipenv install --dev

With the virtual environment set up, activate it with ``pipenv shell``. To exit the virtual environment, run ``exit``.

+++++++
Testing
+++++++

Tests are located in the ``/tests`` directory. To add a new test, add a file with a function name starting with ``test``. To run the full test suite, run:

.. code-block:: bash

    $ invoke test
    $ pytest # can also directly run pytest with custom options

This will run tests serially and fail on the first error, which is preferred for local debugging.
However, this behaviour can also be changed with ``invoke test --[no-]fail-fast --parallel``.

++++++++++
Code style
++++++++++

``tdsr`` provides automatic code formatting and checks code for conformance.

.. code-block:: bash

    $ invoke format   # format python files
    $ invoke lint     # lint python files
    $ invoke type-check     # check for mypy typing errors
    $ invoke pre-commit     # run all of above checks

To pass CI, ensure that ``invoke pre-commit`` passes all checks.

+++++++++++++
Documentation
+++++++++++++

`Sphinx <https://www.sphinx-doc.org/en/master/>`_ is used to build the documentation website, which uses reStructuredText (rst).
See `this website <https://sublime-and-sphinx-guide.readthedocs.io/en/latest/>`_ for examples of how to use rst.

To build and preview the documentation locally, you can use

.. code-block:: bash

    $ invoke docs

This will start a local webserver and open the documentation in the browser for you.
If this is not needed, use the ``--no-serve`` flag.

++++++++++
Releases
++++++++++

**Note:** This is only relevant for maintainers.

Given a clean working tree on the master branch, a release of a new version of ``tdsr`` to `pypi <https://pypi.org/>`_
can be triggered by creating and pushing a new git tag.

The github action CICD pipeline will test, package, and publish the new version automatically.
Before attempting a new release, please make sure all code checks pass.

We use semantic version strings of the form ``{major}.{minor}.{patch}``.
Depending on the type of release, choose the part of the version to be incremented in the command below.

.. code-block:: bash

    $ bump2version patch
    $ git push --follow-tags

Remember to activate the virtual environment as previously outlined if ``bump2version`` is not found.
