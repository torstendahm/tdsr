#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import find_packages, setup
import os

short_description = "TODO"

version = "0.0.1"

PROJECT_ROOT = os.path.dirname(os.path.realpath(__file__))

try:
    readme_rst = os.path.join(PROJECT_ROOT, "README.rst")
    readme_md = os.path.join(PROJECT_ROOT, "README.md")
    if os.path.isfile(readme_rst):
        with open(readme_rst) as readme_file:
            long_description = readme_file.read()
    elif os.path.isfile(readme_md):
        import m2r

        long_description = m2r.parse_from_file(readme_md)
    else:
        raise AssertionError("No readme file")
except (ImportError, AssertionError):
    long_description = short_description

requirements = ["Click>=6.0"]
test_requirements = [
    "tox",
    "pytest",
    "pytest-cov",
    "pytest-xdist",
    "pytest-sugar",
    "mypy",
    "pyfakefs",
    "pytest-subtests",
    "types-setuptools",
]
coverage_requirements = ["coverage"]
formatting_requirements = ["flake8", "black==21.12b0", "isort"]
tool_requirements = [
    "m2r",
    "invoke",
    "pre-commit",
    "cookiecutter",
    "bump2version",
]
dev_requirements = (
    requirements
    + test_requirements
    + coverage_requirements
    + formatting_requirements
    + tool_requirements
)

setup(
    author="torstendahm",
    author_email="contact@romnn.com",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Environment :: Console",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3 :: Only",
        "Operating System :: POSIX",
        "Operating System :: MacOS",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    entry_points={"console_scripts": ["tdsm=tdsm.cli:tdsm"]},
    python_requires=">=3.6",
    install_requires=requirements,
    setup_requires=tool_requirements,
    tests_require=test_requirements,
    extras_require=dict(dev=dev_requirements, test=test_requirements),
    license="GPLv3",
    description=short_description,
    long_description=long_description,
    include_package_data=True,
    package_data={"tdsm": []},
    keywords="tdsm",
    name="tdsm",
    packages=find_packages(include=["tdsm"]),
    test_suite="tests",
    url="https://github.com/torstendahm/tdsm",
    version=version,
    zip_safe=False,
)
