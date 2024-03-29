#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import find_packages, setup
from pathlib import Path
import importlib

spec = importlib.util.spec_from_file_location("version", "tdsr/version.py")
version = importlib.util.module_from_spec(spec)
spec.loader.exec_module(version)

Path().expanduser()

REPO_ROOT = Path(__file__).parent

short_description = "tdsr"

try:
    readme_rst = REPO_ROOT / "README.rst"
    with open(readme_rst) as readme_file:
        long_description = readme_file.read()
except (ImportError, AssertionError):
    long_description = short_description

requirements = ["click", "numpy", "toml"]
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
coverage_requirements = []
formatting_requirements = ["flake8", "black", "isort"]
tool_requirements = [
    "invoke",
    "pre-commit",
    "bump2version",
]
dev_requirements = (
    requirements + test_requirements + coverage_requirements + formatting_requirements + tool_requirements
)

setup(
    author="torstendahm",
    author_email="torsten.dahm@gfz-potsdam.de",
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
    entry_points={"console_scripts": ["tdsr=tdsr.cli:tdsr"]},
    python_requires=">=3.7",
    install_requires=requirements,
    setup_requires=tool_requirements,
    tests_require=test_requirements,
    extras_require=dict(dev=dev_requirements, test=test_requirements),
    license="GPLv3",
    description=short_description,
    long_description=long_description,
    include_package_data=True,
    package_data={"tdsr": []},
    keywords="tdsr",
    name="tdsr",
    packages=find_packages(exclude=["tests"]),
    test_suite="tests",
    url="https://github.com/torstendahm/tdsr",
    version=version.version,
    zip_safe=False,
)
