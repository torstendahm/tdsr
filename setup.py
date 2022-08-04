#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import find_packages, setup
from pathlib import Path
import tdsr

Path().expanduser()

REPO_ROOT = Path(__file__).parent

short_description = "tdsr"

# version = "0.0.1"

try:
    readme_rst = REPO_ROOT / "README.rst"
    with open(readme_rst) as readme_file:
        long_description = readme_file.read()
except (ImportError, AssertionError):
    long_description = short_description

requirements = ["click"]
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
    entry_points={"console_scripts": ["tdsr=tdsr.cli:tdsr"]},
    python_requires=">=3.6",
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
    packages=find_packages(include=["tdsr"]),
    test_suite="tests",
    url="https://github.com/torstendahm/tdsr",
    version=tdsr.__version__,
    zip_safe=False,
)
