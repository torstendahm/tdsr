#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for tdsm package."""

import os
import tempfile
import typing

import pytest
from click.testing import CliRunner

from tdsm import LCM, TDSM, Traditional, cli

TEST_DIR = os.path.dirname(os.path.realpath(__file__))


def test_tdsm_cli_tdsm_method() -> None:
    """Test the CLI tdsm method."""
    runner = CliRunner()
    result = runner.invoke(cli.tdsm)
    assert result.exit_code == 0
    help_result = runner.invoke(cli.tdsm, ["--help"])
    print(help_result.output)
    assert help_result.exit_code == 0
    assert "TDSM method" in help_result.output


def test_tdsm_cli_lcm_method() -> None:
    """Test the CLI lcm method."""
    runner = CliRunner()
    result = runner.invoke(cli.tdsm)
    assert result.exit_code == 0
    help_result = runner.invoke(cli.tdsm, ["lcm", "--help"])
    print(help_result.output)
    assert help_result.exit_code == 0
    assert "LCM method" in help_result.output


def test_tdsm_cli_traditional_method() -> None:
    """Test the CLI traditional method."""
    runner = CliRunner()
    result = runner.invoke(cli.tdsm)
    assert result.exit_code == 0
    help_result = runner.invoke(cli.tdsm, ["traditional", "--help"])
    assert help_result.exit_code == 0
    assert "Traditional method" in help_result.output


def test_tdsm() -> None:
    tdsm = TDSM()
    config, t, cf, ratez, neqz = tdsm()
    # todo: use a results store of pickled data here?
    pass
