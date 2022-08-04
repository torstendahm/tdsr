#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for tdsr package."""

import os
import tempfile
import typing
from pathlib import Path

import pytest
from click.testing import CliRunner

from tdsr import LCM, TDSR, Traditional, cli

TEST_DIR = Path(__file__).parent


def test_tdsr_cli_tdsr_method() -> None:
    """Test the CLI tdsr method."""
    runner = CliRunner()
    result = runner.invoke(cli.tdsr)
    assert result.exit_code == 0
    help_result = runner.invoke(cli.tdsr, ["--help"])
    print(help_result.output)
    assert help_result.exit_code == 0
    assert "tdsr method" in help_result.output


def test_tdsr_cli_lcm_method() -> None:
    """Test the CLI lcm method."""
    runner = CliRunner()
    result = runner.invoke(cli.tdsr)
    assert result.exit_code == 0
    help_result = runner.invoke(cli.tdsr, ["lcm", "--help"])
    print(help_result.output)
    assert help_result.exit_code == 0
    assert "LCM method" in help_result.output


def test_tdsr_cli_traditional_method() -> None:
    """Test the CLI traditional method."""
    runner = CliRunner()
    result = runner.invoke(cli.tdsr)
    assert result.exit_code == 0
    help_result = runner.invoke(cli.tdsr, ["traditional", "--help"])
    assert help_result.exit_code == 0
    assert "Traditional method" in help_result.output
