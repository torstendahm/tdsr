#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for tdsm package."""

import os
import tempfile
import typing

import pytest
from click.testing import CliRunner

from tdsm import cli

TEST_DIR = os.path.dirname(os.path.realpath(__file__))


def test_tdsm_cli() -> None:
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.tdsm)
    assert result.exit_code == 0
    # assert "Usage: proto-compile" in result.output
    # help_result = runner.invoke(cli.proto_compile, ["--help"])
    # assert help_result.exit_code == 0
    # assert "Show this message and exit." in help_result.output


def test_tdsm() -> None:
    pass
