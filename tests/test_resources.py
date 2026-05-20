"""Tests for resource policy helpers."""

import argparse

import pytest

from rectify.core.commands.run_command import create_run_parser
from rectify.core.utils.resources import (
    detect_machine_class,
    resolve_aligner_concurrency,
)


def test_detect_machine_class_small_laptop():
    assert detect_machine_class(total_mem_gb=8, cpu_count=4) == "m1_laptop"


def test_detect_machine_class_cluster_when_memory_or_cpu_large():
    assert detect_machine_class(total_mem_gb=64, cpu_count=4) == "cluster_node"
    assert detect_machine_class(total_mem_gb=8, cpu_count=16) == "cluster_node"


def test_resolve_aligner_concurrency_auto_laptop_disables_axis_b():
    assert resolve_aligner_concurrency("auto", 8, "m1_laptop") == 1


def test_resolve_aligner_concurrency_auto_cluster_reserves_two_threads():
    assert resolve_aligner_concurrency("auto", 16, "cluster_node") == 14
    assert resolve_aligner_concurrency("auto", 1, "cluster_node") == 1


def test_resolve_aligner_concurrency_numeric_clamps_to_total_threads():
    assert resolve_aligner_concurrency("4", 16, "cluster_node") == 4
    assert resolve_aligner_concurrency("99", 16, "cluster_node") == 16


@pytest.mark.parametrize("value", ["0", "-1", "many"])
def test_resolve_aligner_concurrency_rejects_invalid_values(value):
    with pytest.raises(ValueError):
        resolve_aligner_concurrency(value, 16, "cluster_node")


def test_run_all_parser_accepts_aligner_concurrency():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="cmd")
    create_run_parser(subparsers)

    args = parser.parse_args(
        ["run-all", "input.bam", "--Scer", "-o", "out", "--aligner-concurrency", "4"]
    )

    assert args.aligner_concurrency == "4"


@pytest.mark.parametrize("value", ["0", "-1", "many"])
def test_run_all_parser_rejects_invalid_aligner_concurrency(value):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="cmd")
    create_run_parser(subparsers)

    with pytest.raises(SystemExit):
        parser.parse_args(["run-all", "input.bam", "--Scer", "-o", "out", "--aligner-concurrency", value])
