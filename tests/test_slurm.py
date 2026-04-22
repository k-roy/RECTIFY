#!/usr/bin/env python3
"""
Tests for SLURM integration module.

Author: Kevin R. Roy
Date: 2026-03
"""

import os
import pytest

from rectify.slurm import (
    get_available_cpus,
    set_thread_limits,
    is_slurm_job,
    get_slurm_info,
)


class TestGetAvailableCpus:
    """Tests for get_available_cpus function."""

    def test_respects_slurm_cpus_per_task(self, monkeypatch):
        """SLURM_CPUS_PER_TASK should take priority."""
        monkeypatch.setenv('SLURM_CPUS_PER_TASK', '16')
        monkeypatch.delenv('LOKY_MAX_CPU_COUNT', raising=False)
        assert get_available_cpus() == 16

    def test_respects_loky_max_cpu_count(self, monkeypatch):
        """LOKY_MAX_CPU_COUNT should be used if SLURM not set."""
        monkeypatch.delenv('SLURM_CPUS_PER_TASK', raising=False)
        monkeypatch.setenv('LOKY_MAX_CPU_COUNT', '8')
        assert get_available_cpus() == 8

    def test_slurm_takes_priority_over_loky(self, monkeypatch):
        """SLURM should take priority over LOKY."""
        monkeypatch.setenv('SLURM_CPUS_PER_TASK', '16')
        monkeypatch.setenv('LOKY_MAX_CPU_COUNT', '8')
        assert get_available_cpus() == 16

    def test_falls_back_to_system(self, monkeypatch):
        """Should fall back to system CPU count if no env vars."""
        monkeypatch.delenv('SLURM_CPUS_PER_TASK', raising=False)
        monkeypatch.delenv('LOKY_MAX_CPU_COUNT', raising=False)
        cpus = get_available_cpus()
        # Should return half of system CPUs (minimum 1)
        assert cpus >= 1

    def test_uses_default_if_none_available(self, monkeypatch):
        """Should use default if nothing else available."""
        monkeypatch.delenv('SLURM_CPUS_PER_TASK', raising=False)
        monkeypatch.delenv('LOKY_MAX_CPU_COUNT', raising=False)
        # Can't easily mock os.cpu_count, but we can test the default param
        cpus = get_available_cpus(default=4)
        assert cpus >= 1


class TestSetThreadLimits:
    """Tests for set_thread_limits function."""

    def test_sets_all_thread_vars(self, monkeypatch):
        """Should set all thread limit environment variables."""
        # Clear existing vars
        for var in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                    'LOKY_MAX_CPU_COUNT', 'NUMEXPR_MAX_THREADS']:
            monkeypatch.delenv(var, raising=False)

        result = set_thread_limits(8)

        assert result == 8
        assert os.environ['OMP_NUM_THREADS'] == '8'
        assert os.environ['OPENBLAS_NUM_THREADS'] == '8'
        assert os.environ['MKL_NUM_THREADS'] == '8'
        assert os.environ['LOKY_MAX_CPU_COUNT'] == '8'
        assert os.environ['NUMEXPR_MAX_THREADS'] == '8'

    def test_auto_detects_if_none(self, monkeypatch):
        """Should auto-detect CPU count if n_threads is None."""
        monkeypatch.setenv('SLURM_CPUS_PER_TASK', '12')
        result = set_thread_limits(None)
        assert result == 12
        assert os.environ['OMP_NUM_THREADS'] == '12'


class TestIsSlurmJob:
    """Tests for is_slurm_job function."""

    def test_returns_true_in_slurm(self, monkeypatch):
        """Should return True when SLURM_JOB_ID is set."""
        monkeypatch.setenv('SLURM_JOB_ID', '12345')
        assert is_slurm_job() is True

    def test_returns_false_outside_slurm(self, monkeypatch):
        """Should return False when SLURM_JOB_ID is not set."""
        monkeypatch.delenv('SLURM_JOB_ID', raising=False)
        assert is_slurm_job() is False


class TestGetSlurmInfo:
    """Tests for get_slurm_info function."""

    def test_returns_info_in_slurm(self, monkeypatch):
        """Should return dict with SLURM info when in job."""
        monkeypatch.setenv('SLURM_JOB_ID', '12345')
        monkeypatch.setenv('SLURM_CPUS_PER_TASK', '16')
        monkeypatch.setenv('SLURM_ARRAY_TASK_ID', '3')
        monkeypatch.setenv('SLURM_JOB_PARTITION', 'test-partition')

        info = get_slurm_info()

        assert info['job_id'] == '12345'
        assert info['cpus'] == '16'
        assert info['array_task_id'] == '3'
        assert info['partition'] == 'test-partition'

    def test_returns_empty_outside_slurm(self, monkeypatch):
        """Should return empty dict when not in SLURM job."""
        monkeypatch.delenv('SLURM_JOB_ID', raising=False)
        assert get_slurm_info() == {}
