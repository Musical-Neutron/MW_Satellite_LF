#!/usr/bin/env python3

import numpy as np
import pytest

from mwsatslf.cosmology import critical_density_200, compute_r200

REF_M200 = 1.3432e12  # Msun/h
REF_R200 = 0.2458  # Mpc/h
EXPECTED_RHO_CRIT200 = 3.0 / (4.0 * np.pi) * REF_M200 / REF_R200**3


class TestComputeR200:
    """Tests for compute_r200."""

    def test_scalar(self):
        """Return correct R200 for known values."""
        rho = EXPECTED_RHO_CRIT200
        r200 = compute_r200(REF_M200, rho)
        # Should be very close to original REF_R200
        assert np.isclose(r200, REF_R200, rtol=1e-12)

    def test_array_mass(self):
        """Works element-wise for an array of masses."""
        m200 = np.array([1e11, 1e12, 1e13])
        rho = 200.0  # arbitrary constant
        r200 = compute_r200(m200, rho)
        assert r200.shape == m200.shape
        # Check one element manually
        expected = (m200 / rho / (4.0 * np.pi / 3.0)) ** (1.0 / 3.0)
        np.testing.assert_allclose(r200, expected)

    def test_zero_mass(self):
        """Zero mass gives zero radius."""
        r200 = compute_r200(0.0, 200.0)
        assert r200 == 0.0

    def test_negative_mass(self):
        """Negative mass should produce a negative cube root
        (as defined by numpy)."""
        r200 = compute_r200(-8.0, 1.0)
        expected = (-8.0 / (4.0 * np.pi / 3.0)) ** (1.0 / 3.0)
        assert np.isclose(r200, expected)

    def test_negative_rho(self):
        """Negative rho_crit200 yields a real negative root for negative
        numbers (defined in numpy)."""
        r200 = compute_r200(8.0, -1.0)
        expected = (-8.0 / (4.0 * np.pi / 3.0)) ** (1.0 / 3.0)
        assert np.isclose(r200, expected)

    def test_rho_zero(self):
        """Division by zero should raise an error."""
        with pytest.raises(ZeroDivisionError):
            compute_r200(1.0, 0.0)


class TestCriticalDensity200:
    """Tests for critical_density_200."""

    def test_scalar(self):
        """Return correct 200*rho_crit for the reference halo."""
        rho = critical_density_200(REF_M200, REF_R200)
        assert np.isclose(rho, EXPECTED_RHO_CRIT200, rtol=1e-12)

    def test_inverse_relationship(self):
        """For any positive m200 and r200, compute_r200 applied to
        the result should recover the original r200."""
        m200 = 5e11
        r200 = 0.3
        rho = critical_density_200(m200, r200)
        r_back = compute_r200(m200, rho)
        assert np.isclose(r_back, r200)

    def test_array_inputs(self):
        """Works element-wise with array inputs."""
        m200 = np.array([1e11, 2e11, 3e11])
        r200 = np.array([0.1, 0.15, 0.2])
        rho = critical_density_200(m200, r200)
        expected = 3.0 / (4.0 * np.pi) * m200 / r200**3
        np.testing.assert_allclose(rho, expected)

    def test_zero_radius(self):
        """Division by zero should raise an error."""
        with pytest.raises(ZeroDivisionError):
            critical_density_200(1.0, 0.0)

    def test_zero_mass(self):
        """Zero mass gives zero density (no division by zero)."""
        rho = critical_density_200(0.0, 0.3)
        assert rho == 0.0
