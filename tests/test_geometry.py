#!/usr/bin/env python3

import unittest
from unittest.mock import patch
import pytest
import numpy as np
from mwsatslf.geometry import (
    cartesian_to_spherical,
    random_points_on_sphere_surface,
    uniform_points_on_sphere_surface,
)


@pytest.fixture
def single_vector():
    """Single (x, y, z) as 1D array."""
    return np.array([1.0, 2.5, 3.0])


@pytest.fixture
def single_list():
    """Same vector as Python list."""
    return [1.0, 2.5, 3.0]


@pytest.fixture
def multiple_vectors():
    """2x3 array of vectors."""
    return np.array([[1.0, 2.5, 3.0], [2.5, 3.0, 3.5]])


@pytest.fixture
def multiple_lists():
    """List of lists."""
    return [[1.0, 2.5, 3.0], [2.5, 3.0, 3.5]]


@pytest.fixture
def expected_single():
    """Expected spherical coordinates for single_vector."""
    r = np.sqrt(1.0**2 + 2.5**2 + 3.0**2)
    theta = np.arccos(3.0 / r)
    phi = np.arctan2(2.5, 1.0)
    return np.array([r, theta, phi])


@pytest.fixture
def expected_multiple():
    """Expected spherical coordinates for multiple_vectors."""
    r1 = np.sqrt(1.0**2 + 2.5**2 + 3.0**2)
    theta1 = np.arccos(3.0 / r1)
    phi1 = np.arctan2(2.5, 1.0)
    r2 = np.sqrt(2.5**2 + 3.0**2 + 3.5**2)
    theta2 = np.arccos(3.5 / r2)
    phi2 = np.arctan2(3.0, 2.5)
    return np.array([[r1, theta1, phi1], [r2, theta2, phi2]])


class TestCartesianToSpherical:
    """Test cartesian_to_spherical with various inputs."""

    def test_single_vector_array(self, single_vector, expected_single):
        """1D array returns 1D array of (r, theta, phi)."""
        result = cartesian_to_spherical(single_vector)
        assert isinstance(result, np.ndarray)
        assert result.shape == (3,)
        np.testing.assert_allclose(result, expected_single)

    def test_single_vector_list(self, single_list, expected_single):
        """Python list should be handled identically."""
        result = cartesian_to_spherical(single_list)
        assert result.shape == (3,)
        np.testing.assert_allclose(result, expected_single)

    def test_multiple_vectors(self, multiple_vectors, expected_multiple):
        """2D array returns shape (N,3)."""
        result = cartesian_to_spherical(multiple_vectors)
        assert result.shape == (2, 3)
        np.testing.assert_allclose(result, expected_multiple)

    def test_multiple_lists(self, multiple_lists, expected_multiple):
        """List of lists returns shape (N,3)."""
        result = cartesian_to_spherical(multiple_lists)
        assert result.shape == (2, 3)
        np.testing.assert_allclose(result, expected_multiple)

    def test_single_2d_array(self):
        """A single vector given as shape (1,3) returns (1,3)."""
        coord = np.array([[1.0, 2.5, 3.0]])
        result = cartesian_to_spherical(coord)
        assert result.shape == (1, 3)
        r = np.sqrt(1**2 + 2.5**2 + 3**2)
        expected = np.array([[r, np.arccos(3 / r), np.arctan2(2.5, 1)]])
        np.testing.assert_allclose(result, expected)

    # Validation / error handling

    @pytest.mark.parametrize(
        "bad_input",
        [
            "s",
            1,
            1.5,
            None,
        ],
    )
    def test_type_errors(self, bad_input):
        """Non-list/array inputs raise TypeError."""
        with pytest.raises(TypeError):
            cartesian_to_spherical(bad_input)

    @pytest.mark.parametrize(
        "bad_shape",
        [
            [1, 2.5],
            [1.0, 2.5, 3.0, 4.5],
            [[1.0, 2.5, 3.0, 4.5], [5.0, 6.5]],
        ],
    )
    def test_value_errors(self, bad_shape):
        """Incorrect shape raises ValueError."""
        with pytest.raises(ValueError):
            cartesian_to_spherical(bad_shape)


class TestRandomPointsOnSphereSurface:
    """Tests for random_points_on_sphere_surface."""

    def test_shape(self):
        n = 50
        points = random_points_on_sphere_surface(n)
        assert points.shape == (n, 3)

    def test_points_on_unit_sphere(self):
        """All points must lie exactly on the unit sphere."""
        points = random_points_on_sphere_surface(1000)
        norms = np.linalg.norm(points, axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-15)

    def test_default_full_sphere(self):
        """With no bounds, z covers [-1, 1]."""
        np.random.seed(0)
        points = random_points_on_sphere_surface(5000)
        z = points[:, 2]
        assert np.min(z) >= -1.0
        assert np.max(z) <= 1.0
        # With many points, extremes should be reached
        assert np.min(z) < -0.9
        assert np.max(z) > 0.9

    def test_cos_theta_min_clips(self):
        """cos_theta_min below -1.0 is clipped to -1.0."""
        points = random_points_on_sphere_surface(500, cos_theta_min=-5.0)
        z = points[:, 2]
        assert np.min(z) >= -1.0

    def test_cos_theta_max_clips(self):
        """cos_theta_max above 1.0 is clipped to 1.0."""
        points = random_points_on_sphere_surface(500, cos_theta_max=5.0)
        z = points[:, 2]
        assert np.max(z) <= 1.0

    def test_bounded_band(self):
        """Given bounds, all z lie inside the interval."""
        n = 2000
        low, high = -0.5, 0.3
        points = random_points_on_sphere_surface(
            n, cos_theta_min=low, cos_theta_max=high
        )
        z = points[:, 2]
        assert np.all((z >= low) & (z <= high))
        # Closeness to boundaries should be plausible
        assert np.min(z) < low + 0.1
        assert np.max(z) > high - 0.1

    def test_swapped_bounds(self):
        """If cos_theta_min > cos_theta_max, they are swapped internally."""
        n = 1000
        # Swap them
        points = random_points_on_sphere_surface(
            n, cos_theta_min=0.5, cos_theta_max=-0.1
        )
        z = points[:, 2]
        # After swapping, the effective bounds are [-0.1, 0.5]
        assert np.all((z >= -0.1) & (z <= 0.5))

    def test_zero_points(self):
        """Zero points returns an empty array."""
        points = random_points_on_sphere_surface(0)
        assert points.shape == (0, 3)

    def test_n_negative(self):
        """Negative n raises ValueError from np.empty (not specifically handled)."""
        with pytest.raises(ValueError):
            random_points_on_sphere_surface(-5)

    def test_single_point(self):
        """Single point is generated correctly."""
        points = random_points_on_sphere_surface(1)
        assert points.shape == (1, 3)
        assert np.isclose(np.linalg.norm(points[0]), 1.0)


class TestUniformPointsOnSphereSurface:
    """Tests for uniform_points_on_sphere_surface."""

    def test_shape(self):
        n = 50
        points = uniform_points_on_sphere_surface(n)
        assert points.shape == (n, 3)

    def test_points_on_unit_sphere(self):
        """All points lie on the unit sphere."""
        points = uniform_points_on_sphere_surface(200)
        norms = np.linalg.norm(points, axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-15)

    def test_uniform_z_distribution(self):
        """The z coordinates should be uniformly spaced between -1 and 1."""
        n = 100
        points = uniform_points_on_sphere_surface(n)
        z = points[:, 2]
        expected_z = 1.0 - (2.0 / n) * np.arange(n)
        np.testing.assert_allclose(z, expected_z, atol=1e-15)

    def test_xy_coordinates(self):
        """Directly check x and y against golden-angle spiral formula."""
        n = 30
        dPhi = np.pi * (3.0 - 5.0**0.5)
        points = uniform_points_on_sphere_surface(n)
        step = np.arange(n)
        z = points[:, 2]
        r_xy = np.sqrt(1.0 - z**2)
        phi_expected = dPhi * step
        x_expected = r_xy * np.cos(phi_expected)
        y_expected = r_xy * np.sin(phi_expected)
        # All x,y must match, even at the pole where r_xy=0 -> x=y=0.
        np.testing.assert_allclose(points[:, 0], x_expected, atol=1e-15)
        np.testing.assert_allclose(points[:, 1], y_expected, atol=1e-15)

    def test_xyz_normalization(self):
        """Quick check that x,y,z correspond to the unit sphere formula."""
        n = 50
        points = uniform_points_on_sphere_surface(n)
        z = points[:, 2]
        rxy = np.sqrt(1.0 - z**2)
        x = points[:, 0]
        y = points[:, 1]
        np.testing.assert_allclose(np.sqrt(x**2 + y**2), rxy, atol=1e-15)

    def test_zero_points_raises(self):
        """Zero points should raise ZeroDivisionError (dz = 2/0)."""
        with pytest.raises(ZeroDivisionError):
            uniform_points_on_sphere_surface(0)

    def test_negative_n_raises(self):
        """Negative n produces empty arange? It raises ValueError in np.empty."""
        with pytest.raises(ValueError):
            uniform_points_on_sphere_surface(-5)
