#!/usr/bin/env python3

import unittest

import numpy as np
from compute_satellites import spherical_to_cartesian


class SphericalToCartesianTests(unittest.TestCase):

    def setUp(self) -> None:
        self.sph_vector = np.asarray([5., 0.25, 0.3])
        self.sph_vector_list = [5., 0.25, 0.3]
        self.sph_vector_set = np.row_stack(([5., 0.25, 0.3], [7., 0.45, 0.9]))
        self.sph_vector_set_list = [[5., 0.25, 0.3], [7., 0.45, 0.9]]

        self.expected_cart_vector = np.asarray(
            [1.181770149, 0.3655643458, 4.844562109])
        self.expected_cart_vector_set = np.row_stack(
            ([1.181770149, 0.3655643458,
              4.844562109], [1.892652383, 2.385041453, 6.303129716]))
        return super().setUp()

    def test_spherical_to_cartesian_validation(self):
        """Check function validation handling.
        """
        self.assertRaises(TypeError, spherical_to_cartesian, 's')
        self.assertRaises(TypeError, spherical_to_cartesian, 1)
        self.assertRaises(TypeError, spherical_to_cartesian, 1.5)
        self.assertRaises(ValueError, spherical_to_cartesian, [1, 2.5])
        self.assertRaises(ValueError, spherical_to_cartesian,
                          [1., 2.5, 3., 4.5])
        self.assertRaises(ValueError, spherical_to_cartesian,
                          [[1., 2.5, 3., 4.5], [5., 6.5]])

        return None

    def test_spherical_to_cartesian_function(self):
        """Check function operates as expected.
        """
        self.assertIsNone(
            np.testing.assert_array_almost_equal(
                spherical_to_cartesian(self.sph_vector),
                self.expected_cart_vector))
        self.assertIsNone(
            np.testing.assert_array_almost_equal(
                spherical_to_cartesian(self.sph_vector_list),
                self.expected_cart_vector))
        self.assertIsNone(
            np.testing.assert_array_almost_equal(
                spherical_to_cartesian(self.sph_vector_set),
                self.expected_cart_vector_set))
        self.assertIsNone(
            np.testing.assert_array_almost_equal(
                spherical_to_cartesian(self.sph_vector_set_list),
                self.expected_cart_vector_set))

        return None
