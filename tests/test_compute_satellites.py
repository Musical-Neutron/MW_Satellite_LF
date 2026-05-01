#!/usr/bin/env python3

import unittest

import numpy as np
from compute_satellites import cartesian_to_spherical


class CartesianToSphericalTests(unittest.TestCase):

    def setUp(self) -> None:
        self.cart_vector = np.asarray([1., 2.5, 3.])
        self.cart_vector_list = [1., 2.5, 3.]
        self.cart_vector_set = np.row_stack(([1., 2.5, 3.], [2.5, 3., 3.5]))
        self.cart_vector_set_list = [[1., 2.5, 3.], [2.5, 3., 3.5]]

        self.expected_sph_vector = np.asarray(
            [4.0311288741493, 0.73144738125492, 1.1902899496825])
        self.expected_sph_vector_set = np.row_stack(
            ([4.0311288741493, 0.73144738125492, 1.1902899496825],
             [5.2440442408508, 0.8400523908062, 0.87605805059819]))
        return super().setUp()

    def test_cartesian_to_spherical_validation(self):
        """Check function validation handling.
        """
        self.assertRaises(TypeError, cartesian_to_spherical, 's')
        self.assertRaises(TypeError, cartesian_to_spherical, 1)
        self.assertRaises(TypeError, cartesian_to_spherical, 1.5)
        self.assertRaises(ValueError, cartesian_to_spherical, [1, 2.5])
        self.assertRaises(ValueError, cartesian_to_spherical,
                          [1., 2.5, 3., 4.5])
        self.assertRaises(ValueError, cartesian_to_spherical,
                          [[1., 2.5, 3., 4.5], [5., 6.5]])

        return None

    def test_cartesian_to_spherical_function(self):
        """Check function operates as expected.
        """
        sph_vec = cartesian_to_spherical(self.cart_vector)
        sph_vec_list = cartesian_to_spherical(self.cart_vector_list)
        sph_vec_set = cartesian_to_spherical(self.cart_vector_set)
        sph_vec_set_list = cartesian_to_spherical(self.cart_vector_set_list)

        self.assertIsNone(
            np.testing.assert_array_almost_equal(sph_vec,
                                                 self.expected_sph_vector))
        self.assertIsNone(
            np.testing.assert_array_almost_equal(sph_vec_set,
                                                 self.expected_sph_vector_set))
        self.assertIsNone(
            np.testing.assert_array_almost_equal(sph_vec_list,
                                                 self.expected_sph_vector))
        self.assertIsNone(
            np.testing.assert_array_almost_equal(sph_vec_set_list,
                                                 self.expected_sph_vector_set))

        return None
