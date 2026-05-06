#!/usr/bin/env python3
from typing import Union

import numpy as np


def cartesian_to_spherical(coordinates: Union[list, np.ndarray]) -> np.ndarray:
    """Converts from Cartesian to Spherical basis.

    Args:
        coordinates (nd.array (N,3)): (x, y, z) values to convert.

    Returns:
        np.ndarray (N,3): coordinates in Spherical basis
            (r, theta, phi).
    """
    # Validation checks
    coordinates_permitted_types = (list, np.ndarray)
    if not isinstance(coordinates, coordinates_permitted_types):
        raise TypeError(
            "coordinates must be one of: {}".format(coordinates_permitted_types)
        )
    else:
        coordinates = np.asarray(coordinates)

    # Check if single vector is supplied or set of vectors
    array_of_arrays = isinstance(coordinates[0], np.ndarray)
    if not array_of_arrays:
        coordinates = coordinates.reshape((1, 3))

    # Function is designed to work in 3D
    if np.size(coordinates, 1) != 3:
        raise ValueError("coordinates should be (N,3) array")

    # Convert to spherical coordinates
    r = np.sqrt((coordinates * coordinates).sum(axis=1))
    theta = np.arccos(coordinates[:, 2] / r)
    phi = np.arctan2(coordinates[:, 1], coordinates[:, 0])

    # Output based on input array
    if array_of_arrays:
        return_array = np.column_stack((r, theta, phi))
    else:
        return_array = np.concatenate((r, theta, phi))

    return return_array


def random_points_on_sphere_surface(n_points, cos_theta_min=None, cos_theta_max=None):
    """Generates an array of random points on the surface of a unit
    sphere.

    Args:
        n_points (int): Number of points to generate.
        cos_theta_min (fl, optional): Lower bound on cos_theta values.
            Defaults to None.
        cos_theta_max (fl, optional): Upper bound on cos_theta values.
            Defaults to None.

    Returns:
        Nx3 arr: Cartesian vector for each point.
    """
    sph = np.empty((n_points, 3), float)

    if cos_theta_min is None or cos_theta_min < -1.0:
        cos_theta_min = -1.0
    if cos_theta_max is None or cos_theta_max > +1.0:
        cos_theta_max = +1.0
    if cos_theta_min > cos_theta_max:
        cos_theta_min, cos_theta_max = cos_theta_max, cos_theta_min

    # z-coordinates
    sph[:, 2] = np.random.uniform(cos_theta_min, cos_theta_max, n_points)
    z2 = np.sqrt(1.0 - sph[:, 2] ** 2)
    phi = (2.0 * np.pi) * np.random.random(n_points)
    sph[:, 0] = z2 * np.cos(phi)  # x
    sph[:, 1] = z2 * np.sin(phi)  # y

    return sph


def uniform_points_on_sphere_surface(n_points):
    """Generates an array of points uniformly distributed on the surface
        of a unit sphere.
        http://bit.ly/2cQClOc
        Courtesy of Marius Cautun.

    Args:
        n_points (int): Number of points to generate.

    Returns:
        Nx3 arr: Cartesian vectors for the generated points.
    """
    dz = 2.0 / n_points
    dPhi = np.pi * (3.0 - 5.0**0.5)
    points = np.empty((n_points, 3), np.float64)

    # Calculate z coordinate
    step = np.arange(n_points)
    points[:, 2] = 1.0 - dz * step

    # Calculate x and y coordinates
    r = (1.0 - points[:, 2] * points[:, 2]) ** 0.5
    phi = dPhi * step
    points[:, 0] = r * np.cos(phi)
    points[:, 1] = r * np.sin(phi)

    return points
