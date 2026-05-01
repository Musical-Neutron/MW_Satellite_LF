#!/usr/bin/env python3
import numpy as np


def second_pointings(pointings, alpha_off=2.0 * np.pi / 3.0):
    """Generates a set of pointing vectors with respect to an initial
        set.

    Args:
        pointings (Nx3 arr): Cartesian vectors with respect to which to
            generate the second set of vectors.
        alpha_off (fl, optional): Offset angle (in radians) with respect
            to the intial vectors provided to the function.
            Defaults to 2.*np.pi/3 (120 degrees).

    Returns:
        Nx3 arr: Cartesian vectors.
    """

    # Define new (u, v, w) orthogonal basis, with w aligned along the
    # initial pointing directions (pointings).
    # Define random rotations around w axis
    new_vec_orientation = np.random.rand(len(pointings)) * 2 * np.pi
    new_vecs = np.empty((len(pointings), 3))

    # Populate new_vecs in transformed (u, v, w) basis
    new_vecs[:, 0] = np.sin(new_vec_orientation) * np.sin(alpha_off)
    new_vecs[:, 1] = -np.cos(new_vec_orientation) * np.sin(alpha_off)
    new_vecs[:, 2] = np.cos(alpha_off)

    # Compute transformation matrix to cartesian basis
    mod_pointings = np.linalg.norm(pointings, axis=1)

    w = pointings / mod_pointings[:, np.newaxis]

    u = np.empty((len(pointings), 3), np.float64)
    v = np.empty((len(pointings), 3), np.float64)

    # wz != 0
    temp_w = w[w[:, 2] != 0.0]
    temp_u = u[w[:, 2] != 0.0]

    temp_u[:, 0] = 1.0
    temp_u[:, 1] = 1.0
    temp_u[:, 2] = -(temp_w[:, 0] + temp_w[:, 1]) / temp_w[:, 2]

    u[w[:, 2] != 0.0] = temp_u

    # wz == 0, wy != 0
    temp_w = w[(w[:, 2] == 0.0) & (w[:, 1] != 0.0)]
    temp_u = u[(w[:, 2] == 0.0) & (w[:, 1] != 0.0)]

    temp_u[:, 0] = 1.0
    temp_u[:, 1] = -(temp_w[:, 0] + temp_w[:, 2]) / temp_w[:, 1]
    temp_u[:, 2] = 1.0

    u[(w[:, 2] == 0.0) & (w[:, 1] != 0.0)] = temp_u

    # wz, wy ==0; wx !=0
    temp_w = w[(w[:, 1] == 0.0) & (w[:, 2] == 0.0)]
    temp_u = u[(w[:, 1] == 0.0) & (w[:, 2] == 0.0)]

    temp_u[:, 0] = -(temp_w[:, 1] + temp_w[:, 2]) / temp_w[:, 0]
    temp_u[:, 1] = 1.0
    temp_u[:, 2] = 1.0

    u[(w[:, 1] == 0.0) & (w[:, 2] == 0)] = temp_u

    v[:, 0] = u[:, 1] * w[:, 2] - w[:, 1] * u[:, 2]
    v[:, 1] = u[:, 2] * w[:, 0] - w[:, 2] * u[:, 0]
    v[:, 2] = u[:, 0] * w[:, 1] - w[:, 0] * u[:, 1]

    # Normalise vectors
    mod_u = np.linalg.norm(u, axis=1)
    u /= mod_u[:, np.newaxis]
    mod_v = np.linalg.norm(v, axis=1)
    v /= mod_v[:, np.newaxis]

    # Compute new pointing vectors in Cartesian basis
    new_pointings = np.empty((len(pointings), 3), np.float64)
    new_pointings[:, 0] = (
        new_vecs[:, 0] * u[:, 0] + new_vecs[:, 1] * v[:, 0] + new_vecs[:, 2] * w[:, 0]
    )
    new_pointings[:, 1] = (
        new_vecs[:, 0] * u[:, 1] + new_vecs[:, 1] * v[:, 1] + new_vecs[:, 2] * w[:, 1]
    )
    new_pointings[:, 2] = (
        new_vecs[:, 0] * u[:, 2] + new_vecs[:, 1] * v[:, 2] + new_vecs[:, 2] * w[:, 2]
    )

    # Compute unit vectors of new pointings
    mod_new_pointings = np.linalg.norm(new_pointings, axis=1)
    new_pointings /= mod_new_pointings[:, np.newaxis]

    # Compute original pointing unit vectors
    pointings /= mod_pointings[:, np.newaxis]

    return new_pointings


def survey_Reff(Mv, a, b, survey_area):
    """Computes the effective radius of detection of satellites of a
        given magnitude in a given survey.

    Args:
        Mv (arr): Absolute V-band mangitude values.
        a (fl): a* parameter from eq. 3 in Newton et al. (2018)
        b (fl): b* parameter from eq. 3 in Newton et al. (2018)
        survey_area (fl): Fractional sky area covered by given survey.

    Returns:
        arr: Reff [kpc].
    """
    return 10 ** ((-a * Mv - b)) * 1.0e3


def survey_cone(survey_area):
    """Computes the fractional sky area coverage and cosine of a mock
        conical survey region.

    Args:
        survey_area (fl). Survey sky area coverage [deg^2].

    Returns:
        tuple:
            [0]: Fractional survey area
            [1]: Cosine(survey angle).
    """
    fullSkyArea = 4 * np.pi * (180.0 / np.pi) ** 2  # in deg^2
    fracArea = survey_area / fullSkyArea
    cos_surveyAngle = 1.0 - 2.0 * fracArea

    return (fracArea, cos_surveyAngle)
