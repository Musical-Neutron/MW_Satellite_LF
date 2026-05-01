#!/usr/bin/env python3
import numpy as np


def combined_satellite_estimate(
    sub_pos,
    sub_dis,
    sub_index,
    sdss_direction,
    sdss_cos_survey_angle,
    sdss_rmax,
    des_direction,
    des_cos_survey_angle,
    des_rmax,
    mv,
    classical_sat_count,
    mv_bins,
):
    """Creates a luminosity function estimate for a combined observation
        of two mock surveys.

    Args:
        sub_pos (Nx3 arr): Cartesian subhalo positions [kpc].
        sub_dis (Nx1 arr): Subhalo distances from halo centre [kpc].
        sub_index (Nx1 arr): Indexes that order the subhalo list.
        sdss_direction (1x3 arr): Vector specifying the direction of the
            SDSS mock survey cone.
        sdss_cos_survey_angle (fl): Cosine of the SDSS survey opening
            angle.
        sdss_rmax (1xM arr): The SDSS Rmax value calculated for the
            associated M_V.
        des_direction (1x3 arr): Vector specifying the direction of the
            DES mock survey cone.
        des_cos_survey_angle (fl): Cosine of the DES survey opening
            angle.
        des_rmax (1xM arr): The DES Rmax value calculated for the
            associated M_V.
        mv (1xM arr): Absolute V-band magnitudes to be used in the
            estimate.
        classical_sat_count: Output from 'update classical satellites'.
        mv_bins (arr): Magnitude bins for estimate.

    Returns:
        arr: Satellite galaxy luminosity function.
    """
    # Find all subhaloes inside the mock surveys
    sdss_cos_theta = (sub_pos * sdss_direction).sum(axis=1) / sub_dis
    des_cos_theta = (sub_pos * des_direction).sum(axis=1) / sub_dis
    # Subhaloes inside the SDSS mock survey volume
    s_sdss = sdss_cos_theta >= sdss_cos_survey_angle
    # Subhaloes inside the DES mock survey volume
    s_des = des_cos_theta >= des_cos_survey_angle
    # Select subhaloes that are in either the first or second mocks
    select = s_sdss + s_des

    # Construct arrays of distances and indices for all subhaloes inside
    # the survey volumes.
    _dis = sub_dis[select]
    _index = sub_index[select]
    _is_sdss = s_sdss[select]  # True if subhalo is inside the SDSS mock

    # Temporary array to store the total number of subhaloes at each
    # observed satellite magnitude.
    N_tot = np.zeros(mv.shape[0], np.float32)
    """ Loop over each observed satellite. The first and last entries of
    the array are not actual observations so skip them. """
    for j in range(1, mv.shape[0]):
        # Select subhaloes inside SDSS and DES, respectively.
        s_sdss = _is_sdss * (_dis <= sdss_rmax[j])
        s_des = ~_is_sdss * (_dis <= des_rmax[j])
        select = s_sdss + s_des
        if select.sum() < 2:
            # Check that we have at least two subhaloes for given satellite
            # (except for the last value of j)
            N_tot[j] = np.nan
        elif j < mv.shape[0] - 1:
            # Typical case, where we have 1 observation
            index_min, index_max = _index[select][[0, 1]]
            N_tot[j] = np.random.randint(index_min + 1, index_max + 1, 1)

            # Only keep subhaloes NOT used for the current estimate
            s2 = _index >= N_tot[j]
            # Update the _dis and _index array to reflect that objects were
            # removed.
            _dis = _dis[s2]
            _index = _index[s2]
            _is_sdss = _is_sdss[s2]
        else:
            """The last value of j, corresponding to the left-most point on
            the graph. No observation so we make an estimate of how many
            satellites could be there."""
            index_max = _index[select][0]
            N_tot[j] = np.random.randint(N_tot[j - 1], index_max + 1, 1)
            # This is the last iteration, so no need to get rid of any
            # subhaloes as the arrays won't be used again.

    # Add the classical satellites
    N_tot += classical_sat_count

    # Interpolate the satellite count on a regular grid in MV
    return 10 ** np.interp(mv_bins, mv, np.log10(N_tot))


def indiv_satellite_estimate(
    sub_pos,
    sub_dis,
    sub_index,
    direction,
    cos_survey_angle,
    rmax,
    mv,
    classical_sat_count,
    mv_bins,
):
    """Creates a luminosity function estimate for a single mock survey.

    Args:
        sub_pos (Nx3 arr): Cartesian subhalo positions [kpc].
        sub_dis (Nx1 arr): Subhalo distances from halo centre [kpc].
        sub_index (Nx1 arr): Indexes that order the subhalo list.
        direction (1x3 arr): Vector specifying the direction of the mock
            survey cone.
        cos_survey_angle (fl): Cosine of the survey opening angle.
        rmax (1xM arr): The Rmax value calculated for the associated
            M_V.
        mv (1xM arr): Absolute V-band magnitudes to be used in the
            estimate.
        classical_sat_count: Output from 'update classical satellites'.
        mv_bins (arr): Magnitude bins for estimate.

    Returns:
        arr: Satellite galaxy luminosity function.
    """
    # Find all subhaloes inside the mock survey
    cos_theta = (sub_pos * direction).sum(axis=1) / sub_dis
    select = cos_theta >= cos_survey_angle

    # Construct arrays of distances and indices for all subhaloes inside
    # the survey.
    _dis = sub_dis[select]
    _index = sub_index[select]

    # Temporary array to store the total number of subhaloes at each
    # observed satellite magnitude.
    N_tot = np.zeros(mv.shape[0], np.float32)
    """ Loop over each observed satellite. The first and last entries of
    the array are not actual observations so skip them. """
    for j in range(1, mv.shape[0]):
        # Select subhaloes inside survey volume.
        select = _dis <= rmax[j]
        if select.sum() < 2:
            # Check that we have at least two subhaloes for given satellite
            # (except for the last value of j)
            N_tot[j] = np.nan
        elif j < mv.shape[0] - 1:
            # Typical case, where we have 1 observation
            index_min, index_max = _index[select][[0, 1]]
            N_tot[j] = np.random.randint(index_min + 1, index_max + 1, 1)

            # Only keep subhaloes NOT used for the current estimate
            s2 = _index >= N_tot[j]
            # Update the _dis and _index array to reflect that objects were
            # removed.
            _dis = _dis[s2]
            _index = _index[s2]
        else:
            """The last value of j, corresponding to the left-most point on
            the graph. No observation so we make an estimate of how many
            satellites could be there."""
            index_max = _index[select][0]
            N_tot[j] = np.random.randint(N_tot[j - 1], index_max + 1, 1)
            # This is the last iteration, so no need to get rid of any
            # subhaloes as the arrays won't be used again.

    # Add the classical satellites
    N_tot += classical_sat_count

    # Interpolate the satellite count on a regular grid in MV
    return 10 ** np.interp(mv_bins, mv, np.log10(N_tot))


def update_classical_satellites(mv_keep, mv_lmc, classical_sats):
    """Updates the count of classical satellites given the selection of
        satellites that are associated to the LMC.

    Args:
        mv_keep (arr): Absolute V-band mangnitudes of the satellites to
            use in the procedure.
        mv_lmc (arr): Absolute V-band magnitudes of the satellites to be
            removed from the procedure.
        classical_sats (arr): Classical satellites to add in.

    Returns:
        1xN arr: Cumulative number of satellites that have been treated
        'classically' in this analysis.
    """
    out = np.zeros(mv_keep.shape[0], np.int32)
    for i in range(mv_keep.shape[0]):
        out[i] = (mv_keep[i] >= mv_lmc).sum() + (mv_keep[i] >= classical_sats).sum()
    return out
