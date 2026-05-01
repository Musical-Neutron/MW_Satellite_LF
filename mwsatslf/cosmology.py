#!/usr/bin/env python3
import numpy as np


def virial_radius(m200, rho_crit200):
    """Calculates R_200 of a halo given M_200 and rho_crit.

    Args:
        m200 (fl): Mass of halo [Msun / h].
        rho_crit200 (fl): Critical density of the Universe
            [Msun/h (Mpc/h)^-3].

    Returns:
        fl: R_200 [Mpc / h].
    """
    return (m200 / rho_crit200 / (4.0 * np.pi / 3.0)) ** (1.0 / 3.0)


def compute_200rho_crit(m200, r200):
    """Returns the 200 * rho_critical value given M200 and R200.
            Courtesy of Marius Cautun.

    Args:
        m200 (fl) : Mass of halo [Msun/h].
        r200 (fl) : R_200 of halo [Mpc/h].

    Returns:
        fl: 200 * rho_critical.
    """
    return 3.0 / (4.0 * np.pi) * m200 / r200**3
