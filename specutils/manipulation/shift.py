# -*- coding: utf-8 -*-

from astropy import constants
from astropy import units as u


def gravitational_redshift(solar=True, obj_mass=None, obj_radius=None,
                           obj_distance=None):
    """
    Calculates the gravitational redshift z of a distant object
    for an observer at Earth. Also includes the Earth's gravitational blueshift.

    Parameters
    ----------
    solar : bool, optional
        If True, will use the Sun as distant object so that mass, radius,
        and distance are taken from constants and not needed as arguments.
    obj_mass : `Quantity` object (number with mass units)
        Mass of the distant object.
    obj_radius: `Quantity` object (number with length units)
        Radius of the distant object.
    obj_distance: `Quantity` object (number with length units)
        Distance from Earth to the distant object. Can be inf if unknown.

    Returns
    -------
    grav_red : `Quantity` (dimensionless)
        Gravitational redshift z. Convert to wavelength shift by multiplying
        by reference wavelength, or to velocity by multiplying by the
        speed of light.
    """
    if solar:
        obj_mass = constants.M_sun
        obj_radius = constants.R_sun
        obj_distance = constants.au
    else:
        qtype = u.Quantity
        assert all([isinstance(obj_mass, qtype),
                    isinstance(obj_radius, qtype),
                    isinstance(obj_distance, qtype)]), ("If object is not 'sun',"
                      " mass, radius and distance must be given as Quantity.")
        obj_mass = obj_mass.si
        obj_radius = obj_radius.si
        obj_distance = obj_distance.si
    # Schwarzschild radii of object and Earth
    r0 = constants.G * obj_mass / (constants.c ** 2)
    r0_earth = constants.G * constants.M_earth / (constants.c ** 2)
    shift_obj = (1 - r0 / obj_distance) / (1 - r0 / obj_radius)
    shift_earth = (1 - r0_earth / obj_distance) / (1 - r0_earth / constants.R_earth)
    return shift_obj - shift_earth
