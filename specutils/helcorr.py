# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst

""" Heliocentric velocity corrections """

from __future__ import division, print_function

# This code has been modified from legacy codes. Credits are chronologically
# attributable to: Sergey Koposov, Kochukho, Kudryavtsev, W. Landsman, 
# Chris McCarthy, and Jeff Valenti.

# Third-party
import numpy as np

import astropy.constants as constants
import astropy.coordinates as coord
import astropy.time
import astropy.units as u

__all__ = ["helcorr", "baryvel"]

def baryvel(dje):
    """
    Calculates the heliocentric and barycentric velocity components of Earth.

    Parameters
    ----------
    dje : `~astropy.time.Time` or float
        The Julian ephemeris date.

    Returns
    -------
    dvelh : `~astropy.units.Quantity`
        The heliocentric velocity components in (X, Y, Z) coordinates.

    dvelb : `~astropy.units.Quantity`
        The barycentric velocity components in (X, Y, Z) coordinates.
    """

    if isinstance(dje, astropy.time.Time):
        dje = dje.jd

    # Prepare for the pain.
    dcto = 2415020.
    dcjul = 36525. # Days in Julian year
    dc1 = 1.

    # Constants dcfel(i,k) of fast changing elements.
    dcfel = np.array([
        1.7400353e00, 6.2833195099091e02,  5.2796e-6, 
        6.2565836e00, 6.2830194572674e02, -2.6180e-6,
        4.7199666e00, 8.3997091449254e03, -1.9780e-5,
        1.9636505e-1, 8.4334662911720e03, -5.6044e-5, 
        4.1547339e00, 5.2993466764997e01,  5.8845e-6, 
        4.6524223e00, 2.1354275911213e01,  5.6797e-6, 
        4.2620486e00, 7.5025342197656e00,  5.5317e-6, 
        1.4740694e00, 3.8377331909193e00,  5.6093e-6]).reshape(8, 3)

    # Constants dceps and ccsel(i,k) of slowly changing elements.
    dceps = np.array([4.093198e-1, -2.271110e-4, -2.860401e-8])
    ccsel = np.array([
        1.675104e-2, -4.179579e-5, -1.260516e-7,
        2.220221e-1,  2.809917e-2,  1.852532e-5,
        1.589963e00,  3.418075e-2,  1.430200e-5,
        2.994089e00,  2.590824e-2,  4.155840e-6,  
        8.155457e-1,  2.486352e-2,  6.836840e-6,
        1.735614e00,  1.763719e-2,  6.370440e-6, 
        1.968564e00,  1.524020e-2, -2.517152e-6, 
        1.282417e00,  8.703393e-3,  2.289292e-5,
        2.280820e00,  1.918010e-2,  4.484520e-6, 
        4.833473e-2,  1.641773e-4, -4.654200e-7,
        5.589232e-2, -3.455092e-4, -7.388560e-7, 
        4.634443e-2, -2.658234e-5,  7.757000e-8,  
        8.997041e-3,  6.329728e-6, -1.939256e-9,
        2.284178e-2, -9.941590e-5,  6.787400e-8, 
        4.350267e-2, -6.839749e-5, -2.714956e-7,
        1.348204e-2,  1.091504e-5,  6.903760e-7, 
        3.106570e-2, -1.665665e-4, -1.590188e-7]).reshape(17, 3)

    # Constants of the arguments of the short-period perturbations.
    dcargs = np.array([
        5.0974222e0, -7.8604195454652e2,
        3.9584962e0, -5.7533848094674e2,
        1.6338070e0, -1.1506769618935e3, 
        2.5487111e0, -3.9302097727326e2, 
        4.9255514e0, -5.8849265665348e2, 
        1.3363463e0, -5.5076098609303e2, 
        1.6072053e0, -5.2237501616674e2, 
        1.3629480e0, -1.1790629318198e3, 
        5.5657014e0, -1.0977134971135e3, 
        5.0708205e0, -1.5774000881978e2, 
        3.9318944e0,  5.2963464780000e1, 
        4.8989497e0,  3.9809289073258e1, 
        1.3097446e0,  7.7540959633708e1, 
        3.5147141e0,  7.9618578146517e1, 
        3.5413158e0, -5.4868336758022e2]).reshape(15, 2)

    # Amplitudes ccamps(n,k) of the short-period perturbations.
    ccamps = np.array([
        -2.279594e-5,  1.407414e-5,  8.273188e-6,  1.340565e-5, -2.490817e-7,
        -3.494537e-5,  2.860401e-7,  1.289448e-7,  1.627237e-5, -1.823138e-7,
         6.593466e-7,  1.322572e-5,  9.258695e-6, -4.674248e-7, -3.646275e-7,
         1.140767e-5, -2.049792e-5, -4.747930e-6, -2.638763e-6, -1.245408e-7,
         9.516893e-6, -2.748894e-6, -1.319381e-6, -4.549908e-6, -1.864821e-7, 
         7.310990e-6, -1.924710e-6, -8.772849e-7, -3.334143e-6, -1.745256e-7, 
        -2.603449e-6,  7.359472e-6,  3.168357e-6,  1.119056e-6, -1.655307e-7, 
         3.228859e-6,  1.308997e-7,  1.013137e-7,  2.403899e-6, -3.736225e-7, 
         3.442177e-7,  2.671323e-6,  1.832858e-6, -2.394688e-7, -3.478444e-7, 
         8.702406e-6, -8.421214e-6, -1.372341e-6, -1.455234e-6, -4.998479e-8, 
        -1.488378e-6, -1.251789e-5,  5.226868e-7, -2.049301e-7, 0,
        -8.043059e-6, -2.991300e-6,  1.473654e-7, -3.154542e-7, 0, 
         3.699128e-6, -3.316126e-6,  2.901257e-7,  3.407826e-7, 0, 
         2.550120e-6, -1.241123e-6,  9.901116e-8,  2.210482e-7, 0, 
        -6.351059e-7,  2.341650e-6,  1.061492e-6,  2.878231e-7, 0]).reshape(15, 5)

    # Constants csec3 and ccsec(n,k) of the secular perturbations in longitude.
    ccsec3 = -7.757020e-8
    ccsec = np.array([
        1.289600e-6, 5.550147e-1, 2.076942e00,
        3.102810e-5, 4.035027e00, 3.525565e-1, 
        9.124190e-6, 9.990265e-1, 2.622706e00, 
        9.793240e-7, 5.508259e00, 1.559103e01]).reshape(4, 3)


    # Sidereal rates.
    dcsld = 1.990987e-7                   #sidereal rate in longitude
    ccsgd = 1.990969e-7                   #sidereal rate in mean anomaly

    # Constants used in the calculation of the lunar contribution.
    cckm = 3.122140e-5
    ccmld = 2.661699e-6
    ccfdi = 2.399485e-7

    # Constants dcargm(i,k) of the arguments of the perturbations of the motion
    # of the moon.
    dcargm = np.array([5.1679830e0, 8.3286911095275e3, 5.4913150e0, 
        -7.2140632838100e3, 5.9598530e0, 1.5542754389685e4]).reshape(3, 2)

    # Amplitudes ccampm(n,k) of the perturbations of the moon.
    ccampm = np.array([
         1.097594e-1, 2.896773e-7, 5.450474e-2,  1.438491e-7,
        -2.223581e-2, 5.083103e-8, 1.002548e-2, -2.291823e-8, 
         1.148966e-2, 5.658888e-8, 8.249439e-3,  4.063015e-8]).reshape(3, 4)

    # ccpamv(k) = a*m*dl,dt (planets), dc1mme = 1-mass(earth+moon)
    ccpamv = np.array([8.326827e-11, 1.843484e-11, 1.988712e-12, 1.881276e-12])
    dc1mme = 0.99999696e0

    # Time arguments.
    dt = (dje - dcto) / dcjul
    tvec = np.array([1e0, dt, dt * dt])

    # Values of all elements for the instant(aneous?) dje.
    temp = np.dot(tvec.T, dcfel.T).T % (2 * np.pi)
    dml = temp[0]
    forbel = temp[1:8]
    g = forbel[0]

    deps = (tvec * dceps).sum() % (2 * np.pi)
    sorbel = np.dot(tvec.T, ccsel.T).T % (2 * np.pi)
    e = sorbel[0]

    # Secular perturbations in longitude.
    sn = np.sin(np.dot(tvec[0:2].T, ccsec[:, 1:3].T).T % (2 * np.pi))

    # Periodic perturbations of the Earth-Moon barycenter.
    pertl = (ccsec[:,0] * sn).sum() + dt * ccsec3 * sn[2]
    pertld, pertr, pertrd = 0, 0, 0
    for k in range(0, 15):
        a = (dcargs[k,0] + dt * dcargs[k,1]) % 2 * np.pi
        cosa, sina = np.cos(a), np.sin(a)
        pertl += ccamps[k,0] * cosa + ccamps[k,1] * sina
        pertr += ccamps[k,2] * cosa + ccamps[k,3] * sina
        if k < 11:   
            pertld += (ccamps[k,1] * cosa - ccamps[k,0] * sina) * ccamps[k,4]
            pertrd += (ccamps[k,3] * cosa - ccamps[k,2] * sina) * ccamps[k,4]

    # Elliptic part of the motion of the Earth-Moon barycenter.
    phi = (e * e / 4e0) * (((8e0 / e) - e) * np.sin(g) + 5 * np.sin(2 * g) \
        + (13 / 3.) * e * np.sin(3 * g))
    f = g + phi
    sinf, cosf = np.sin(f), np.cos(f)
    dpsi = (dc1 - e * e) / (dc1 + e * cosf)
    phid = 2 * e * ccsgd * ((1 + 1.5 * e**2) * cosf + e * (1.25 - 0.5 * sinf**2))
    psid = ccsgd * e * sinf / np.sqrt(dc1 - e * e)

    # Perturbed heliocentric motion of the Earth-Moon barycenter.
    d1pdro = dc1 + pertr
    drd = d1pdro * (psid + dpsi * pertrd)
    drld = d1pdro * dpsi * (dcsld + phid + pertld)
    dtl = (dml + phi + pertl) % (2 * np.pi)
    dsinls = np.sin(dtl)
    dcosls = np.cos(dtl)
    dxhd = drd * dcosls - drld * dsinls
    dyhd = drd * dsinls + drld * dcosls

    # Influence of eccentricity, evection and variation on the geocentric
    # motion of the moon.
    pertl, pertld, pertp, pertpd = 0, 0, 0, 0
    for k in range(0, 3):
        a = (dcargm[k,0] + dt * dcargm[k,1]) % (2 * np.pi)
        sina = np.sin(a)
        cosa = np.cos(a)
        pertl += ccampm[k,0] * sina
        pertld += ccampm[k,1] * cosa
        pertp += ccampm[k,2] * cosa
        pertpd -= ccampm[k,3] * sina

    # Heliocentric motion of the Earth.
    tl = forbel[1] + pertl
    sinlm = np.sin(tl)
    coslm = np.cos(tl)
    sigma = cckm / (1.0 + pertp)
    a = sigma * (ccmld + pertld)
    b = sigma * pertpd
    dxhd = dxhd + a * sinlm + b * coslm
    dyhd = dyhd - a * coslm + b * sinlm
    dzhd = -sigma * ccfdi * np.cos(forbel[2])

    # Barycentric motion of the Earth.
    dxbd = dxhd * dc1mme
    dybd = dyhd * dc1mme
    dzbd = dzhd * dc1mme
    for k in range(0, 4):
        plon = forbel[k + 3]
        pomg = sorbel[k + 1]
        pecc = sorbel[k + 9]
        tl = (plon + 2.0 * pecc * np.sin(plon - pomg)) % (2 * np.pi)
        dxbd += ccpamv[k] * (np.sin(tl) + pecc * np.sin(pomg))
        dybd -= ccpamv[k] * (np.cos(tl) + pecc * np.cos(pomg))
        dzbd -= ccpamv[k] * sorbel[k + 13] * np.cos(plon - sorbel[k + 5])

    # Transition to mean equator of date.
    dcosep = np.cos(deps)
    dsinep = np.sin(deps)
    dyahd = dcosep * dyhd - dsinep * dzhd
    dzahd = dsinep * dyhd + dcosep * dzhd
    dyabd = dcosep * dybd - dsinep * dzbd
    dzabd = dsinep * dybd + dcosep * dzbd

    dvelh = constants.au * (np.array([dxhd, dyahd, dzahd])) / u.second
    dvelb = constants.au * (np.array([dxbd, dyabd, dzabd])) / u.second
    return (dvelh, dvelb)


# NOTE:
# We may want to change the syntax input for helcorr so that it accepts a single
# sky coordinate instead of ra/dec.
# Similarly lon/lat/alt/jd could be replaced with a single astropy.units.Time
# class.

def helcorr(lon, lat, alt, ra, dec, mjd):
    """
    Calculate the heliocentric radial velocity corrections for an astronomical 
    source.

    Parameters
    ----------
    lon : `~astropy.coordinates.Longitude` or float
        Earth longitude of the observatory (western direction is positive). Can
        be anything that initialises an `~astropy.coordinates.Angle` object
        (if float, in degrees).
    lat : `~astropy.coordinates.Latitude` or float
        Earth latitude of observatory. Can be anything that initialises an
        `~astropy.coordinates.Latitude` object (if float, in degrees).
    alt : `~astropy.units.Quantity` or float
        Altitude of the observatory (if float, in meters).
    ra : `~astropy.coordinates.Angle` or float
        Right ascension of the object for epoch J2000 (if float, in degrees).
    dec : `~astropy.coordinates.Angle` or float
        Declination of the object for epoch J2000 (if float, in degrees).
    mjd : float
        The modified Julian date for the middle of exposure.

    Returns
    -------
    helcorr : `~astropy.units.Quantity`
        The heliocentric velocity correction.
    """

    if not isinstance(lon, coord.Longitude):
        lon = coord.Longitude(lon * u.deg)

    if not isinstance(lat, coord.Latitude):
        lat = coord.Latitude(lat * u.deg)

    if not isinstance(alt, u.Quantity):
        alt *= u.m

    if not isinstance(ra, u.Quantity):
        ra *= u.deg

    if not isinstance(dec, u.Quantity):
        dec *= u.deg

    # Here we specify the location so that we can easily calculate the mean
    # local siderial time later on
    time = astropy.time.Time(2.4e6 + mjd, format="jd", location=(lon, lat, alt))
    epoch = time.datetime.year + time.datetime.month/12. \
        + time.datetime.day/365.

    # Precess the coordinates to the current epoch
    coordinate = coord.SkyCoord(ra, dec, frame="fk5").transform_to(
        coord.FK5(equinox="J{}".format(epoch)))

    # Convert geodetic latitude into geocentric latitude to correct for rotation
    # of the Earth
    dlat = ((-11. * 60. + 32.743) * np.sin(2 * lat) + 1.1633 * np.sin(4 * lat) \
        - 0.0026 * np.sin(6 * lat)) * u.degree
    geocentric_lat = lat + dlat / 3600.

    # Calculate distance of observer from Earth center
    r = alt + 6378160.0 * u.m * (0.998327073 \
        + 0.001676438 * np.cos(2 * geocentric_lat) \
        - 0.000003510 * np.cos(4 * geocentric_lat) \
        + 0.000000008 * np.cos(6 * geocentric_lat))

    # Calculate rotational velocity perpendicular to the radius vector
    # Note: 23.934469591229 is the siderial day in hours for 1986
    v = 2 * np.pi * r / (23.934469591229 * 3600 * u.second)

    # Calculate vdiurnal velocity
    vdiurnal = v * np.cos(lat) * np.cos(coordinate.dec) \
      * np.sin(coordinate.ra - time.sidereal_time("mean"))

    # Calculate baricentric and heliocentric velocities
    vh, vb = baryvel(time)

    # Project along the line of sight
    projection = np.array([
        np.cos(coordinate.dec) * np.cos(coordinate.ra),
        np.cos(coordinate.dec) * np.sin(coordinate.ra),
        np.sin(coordinate.dec)])
    vbar = (vb * projection).sum()
    vhel = (vh * projection).sum()

    # Using baricentric velocity for correction
    correction = vdiurnal + vbar

    # [TODO] it may be useful to return other components of velocity or extra
    # information about the transforms (e.g., gmst, ut, lmst, dlat, lat, vbar,
    # vhel, etc)
    return correction
