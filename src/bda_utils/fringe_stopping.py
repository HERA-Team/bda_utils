# -*- coding: utf-8 -*-
# Copyright (c) 2021 The HERA Collaboration
# Licensed under the 2-clause BSD License

import numpy as np
from astropy import coordinates
from astropy import constants
from astropy import units
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, FK5, Angle
import redis
from pyuvdata import utils as uvutils


def compute_phasors(antenna_positions, cofa_location, t_int, t_fs):
    """
    Compute the phasors needed for fringe stopping.

    Parameters
    ----------

    Returns
    -------
    None
    """
    # make sure that antenna_positions is the right size
    antpos_shape = antenna_positions.shape
    if len(antpos_shape) != 2:
        raise ValueError("`antenna_positions` must be a 2d array")
    if antpos_shape[1] != 3:
        raise ValueError("`antenna_positions` must have a shape of (Nants, 3)")

    # unpack cofa location
    (cofa_lat, cofa_lon, cofa_alt) = cofa_location

    cofa_xyz = uvutils.XYZ_from_LatLonAlt(cofa_lat, cofa_lon, cofa_alt)
    cofa_enu = uvutils.ENU_from_ECEF(cofa_xyz, cofa_lat, cofa_lon, cofa_alt)
    antpos_enu = uvutils.ENU_from_ECEF(antenna_positions + cofa_xyz, cofa_lat, cofa_lon, cofa_alt)

    # lay down time grid
    ntimes = int(np.round(t_int / t_fs))
    t_mid = t_int / 2.0
    times = np.zeros(ntimes, dtype=np.float64)
    for i in range(ntimes):
        times[i] = t_mid - t_int + (i + 0.5) * t_fs

    # define telescope location in astropy
    telescope_location = EarthLocation.from_geocentric(*cofa_xyz, unit=units.m)
    phase_frame = "icrs"
    time0 = Time.now()
    time0 = Time(time0.jd, format="jd", location=telescope_location)
    times_jd = time0 + times * units.s

    # define phase center for center of integration
    lst_rad = time0.sidereal_time("apparent").to("rad").value
    itrs_telescope_locations = telescope_location.get_itrs(obstime=times_jd)
    phase_center_coord = SkyCoord(ra=lst_rad, dec=cofa_lat, unit="radian", equinox=time0, frame=FK5)
    icrs_coord = phase_center_coord.transform_to("icrs")

    icrs = coordinates.ICRS()
    frame_telescope_locations = itrs_telescope_locations.transform_to(icrs)
    frame_telescope_locations.representation_type = "cartesian"

    # set up nbl-sized arrays
    nants = antpos_shape[0]
    nbls = int(nants * (nants + 1) / 2)
    antenna_numbers = np.arange(nants)
    ant1_array = np.zeros((nbls,), dtype=np.int32)
    ant2_array = np.zeros((nbls,), dtype=np.int32)
    idx = 0
    for i in range(nants):
        for j in range(i, nants):
            ant1_array[idx] = i
            ant2_array[idx] = j
            idx += 1
    ant1_array_nblt = np.tile(ant1_array, ntimes)
    ant2_array_nblt = np.tile(ant2_array, ntimes)

    antpos_xyz = np.float64(antenna_positions + cofa_xyz)
    frame_phase_center = icrs_coord
    uvw_array = np.zeros((nbls * ntimes, 3), dtype=np.float64)

    tau_astropy = np.zeros((nants, ntimes), dtype=np.float64)

    # compute uvw coordinates
    ant_sort = np.argsort(antenna_numbers)
    for i in range(len(times)):
        i1 = i * nbls
        i2 = (i + 1) * nbls
        inds = np.ix_(np.arange(i1, i2))
        time = times_jd[i]
        frame_telescope_location = frame_telescope_locations[i]
        antpos_itrs = SkyCoord(
            x=antpos_xyz[:, 0] * units.m,
            y=antpos_xyz[:, 1] * units.m,
            z=antpos_xyz[:, 2] * units.m,
            frame="itrs",
            obstime=time,
        )
        frame_ant_coord = antpos_itrs.transform_to(phase_frame)
        frame_ant_rel = (
            (frame_ant_coord.cartesian - frame_telescope_location.cartesian)
            .get_xyz().T.value
        )
        frame_ant_uvw = uvutils.phase_uvw(
            frame_phase_center.ra.rad, frame_phase_center.dec.rad, frame_ant_rel
        )
        ant1_index = np.searchsorted(antenna_numbers[ant_sort], ant1_array_nblt[inds])
        ant2_index = np.searchsorted(antenna_numbers[ant_sort], ant2_array_nblt[inds])
        uvw_array[inds] = (
            frame_ant_uvw[ant_sort][ant2_index, :]
            - frame_ant_uvw[ant_sort][ant1_index, :]
        )
        tau_astropy[:, 1] = frame_ant_uvw[:, 2] / constants.c.to("m/s").value

    delays_astropy = uvw_array[:, 2] / constants.c.to("m/s").value

    return delays_astropy
