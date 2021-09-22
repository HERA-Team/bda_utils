# -*- coding: utf-8 -*-
# Copyright (c) 2021 The HERA Collaboration
# Licensed under the 2-clause BSD License

"""Tests for fringe_stopping.py."""

import pytest
import numpy as np
from astropy import constants
from astropy import units

from bda_utils import fringe_stopping


def test_compute_phasors():
# from a recent-ish data file
    antenna_positions = np.array(
        [
            [  11.36289375,  -77.30163178,  -29.87911092],
            [  41.19076854, -136.26376198,  -19.14741053],
            [  38.00914182, -145.30681401,  -30.17954675],
            [   7.33295749,  -31.90930428,   -7.87568616],
            [  39.15360811, -113.57628395,   -8.13735459],
        ]
    )

    cofa_lat = -30.72152612068925 * np.pi / 180.0
    cofa_lon = 21.42830382686301 * np.pi / 180.0
    cofa_alt = 1051.49

    cofa_location = (cofa_lat, cofa_lon, cofa_alt)

    t_int = 10.0  # seconds
    t_fs = 0.1  # seconds
    ntimes = int(t_int / t_fs)
    nants = antenna_positions.shape[0]
    nbls = int(nants * (nants + 1) / 2)
    nblts = nbls * ntimes

    delays = fringe_stopping.compute_phasors(
        antenna_positions, cofa_location, t_int, t_fs
    )

    assert delays.shape == (nblts,)
