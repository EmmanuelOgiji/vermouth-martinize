# -*- coding: utf-8 -*-
# Copyright 2018 University of Groningen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Tests for the geometry module.
"""

import pytest
import itertools
import numpy as np
from vermouth import geometry


def _generate_test_angles(n_angles):
    """
    Gererate a series of coordinates of 3 points at a different angles.

    Generate 'n_angles' structures with angles between 0 and pi with regular
    angle spacing.

    Parameters
    ----------
    n_angles: int
        Number of structures to generate.

    Yields
    ------
    coordinates: np.ndrarray
        The coordinates of the points. Each row corresponds to a pointm and
        each column corresponds to a dimension.
    angle: float
        The angle between the 3 points.
    """
    # The shift is arbitrary. Its purpose is to avoid having the angle centered
    # at the origin as it may be a special case.
    shift = np.array([1.1, -7.2, 9.1])
    coordinates = np.array([
        [2, 0, 0],
        [0, 0, 0],
        [0, 0, 0],  # Will be redifined
    ]) + shift
    # The distance of the 3rd point to the second one, also the radius in polar
    # coordinates. The value is totally arbitrary, and it should not change the
    # angle; but avoid setting the radius to 1 as it may hide normalization
    # errors.
    radius = 4.2
    for angle in np.linspace(0, np.pi, num=n_angles):
        coordinates[-1, 0] = radius * np.cos(angle) + shift[0]
        coordinates[-1, 1] = radius * np.sin(angle) + shift[1]
        yield coordinates.copy(), angle


def _generate_test_dihedrals(n_angles):
    """
    Generate a series of coordinates for 4 points at different torsion angles.

    Generate 'n_angles' structures with torsion angles between -pi and +pi
    with regular angle spacing.

    Parameters
    ----------
    n_angles: int
        Number of structures to generate.

    Yields
    ------
    coordinates: np.ndrarray
        The coordinates of the points. Each row corresponds to a pointm and
        each column corresponds to a dimension.
    angle: float
        The angle around the ais between the middle points.
    """
    # The shift is arbitrary. Its purpose is to avoid having the angle centered
    # at the origin as it may be a special case.
    shift = np.array([1.1, -7.2, 9.1])
    coordinates = np.array([
        [2, 0, 0],
        [0, 0, 0],
        [0, 0, 5],
        [0, 0, 5],  # Will be redifined
    ]) + shift
    # The distance of the 3rd point to the second one, also the radius in polar
    # coordinates. The value is totally arbitrary, and it should not change the
    # angle; but avoid setting the radius to 1 as it may hide normalization
    # errors.
    radius = 4.2
    for angle in np.linspace(-np.pi, +np.pi, num=n_angles):
        coordinates[-1, 0] = radius * np.cos(angle) + shift[0]
        coordinates[-1, 1] = radius * np.sin(angle) + shift[1]
        yield coordinates.copy(), angle


@pytest.mark.parametrize(
    'points, angle',
    itertools.chain(
        _generate_test_angles(10),
        ((np.array([[0,  3, 0], [0, 0, 0], [0, 6, 0]]), 0), ),
        ((np.array([[0, -9, 0], [0, 0, 0], [0, 2, 0]]), np.pi), ),
    )
)
def test_angle(points, angle):
    vectorBA = points[0, :] - points[1, :]
    vectorBC = points[2, :] - points[1, :]
    assert np.allclose(geometry.angle(vectorBA, vectorBC), angle)


@pytest.mark.parametrize(
    'points, angle',
    itertools.chain(
        _generate_test_dihedrals(10),
        ((np.array([[0,  3, 0], [0, 0, 0], [4, 0, 0], [7, 6, 0]]), 0), ),
        ((np.array([[0, -9, 0], [0, 0, 0], [4, 0, 0], [5, 6, 0]]), np.pi), ),
    )
)
def test_dihedral(points, angle):
    calc_angle = geometry.dihedral(points)
    # +pi and -pi are the same angle; we normalize them to pi
    if np.allclose(calc_angle, -np.pi):
        calc_angle *= -1
    if np.allclose(angle, -np.pi):
        angle *= -1
    assert np.allclose(calc_angle, angle)


@pytest.mark.parametrize(
    'points, angle',
    itertools.chain(
        _generate_test_dihedrals(10),
        ((np.array([[0,  3, 0], [0, 0, 0], [4, 0, 0], [7, 6, 0]]), 0), ),
        ((np.array([[0, -9, 0], [0, 0, 0], [4, 0, 0], [5, 6, 0]]), np.pi), ),
    )
)
def test_dihedral_phase(points, angle):
    angle_phase = angle + np.pi
    if angle_phase > np.pi:
        angle_phase -= 2 * np.pi
    if angle_phase < -np.pi:
        angle_phase += 2 * np.pi
    calc_angle = geometry.dihedral_phase(points)
    # +pi and -pi are the same angle; we normalize them to pi
    if np.allclose(calc_angle, -np.pi):
        calc_angle *= -1
    if np.allclose(angle_phase, -np.pi):
        angle_phase *= -1
    assert np.allclose(calc_angle, angle_phase)
