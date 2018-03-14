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

import itertools
import numpy as np
import networkx as nx
from .processor import Processor
from .. import selectors

DEFAULT_BOND_TYPE = 6


def self_distance_matrix(coordinates):
    return np.sqrt(((coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]) ** 2).sum(axis=-1))


def compute_decay(distance, shift, rate, power):
    return np.exp(-rate * ((distance - shift) **  power))


def compute_force_constants(distance_matrix, lower_bound, upper_bound,
                            decay_factor, decay_power, base_constant,
                            minimum_force):
    constants = compute_decay(distance_matrix, lower_bound, decay_factor, decay_power)
    np.fill_diagonal(constants, 0)
    constants *= base_constant
    constants[constants < minimum_force] = 0
    constants[distance_matrix > upper_bound] = 0
    return constants


def build_connectivity_matrix(graph, separation, selection=None):
    if separation <= 0:
        raise ValueError('Separation has to be strictly positive.')
    if separation == 1:
        # The connectivity matrix with a separation of 1 is the adjacency
        # matrix. Thanksfully, networkx can directly give it to us a a numpy
        # array.
        return nx.to_numpy_matrix(graph, nodelist=selection).astype(bool)
    subgraph = graph.subgraph(selection)
    connectivity = np.zeros((len(subgraph), len(subgraph)), dtype=bool)
    for (i, key_i),  (j, key_j) in itertools.combinations(enumerate(subgraph.nodes), 2):
        shortest_path = len(nx.shortest_path(subgraph, key_i, key_j))
        # The source and the target are counted in the shortest path
        connectivity[i, j] = shortest_path <= separation + 2
        connectivity[j, i] = connectivity[i, j]
    return connectivity



def apply_rubber_band(molecule, selector,
                      lower_bound, upper_bound,
                      decay_factor, decay_power,
                      base_constant, minimum_force,
                      bond_type):
    selection = []
    coordinates = []
    missing = []
    for node_key, attributes in molecule.nodes.items():
        if selector(attributes):
            selection.append(node_key)
            coordinates.append(attributes.get('position'))
            if coordinates[-1] is None:
                missing.append(node_key)
    if missing:
        raise ValueError('All atoms from the selection must have coordinates. '
                         'The following atoms do not have some: {}.'
                         .format(' '.join(missing)))
    coordinates = np.stack(coordinates)
    distance_matrix = self_distance_matrix(coordinates)
    constants = compute_force_constants(distance_matrix, lower_bound,
                                        upper_bound, decay_factor, decay_power,
                                        base_constant, minimum_force)
    connectivity = build_connectivity_matrix(molecule, 2, selection=selection)
    # Set the force constant to 0 for pairs that are connected. `connectivity`
    # is a matrix of booleans that is True when a pair is connected. Because
    # booleans acts as 0 or 1 in operation, we multiply the force constant
    # matrix by the oposite (OR) of the connectivity matrix.
    constants *= ~connectivity
    distance_matrix = distance_matrix.round(5)  # For compatibility with legacy
    for from_idx, to_idx in zip(*np.triu_indices_from(constants)):
        from_key = selection[from_idx]
        to_key = selection[to_idx]
        force_constant = constants[from_idx, to_idx]
        length = distance_matrix[from_idx, to_idx]
        if force_constant > minimum_force:
            molecule.add_interaction(
                type_='bonds',
                atoms=(from_key, to_key),
                parameters=[bond_type, length, force_constant],
                meta={'group': 'Rubber band'},
            )


class ApplyRubberBand(Processor):
    def __init__(self, lower_bound, upper_bound, decay_factor, decay_power,
                 base_constant, minimum_force,
                 bond_type=None,
                 selector=selectors.select_backbone,
                 bond_type_variable='elastic_network_bond_type'):
        super().__init__()
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.decay_factor = decay_factor
        self.decay_power = decay_power
        self.base_constant = base_constant
        self.minimum_force = minimum_force
        self.bond_type = bond_type
        self.selector = selector
        self.bond_type_variable = bond_type_variable

    def run_molecule(self, molecule):
        # Choose the bond type. From high to low, the priority order is:
        # * what is set as an argument to the processor
        # * what is written in the force field variables
        #   under the key `self.bond_type_variable`
        # * the default value set in DEFAULT_BOND_TYPE
        bond_type = self.bond_type
        if self.bond_type is None:
            bond_type = molecule.force_field.variables.get(self.bond_type_variable,
                                                           DEFAULT_BOND_TYPE)

        apply_rubber_band(molecule, self.selector,
                          lower_bound=self.lower_bound,
                          upper_bound=self.upper_bound,
                          decay_factor=self.decay_factor,
                          decay_power=self.decay_power,
                          base_constant=self.base_constant,
                          minimum_force=self.minimum_force,
                          bond_type=bond_type)
        return molecule
