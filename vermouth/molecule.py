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
Created on Thu Sep 14 10:58:04 2017

@author: Peter Kroon
"""

from collections import defaultdict, OrderedDict, namedtuple
import copy
from functools import partial

import networkx as nx
import numpy as np

from . import graph_utils
from . import geometry


Interaction = namedtuple('Interaction', 'atoms parameters meta')
DeleteInteraction = namedtuple('DeleteInteraction',
                               'atoms atom_attrs parameters meta')


class LinkPredicate:
    def __init__(self, value):
        self.value = value

    def match(self, node, key):
        raise NotImplementedError

    def __repr__(self):
        return '<{} at {:x} value={}>'.format(self.__class__.__name__, id(self), self.value)


class Choice(LinkPredicate):
    def match(self, node, key):
        return node.get(key) in self.value


class NotDefinedOrNot(LinkPredicate):
    def match(self, node, key):
        return key not in node or node[key] != self.value


class LinkParameterEffector:
    n_keys_asked = None

    def __init__(self, keys, format=None):
        self.keys = keys
        if self.n_keys_asked is not None and len(self.keys) != self.n_keys_asked:
            raise ValueError(
                'Unexpected number of keys provided in {}: '
                '{} were expected, but {} were rovided.'
                .format(self.__class__.name, self.n_keys_asked, len(keys))
            )
        self.format = format

    def __call__(self, molecule, match):
        keys = [match[key] for key in self.keys]
        result = self.apply(molecule, keys)
        if self.format is not None:
            result = '{value:{format}}'.format(value=result, format=self.format)
        return result

    def apply(self, molecule, keys):
        msg = 'The method need to be implemented by the children class.'
        raise NotImplementedError(msg)


class ParamDistance(LinkParameterEffector):
    n_keys_asked = 2

    def apply(self, molecule, keys):
        # This will raise a ValueError if an atom is missing, or if an
        # atom does not have position.
        positions = np.stack([molecule.nodes[key]['position'] for key in keys])
        # We assume there are two rows; which we can since we checked earlier
        # that exactly two atom keys were passed.
        distance = np.sqrt(np.sum(np.diff(positions, axis=0)**2))
        return distance


class ParamAngle(LinkParameterEffector):
    n_keys_asked = 3

    def apply(self, molecule, keys):
        # This will raise a ValueError if an atom is missing, or if an
        # atom does not have position.
        positions = np.stack([molecule.nodes[key]['position'] for key in keys])
        vectorBA = positions[0, :] - positions[1, :]
        vectorBC = positions[2, :] - positions[1, :]
        angle = geometry.angle(vectorBA, vectorBC)
        return np.degrees(angle)


class ParamDihedral(LinkParameterEffector):
    n_keys_asked = 4

    def apply(self, molecule, keys):
        # This will raise a ValueError if an atom is missing, or if an
        # atom does not have position.
        positions = np.stack([molecule.nodes[key]['position'] for key in keys])
        vectorAB = positions[1, :] - positions[0, :]
        vectorBC = positions[2, :] - positions[1, :]
        vectorCD = positions[3, :] - positions[2, :]
        angle = geometry.dihedral(vectorAB, vectorBC, vectorCD)
        return np.degrees(angle)


class ParamDihedralLeft(LinkParameterEffector):
    n_keys_asked = 4

    def apply(self, molecule, keys):
        # This will raise a ValueError if an atom is missing, or if an
        # atom does not have position.
        positions = np.stack([molecule.nodes[key]['position'] for key in keys])
        vectorAB = positions[1, :] - positions[0, :]
        vectorBC = positions[2, :] - positions[1, :]
        vectorCD = positions[3, :] - positions[2, :]
        angle = geometry.dihedral_left(vectorAB, vectorBC, vectorCD)
        return np.degrees(angle)


class Molecule(nx.Graph):
    # As the particles are stored as nodes, we want the nodes to stay
    # ordered.
    node_dict_factory = OrderedDict

    def __init__(self, *args, **kwargs):
        self.meta = kwargs.pop('meta', {})
        self._force_field = kwargs.pop('force_field', None)
        super().__init__(*args, **kwargs)
        self.interactions = defaultdict(list)
        self.nrexcl = None

    @property
    def force_field(self):
        """
        The force field the molecule is described for.

        The force field is assumed to be consistent for all the molecules of
        a system. While it is possible to reassign
        :attr:`Molecule._force_field`, it is recommended to assign the force
        field at the system level as reassigning :attr:`System.force_field`
        will propagate the change to all the molecules in that system.
        """
        return self._force_field

    @property
    def atoms(self):
        for node in self.nodes():
            node_attr = self.node[node]
            yield node, node_attr

    def copy(self, as_view=False):
        copy = super().copy(as_view)
        if not as_view:
            copy = self.__class__(copy)
        copy._force_field = self.force_field
        copy.meta = self.meta.copy()
        return copy

    def subgraph(self, *args, **kwargs):
        return self.__class__(super().subgraph(*args, **kwargs))

    def add_interaction(self, type_, atoms, parameters, meta=None):
        if meta is None:
            meta = {}
        for atom in atoms:
            if atom not in self:
                # KeyError?
                raise ValueError('Unknown atom {}'.format(atom))
        self.interactions[type_].append(
            Interaction(atoms=tuple(atoms), parameters=parameters, meta=meta)
        )

    def add_or_replace_interaction(self, type_, atoms, parameters, meta=None):
        if meta is None:
            meta = {}
        for idx, interaction in enumerate(self.interactions[type_]):
            if (interaction.atoms == tuple(atoms)
                    and interaction.meta.get('version', 0) == meta.get('version', 0)):
                new_interaction = Interaction(
                    atoms=tuple(atoms), parameters=parameters, meta=meta,
                )
                self.interactions[type_][idx] = new_interaction
                break
        else:  # no break
            self.add_interaction(type_, atoms, parameters, meta)

    def get_interaction(self, type_):
        return self.interactions[type_]

    def remove_interaction(self, type_, atoms, version=0):
        for idx, interaction in enumerate(self.interactions[type_]):
            if interaction.atoms == atoms and interaction.meta.get('version', 0):
                break
        else:  # no break
            msg = ("Can't find interaction of type {} between atoms {} "
                   "and with version {}")
            raise KeyError(msg.format(type_, atoms, version))
        del self.interactions[type_][idx]

    def remove_matching_interaction(self, type_, template_interaction):
        for idx, interaction in enumerate(self.interactions[type_]):
            if interaction_match(self, interaction, template_interaction):
                del self.interactions[type_][idx]
                break
        else:  # no break
            raise ValueError('Cannot find a matching interaction.')

    def find_atoms(self, **attrs):
        for node_idx in self:
            node = self.nodes[node_idx]
            if all(node.get(attr, None) == val for attr, val in attrs.items()):
                yield node_idx

    def __getattr__(self, name):
        # TODO: DRY
        if name.startswith('get_') and name.endswith('s'):
            type_ = name[len('get_'):-len('s')]
            return partial(self.get_interaction, type_)
        elif name.startswith('add_'):
            type_ = name[len('add_'):]
            return partial(self.add_interaction, type_)
        elif name.startswith('remove_'):
            type_ = name[len('remove_'):]
            return partial(self.remove_interaction, type_)
        else:
            raise AttributeError

    def merge_molecule(self, molecule):
        """
        Add the atoms and the interactions of a molecule at the end of this one.

        Atom and residue index of the new atoms are offset to follow the last
        atom of this molecule.
        """
        if self.force_field != molecule.force_field:
            raise ValueError(
                'Cannot merge molecules with different force fields.'
            )
        if self.nrexcl != molecule.nrexcl:
            raise ValueError(
                'Cannot merge molecules with different nrexcl. '
                'This molecule has nrexcl={}, while the other has nrexcl={}.'
                .format(self.nrexcl, molecule.nrexcl)
            )
        if len(self.nodes()):
            # We assume that the last id is always the largest.
            last_node_idx = max(self) 
            offset = last_node_idx + 1
            residue_offset = self.node[last_node_idx]['resid'] + 1
            offset_charge_group = self.node[last_node_idx].get('charge_group', -1) + 1
        else:
            offset = 0
            residue_offset = 0
            offset_charge_group = 0

        correspondence = {}
        for idx, node in enumerate(molecule.nodes(), start=offset):
            correspondence[node] = idx
            new_atom = copy.copy(molecule.node[node])
            new_atom['resid'] += residue_offset
            new_atom['charge_group'] = (new_atom.get('charge_group', 0)
                                        + offset_charge_group)
            self.add_node(idx, **new_atom)

        for name, interactions in molecule.interactions.items():
            for interaction in interactions:
                atoms = tuple(correspondence[atom] for atom in interaction.atoms)
                self.add_interaction(name, atoms, interaction.parameters, interaction.meta)

        for edge in molecule.edges:
            self.add_edge(*(correspondence[node] for node in edge))

    def share_moltype_with(self, other):
        # TODO: Test the node attributes, the molecule attributes, and
        # the interactions.
        return nx.is_isomorphic(self, other)

    def iter_residues(self):
        residue_graph = graph_utils.make_residue_graph(self)
        return (tuple(residue_graph.nodes[res]['graph'].nodes) for res in residue_graph.nodes)


class Block(nx.Graph):
    """
    Residue topology template

    Attributes
    ----------
    name: str or None
        The name of the residue. Set to `None` if undefined.
    atoms: iterator of dict
        The atoms in the residue. Each atom is a dict with *a minima* a key
        'name' for the name of the atom, and a key 'atype' for the atom type.
        An atom can also have a key 'charge', 'charge_group', 'comment', or any
        arbitrary key. 
    interactions: dict
        All the known interactions. Each item of the dictionary is a type of
        interaction, with the key being the name of the kind of interaction
        using Gromacs itp/rtp conventions ('bonds', 'angles', ...) and the
        value being a list of all the interactions of that type in the residue.
        An interaction is a dict with a key 'atoms' under which is stored the
        list of the atoms involved (referred by their name), a key 'parameters'
        under which is stored an arbitrary list of non-atom parameters as
        written in a RTP file, and arbitrary keys to store custom metadata. A
        given interaction can have a comment under the key 'comment'.
    """
    # As the particles are stored as nodes, we want the nodes to stay
    # ordered.
    node_dict_factory = OrderedDict

    def __init__(self):
        super(Block, self).__init__(self)
        self.name = None
        self.interactions = {}
        self._apply_to_all_interactions = defaultdict(dict)

    def __repr__(self):
        name = self.name
        if name is None:
            name = 'Unnamed'
        return '<{} "{}" at 0x{:x}>'.format(self.__class__.__name__,
                                          name, id(self))

    def add_atom(self, atom):
        try:
            name = atom['atomname']
        except KeyError:
            raise ValueError('Atom has no atomname: "{}".'.format(atom))
        self.add_node(name, **atom)

    @property
    def atoms(self):
        for node in self.nodes():
            node_attr = self.node[node]
            # In pre-blocks, some nodes correspond to particles in neighboring
            # residues. These node do not carry particle information and should
            # not appear as particles.
            if node_attr:
                yield node_attr

    def make_edges_from_interaction_type(self, type_):
        """
        Create edges from the interactions of a given type.

        The interactions must be described so that two consecutive atoms in an
        interaction should be linked by an edge. This is the case for bonds,
        angles, proper dihedral angles, and cmap torsions. It is not always
        true for improper torsions.

        Cmap are described as two consecutive proper dihedral angles. The
        atoms for the interaction are the 4 atoms of the first dihedral angle
        followed by the next atom forming the second dihedral angle with the
        3 previous ones. Each pair of consecutive atoms generate an edge.

        .. warning::

            If there is no interaction of the required type, it will be
            silently ignored.

        Parameters
        ----------
        type_: str
            The name of the interaction type the edges should be built from.
        """
        for interaction in self.interactions.get(type_, []):
            if interaction.meta.get('edge', True):
                atoms = interaction.atoms
                self.add_edges_from(zip(atoms[:-1], atoms[1:]))

    def make_edges_from_interactions(self):
        """
        Create edges from the interactions we know how to convert to edges.

        The known interactions are bonds, angles, proper dihedral angles, and
        cmap torsions.
        """
        known_types = ('bonds', 'angles', 'dihedrals', 'cmap', 'constraints')
        for type_ in known_types:
            self.make_edges_from_interaction_type(type_)

    def guess_angles(self):
        for a in self.nodes():
            for b in self.neighbors(a):
                for c in self.neighbors(b):
                    if c == a:
                        continue
                    yield (a, b, c)

    def guess_dihedrals(self, angles=None):
        angles = angles if angles is not None else self.guess_angles()
        for a, b, c in angles:
            for d in self.neighbors(c):
                if d not in (a, b):
                    yield (a, b, c, d)

    def has_dihedral_around(self, center):
        """
        Returns True if the block has a dihedral centered around the given bond.

        Parameters
        ----------
        center: tuple
            The name of the two central atoms of the dihedral angle. The
            method is sensitive to the order.

        Returns
        -------
        bool
        """
        all_centers = [tuple(dih['atoms'][1:-1])
                       for dih in self.interactions.get('dihedrals', [])]
        return tuple(center) in all_centers

    def has_improper_around(self, center):
        """
        Returns True if the block has an improper centered around the given bond.

        Parameters
        ----------
        center: tuple
            The name of the two central atoms of the improper torsion. The
            method is sensitive to the order.

        Returns
        -------
        bool
        """
        all_centers = [tuple(dih.atoms[1:-1])
                       for dih in self.interactions.get('impropers', [])]
        return tuple(center) in all_centers

    def to_molecule(self, atom_offset=0, resid=1, offset_charge_group=1):
        name_to_idx = {}
        mol = Molecule()
        for idx, atom in enumerate(self.atoms, start=atom_offset):
            name_to_idx[atom['atomname']] = idx
            new_atom = copy.copy(atom)
            new_atom['resid'] = resid
            new_atom['resname'] = self.name
            new_atom['charge_group'] = (new_atom.get('charge_group', 0)
                                        + offset_charge_group)
            mol.add_node(idx, **new_atom)
        for name, interactions in self.interactions.items():
            for interaction in interactions:
                atoms = tuple(
                    name_to_idx[atom] for atom in interaction.atoms
                )
                mol.add_interaction(
                    name, atoms,
                    interaction.parameters
                )
        for edge in self.edges:
            mol.add_edge(*(name_to_idx[node] for node in edge))

        try:
            mol.nrexcl = self.nrexcl
        except AttributeError:
            pass

        return mol


class Link(Block):
    """
    Template link between two residues.
    """
    node_dict_factory = OrderedDict

    def __init__(self):
        super().__init__()
        self.non_edges = []
        self.removed_interactions = {}
        self._apply_to_all_nodes = {}
        self.molecule_meta = {}
        self.patterns = []


def attributes_match(attributes, template_attributes, ignore_keys=()):
    for attr, value in template_attributes.items():
        if attr in ignore_keys:
            continue
        if isinstance(value, LinkPredicate):
            match = value.match(attributes, attr)
            if not value.match(attributes, attr):
                return False
        elif attributes.get(attr) != value:
            return False
    return True


def interaction_match(molecule, interaction, template_interaction):
    atoms_match = tuple(template_interaction.atoms) == tuple(interaction.atoms)
    parameters_match = (
        not template_interaction.parameters
        or tuple(template_interaction.parameters) == tuple(interaction.parameters)
    )
    if atoms_match and parameters_match:
        try:
            atom_attrs = template_interaction.atom_attrs
        except AttributeError:
            atom_attrs = [{}, ] * len(template_interaction.atoms)
        nodes = [molecule.nodes[atom] for atom in interaction.atoms]
        for atom, template_atom in zip(nodes, atom_attrs):
            if not attributes_match(atom, template_atom):
                return False
        return attributes_match(interaction.meta, template_interaction.meta)
    return False


if __name__ == '__main__':
    mol = Molecule()
    mol.add_edge(0, 1)
    mol.add_edge(1, 2)
    nx.subgraph(mol, (0, 1))

    mol.add_interaction('bond', (0, 1), tuple((1, 2)))
    mol.add_interaction('bond', (1, 2), tuple((10, 20)))
    mol.add_angle((0, 1, 2), tuple([10, 2, 3]))

    print(mol.interactions)
    print(mol.get_interaction('bond'))
    print(mol.get_bonds())
    print(mol.get_angles())

    mol.remove_interaction('bond', (0, 3))
    print(mol.get_bonds())
