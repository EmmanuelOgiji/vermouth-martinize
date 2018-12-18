#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 14:54:23 2018

@author: peterkroon
"""

from .log_helpers import StyleAdapter, get_logger
from .forcefield import FORCE_FIELDS
from .ffinput import _tokenize, _parse_atom_attributes
from .graph_utils import MappingGraphMatcher

from collections import defaultdict
from functools import partial

import networkx as nx

LOGGER = StyleAdapter(get_logger(__name__))


def get_block(ff_name, type, resname):
    return getattr(FORCE_FIELDS[ff_name], type+'s')[resname]


class Mapping:
    def __init__(self, block_from, block_to, mapping, references,
                 ff_from=None, ff_to=None, extra=(), normalize_weights=False,
                 type='block', names=tuple()):
        self.blocks_from = block_from
        self.blocks_to = block_to
        self.blocks_to.extra = extra
        self.references = references
        self.ff_from = ff_from
        self.ff_to = ff_to
        self.type = type
        self.names = names
        self.mapping = mapping
        # Remove nodes not mapped from blocks_from
        unmapped = set(self.blocks_from.nodes.keys()) - set(self.mapping.keys())
        self.blocks_from.remove_nodes_from(unmapped)

        # Normalize the weights
        if normalize_weights:
            self._normalize_weights()

    @property
    def reverse_mapping(self):
        rev_mapping = defaultdict(dict)  # {idx to: {idx from: weight}}
        for idx_from in self.mapping:
            for idx_to, weight in self.mapping[idx_from].items():
                rev_mapping[idx_to][idx_from] = weight
        return dict(rev_mapping)

    def map(self, graph, node_match=None, edge_match=None):
        if node_match is None:
            def node_match(node1, node2):
                return True

        if edge_match is None:
            def edge_match(node11, node12, node21, node22):
                return True
        else:
            edge_match = partial(edge_match, graph, self.blocks_from)

        return self._graph_map(graph, node_match, edge_match)

    def _graph_map(self, graph, node_match, edge_match):
        # 1 Find subgraph isomorphism between blocks_from and graph
        # 2 Translate found match ({graph idx: blocks from idx}) to indices in
        #   blocks_to using self.mapping
        # 3 Return found matches and blocks_to?
        # PS. I don't really like this, because this object is becoming too
        # intelligent by also having to do the isomorphism. On the other hand,
        # it makes sense from the maths point of view.
        matcher = MappingGraphMatcher(graph, self.blocks_from,
                                      node_match=node_match,
                                      edge_match=edge_match)
        for match in matcher.subgraph_isomorphisms_iter():
            new_match = defaultdict(dict)
            for graph_idx, from_idx in match.items():
                new_match[graph_idx].update(self.mapping[from_idx])
            yield new_match, self.blocks_to

    def _normalize_weights(self):
        # Normalize weights such that the the sum of the weights of nodes
        # mapping to something is one
        rev_mapping = self.reverse_mapping
        for idx_from in self.mapping:
            for idx_to, weight in self.mapping[idx_from].items():
                self.mapping[idx_from][idx_to] = weight/sum(rev_mapping[idx_to].values())


class MappingBuilder:
    def __init__(self):
        self.reset()

    def reset(self):
        self.mapping = defaultdict(dict)
        self.blocks_from = None
        self.blocks_to = None
        self.ff_from = None
        self.ff_to = None
        self.names = []
        self.references = {}

    def to_ff(self, ff_name):
        self.ff_to = ff_name

    def from_ff(self, ff_name):
        self.ff_from = ff_name

    def add_block_from(self, block):
        self.names.append(block.name)
        block = block.copy()
        for node in block.nodes.values():
            if 'replace' in node:
                node.update(node['replace'])
                del node['replace']
        if self.blocks_from is None:
            self.blocks_from = nx.Graph(block)
            self.blocks_from = block.to_molecule(force_field=FORCE_FIELDS[self.ff_from])
        else:
            block = block.to_molecule(force_field=FORCE_FIELDS[self.ff_from])
            self.blocks_from.merge_molecule(block)

    def add_block_to(self, block):
        block = block.copy()
        for node in block.nodes.values():
            if 'replace' in node:
                node.update(node['replace'])
                del node['replace']

        if self.blocks_to is None:
            self.blocks_to = block.to_molecule(force_field=FORCE_FIELDS[self.ff_to])
        else:
            block = block.to_molecule(force_field=FORCE_FIELDS[self.ff_to])
            self.blocks_to.merge_molecule(block)

    def add_node_from(self, attrs):
        idx = max(self.blocks_from.nodes) + 1
        self.blocks_from.add_node(idx, **attrs)
        print(self.blocks_from.nodes(data=True))
    
    def add_node_to(self, attrs):
        idx = max(self.blocks_to.nodes) + 1
        self.blocks_to.add_node(idx, **attrs)

    def add_edge_from(self, attrs1, attrs2):
        nodes1 = list(self.blocks_from.find_atoms(**attrs1))
        nodes2 = list(self.blocks_from.find_atoms(**attrs2))
        print(attrs1, nodes1)
        print(attrs2, nodes2)
        print(self.blocks_from.nodes(data=True))
        assert len(nodes1) == len(nodes2) == 1
        assert nodes1 != nodes2
        self.blocks_from.add_edge(nodes1[0], nodes2[0])

    def add_edge_to(self, attrs1, attrs2):
        nodes1 = list(self.blocks_to.find_atoms(**attrs1))
        nodes2 = list(self.blocks_to.find_atoms(**attrs2))
        assert len(nodes1) == len(nodes2) == 1
        assert nodes1 != nodes2
        self.blocks_to.add_edge(nodes1[0], nodes2[0])

    def add_mapping(self, attrs_from, attrs_to, weight):
        nodes_from = list(self.blocks_from.find_atoms(**attrs_from))
        nodes_to = list(self.blocks_to.find_atoms(**attrs_to))
        assert len(nodes_from) == len(nodes_to) == 1
        self.mapping[nodes_from[0]][nodes_to[0]] = weight

    def add_reference(self, attrs_to, attrs_from):
        nodes_to = list(self.blocks_to.find_atoms(**attrs_to))
        assert len(nodes_to) == 1
        node_to = nodes_to[0]
        nodes_from = set(self.blocks_from.find_atoms(**attrs_from))
        mapped_nodes = {from_ for from_ in self.mapping if node_to in self.mapping[from_]}
        nodes_from = nodes_from.intersection(mapped_nodes)
        assert len(nodes_from) == 1
        self.references[node_to] = next(iter(nodes_from))

    def get_mapping(self, map_type):
        if self.blocks_from is None:
            return None
        mapping = Mapping(self.blocks_from, self.blocks_to, dict(self.mapping),
                          self.references, ff_from=self.ff_from, ff_to=self.ff_to,
                          type=map_type, names=tuple(self.names))
        return mapping


class MappingDirector:
    RESNAME_NUM_SEP = '#'
    RESIDUE_ATOM_SEP = ':'
    MAP_TYPES = ['block', 'modification']

    def __init__(self, builder=None):
        if builder is None:
            self.builder = MappingBuilder()
        else:
            self.builder = builder
        self.identifiers = {}
        self.section = None
        self.map_type = None
        self.from_ff = None
        self.to_ff = None
        self._current_id = {'from_': None, 'to_': None}

    def parse(self, file_handle):
        for lineno, line in enumerate(file_handle, 1):
            # TODO split off comments
            line = line.strip()
            if not line:
                continue
            if self._is_section_header(line):
                map_type = self.map_type
                new_mapping = self._header(line)
                if new_mapping:
                    self._current_id = {'from_': None, 'to_': None}
                    self.identifiers = {}
                    mapping = self.builder.get_mapping(map_type)
                    if mapping is not None:
                        yield mapping
                        self.builder.reset()
            else:
                self._parse_line(line, lineno)
        mapping = self.builder.get_mapping(self.map_type)
        if mapping is not None:
            yield mapping

    @staticmethod
    def _is_section_header(line):
        return line.startswith('[') and line.endswith(']')

    def _parse_line(self, line, lineno):
        try:
            getattr(self, self.section)(line)
        except Exception:
            LOGGER.error("Problems parsing line {}. I think it should be a '{}'"
                         " line, but I can't parse it as such.",
                         lineno, self.section)
            raise

    def __getattr__(self, key):
        key = key.replace(' ', '_')
        key = '_' + key
        return self.__getattribute__(key)

    def _header(self, line):
        section = line.strip('[ ]').casefold()
        try:
            if section not in self.MAP_TYPES:
                getattr(self, section)
        except AttributeError as err:
            raise KeyError('Section "{}" is unknown'.format(section)) from err
        else:
            if section in self.MAP_TYPES:
                self.map_type = section
                return True
            else:
                self.section = section
                return False

    def _parse_blocks(self, line):
        tokens = list(_tokenize(line))
        if len(tokens) == 2 and tokens[1].startswith('{') and tokens[1].endswith('}'):
            # It's definitely full spec.
            identifier = tokens[0]
            attrs = _parse_atom_attributes(tokens[1])
            yield identifier, attrs
        else:  # It must be shorthand
            for identifier in tokens:
                resname, resid = self._parse_block_shorthand(identifier)
                attrs = {'resname': resname, 'resid': resid}
                yield identifier, attrs

    def _parse_block_shorthand(self, token):
        if self.RESNAME_NUM_SEP in token:
            # PO4#3
            resname, resid = token.split(self.RESNAME_NUM_SEP)
        elif token[-1].isdigit():
            # ALA2
            for idx, char in enumerate(reversed(token)):
                if not char.isdigit():
                    idx = len(token) - idx
                    resname = token[:idx]
                    resid = int(token[idx:])
                    break
        else:
            # ALA
            resname = token
            resid = 1
        return resname, resid

    def _resolve_atom_spec(self, atom_str, prefix=''):
        if self.RESIDUE_ATOM_SEP in atom_str:
            id_, name = atom_str.split(self.RESIDUE_ATOM_SEP)
        else:
            id_, name = None, atom_str

        if id_ is None:
            options = {name for name in self.identifiers if name.startswith(prefix)}
            if len(options) == 1:
                id_ = next(iter(options))
                id_ = id_[len(prefix):]

        if id_ is None:
            attrs = self._current_id[prefix][0].copy()
        else:
            attrs = self.identifiers[prefix + id_][0].copy()
            self._current_id[prefix] = self.identifiers[prefix + id_]

        attrs['atomname'] = name
        return attrs

    def _to(self, line):
        self.to_ff = line
        self.builder.to_ff(self.to_ff)

    def _from(self, line):
        self.from_ff = line
        self.builder.from_ff(self.from_ff)

    def _from_blocks(self, line):
        for identifier, attrs in self._parse_blocks(line):
            block = get_block(self.from_ff, self.map_type, attrs['resname'])
            self.builder.add_block_from(block)
            self.identifiers['from_' + identifier] = attrs, block

    def _to_blocks(self, line):
        for identifier, attrs in self._parse_blocks(line):
            block = get_block(self.to_ff, self.map_type, attrs['resname'])
            self.builder.add_block_to(block)
            self.identifiers['to_' + identifier] = attrs, block

    def _from_nodes(self, line):
        name, *attrs = _tokenize(line)
        if attrs:
            attrs = _parse_atom_attributes(*attrs)
        else:
            attrs = {}
        if 'atomname' not in attrs:
            attrs['atomname'] = name
        print(0, attrs)
        self.builder.add_node_from(attrs)
    
    def _to_nodes(self, line):
        name, *attrs = _tokenize(line)
        
        if attrs:
            attrs = _parse_atom_attributes(*attrs)
        else:
            attrs = {}
        if 'atomname' not in attrs:
            attrs['atomname'] = name
        self.builder.add_node_to(attrs)

    def _from_edges(self, line):
        at1, at2 = line.split()
        attrs1 = self._resolve_atom_spec(at1, 'from_')
        attrs2 = self._resolve_atom_spec(at2, 'from_')
        # ..? For some reason modifications get a resname=''.
        if self.map_type == 'modification':
            attrs1['resname'] = ''
            attrs2['resname'] = ''
        self.builder.add_edge_from(attrs1, attrs2)

    def _to_edges(self, line):
        at1, at2 = line.split()
        attrs1 = self._resolve_atom_spec(at1, 'to_')
        attrs2 = self._resolve_atom_spec(at2, 'to_')
        # ..? For some reason modifications get a resname=''.
        if self.map_type == 'modification':
            attrs1['resname'] = ''
            attrs2['resname'] = ''
        self.builder.add_edge_to(attrs1, attrs2)

    def _mapping(self, line):
        from_, to_, *weight = line.split()
        if weight:
            weight = int(weight[0])
        else:
            weight = 1

        attrs_from = self._resolve_atom_spec(from_, 'from_')
        attrs_to = self._resolve_atom_spec(to_, 'to_')
        # ..? For some reason modifications get a resname=''.
        if self.map_type == 'modification':
            attrs_from['resname'] = ''
            attrs_to['resname'] = ''

        self.builder.add_mapping(attrs_from, attrs_to, weight)

    def _reference_atoms(self, line):
        to_, from_ = line.split()
        attrs_to = self._resolve_atom_spec(to_, 'to_')
        attrs_from = self._resolve_atom_spec(from_, 'to_')
        self.builder.add_reference(attrs_to, attrs_from)


def parse_mapping_file(filepath):
    with open(filepath) as map_in:
        director = MappingDirector()
        mappings = list(director.parse(map_in))
    return mappings
