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
Helper functions for parsers
"""
from .ffinput import _tokenize, _parse_macro, _substitute_macros
from collections import deque


class SectionParser(type):
    def __new__(cls, name, bases, attrs, **kwargs):
        obj = super().__new__(cls, name, bases, attrs, **kwargs)
        if not hasattr(obj, 'METH_DICT'):
            obj.METH_DICT = {}
        mapping = obj.METH_DICT
        for attribute_name in dir(obj):
            attribute = getattr(obj, attribute_name)
            if hasattr(attribute, '_section_names'):
                for names, kwargs in attribute._section_names.items():
                    mapping[names] = (attribute, kwargs)
        return obj

    @staticmethod
    def section_parser(*names, **kwargs):
        def wrapper(method):
            if not hasattr(method, '_section_names'):
                method._section_names = {}
            method._section_names[names] = kwargs
            return method
        return wrapper


class LineParser:
    def parse(self, file_handle):
        """
        Reads lines from `file_handle`, and calls :meth:`dispatch` to find
        which method to call to do the actual parsing. Yields the result of
        that call, if it's not `None`.
        At the end, calls :meth:`finalize`, and yields its results, iff
        it's not None.

        file_handle: collections.abc.Iterable[str]
            The data to parse. Should produce lines of data.

        Yields
        ------
        object
            The results of dispatching to parsing methods, and of
            :meth:`finalize`.
        """
        for lineno, line in enumerate(file_handle, 1):
            # TODO split off comments
            line, _ = split_comments(line, self.COMMENT_CHAR)
            if not line:
                continue
            result = self.dispatch(line)(line, lineno)
            if result is not None:
                yield result

        result = self.finalize(lineno)
        if result is not None:
            yield result

    def dispatch(self, line):
        """
        Finds the correct method to parse `line`. Currently always returns
        :meth:`parse_line`
        """
        return self.parse_line

    def parse_line(self, line, lineno):
        """
        Does nothing and should be overridden by subclasses.
        """
        return


class SectionLineParser(LineParser, metaclass=SectionParser):
    def __init__(self, *args, **kwargs):
        self.macros = {}
        self.section = None
        # TODO: get rid of map_type
        self.map_type = None
        super().__init__(*args, **kwargs)

    def dispatch(self, line):
        """
        Looks at `line` to see what kind of line it is, and returns either
        :meth:`parse_header` if `line` is a section header or
        :meth:`parse_section` otherwise. Calls :meth:`is_section_header` to see
        whether `line` is a section header or not.

        Parameters
        ----------
        line: str

        Returns
        -------
        callable
        """
        if self.is_section_header(line):
            return self.parse_header
        else:
            return self.parse_section

    def finalize(self, lineno=0):
        """
        Called after the last line has been parsed to wrap up. Resets
        the instance and calls :meth:`finalize_section`.
        """
        result = self.finalize_section(self.section)
        self.macros = {}
        self.section = None
        # TODO: get rid of map_type
        self.map_type = None
        return result

    def finalize_section(self, section_name):
        """
        Called once a section is finished. Currently does nothing.
        """
        return

    def parse_section(self, line, lineno):
        """
        Parse `line` with line number `lineno` by looking up the section in
        :attr:`METH_DICT` and calling that method.

        Parameters
        ----------
        line: str
        lineno: int

        Returns
        -------
        object
            The result returned by calling the registered method.
        """
        line = _substitute_macros(line, self.macros)
        if self.section not in self.METH_DICT:
            raise IOError("Can't parse line {} in section '{}' because the "
                          "section is unknown".format(lineno, self.section))
        try:
            method, kwargs = self.METH_DICT[self.section]
            return method(self, line, lineno, **kwargs)
        except Exception as error:
            raise IOError("Problems parsing line {}. I think it should be a "
                          "'{}' line, but I can't parse it as such."
                          "".format(lineno, self.section)) from error

    def parse_header(self, line, lineno=0):
        """
        Parses a section header. Sets :attr:`section` and :attr:`map_type` when
        applicable. Does not check whether `line` is a valid section header.

        Parameters
        ----------
        line: str

        Returns
        -------
        object
            The result of calling :meth:`finalize_section`, which is called
            if the header is specified in :attr:`MAP_TYPES`

        Raises
        ------
        KeyError
            If the section header is unknown.
        """
        # TODO: This method needs to be purged of references to self.MAP_TYPES
        #       and self.map_type
        prev_section = self.section

        section = (line.strip('[ ]').casefold(), )
        self.section = section
        if section not in self.METH_DICT and section not in self.MAP_TYPES:
            raise IOError("Section '{}' on line {} is unknown. The following "
                          "sections are known: {}."
                          "".format(section, lineno,
                                    list(self.METH_DICT.keys()) + self.MAP_TYPES))

        section_end = self.section in self.MAP_TYPES
        if section_end:
            self.map_type = self.section
            result = self.finalize_section(prev_section)
            return result

    @staticmethod
    def is_section_header(line):
        """
        Parameters
        ----------
        line: str
            A line of text.

        Returns
        -------
        bool
            ``True`` iff `line` is a section header.
        """
        return line.startswith('[') and line.endswith(']')

    @SectionParser.section_parser('macros')
    def _macros(self, line, lineno=0):
        """
        Parses a "macros" section. Adds to :attr:`macros`.

        Parameters
        ----------
        line: str
        """
        line = deque(_tokenize(line))
        _parse_macro(line, self.macros)


def split_comments(line, comment_char=';'):
    split = line.split(comment_char, 1)
    data = split[0].strip()
    if len(split) == 1:
        comment = ''
    else:
        comment = split[1].strip()
    return data, comment
