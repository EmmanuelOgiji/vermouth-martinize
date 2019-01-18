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


def parser_class(dict_name='METH_DICT', attr_name='_name'):
    def _parser_class(cls):
        if not hasattr(cls, dict_name):
            setattr(cls, dict_name, {})
        mapping = getattr(cls, dict_name)
        for attribute_name in dir(cls):
            attribute = getattr(cls, attribute_name)
            if hasattr(attribute, attr_name):
                mapping[getattr(attribute, attr_name)] = attribute
        return cls
    return _parser_class


def section_parser(name):
    def wrapper(method):
        method._section_name = name
        return method
    return wrapper
