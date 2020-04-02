###############################################################################
#
# keldy
#
# Copyright (c) 2019-2020 The Simons Foundation
#   authors: Philipp Dumitrescu
#
# keldy is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# keldy is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# keldy. If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################/
r"""
DOC

"""

from cpp2py import Cpp2pyInfoBase


class Cpp2pyInfo(Cpp2pyInfoBase):

    table_imports = {
        'keldy': 'keldy.common',
        'keldy::warpers': 'keldy.warpers',
        'keldy::impurity_oneband': 'keldy.impurity_oneband',
    }

    # _table_converters = {}

    table_converters = dict(
    )  # (k, "triqs/cpp2py_converters/%s.hpp"%v) for (k,v) in _table_converters.items())


# def _get_cpp2py_wrapped_class_enums():
#     return {'module_name' : 'UNUSED', 'includes' : "['<triqs/cpp2py_converters.hpp>']"}

__all__ = [
    'Cpp2pyInfo', 'common', 'warpers', 'impurity_oneband', 'visualization'
]
