################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2016-2018, N. Wentzell
# Copyright (C) 2018-2019, The Simons Foundation
#   author: N. Wentzell
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

r"""
DOC

"""

from cpp2py import Cpp2pyInfoBase

class Cpp2pyInfo(Cpp2pyInfoBase):

    table_imports = {
        'keldy::warpers' : 'keldy.warpers',
        'keldy::impurity_oneband' : 'keldy.impurity_oneband',
        }

    # _table_converters = {
    #     'triqs::arrays::array' : 'arrays',
    #     'triqs::arrays::matrix' : 'arrays',
    #     'triqs::arrays::vector' : 'arrays',
    #     'triqs::gfs::gf*' : 'gf',
    #     'triqs::gfs::block_gf*' : 'gf',
    #     'triqs::gfs::block2_gf*' : 'gf',
    #     'triqs::operators::many_body_operator*' : 'operators_real_complex',
    #     'triqs::hilbert_space::fundamental_operator_set' : 'fundamental_operator_set',
    #     'triqs::utility::real_or_complex' : 'real_or_complex',
    #     'triqs::h5::group' : 'h5'
    #     }

    table_converters = dict ()# (k, "triqs/cpp2py_converters/%s.hpp"%v) for (k,v) in _table_converters.items())

# def _get_cpp2py_wrapped_class_enums():
#     return {'module_name' : 'UNUSED', 'includes' : "['<triqs/cpp2py_converters.hpp>']"}

# import visualization

__all__ = ['Cpp2pyInfo', 'warpers', 'impurity_oneband', 'visualization']



# import impurity_oneband_module as impurity_oneband
# import gsl_vegas_wrap

# __all__ = ['Model', 'WarperPlasmaSimpleT']
