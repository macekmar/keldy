# NOT generated automatically
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "binner", doc = "Generic multidimentional binner", app_name = "keldy")

module.add_include("keldy/binner.hpp")
module.add_include("keldy/common.hpp")
module.add_include("<triqs/arrays.hpp>")

module.add_include("<cpp2py/converters/pair.hpp>")
module.add_include("<cpp2py/converters/vector.hpp>")
module.add_include("<cpp2py/converters/std_array.hpp>")
module.add_include("<cpp2py/converters/tuple.hpp>")
module.add_include("<cpp2py/converters/variant.hpp>")
module.add_include("<triqs/cpp2py_converters.hpp>")

module.add_using("namespace keldy::binner")

for N, M in [(1, 0), (1, 1), (2, 0), (2, 1)]:
    # The class sparse_binner_t<N, M>
    c = class_(
               py_type = "SparseBinner_{}_{}".format(N, M),
               c_type = "sparse_binner_t<{}, {}>".format(N, M),
               c_type_absolute = "keldy::binner::sparse_binner_t<{}, {}>".format(N, M),
               hdf5 = False,
               is_printable = False
    )

    c.add_constructor("""()""",
                    doc = """Construct empty sparse binner""")

    c.add_member(c_name = "data",
            c_type = "std::vector<std::pair<sparse_binner_t<{}, {}>::coord_arr_t, dcomplex>>".format(N, M),
                 read_only = False,
    )

    c.add_method_copy()

    args = ""
    for n in range(N):
        args += "double a{}, ".format(n)
    for m in range(M):
        args += "long b{}, ".format(m)
    c.add_method("void accumulate(dcomplex value, " + args[:-2] + ")")

    c.add_method("double sum_moduli()")

    module.add_class(c)

    ###########################
    # The class binner_t<N, M>
    for data_t in ['double', 'dcomplex']:
        c = class_(
                py_type = "Binner_{}_{}_{}".format(N, M, data_t),
                c_type = "binner_t<{}, {}, {}>".format(N, M, data_t),
                c_type_absolute = "keldy::binner::binner_t<{}, {}, {}>".format(N, M, data_t),
                hdf5 = False,
                is_printable = False
        )

        c.add_constructor("""(std::array<std::tuple<double, double, long>, {}> _continuous_axes, std::array<long, {}> _discreet_axes = {{}})""".format(N, M),
                        doc = "")

        c.add_method_copy()

        c.add_method("mda::array<{}, {}> const& get_data()".format(data_t, N + M))
        c.add_method("mda::array<long, {}> const& get_nr_values_added()".format(N + M))
        c.add_method("auto const& get_nr_values_dropped()")
        c.add_method("auto const& get_continuous_axes()")
        c.add_method("auto const& get_discreet_axes()")
        c.add_method("int get_nr_bins(int axis = 0)")

        c.add_method("auto get_bin_coord(int axis = 0)")
        c.add_method("auto get_bin_size(int axis = 0)")

        module.add_class(c)


module.generate_code()
