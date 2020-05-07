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
                 c_type = "std::vector<std::pair<std::array<double_or_int, {}>, dcomplex>>".format(N + M),
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

    # module.add_function("friend double sum_moduli(keldy::binner::sparse_binner_t const& in)")
    # module.add_function("friend double sum_moduli(keldy::binner::sparse_binner_t<{}, {}> const& in)".format(N, M))

    ###########################
    # The class binner_t<N, M>
    c = class_(
               py_type = "Binner_{}_{}".format(N, M),
               c_type = "binner_t<{}, {}>".format(N, M),
               c_type_absolute = "keldy::binner::binner_t<{}, {}>".format(N, M),
               hdf5 = False,
               is_printable = False
    )

    c.add_constructor("""(std::array<std::tuple<double, double, size_t>, {}> _continuous_axes, std::array<size_t, {}> _discreet_axes = {{}})""".format(N, M),
                      doc = "")

    c.add_method_copy()

    args = ""
    for n in range(N):
        args += "double a{}, ".format(n)
    for m in range(M):
        args += "long b{}, ".format(m)
    c.add_method("void accumulate(dcomplex value, " + args[:-2] + ")")

    c.add_method("array<dcomplex, {}> get_data()".format(N + M))
    c.add_method("array<unsigned long, {}> get_nr_values_added()".format(N + M))
    c.add_method("unsigned long  get_nr_values_dropped()")
    c.add_method("auto get_continuous_axes()")
    c.add_method("auto get_discreet_axes()")

    c.add_method("auto get_bin_coord(size_t axis = 0)")
    c.add_method("auto get_bin_size(size_t axis = 0)")

    module.add_class(c)

module.generate_code()
