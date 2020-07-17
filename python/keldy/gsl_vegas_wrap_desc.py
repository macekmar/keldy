# Generated automatically using the command :
# c++2py ../../c++/keldy/interfaces/gsl_vegas_wrap.hpp -C triqs --cxxflags="-std=c++17 " --includes ../../c++ -N keldy
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "gsl_vegas_wrap", doc = r"", app_name = "gsl_vegas_wrap")

# Imports

# Add here all includes
module.add_include("keldy/interfaces/gsl_vegas_wrap.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/function.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>

using namespace keldy;
""")


# The class gsl_vegas_wrapper_t
c = class_(
        py_type = "GslVegasWrapperT",  # name of the python class
        c_type = "keldy::gsl_vegas_wrapper_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(std::function<double (std::vector<double>)> f, int dim, double hypercube_extent, std::string gsl_rng_name)""", doc = r"""""")

c.add_method("""double operator() (std::vector<double> x)""",
             name = "__call__",
             doc = r"""""")

c.add_method("""double get_result ()""",
             doc = r"""""")

c.add_method("""double get_error ()""",
             doc = r"""""")

c.add_method("""uint64_t get_nr_points_run ()""",
             doc = r"""""")

c.add_method("""double chisq ()""",
             doc = r"""""")

c.add_method("""void run (int nr_steps)""",
             doc = r"""""")

c.add_method("""void set_params (**keldy::gsl_monte_vegas_params_wrap)""",
             doc = r"""



+----------------+--------+---------+---------------+
| Parameter Name | Type   | Default | Documentation |
+================+========+=========+===============+
| alpha          | double | --      |               |
+----------------+--------+---------+---------------+
| iterations     | size_t | --      |               |
+----------------+--------+---------+---------------+
| stage          | int    | --      |               |
+----------------+--------+---------+---------------+
| mode           | int    | --      |               |
+----------------+--------+---------+---------------+
| verbose        | int    | --      |               |
+----------------+--------+---------+---------------+
""")

c.add_method("""keldy::gsl_monte_vegas_params_wrap get_params ()""",
             doc = r"""""")

module.add_class(c)


# Converter for gsl_monte_vegas_params_wrap
c = converter_(
        c_type = "keldy::gsl_monte_vegas_params_wrap",
        doc = r"""""",
)
c.add_member(c_name = "alpha",
             c_type = "double",
             initializer = """  """,
             doc = r"""""")

c.add_member(c_name = "iterations",
             c_type = "size_t",
             initializer = """  """,
             doc = r"""""")

c.add_member(c_name = "stage",
             c_type = "int",
             initializer = """  """,
             doc = r"""""")

c.add_member(c_name = "mode",
             c_type = "int",
             initializer = """  """,
             doc = r"""""")

c.add_member(c_name = "verbose",
             c_type = "int",
             initializer = """  """,
             doc = r"""""")

module.add_converter(c)


module.generate_code()