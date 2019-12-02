# Generated automatically using the command :
# c++2py /Users/pdumitrescu/Documents/coding/repos/keldy/c++/keldy/interfaces/cuba_interface.hpp -C pytriqs --cxxflags="-std=c++17 " --includes ../../c++ --includes /usr/local/Cellar/gsl/2.6/include/gsl --target_file_only
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "cuba_interface", doc = r"", app_name = "cuba_interface")

# Imports

# Add here all includes
module.add_include("keldy/interfaces/cuba_interface.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/function.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>

""")


# The class cuba_wrapper
c = class_(
        py_type = "CubaWrapper",  # name of the python class
        c_type = "keldy::cuba_wrapper",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "f",
             c_type = "std::function<double (std::vector<double>)>",
             read_only= False,
             doc = r"""""")

c.add_constructor("""(std::function<double(std::vector<double>)> f_, int dim, keldy::cuba_common_param in_)""", doc = r"""""")

c.add_method("""void run_vegas (keldy::cuba_vegas_param in_v)""",
             doc = r"""""")

c.add_method("""void run_suave (keldy::cuba_suave_param in_s)""",
             doc = r"""""")

c.add_method("""void run_cuhre (int key_integration_order)""",
             doc = r"""""")

c.add_method("""keldy::cuba_output get_output ()""",
             doc = r"""""")

module.add_class(c)

module.add_function ("void keldy::fake_output (**keldy::cuba_output)", doc = r"""



+----------------+--------+---------+---------------+
| Parameter Name | Type   | Default | Documentation |
+================+========+=========+===============+
| n_regions      | int    | 0       |               |
+----------------+--------+---------+---------------+
| n_eval         | int    | 0       |               |
+----------------+--------+---------+---------------+
| error_flag     | int    | 0       |               |
+----------------+--------+---------+---------------+
| result         | double | 0.0     |               |
+----------------+--------+---------+---------------+
| error          | double | 0.0     |               |
+----------------+--------+---------+---------------+
| chi_sq_prob    | double | 0.0     |               |
+----------------+--------+---------+---------------+
""")

module.add_function ("void keldy::fake_common (**keldy::cuba_common_param)", doc = r"""



+-------------------------------+-------------+-----------+---------------+
| Parameter Name                | Type        | Default   | Documentation |
+===============================+=============+===========+===============+
| n_dim                         | int         | --        |               |
+-------------------------------+-------------+-----------+---------------+
| n_components                  | int         | 1         |               |
+-------------------------------+-------------+-----------+---------------+
| n_points_vectorization        | int         | 1         |               |
+-------------------------------+-------------+-----------+---------------+
| error_eps_rel                 | double      | 1e-12     |               |
+-------------------------------+-------------+-----------+---------------+
| error_eps_abs                 | double      | 1e-12     |               |
+-------------------------------+-------------+-----------+---------------+
| flags                         | int         | 0         |               |
+-------------------------------+-------------+-----------+---------------+
| verbosity                     | int         | 0         |               |
+-------------------------------+-------------+-----------+---------------+
| use_last_sampleset_only       | bool        | false     |               |
+-------------------------------+-------------+-----------+---------------+
| sample_function_smoothing_off | bool        | false     |               |
+-------------------------------+-------------+-----------+---------------+
| store_state_after_run         | bool        | false     |               |
+-------------------------------+-------------+-----------+---------------+
| rng_type                      | std::string | "sobol"   |               |
+-------------------------------+-------------+-----------+---------------+
| seed                          | int         | 1         |               |
+-------------------------------+-------------+-----------+---------------+
| randlux_level                 | int         | 0         |               |
+-------------------------------+-------------+-----------+---------------+
| min_number_evaluations        | int         | 1000      |               |
+-------------------------------+-------------+-----------+---------------+
| max_number_evaluations        | int         | int(1e10) |               |
+-------------------------------+-------------+-----------+---------------+
""")

module.add_function ("void keldy::fake_vegas (**keldy::cuba_vegas_param)", doc = r"""



+--------------------------------+------+---------+---------------+
| Parameter Name                 | Type | Default | Documentation |
+================================+======+=========+===============+
| n_evals_per_iteration_start    | int  | 1000    |               |
+--------------------------------+------+---------+---------------+
| n_evals_per_iteration_increase | int  | 1000    |               |
+--------------------------------+------+---------+---------------+
| n_samples_per_batch            | int  | 500     |               |
+--------------------------------+------+---------+---------------+
| internal_store_grid_nr         | int  | 0       |               |
+--------------------------------+------+---------+---------------+
""")

module.add_function ("void keldy::fake_suave (**keldy::cuba_suave_param)", doc = r"""



+---------------------------------+--------+---------+---------------+
| Parameter Name                  | Type   | Default | Documentation |
+=================================+========+=========+===============+
| n_new_evals_each_subdivision    | int    | 1000    |               |
+---------------------------------+--------+---------+---------------+
| n_min_samples_region_threashold | int    | 1000    |               |
+---------------------------------+--------+---------+---------------+
| flatness_parameter_p            | double | 10.0    |               |
+---------------------------------+--------+---------+---------------+
""")


# Converter for cuba_output
c = converter_(
        c_type = "keldy::cuba_output",
        doc = r"""""",
)
c.add_member(c_name = "n_regions",
             c_type = "int",
             initializer = """ 0 """,
             doc = r"""""")

c.add_member(c_name = "n_eval",
             c_type = "int",
             initializer = """ 0 """,
             doc = r"""""")

c.add_member(c_name = "error_flag",
             c_type = "int",
             initializer = """ 0 """,
             doc = r"""""")

c.add_member(c_name = "result",
             c_type = "double",
             initializer = """ 0.0 """,
             doc = r"""""")

c.add_member(c_name = "error",
             c_type = "double",
             initializer = """ 0.0 """,
             doc = r"""""")

c.add_member(c_name = "chi_sq_prob",
             c_type = "double",
             initializer = """ 0.0 """,
             doc = r"""""")

module.add_converter(c)

# Converter for cuba_common_param
c = converter_(
        c_type = "keldy::cuba_common_param",
        doc = r"""""",
)
c.add_member(c_name = "n_dim",
             c_type = "int",
             initializer = """  """,
             doc = r"""""")

c.add_member(c_name = "n_components",
             c_type = "int",
             initializer = """ 1 """,
             doc = r"""""")

c.add_member(c_name = "n_points_vectorization",
             c_type = "int",
             initializer = """ 1 """,
             doc = r"""""")

c.add_member(c_name = "error_eps_rel",
             c_type = "double",
             initializer = """ 1e-12 """,
             doc = r"""""")

c.add_member(c_name = "error_eps_abs",
             c_type = "double",
             initializer = """ 1e-12 """,
             doc = r"""""")

c.add_member(c_name = "flags",
             c_type = "int",
             initializer = """ 0 """,
             doc = r"""""")

c.add_member(c_name = "verbosity",
             c_type = "int",
             initializer = """ 0 """,
             doc = r"""""")

c.add_member(c_name = "use_last_sampleset_only",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""""")

c.add_member(c_name = "sample_function_smoothing_off",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""""")

c.add_member(c_name = "store_state_after_run",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""""")

c.add_member(c_name = "rng_type",
             c_type = "std::string",
             initializer = """ "sobol" """,
             doc = r"""""")

c.add_member(c_name = "seed",
             c_type = "int",
             initializer = """ 1 """,
             doc = r"""""")

c.add_member(c_name = "randlux_level",
             c_type = "int",
             initializer = """ 0 """,
             doc = r"""""")

c.add_member(c_name = "min_number_evaluations",
             c_type = "int",
             initializer = """ 1000 """,
             doc = r"""""")

c.add_member(c_name = "max_number_evaluations",
             c_type = "int",
             initializer = """ int(1e10) """,
             doc = r"""""")

module.add_converter(c)

# Converter for cuba_vegas_param
c = converter_(
        c_type = "keldy::cuba_vegas_param",
        doc = r"""""",
)
c.add_member(c_name = "n_evals_per_iteration_start",
             c_type = "int",
             initializer = """ 1000 """,
             doc = r"""""")

c.add_member(c_name = "n_evals_per_iteration_increase",
             c_type = "int",
             initializer = """ 1000 """,
             doc = r"""""")

c.add_member(c_name = "n_samples_per_batch",
             c_type = "int",
             initializer = """ 500 """,
             doc = r"""""")

c.add_member(c_name = "internal_store_grid_nr",
             c_type = "int",
             initializer = """ 0 """,
             doc = r"""""")

module.add_converter(c)

# Converter for cuba_suave_param
c = converter_(
        c_type = "keldy::cuba_suave_param",
        doc = r"""""",
)
c.add_member(c_name = "n_new_evals_each_subdivision",
             c_type = "int",
             initializer = """ 1000 """,
             doc = r"""""")

c.add_member(c_name = "n_min_samples_region_threashold",
             c_type = "int",
             initializer = """ 1000 """,
             doc = r"""""")

c.add_member(c_name = "flatness_parameter_p",
             c_type = "double",
             initializer = """ 10.0 """,
             doc = r"""""")

module.add_converter(c)


module.generate_code()