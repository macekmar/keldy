# Generated automatically using the command :
# c++2py ../../c++/keldy/qrng.hpp -a keldy -m common -o common -C pytriqs --cxxflags="-std=c++17 " --includes ../../c++ --members_read_only -N keldy
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "common", doc = r"", app_name = "keldy")

# Imports

# Add here all includes
module.add_include("keldy/qrng.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/vector.hpp>

using namespace keldy;
""")

module.add_enum("spin_t", ['spin_t::up', 'spin_t::down'], "keldy", doc = r"""""")
module.add_enum("keldysh_idx_t", ['keldysh_idx_t::forward', 'keldysh_idx_t::backward'], "keldy", doc = r"""""")

# The class contour_pt_t
c = class_(
        py_type = "ContourPtT",  # name of the python class
        c_type = "keldy::contour_pt_t",   # name of the C++ class
        doc = r"""Point of the Keldysh Contour (time, keldysh_idx, timesplit)""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "time",
             c_type = "keldy::time_real_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "k_idx",
             c_type = "keldy::keldysh_idx_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "timesplit_n",
             c_type = "int",
             read_only= True,
             doc = r"""""")

module.add_class(c)

# The class sobol
c = class_(
        py_type = "Sobol",  # name of the python class
        c_type = "keldy::sobol",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(int dim, int rng_state_seed, int log_max_points_ = 31)""", doc = r"""""")

c.add_constructor("""(int dim, int rng_state_seed, bool do_shift, bool do_scramble, int rng_seed_shift, int log_max_points_ = 31)""", doc = r"""""")

c.add_method("""std::vector<double> operator() ()""",
             name = "__call__",
             doc = r"""""")

c.add_method("""void seed (int k)""",
             doc = r"""""")

c.add_method("""void discard (int nr_discard)""",
             doc = r"""""")

module.add_class(c)

module.add_function ("int keldy::compare_3way (keldy::contour_pt_t a, keldy::contour_pt_t b)", doc = r"""""")



module.generate_code()