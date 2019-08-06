# Generated automatically using the command :
# c++2py ../../c++/keldy/impurity_oneband/compute_obs.hpp --members_read_only -N keldy -N keldy::impurity_oneband -a keldy -m impurity_oneband_module -o impurity_oneband_module -C pytriqs --cxxflags="-std=c++17 " --includes ../../c++
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "impurity_oneband_module", doc = r"", app_name = "keldy")

# Imports
module.add_imports(*['pytriqs.gf'])

# Add here all includes
module.add_include("keldy/impurity_oneband/compute_obs.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <cpp2py/converters/function.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/gf.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;
""")

module.add_enum("spin_t", ['spin_t::up', 'spin_t::down'], "keldy", doc = r"""""")
module.add_enum("keldysh_idx_t", ['keldysh_idx_t::forward', 'keldysh_idx_t::backward'], "keldy", doc = r"""""")

# The class gf_index_t
c = class_(
        py_type = "GfIndexT",  # name of the python class
        c_type = "keldy::impurity_oneband::gf_index_t",   # name of the C++ class
        doc = r"""Point of the Contour Keldysh Green Function (time, spin, keldysh_idx)""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "time",
             c_type = "keldy::time_real_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "spin",
             c_type = "keldy::spin_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "k_idx",
             c_type = "keldy::keldysh_idx_t",
             read_only= True,
             doc = r"""""")

c.add_constructor("""()""", doc = r"""Constructor: (time, spin, keldysh_idx)""")

c.add_constructor("""(keldy::time_real_t time_, int spin_, int k_idx_)""", doc = r"""""")

module.add_class(c)

# The class g0_model
c = class_(
        py_type = "G0Model",  # name of the python class
        c_type = "keldy::impurity_oneband::g0_model",   # name of the C++ class
        doc = r"""Defines model throuh non-interacting Green function g_lesser / g_greater""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "g0_lesser",
             c_type = "block_gf<triqs::gfs::retime, triqs::gfs::scalar_valued>",
             read_only= True,
             doc = r"""Lesser Green function $G^{<}_{\sigma}(t)$; block spin $\sigma$ {up, down}""")

c.add_member(c_name = "g0_greater",
             c_type = "block_gf<triqs::gfs::retime, triqs::gfs::scalar_valued>",
             read_only= True,
             doc = r"""Greater Green function $G^{>}_{\sigma}(t)$; block spin $\sigma$ {up, down}""")

c.add_member(c_name = "param_",
             c_type = "keldy::impurity_oneband::model_param_t",
             read_only= True,
             doc = r"""""")

c.add_constructor("""(**keldy::impurity_oneband::model_param_t)""", doc = r"""



+-------------------+-------------+------------+---------------+
| Parameter Name    | Type        | Default    | Documentation |
+===================+=============+============+===============+
| beta              | double      | 1.0        |               |
+-------------------+-------------+------------+---------------+
| bias_V            | double      | 0.0        |               |
+-------------------+-------------+------------+---------------+
| eps_d             | double      | 0.0        |               |
+-------------------+-------------+------------+---------------+
| Gamma             | double      | 1.0        |               |
+-------------------+-------------+------------+---------------+
| time_max          | double      | +100.0     |               |
+-------------------+-------------+------------+---------------+
| nr_time_points_gf | int         | 1000       |               |
+-------------------+-------------+------------+---------------+
| bath_type         | std::string | "flatband" |               |
+-------------------+-------------+------------+---------------+
""")

c.add_method("""void make_semicircular_model ()""",
             doc = r"""""")

c.add_method("""void make_flat_band ()""",
             doc = r"""""")

module.add_class(c)

# The class g0_keldysh_contour_t
c = class_(
        py_type = "G0KeldyshContourT",  # name of the python class
        c_type = "keldy::impurity_oneband::g0_keldysh_contour_t",   # name of the C++ class
        doc = r"""Adapt g0_lesser and g0_greater into Green function on Keldysh contour""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "model",
             c_type = "keldy::impurity_oneband::g0_model",
             read_only= True,
             doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_model model_)""", doc = r"""""")

c.add_method("""keldy::dcomplex operator() (keldy::impurity_oneband::gf_index_t a, keldy::impurity_oneband::gf_index_t b, bool use_lesser_at_eq_points)""",
               name = "__call__",
             doc = r"""Evalutate G, passing two Keldysh contour points""")

module.add_class(c)

# The class integrand_g_t1t2_direct
c = class_(
        py_type = "IntegrandGT1t2Direct",  # name of the python class
        c_type = "keldy::impurity_oneband::integrand_g_t1t2_direct",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "g0",
             c_type = "keldy::impurity_oneband::g0_keldysh_contour_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "external_A",
             c_type = "keldy::impurity_oneband::gf_index_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "external_B",
             c_type = "keldy::impurity_oneband::gf_index_t",
             read_only= True,
             doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_keldysh_contour_t g0_, keldy::impurity_oneband::gf_index_t external_A_, keldy::impurity_oneband::gf_index_t external_B_)""", doc = r"""""")

c.add_method("""keldy::dcomplex operator() (std::vector<double> times)""",
               name = "__call__",
             doc = r"""Returns integrand for the specified times""")

module.add_class(c)

# The class new_class
c = class_(
        py_type = "NewClass",  # name of the python class
        c_type = "keldy::new_class",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

module.add_class(c)

# The class sobol
c = class_(
        py_type = "Sobol",  # name of the python class
        c_type = "keldy::sobol",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(int dim, int log_max_points_ = 31)""", doc = r"""""")

c.add_method("""std::vector<double> operator() ()""",
               name = "__call__",
             doc = r"""""")

c.add_method("""void seed (int k)""",
             doc = r"""""")

c.add_method("""void discard (int nr_discard)""",
             doc = r"""""")

module.add_class(c)

# The class warper_plasma_simple_t
c = class_(
        py_type = "WarperPlasmaSimpleT",  # name of the python class
        c_type = "keldy::warper_plasma_simple_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""()""", doc = r"""""")

c.add_constructor("""(std::function<double(double)> f1_, double t_max_, int nr_function_sample_points)""", doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double evaluate_warping_function (std::vector<double> ui_vec)""",
             doc = r"""""")

module.add_class(c)

# The class compute_charge_Q
c = class_(
        py_type = "ComputeChargeQ",  # name of the python class
        c_type = "keldy::impurity_oneband::compute_charge_Q",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(int order, double time, keldy::impurity_oneband::model_param_t params, int nr_sample_points_ansatz)""", doc = r"""""")

c.add_method("""void run (int nr_steps)""",
             doc = r"""""")

c.add_method("""keldy::dcomplex reduce_result ()""",
             doc = r"""""")

c.add_method("""int get_nr_points_run ()""",
             doc = r"""""")

module.add_class(c)

module.add_function ("std::vector<double> keldy::vi_from_ui (double t_max, std::vector<double> u_times)", doc = r"""""")

module.add_function ("std::vector<double> keldy::ui_from_vi (double t_max, std::vector<double> v_times)", doc = r"""""")


# Converter for model_param_t
c = converter_(
        c_type = "keldy::impurity_oneband::model_param_t",
        doc = r"""""",
)
c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """ 1.0 """,
             doc = r"""""")

c.add_member(c_name = "bias_V",
             c_type = "double",
             initializer = """ 0.0 """,
             doc = r"""""")

c.add_member(c_name = "eps_d",
             c_type = "double",
             initializer = """ 0.0 """,
             doc = r"""""")

c.add_member(c_name = "Gamma",
             c_type = "double",
             initializer = """ 1.0 """,
             doc = r"""""")

c.add_member(c_name = "time_max",
             c_type = "double",
             initializer = """ +100.0 """,
             doc = r"""""")

c.add_member(c_name = "nr_time_points_gf",
             c_type = "int",
             initializer = """ 1000 """,
             doc = r"""""")

c.add_member(c_name = "bath_type",
             c_type = "std::string",
             initializer = """ "flatband" """,
             doc = r"""""")

module.add_converter(c)


module.generate_code()