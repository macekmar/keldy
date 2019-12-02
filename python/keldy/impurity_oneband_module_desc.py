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
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
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

c.add_member(c_name = "k_idx",
             c_type = "keldy::keldysh_idx_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "timesplit_n",
             c_type = "int",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "spin",
             c_type = "keldy::spin_t",
             read_only= True,
             doc = r"""""")

c.add_constructor("""()""", doc = r"""Constructor: (time, spin, keldysh_idx)""")

c.add_constructor("""(keldy::time_real_t time_, int spin_, int k_idx_)""", doc = r"""""")

c.add_constructor("""(keldy::time_real_t time_, int spin_, int k_idx_, int timesplit_n_)""", doc = r"""""")

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

c.add_constructor("""(keldy::impurity_oneband::model_param_t parameters)""", doc = r"""""")

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

c.add_method("""keldy::dcomplex operator() (keldy::impurity_oneband::gf_index_t a, keldy::impurity_oneband::gf_index_t b, bool internal_point = true)""",
             name = "__call__",
             doc = r"""Evalutate G, passing two Keldysh contour points""")

module.add_class(c)

# The class integrand_g_direct
c = class_(
        py_type = "IntegrandGT1t2Direct",  # name of the python class
        c_type = "keldy::impurity_oneband::integrand_g_direct",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(keldy::impurity_oneband::g0_keldysh_contour_t g0_, keldy::impurity_oneband::gf_index_t external_A_, keldy::impurity_oneband::gf_index_t external_B_)""", doc = r"""""")

c.add_method("""keldy::impurity_oneband::integrand_g_direct::result_t operator() (std::vector<double> times)""",
             name = "__call__",
             doc = r"""Returns integrand for the specified times""")

module.add_class(c)

# The class sparse_kernel_binner
c = class_(
        py_type = "SparseKernelBinner",  # name of the python class
        c_type = "keldy::impurity_oneband::sparse_kernel_binner",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "data",
             c_type = "std::vector<std::pair<gf_index_t, dcomplex> >",
             read_only= True,
             doc = r"""""")

c.add_method("""double sum_weights ()""",
             doc = r"""""")

module.add_class(c)

# The class kernel_binner
c = class_(
        py_type = "KernelBinner",  # name of the python class
        c_type = "keldy::impurity_oneband::kernel_binner",   # name of the C++ class
        doc = r"""Kernel binner for Green Function :math:`K(Y, X')` with binning happening over :math:`Y`
 and :math:`X'` fixed by boundary conditions.

 TODO: How to include spin up / down separatley.""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""()""", doc = r"""""")

c.add_constructor("""(double t_min_, double t_max_, int n_bins_)""", doc = r"""""")

c.add_method("""triqs::arrays::array<std::complex<double>,2> get_values ()""",
             doc = r"""""")

c.add_method("""triqs::arrays::array<uint64_t,2> get_nr_values ()""",
             doc = r"""""")

c.add_method("""triqs::arrays::array<double,1> get_bin_times ()""",
             doc = r"""""")

c.add_method("""double get_bin_size ()""",
             doc = r"""""")

c.add_method("""int get_nr_point_dropped ()""",
             doc = r"""""")

c.add_method("""void accumulate (keldy::impurity_oneband::gf_index_t a, keldy::dcomplex value)""",
             doc = r"""Includes boundary points, so t_min <= t <= t_max. t_max gets put in last bin""")

module.add_class(c)

# The class integrand_g_kernel
c = class_(
        py_type = "IntegrandGKernel",  # name of the python class
        c_type = "keldy::impurity_oneband::integrand_g_kernel",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(keldy::impurity_oneband::g0_keldysh_contour_t g0_, keldy::impurity_oneband::gf_index_t g_idx_X_)""", doc = r"""""")

c.add_method("""keldy::impurity_oneband::integrand_g_kernel::result_t operator() (std::vector<double> times)""",
             name = "__call__",
             doc = r"""""")

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

c.add_constructor("""(int dim, int rng_seed, int log_max_points_ = 31)""", doc = r"""""")

c.add_method("""std::vector<double> operator() ()""",
             name = "__call__",
             doc = r"""""")

c.add_method("""void seed (int k)""",
             doc = r"""""")

c.add_method("""void discard (int nr_discard)""",
             doc = r"""""")

module.add_class(c)

# The class idenity_function
c = class_(
        py_type = "IdenityFunction",  # name of the python class
        c_type = "keldy::idenity_function",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_method("""double operator() (double t)""",
             name = "__call__",
             doc = r"""""")

module.add_class(c)

# The class warper_plasma_simple_t
c = class_(
        py_type = "WarperPlasmaSimpleT",  # name of the python class
        c_type = "keldy::warper_plasma_simple_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(double t_max_)""", doc = r"""""")

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

# The class compute_charge_Q_direct
c = class_(
        py_type = "ComputeChargeQDirect",  # name of the python class
        c_type = "keldy::impurity_oneband::compute_charge_Q_direct",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, int order, std::string warper_function_name, int nr_sample_points_warper)""", doc = r"""""")

c.add_method("""void run (int nr_steps)""",
             doc = r"""""")

c.add_method("""keldy::dcomplex reduce_result ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_run ()""",
             doc = r"""""")

c.add_method("""keldy::warper_plasma_simple_t get_warper ()""",
             doc = r"""""")

c.add_method("""keldy::impurity_oneband::integrand_g_direct get_integrand ()""",
             doc = r"""""")

module.add_class(c)


# The class compute_charge_Q_direct_gsl_vegas
c = class_(
        py_type = "ComputeChargeQDirectGSLVegas",  # name of the python class
        c_type = "keldy::impurity_oneband::compute_charge_Q_direct_gsl_vegas",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, int order, std::string gsl_rng_name)""", doc = r"""""")

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

c = class_(
        py_type = "ComputeChargeQDirectCuba",  # name of the python class
        c_type = "keldy::impurity_oneband::compute_charge_Q_direct_cuba",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "f",
             c_type = "std::function<double (std::vector<double>)>",
             read_only= False,
             doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, int order, keldy::cuba_common_param in)""", doc = r"""""")

c.add_method("""void run_vegas (keldy::cuba_vegas_param in_v)""",
             doc = r"""""")

c.add_method("""void run_suave (keldy::cuba_suave_param in_s)""",
             doc = r"""""")

c.add_method("""void run_cuhre (int key_integration_order)""",
             doc = r"""""")

c.add_method("""keldy::cuba_output get_output ()""",
             doc = r"""""")

module.add_class(c)

# The class compute_gf_kernel
c = class_(
        py_type = "ComputeGfKernel",  # name of the python class
        c_type = "keldy::impurity_oneband::compute_gf_kernel",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, int order, std::string warper_function_name, int nr_sample_points_warper)""", doc = r"""""")

c.add_method("""void run (int nr_steps)""",
             doc = r"""""")

c.add_method("""keldy::impurity_oneband::kernel_binner reduce_result ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_run ()""",
             doc = r"""""")

c.add_method("""keldy::warper_plasma_simple_t get_warper ()""",
             doc = r"""""")

c.add_method("""keldy::impurity_oneband::integrand_g_kernel get_integrand ()""",
             doc = r"""""")

module.add_class(c)

module.add_function ("void keldy::impurity_oneband::fake (**keldy::impurity_oneband::model_param_t)", doc = r"""



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
| alpha             | double      | 0.0        |               |
+-------------------+-------------+------------+---------------+
| bath_type         | std::string | "flatband" |               |
+-------------------+-------------+------------+---------------+
""")

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

c.add_member(c_name = "alpha",
             c_type = "double",
             initializer = """ 0.0 """,
             doc = r"""""")

c.add_member(c_name = "bath_type",
             c_type = "std::string",
             initializer = """ "flatband" """,
             doc = r"""""")

module.add_converter(c)

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