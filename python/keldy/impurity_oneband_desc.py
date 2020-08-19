# Generated automatically using the command :
# c++2py ../../c++/keldy/impurity_oneband/compute_obs.hpp -N keldy::impurity_oneband -a keldy -m impurity_oneband -o impurity_oneband -C triqs -C keldy --cxxflags="-std=c++17 " --includes ../../c++
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "impurity_oneband", doc = r"", app_name = "keldy")

# Imports
module.add_imports(*['keldy.common', 'keldy.warpers', 'triqs.gf'])

# Add here all includes
module.add_include("keldy/impurity_oneband/compute_obs.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>

using namespace keldy::impurity_oneband;
""")


# The class gf_index_t
c = class_(
        py_type = "GfIndexT",  # name of the python class
        c_type = "keldy::impurity_oneband::gf_index_t",   # name of the C++ class
        doc = r"""Index of the Keldysh Green Function""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "contour",
             c_type = "keldy::contour_pt_t",
             read_only= False,
             doc = r"""""")

c.add_member(c_name = "spin",
             c_type = "keldy::spin_t",
             read_only= False,
             doc = r"""""")

c.add_member(c_name = "orbital",
             c_type = "keldy::orbital_t",
             read_only= False,
             doc = r"""""")

c.add_constructor("""()""", doc = r"""Constructor: (time, spin, keldysh_idx)""")

c.add_constructor("""(keldy::time_real_t time_, int spin_, int k_idx_)""", doc = r"""""")

c.add_constructor("""(keldy::time_real_t time_, int spin_, int k_idx_, int timesplit_n_)""", doc = r"""""")

c.add_constructor("""(keldy::time_real_t time_, int spin_, int k_idx_, int timesplit_n_, int orbital_)""", doc = r"""""")

module.add_class(c)

# The class g0_model_omega
c = class_(
        py_type = "G0ModelOmega",  # name of the python class
        c_type = "keldy::impurity_oneband::g0_model_omega",   # name of the C++ class
        doc = r"""Defines model throuh non-interacting Green function g_lesser / g_greater""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""()""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::model_param_t parameters)""", doc = r"""""")

c.add_method("""keldy::impurity_oneband::model_param_t get_param ()""",
             doc = r"""""")

c.add_method("""double get_param_alpha ()""",
             doc = r"""""")

c.add_method("""double mu_left ()""",
             doc = r"""""")

c.add_method("""double mu_right ()""",
             doc = r"""""")

c.add_method("""keldy::dcomplex g0_dot_R (keldy::dcomplex omega)""",
             doc = r"""""")

c.add_method("""keldy::dcomplex g0_dot_A (keldy::dcomplex omega)""",
             doc = r"""""")

c.add_method("""keldy::dcomplex g0_dot_lesser (keldy::dcomplex omega)""",
             doc = r"""""")

c.add_method("""keldy::dcomplex g0_dot_greater (keldy::dcomplex omega)""",
             doc = r"""""")

c.add_method("""keldy::dcomplex bath_hybrid_R_left (keldy::dcomplex omega)""",
             doc = r"""""")

c.add_method("""keldy::dcomplex bath_hybrid_A_left (keldy::dcomplex omega)""",
             doc = r"""""")

c.add_method("""keldy::dcomplex bath_hybrid_K_left (keldy::dcomplex omega)""",
             doc = r"""""")

c.add_method("""keldy::dcomplex bath_hybrid_R_right (keldy::dcomplex omega)""",
             doc = r"""""")

c.add_method("""keldy::dcomplex bath_hybrid_A_right (keldy::dcomplex omega)""",
             doc = r"""""")

c.add_method("""keldy::dcomplex bath_hybrid_K_right (keldy::dcomplex omega)""",
             doc = r"""""")

c.add_method("""keldy::dcomplex g0_rightlead_dot_lesser (keldy::dcomplex omega)""",
             doc = r"""""")

c.add_method("""keldy::dcomplex g0_rightlead_dot_greater (keldy::dcomplex omega)""",
             doc = r"""""")

module.add_class(c)

# The class g0_model
c = class_(
        py_type = "G0Model",  # name of the python class
        c_type = "keldy::impurity_oneband::g0_model",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = True,
)

c.add_member(c_name = "model_omega",
             c_type = "keldy::impurity_oneband::g0_model_omega",
             read_only= False,
             doc = r"""""")

c.add_member(c_name = "make_dot_lead",
             c_type = "bool",
             read_only= False,
             doc = r"""""")

c.add_member(c_name = "g0_lesser",
             c_type = "block_gf<triqs::gfs::retime, triqs::gfs::matrix_valued>",
             read_only= False,
             doc = r"""Lesser Green function $G^{<}_{\sigma}(t)$; block spin $\sigma$ {up, down}""")

c.add_member(c_name = "lesser_ft_error",
             c_type = "gf<triqs::gfs::retime, triqs::gfs::matrix_valued>",
             read_only= False,
             doc = r"""""")

c.add_member(c_name = "g0_greater",
             c_type = "block_gf<triqs::gfs::retime, triqs::gfs::matrix_valued>",
             read_only= False,
             doc = r"""Greater Green function $G^{>}_{\sigma}(t)$; block spin $\sigma$ {up, down}""")

c.add_member(c_name = "greater_ft_error",
             c_type = "gf<triqs::gfs::retime, triqs::gfs::matrix_valued>",
             read_only= False,
             doc = r"""""")

c.add_constructor("""()""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_model_omega model_omega_, bool make_dot_lead_)""", doc = r"""""")

c.add_method("""std::string hdf5_format ()""",
             is_static = True,
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
             read_only= False,
             doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_model model_)""", doc = r"""""")

c.add_method("""keldy::dcomplex operator() (keldy::impurity_oneband::gf_index_t a, keldy::impurity_oneband::gf_index_t b, bool internal_point = true)""",
             name = "__call__",
             doc = r"""return :math:`g^{ab}(t,t')` in contour basis from :math:`g^<, g^>` functions""")

c.add_method("""double get_time_max ()""",
             doc = r"""""")

module.add_class(c)

# The class integrand_g_direct
c = class_(
        py_type = "IntegrandGDirect",  # name of the python class
        c_type = "keldy::impurity_oneband::integrand_g_direct",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(keldy::impurity_oneband::g0_keldysh_contour_t g0_, keldy::impurity_oneband::gf_index_t external_A_, keldy::impurity_oneband::gf_index_t external_B_, double cutoff_ = 0.)""", doc = r"""""")

c.add_method("""std::pair<keldy::impurity_oneband::integrand_g_direct::result_t, int> operator() (std::vector<double> times, bool keep_u_hypercube = true)""",
             name = "__call__",
             doc = r"""Returns integrand for the specified times""")

module.add_class(c)

# The class integrand_g_direct_time
c = class_(
        py_type = "IntegrandGDirectTime",  # name of the python class
        c_type = "keldy::impurity_oneband::integrand_g_direct_time",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(keldy::impurity_oneband::g0_keldysh_contour_t g0_, keldy::impurity_oneband::gf_index_t external_A_, keldy::impurity_oneband::gf_index_t external_B_, double cutoff_ = 0.)""", doc = r"""""")

c.add_method("""std::pair<keldy::impurity_oneband::integrand_g_direct_time::result_t, int> operator() (std::vector<double> times, bool keep_u_hypercube = true)""",
             name = "__call__",
             doc = r"""Returns integrand for the specified times""")

module.add_class(c)

# The class integrand_g_kernel
c = class_(
        py_type = "IntegrandGKernel",  # name of the python class
        c_type = "keldy::impurity_oneband::integrand_g_kernel",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(keldy::impurity_oneband::g0_keldysh_contour_t g0_, keldy::impurity_oneband::gf_index_t g_idx_X_)""", doc = r"""""")

c.add_method("""std::pair<keldy::impurity_oneband::integrand_g_kernel::result_t, int> operator() (std::vector<double> times, bool keep_u_hypercube = true)""",
             name = "__call__",
             doc = r"""""")

# The class integrand_g_kernel_single_omega
c = class_(
        py_type = "IntegrandGKernelSingleOmega",  # name of the python class
        c_type = "keldy::impurity_oneband::integrand_g_kernel_single_omega",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(keldy::impurity_oneband::g0_keldysh_contour_t g0_, keldy::impurity_oneband::gf_index_t g_idx_X_, double omega_)""", doc = r"""""")

c.add_method("""std::pair<keldy::impurity_oneband::integrand_g_kernel_single_omega::result_t, int> operator() (std::vector<double> times, bool keep_u_hypercube = true)""",
             name = "__call__",
             doc = r"""""")

module.add_class(c)

# The class compute_charge_Q_direct
c = class_(
        py_type = "ComputeChargeQDirect",  # name of the python class
        c_type = "keldy::impurity_oneband::compute_charge_Q_direct",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "warper",
             c_type = "keldy::warpers::warper_train_t",
             read_only= False,
             doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_model model, double time, int order, double cutoff_integrand)""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, int order, double cutoff_integrand)""", doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_direct::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec, int start_domain_nr, bool keep_u_hypercube = true)""",
             doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_direct::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""void run (int nr_steps)""",
             doc = r"""""")

c.add_method("""void reset_rng (std::string rng_name, int rng_state_seed, bool do_shift = false, bool do_scramble = false, int rng_seed_shift = 0)""",
             doc = r"""""")

c.add_method("""dcomplex reduce_result ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_run ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_in_domain ()""",
             doc = r"""""")

c.add_method("""keldy::impurity_oneband::integrand_g_direct get_integrand ()""",
             doc = r"""""")

c.add_method("""keldy::warpers::warper_train_t get_warper ()""",
             doc = r"""""")

module.add_class(c)

# The class compute_current_J_direct
c = class_(
        py_type = "ComputeCurrentJDirect",  # name of the python class
        c_type = "keldy::impurity_oneband::compute_current_J_direct",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "warper",
             c_type = "keldy::warpers::warper_train_t",
             read_only= False,
             doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_model model, double time, int order, double cutoff_integrand)""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, int order, double cutoff_integrand)""", doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_direct::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec, int start_domain_nr, bool keep_u_hypercube = true)""",
             doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_direct::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""void run (int nr_steps)""",
             doc = r"""""")

c.add_method("""void reset_rng (std::string rng_name, int rng_state_seed, bool do_shift = false, bool do_scramble = false, int rng_seed_shift = 0)""",
             doc = r"""""")

c.add_method("""dcomplex reduce_result ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_run ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_in_domain ()""",
             doc = r"""""")

c.add_method("""keldy::impurity_oneband::integrand_g_direct get_integrand ()""",
             doc = r"""""")

c.add_method("""keldy::warpers::warper_train_t get_warper ()""",
             doc = r"""""")

module.add_class(c)

# The class compute_direct
c = class_(
        py_type = "ComputeDirect",  # name of the python class
        c_type = "keldy::impurity_oneband::compute_direct",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "warper",
             c_type = "keldy::warpers::warper_train_t",
             read_only= False,
             doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_model model, gf_index_t id_a, gf_index_t id_b, double time, int order, double cutoff_integrand)""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, gf_index_t id_a, gf_index_t id_b, double time, int order, double cutoff_integrand)""", doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_direct::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec, int start_domain_nr, bool keep_u_hypercube = true)""",
             doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_direct::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""void run (int nr_steps)""",
             doc = r"""""")

c.add_method("""void reset_rng (std::string rng_name, int rng_state_seed, bool do_shift = false, bool do_scramble = false, int rng_seed_shift = 0)""",
             doc = r"""""")

c.add_method("""dcomplex reduce_result ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_run ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_in_domain ()""",
             doc = r"""""")

c.add_method("""keldy::impurity_oneband::integrand_g_direct get_integrand ()""",
             doc = r"""""")

c.add_method("""keldy::warpers::warper_train_t get_warper ()""",
             doc = r"""""")

module.add_class(c)

# The class compute_charge_Q_direct_time
c = class_(
        py_type = "ComputeChargeQDirectTime",  # name of the python class
        c_type = "keldy::impurity_oneband::compute_charge_Q_direct_time",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "warper",
             c_type = "keldy::warpers::warper_train_t",
             read_only= False,
             doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_model model, double time, int order, int nr_time_slices, double cutoff_integrand)""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, int order, int nr_time_slices, double cutoff_integrand)""", doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_direct_time::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec, int start_domain_nr, bool keep_u_hypercube = true)""",
             doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_direct_time::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""void run (int nr_steps)""",
             doc = r"""""")

c.add_method("""void reset_rng (std::string rng_name, int rng_state_seed, bool do_shift = false, bool do_scramble = false, int rng_seed_shift = 0)""",
             doc = r"""""")

c.add_method("""keldy::binner::binner_t<1,0> reduce_result ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_run ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_in_domain ()""",
             doc = r"""""")

c.add_method("""keldy::warpers::warper_train_t get_warper ()""",
             doc = r"""""")

c.add_method("""keldy::impurity_oneband::integrand_g_direct_time get_integrand ()""",
             doc = r"""""")

module.add_class(c)

# The class compute_current_J_direct_time
c = class_(
        py_type = "ComputeCurrentJDirectTime",  # name of the python class
        c_type = "keldy::impurity_oneband::compute_current_J_direct_time",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "warper",
             c_type = "keldy::warpers::warper_train_t",
             read_only= False,
             doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_model model, double time, int order, int nr_time_slices, double cutoff_integrand)""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, int order, int nr_time_slices, double cutoff_integrand)""", doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_direct_time::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec, int start_domain_nr, bool keep_u_hypercube = true)""",
             doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_direct_time::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""void run (int nr_steps)""",
             doc = r"""""")

c.add_method("""void reset_rng (std::string rng_name, int rng_state_seed, bool do_shift = false, bool do_scramble = false, int rng_seed_shift = 0)""",
             doc = r"""""")

c.add_method("""keldy::binner::binner_t<1,0> reduce_result ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_run ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_in_domain ()""",
             doc = r"""""")

c.add_method("""keldy::warpers::warper_train_t get_warper ()""",
             doc = r"""""")

c.add_method("""keldy::impurity_oneband::integrand_g_direct_time get_integrand ()""",
             doc = r"""""")

module.add_class(c)

# The class compute_direct_time
c = class_(
        py_type = "ComputeDirectTime",  # name of the python class
        c_type = "keldy::impurity_oneband::compute_direct_time",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "warper",
             c_type = "keldy::warpers::warper_train_t",
             read_only= False,
             doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_model model, gf_index_t id_a, gf_index_t id_b, double time, int order, int nr_time_slices, double cutoff_integrand)""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, gf_index_t id_a, gf_index_t id_b, double time, int order, int nr_time_slices, double cutoff_integrand)""", doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_direct_time::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec, int start_domain_nr, bool keep_u_hypercube = true)""",
             doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_direct_time::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""void run (int nr_steps)""",
             doc = r"""""")

c.add_method("""void reset_rng (std::string rng_name, int rng_state_seed, bool do_shift = false, bool do_scramble = false, int rng_seed_shift = 0)""",
             doc = r"""""")

c.add_method("""keldy::binner::binner_t<1,0> reduce_result ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_run ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_in_domain ()""",
             doc = r"""""")

c.add_method("""keldy::warpers::warper_train_t get_warper ()""",
             doc = r"""""")

c.add_method("""keldy::impurity_oneband::integrand_g_direct_time get_integrand ()""",
             doc = r"""""")

module.add_class(c)

# The class compute_gf_kernel
c = class_(
        py_type = "ComputeGfKernel",  # name of the python class
        c_type = "keldy::impurity_oneband::compute_gf_kernel",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "warper",
             c_type = "keldy::warpers::warper_train_t",
             read_only= False,
             doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_model model, double time, int order, int nr_bins = 100)""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, int order, int nr_bins = 100)""", doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_kernel::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec, int start_domain_nr, bool keep_u_hypercube = true)""",
             doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_kernel::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""void run (int nr_steps)""",
             doc = r"""""")

c.add_method("""void reset_rng (std::string rng_name, int rng_state_seed, bool do_shift = false, bool do_scramble = false, int rng_seed_shift = 0)""",
             doc = r"""""")

c.add_method("""keldy::binner::binner_t<1, 1> reduce_result ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_run ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_in_domain ()""",
             doc = r"""""")

c.add_method("""keldy::impurity_oneband::integrand_g_kernel get_integrand ()""",
             doc = r"""""")

c.add_method("""keldy::warpers::warper_train_t get_warper ()""",
             doc = r"""""")

module.add_class(c)

# The class compute_gf_kernel_single_omega
c = class_(
        py_type = "ComputeGfKernelSingleOmega",  # name of the python class
        c_type = "keldy::impurity_oneband::compute_gf_kernel_single_omega",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "warper",
             c_type = "keldy::warpers::warper_train_t",
             read_only= False,
             doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_model model, double time, double omega, int order)""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, double omega, int order)""", doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_kernel_single_omega::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec, int start_domain_nr, bool keep_u_hypercube = true)""",
             doc = r"""""")

c.add_method("""std::pair<typename keldy::impurity_oneband::integrand_g_kernel_single_omega::result_t, double> evaluate_warped_integrand (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""void run (int nr_steps)""",
             doc = r"""""")

c.add_method("""void reset_rng (std::string rng_name, int rng_state_seed, bool do_shift = false, bool do_scramble = false, int rng_seed_shift = 0)""",
             doc = r"""""")

c.add_method("""dcomplex reduce_result ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_run ()""",
             doc = r"""""")

c.add_method("""uint64_t reduce_nr_points_in_domain ()""",
             doc = r"""""")

c.add_method("""keldy::impurity_oneband::integrand_g_kernel_single_omega get_integrand ()""",
             doc = r"""""")

c.add_method("""keldy::warpers::warper_train_t get_warper ()""",
             doc = r"""""")

module.add_class(c)

module.add_function ("void keldy::impurity_oneband::fake (**keldy::impurity_oneband::model_param_t)", doc = r"""



+-------------------+-------------+--------------+---------------+
| Parameter Name    | Type        | Default      | Documentation |
+===================+=============+==============+===============+
| beta              | double      | 1.0          |               |
+-------------------+-------------+--------------+---------------+
| bias_V            | double      | 0.0          |               |
+-------------------+-------------+--------------+---------------+
| eps_d             | double      | 0.0          |               |
+-------------------+-------------+--------------+---------------+
| Gamma             | double      | 1.0          |               |
+-------------------+-------------+--------------+---------------+
| alpha             | double      | 0.0          |               |
+-------------------+-------------+--------------+---------------+
| half_bandwidth    | double      | 2.           |               |
+-------------------+-------------+--------------+---------------+
| bath_type         | std::string | "semicircle" |               |
+-------------------+-------------+--------------+---------------+
| time_max          | double      | +100.0       |               |
+-------------------+-------------+--------------+---------------+
| nr_time_points_gf | int         | 1000         |               |
+-------------------+-------------+--------------+---------------+
| ft_method         | std::string | "fft"        |               |
+-------------------+-------------+--------------+---------------+
""")

module.add_function ("void keldy::impurity_oneband::h5_write (h5::group h5group, std::string subgroup_name, keldy::impurity_oneband::model_param_t c)", doc = r"""""")

module.add_function ("void keldy::impurity_oneband::h5_read (h5::group h5group, std::string subgroup_name, keldy::impurity_oneband::model_param_t c)", doc = r"""""")

module.add_function ("bool keldy::impurity_oneband::equivalent_without_timesplit (keldy::impurity_oneband::gf_index_t lhs, keldy::impurity_oneband::gf_index_t rhs)", doc = r"""""")


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

c.add_member(c_name = "alpha",
             c_type = "double",
             initializer = """ 0.0 """,
             doc = r"""""")

c.add_member(c_name = "half_bandwidth",
             c_type = "double",
             initializer = """ 2. """,
             doc = r"""""")

c.add_member(c_name = "bath_type",
             c_type = "std::string",
             initializer = """ "semicircle" """,
             doc = r"""""")

c.add_member(c_name = "time_max",
             c_type = "double",
             initializer = """ +100.0 """,
             doc = r"""""")

c.add_member(c_name = "nr_time_points_gf",
             c_type = "int",
             initializer = """ 1000 """,
             doc = r"""""")

c.add_member(c_name = "ft_method",
             c_type = "std::string",
             initializer = """ "fft" """,
             doc = r"""""")

module.add_converter(c)


module.generate_code()
