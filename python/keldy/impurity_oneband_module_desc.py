# Generated automatically using the command :
# c++2py ../../c++/keldy/impurity_oneband/compute_obs.hpp --members_read_only -N keldy -N keldy::warpers -N keldy::impurity_oneband -a keldy -m impurity_oneband_module -o impurity_oneband_module -C pytriqs --cxxflags="-std=c++17 " --includes ../../c++
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
#include <cpp2py/converters/variant.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/h5.hpp>
#include <triqs/cpp2py_converters/h5.hpp>

using namespace keldy;
using namespace keldy::warpers;
using namespace keldy::impurity_oneband;
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

# The class warper_plasma_uv_t
c = class_(
        py_type = "WarperPlasmaUvT",  # name of the python class
        c_type = "keldy::warpers::warper_plasma_uv_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(double t_max_)""", doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double operator() (std::vector<double> ui_vec)""",
             name = "__call__",
             doc = r"""""")

module.add_class(c)

# The class warper_product_1d_simple_t
c = class_(
        py_type = "WarperProduct1dSimpleT",  # name of the python class
        c_type = "keldy::warpers::warper_product_1d_simple_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(double t_max_)""", doc = r"""""")

c.add_constructor("""(std::function<double (double)> f1_, double t_max_, int nr_function_sample_points)""", doc = r"""""")

c.add_constructor("""(std::function<double (double)> f1_, std::function<double (double)> f1_integrated_, std::function<double (double)> f1_integrated_inverse_, double t_max_, int nr_function_sample_points)""", doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double operator() (std::vector<double> ui_vec)""",
             name = "__call__",
             doc = r"""""")

module.add_class(c)

# The class gf_index_t
c = class_(
        py_type = "GfIndexT",  # name of the python class
        c_type = "keldy::impurity_oneband::gf_index_t",   # name of the C++ class
        doc = r"""Index of the Keldysh Green Function""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "contour",
             c_type = "keldy::contour_pt_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "spin",
             c_type = "keldy::spin_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "orbital",
             c_type = "keldy::orbital_t",
             read_only= True,
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
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "make_dot_lead",
             c_type = "bool",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "g0_lesser",
             c_type = "block_gf<triqs::gfs::retime, triqs::gfs::matrix_valued>",
             read_only= True,
             doc = r"""Lesser Green function $G^{<}_{\sigma}(t)$; block spin $\sigma$ {up, down}""")

c.add_member(c_name = "lesser_ft_error",
             c_type = "gf<triqs::gfs::retime, triqs::gfs::matrix_valued>",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "g0_greater",
             c_type = "block_gf<triqs::gfs::retime, triqs::gfs::matrix_valued>",
             read_only= True,
             doc = r"""Greater Green function $G^{>}_{\sigma}(t)$; block spin $\sigma$ {up, down}""")

c.add_member(c_name = "greater_ft_error",
             c_type = "gf<triqs::gfs::retime, triqs::gfs::matrix_valued>",
             read_only= True,
             doc = r"""""")

c.add_constructor("""()""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_model_omega model_omega_, bool make_dot_lead_)""", doc = r"""""")

c.add_method("""std::string hdf5_scheme ()""",
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
             read_only= True,
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

# The class singleton_binner
c = class_(
        py_type = "SingletonBinner",  # name of the python class
        c_type = "keldy::impurity_oneband::singleton_binner",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "time",
             c_type = "double",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "value",
             c_type = "keldy::dcomplex",
             read_only= True,
             doc = r"""""")

module.add_class(c)

# The class binner_1d
c = class_(
        py_type = "Binner1d",  # name of the python class
        c_type = "keldy::impurity_oneband::binner_1d",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""()""", doc = r"""""")

c.add_constructor("""(double t_min_, double t_max_, int n_bins_)""", doc = r"""""")

c.add_method("""triqs::arrays::array<std::__1::complex<double>, 1> get_values ()""",
             doc = r"""""")

c.add_method("""triqs::arrays::array<unsigned long long, 1> get_nr_values ()""",
             doc = r"""""")

c.add_method("""triqs::arrays::array<double, 1> get_bin_times ()""",
             doc = r"""""")

c.add_method("""double get_bin_size ()""",
             doc = r"""""")

c.add_method("""int get_nr_point_dropped ()""",
             doc = r"""""")

c.add_method("""void accumulate (double time, keldy::dcomplex value)""",
             doc = r"""Includes boundary points, so t_min <= t <= t_max. t_max gets put in last bin""")

module.add_class(c)

# The class integrand_g_direct_time
c = class_(
        py_type = "IntegrandGDirectTime",  # name of the python class
        c_type = "keldy::impurity_oneband::integrand_g_direct_time",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(keldy::impurity_oneband::g0_keldysh_contour_t g0_, keldy::impurity_oneband::gf_index_t external_A_, keldy::impurity_oneband::gf_index_t external_B_, double cutoff_ = 0.)""", doc = r"""""")

c.add_method("""std::pair<keldy::impurity_oneband::integrand_g_direct_time::result_t, int> operator() (std::vector<double> times)""",
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

c.add_method("""void bin_data (std::pair<gf_index_t, dcomplex> in)""",
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

c.add_method("""triqs::arrays::array<std::complex<double>, 2> get_values ()""",
             doc = r"""""")

c.add_method("""triqs::arrays::array<unsigned long long, 2> get_nr_values ()""",
             doc = r"""""")

c.add_method("""triqs::arrays::array<double, 1> get_bin_times ()""",
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

c.add_method("""std::pair<keldy::impurity_oneband::integrand_g_kernel::result_t, int> operator() (std::vector<double> times)""",
             name = "__call__",
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

# The class warper_identity_t
c = class_(
        py_type = "WarperIdentityT",  # name of the python class
        c_type = "keldy::warpers::warper_identity_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""()""", doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double operator() (std::vector<double> ui_vec)""",
             name = "__call__",
             doc = r"""""")

module.add_class(c)

# The class warper_product_1d_t
c = class_(
        py_type = "WarperProduct1dT",  # name of the python class
        c_type = "keldy::warpers::warper_product_1d_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(double t_max_)""", doc = r"""""")

c.add_constructor("""(std::vector<std::function<double (double)> > fn_, double t_max_, int nr_function_sample_points)""", doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double operator() (std::vector<double> ui_vec)""",
             name = "__call__",
             doc = r"""""")

module.add_class(c)

# The class hist_xi
c = class_(
        py_type = "HistXi",  # name of the python class
        c_type = "keldy::warpers::hist_xi",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "bins",
             c_type = "std::vector<double>",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "values",
             c_type = "std::vector<double>",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "counts",
             c_type = "std::vector<int>",
             read_only= True,
             doc = r"""""")

module.add_class(c)

# The class warper_plasma_projection_t
c = class_(
        py_type = "WarperPlasmaProjectionT",  # name of the python class
        c_type = "keldy::warpers::warper_plasma_projection_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(std::function<double (std::vector<double>)> integrand_, keldy::warpers::warper_product_1d_t w_warper_, int order, int nr_function_sample_points)""", doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double evaluate_warping_function (std::vector<double> ui_vec)""",
             doc = r"""""")

module.add_class(c)

# The class warper_train_t
c = class_(
        py_type = "WarperTrainT",  # name of the python class
        c_type = "keldy::warpers::warper_train_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_member(c_name = "warpers",
             c_type = "std::vector<warper_variant>",
             read_only= True,
             doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double operator() (std::vector<double> ui_vec)""",
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

c.add_constructor("""(keldy::impurity_oneband::g0_model model, double time, int order, double cutoff_integrand, std::string warper_function_name, int nr_sample_points_warper, double warper_scale = 1)""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, int order, double cutoff_integrand, std::string warper_function_name, int nr_sample_points_warper, double warper_scale = 1)""", doc = r"""""")

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

# The class compute_charge_Q_direct_plasma_1D
c = class_(
        py_type = "ComputeChargeQDirectPlasma1d",  # name of the python class
        c_type = "keldy::impurity_oneband::compute_charge_Q_direct_plasma_1D",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, int order, std::vector<std::function<double (double)> > fn_, int nr_sample_points_warper)""", doc = r"""""")

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

c.add_constructor("""(keldy::impurity_oneband::g0_model model, double time, int order, double cutoff_integrand, std::string warper_function_name, int nr_sample_points_warper, double warper_scale = 1)""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, int order, double cutoff_integrand, std::string warper_function_name, int nr_sample_points_warper, double warper_scale = 1)""", doc = r"""""")

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

c.add_constructor("""(keldy::impurity_oneband::g0_model model, double time, int order, int nr_time_slices, double cutoff_integrand, std::string warper_function_name, int nr_sample_points_warper, double warper_scale = 1)""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, int order, int nr_time_slices, double cutoff_integrand, std::string warper_function_name, int nr_sample_points_warper, double warper_scale = 1)""", doc = r"""""")

c.add_method("""void run (int nr_steps)""",
             doc = r"""""")

c.add_method("""void reset_rng (std::string rng_name, int rng_state_seed, bool do_shift = false, bool do_scramble = false, int rng_seed_shift = 0)""",
             doc = r"""""")

c.add_method("""keldy::impurity_oneband::binner_1d reduce_result ()""",
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

c.add_constructor("""(keldy::impurity_oneband::model_param_t params, double time, int order, std::string warper_function_name, int nr_sample_points_warper, int nr_bins = 100)""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_model model, double time, int order, std::string warper_function_name, int nr_sample_points_warper, double warper_scale, int nr_bins = 100)""", doc = r"""""")

c.add_constructor("""(keldy::impurity_oneband::g0_model model, double time, int order, std::string warper_function_name, bool alternate, int nr_sample_points_warper, double warper_scale, int nb_bins = 100)""", doc = r"""""")

c.add_method("""void run (int nr_steps)""",
             doc = r"""""")

c.add_method("""void reset_rng (std::string rng_name, int rng_state_seed, bool do_shift = false, bool do_scramble = false, int rng_seed_shift = 0)""",
             doc = r"""""")

c.add_method("""kernel_binner reduce_result ()""",
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

module.add_function ("int keldy::compare_3way (keldy::contour_pt_t a, keldy::contour_pt_t b)", doc = r"""""")

module.add_function ("std::vector<double> keldy::warpers::vi_from_ui (double t_max, std::vector<double> u_times)", doc = r"""""")

module.add_function ("std::vector<double> keldy::warpers::ui_from_vi (double t_max, std::vector<double> v_times)", doc = r"""""")

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

module.add_function ("void keldy::impurity_oneband::h5_write (triqs::h5::group h5group, std::string subgroup_name, keldy::impurity_oneband::model_param_t c)", doc = r"""""")

module.add_function ("void keldy::impurity_oneband::h5_read (triqs::h5::group h5group, std::string subgroup_name, keldy::impurity_oneband::model_param_t c)", doc = r"""""")

module.add_function ("bool keldy::impurity_oneband::equivalent_without_timesplit (keldy::impurity_oneband::gf_index_t lhs, keldy::impurity_oneband::gf_index_t rhs)", doc = r"""""")

module.add_function ("void keldy::warpers::bin_values (keldy::warpers::hist_xi xi, int axis, std::vector<std::vector<double> > points, std::vector<double> values)", doc = r"""""")

module.add_function ("void keldy::warpers::convolve (std::vector<double> signal, std::vector<double> window)", doc = r"""""")

module.add_function ("warpers::warper_product_1d_t keldy::impurity_oneband::alternate_product_plasma_warper_factory (std::string label, int order, double time, int nr_sample_points_warper, double warper_scale)", doc = r"""""")

module.add_function ("warpers::warper_product_1d_simple_t keldy::impurity_oneband::simple_plasma_warper_factory (std::string label, keldy::impurity_oneband::integrand_g_direct f, double time, int nr_sample_points_warper, double warper_scale)", doc = r"""""")

module.add_function ("warpers::warper_product_1d_simple_t keldy::impurity_oneband::simple_plasma_warper_factory (std::string label, double time, int nr_sample_points_warper, double warper_scale)", doc = r"""""")

module.add_function ("warpers::warper_product_1d_simple_t keldy::impurity_oneband::simple_plasma_warper_factory_kernel (std::string label, keldy::impurity_oneband::integrand_g_kernel f, double time, int nr_sample_points_warper, double warper_scale)", doc = r"""""")


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