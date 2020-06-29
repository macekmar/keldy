# Generated automatically using the command :
# c++2py ../../c++/keldy/warpers/warpers.hpp -a keldy -m warpers -o warpers  -C pytriqs -C keldy --cxxflags="-std=c++17 " --includes ../../c++ -N keldy::warpers
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "warpers", doc = r"", app_name = "keldy")

# Imports
module.add_imports(*['keldy.binner', 'keldy.common', 'pytriqs.gf'])

# Add here all includes
module.add_include("keldy/warpers/warpers.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <cpp2py/converters/function.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/gf.hpp>

using namespace keldy::warpers;
""")


# The class warper_identity_t
c = class_(
        py_type = "WarperIdentityT",  # name of the python class
        c_type = "keldy::warpers::warper_identity_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""()""", doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_forward (std::vector<double> ui_vec)""",
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

c.add_method("""std::pair<std::vector<double>,double> map_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

module.add_class(c)

# The class warper_product_1d_simple_t
c = class_(
        py_type = "WarperProduct1dSimpleT",  # name of the python class
        c_type = "keldy::warpers::warper_product_1d_simple_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""()""", doc = r"""""")

c.add_constructor("""(std::function<double(double)> f1_, std::function<double(double)> f1_integrated_, std::function<double(double)> f1_integrated_inverse_, double domain_u_max_, bool do_domain_checks = true)""", doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double operator() (std::vector<double> ui_vec)""",
             name = "__call__",
             doc = r"""""")

module.add_class(c)

# The class warper_product_1d_simple_interp_nearest_t
c = class_(
        py_type = "WarperProduct1dSimpleInterpNearestT",  # name of the python class
        c_type = "keldy::warpers::warper_product_1d_simple_interp_nearest_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""()""", doc = r"""""")

c.add_constructor("""(std::function<double(double)> f1_, double domain_u_max_, int nr_sample_points_)""", doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double operator() (std::vector<double> ui_vec)""",
             name = "__call__",
             doc = r"""""")

module.add_class(c)

# The class warper_product_1d_simple_interp_hybrid_t
c = class_(
        py_type = "WarperProduct1dSimpleInterpHybridT",  # name of the python class
        c_type = "keldy::warpers::warper_product_1d_simple_interp_hybrid_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""()""", doc = r"""""")

c.add_constructor("""(std::function<double(double)> f1_, double domain_u_max_, int nr_sample_points_)""", doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_forward (std::vector<double> ui_vec)""",
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

c.add_constructor("""()""", doc = r"""""")

c.add_method("""int size ()""",
             doc = r"""""")

c.add_method("""void emplace_back (keldy::warpers::warper_product_1d_simple_t w)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double operator() (std::vector<double> ui_vec)""",
             name = "__call__",
             doc = r"""""")

module.add_class(c)

# The class warper_product_1d_interp_nearest_t
c = class_(
        py_type = "WarperProduct1dInterpNearestT",  # name of the python class
        c_type = "keldy::warpers::warper_product_1d_interp_nearest_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""()""", doc = r"""""")

c.add_method("""int size ()""",
             doc = r"""""")

c.add_method("""void emplace_back (keldy::warpers::warper_product_1d_simple_interp_nearest_t w)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double operator() (std::vector<double> ui_vec)""",
             name = "__call__",
             doc = r"""""")

module.add_class(c)

# The class warper_product_1d_interp_hybrid_t
c = class_(
        py_type = "WarperProduct1dInterpHybridT",  # name of the python class
        c_type = "keldy::warpers::warper_product_1d_interp_hybrid_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""()""", doc = r"""""")

c.add_method("""int size ()""",
             doc = r"""""")

c.add_method("""void emplace_back (keldy::warpers::warper_product_1d_simple_interp_hybrid_t w)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double operator() (std::vector<double> ui_vec)""",
             name = "__call__",
             doc = r"""""")

module.add_class(c)

# The class warper_projection_t
c = class_(
        py_type = "WarperProjectionT",  # name of the python class
        c_type = "keldy::warpers::warper_projection_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(std::function<dcomplex(std::vector<double>)> warped_integrand, int order, int num_bins, int nr_samples, double sigma, bool optimize_sigma = true)""", doc = r"""""")

c.add_method("""keldy::binner::binner_t<1,0,double> get_xi (int axis)""",
             doc = r"""""")

c.add_method("""std::vector<double> get_sigmas ()""",
             doc = r"""""")

c.add_method("""triqs::gfs::gf<triqs::gfs::retime,triqs::gfs::scalar_real_valued> get_fi (int axis)""",
             doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

module.add_class(c)

# The class warper_train_t
c = class_(
        py_type = "WarperTrainT",  # name of the python class
        c_type = "keldy::warpers::warper_train_t",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""()""", doc = r"""""")

c.add_method("""void emplace_back (keldy::warpers::warper_identity_t w)""",
             doc = r"""""")

c.add_method("""void emplace_back (keldy::warpers::warper_plasma_uv_t w)""",
             doc = r"""""")

c.add_method("""void emplace_back (keldy::warpers::warper_product_1d_simple_t w)""",
             doc = r"""""")

c.add_method("""void emplace_back (keldy::warpers::warper_product_1d_simple_interp_nearest_t w)""",
             doc = r"""""")

c.add_method("""void emplace_back (keldy::warpers::warper_product_1d_t w)""",
             doc = r"""""")

c.add_method("""void emplace_back (keldy::warpers::warper_product_1d_interp_nearest_t w)""",
             doc = r"""""")

c.add_method("""void emplace_back (keldy::warpers::warper_product_1d_interp_hybrid_t w)""",
             doc = r"""""")

c.add_method("""void emplace_back (keldy::warpers::warper_projection_t w)""",
             doc = r"""""")

c.add_method("""void clear ()""",
             doc = r"""""")

c.add_method("""size_t size ()""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_reverse (std::vector<double> li_vec, int start_domain_nr, int end_domain_nr)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_forward (std::vector<double> ui_vec, int start_domain_nr, int end_domain_nr)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::pair<std::vector<double>,double> map_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_forward (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

module.add_class(c)

module.add_function ("std::vector<double> keldy::warpers::vi_from_ui (double t_max, std::vector<double> u_times)", doc = r"""""")

module.add_function ("std::vector<double> keldy::warpers::ui_from_vi (double t_max, std::vector<double> v_times)", doc = r"""""")

module.add_function ("keldy::warpers::warper_product_1d_simple_t keldy::warpers::make_product_1d_simple_exponential_nointerp (double domain_u_max, double w_scale)", doc = r"""""")

module.add_function ("keldy::warpers::warper_product_1d_simple_t keldy::warpers::make_product_1d_simple_inverse_nointerp (double domain_u_max, double w_scale)", doc = r"""""")

module.add_function ("keldy::warpers::warper_product_1d_simple_t keldy::warpers::make_product_1d_simple_inverse_square_nointerp (double domain_u_max, double w_scale)", doc = r"""""")

module.add_function ("keldy::warpers::warper_product_1d_t keldy::warpers::make_product_1d_inverse_cube_alternate (int order, double time, double warper_scale)", doc = r"""""")

module.add_function ("keldy::warpers::warper_product_1d_interp_nearest_t keldy::warpers::make_product_1d_inverse_cube_alternate_interp_nearest (int order, double time, double warper_scale, int nr_sample_points_warper)", doc = r"""""")

module.add_function ("keldy::warpers::warper_product_1d_interp_hybrid_t keldy::warpers::make_product_1d_inverse_cube_alternate_interp_hybrid (int order, double time, double warper_scale, int nr_sample_points_warper)", doc = r"""""")



module.generate_code()