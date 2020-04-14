# Generated automatically using the command :
# c++2py ../../c++/keldy/warpers/warpers.hpp -a keldy -m warpers -o warpers  -C pytriqs -C keldy --cxxflags="-std=c++17 " --includes ../../c++ --members_read_only -N keldy::warpers
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "warpers", doc = r"", app_name = "keldy")

# Imports
module.add_imports(*['keldy.common'])

# Add here all includes
module.add_include("keldy/warpers/warpers.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/function.hpp>
#include <cpp2py/converters/variant.hpp>
#include <cpp2py/converters/vector.hpp>

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

c.add_constructor("""(double t_max_)""", doc = r"""""")

c.add_constructor("""(std::function<double (double)> f1_, double t_max_, int nr_function_sample_points)""", doc = r"""""")

c.add_constructor("""(std::function<double (double)> f1_, std::function<double (double)> f1_integrated_, std::function<double (double)> f1_integrated_inverse_, double t_max_, int nr_function_sample_points)""", doc = r"""""")

c.add_method("""std::vector<double> ui_from_li (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""std::vector<double> li_from_ui (std::vector<double> ui_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_forward (std::vector<double> ui_vec)""",
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

c.add_method("""double jacobian_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_forward (std::vector<double> ui_vec)""",
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

c.add_method("""double jacobian_reverse (std::vector<double> li_vec)""",
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

c.add_method("""double jacobian_reverse (std::vector<double> li_vec)""",
             doc = r"""""")

c.add_method("""double jacobian_forward (std::vector<double> ui_vec)""",
             doc = r"""""")
             doc = r"""""")

module.add_class(c)

module.add_function ("std::vector<double> keldy::warpers::vi_from_ui (double t_max, std::vector<double> u_times)", doc = r"""""")

module.add_function ("std::vector<double> keldy::warpers::ui_from_vi (double t_max, std::vector<double> v_times)", doc = r"""""")

module.add_function ("void keldy::warpers::bin_values (keldy::warpers::hist_xi xi, int axis, std::vector<std::vector<double> > points, std::vector<double> values)", doc = r"""""")

module.add_function ("void keldy::warpers::convolve (std::vector<double> signal, std::vector<double> window)", doc = r"""""")



module.generate_code()