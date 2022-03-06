#pragma once

#include "crosshatch.hpp"
#include <functional>
#include <variant>
#include <vector>

namespace ch {

    using brush_fn = std::function<ch::hatching_range(double t)>;
    using param_adapter_fn = std::function<double (double t)>;
    using brush_adapter_fn = std::function<ch::hatching_range(ch::hatching_range, double t)>;
    using brush_pipeline_item = std::variant<brush_fn, param_adapter_fn, brush_adapter_fn>;
    using brush_pipeline = std::vector<brush_pipeline_item>;
    using param_unit_of_hatching_fn = std::function<hatching_range(double, double, double, double, double)>;

    ch::hatching_range run_brush_pipeline(const brush_pipeline& pipeline, double t);

    param_adapter_fn make_constant_fn(double k);
    param_adapter_fn make_lerp_fn(double v1, double v2);
    param_adapter_fn make_lerped_normal_dist_fn(double mu1, double sigma1, double mu2, double sigma2);

    param_unit_of_hatching_fn make_default_hatching_unit();
    brush_adapter_fn make_one_param_brush_adaptor_fn(std::function< ch::hatching_range(ch::hatching_range, double)> fn, param_adapter_fn param);

    brush_fn make_linear_hatching_brush_fn(const param_adapter_fn& run_length, const param_adapter_fn& space_length, const param_adapter_fn& vert_space,
        const param_unit_of_hatching_fn& h_fn);

}