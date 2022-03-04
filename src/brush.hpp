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

    ch::hatching_range run_brush_pipeline(const brush_pipeline& pipeline, double t);

}