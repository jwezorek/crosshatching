#pragma once

#include <map>
#include "geometry.hpp"
#include "util.hpp"
#include <opencv2/core.hpp>
#include <range/v3/all.hpp>
#include <tuple>
#include <functional>
#include <variant>
#include <vector>
#include <memory>
#include <unordered_map>
#include <string>

namespace ch {

    struct stroke {
        ranges::any_view<ch::point> polyline;
        double pen_thickness;
    };
    using strokes = ranges::any_view<stroke>;

    struct drawing2 {
        strokes content;
        dimensions<int> size;
    };

    struct brush_context {
        polygon poly;
        std::unordered_map<std::string, double> variables;
    };

    using brush_expr_nil = std::monostate;
    using brush_expr_value = std::variant<brush_expr_nil, double, strokes>;

    class brush_expr2 {
    public:
        virtual brush_expr_value eval(brush_context& ctxt) = 0;
    };

    void debug_brushes();
}