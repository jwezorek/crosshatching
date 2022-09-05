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
#include <optional>
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

    using nil_value = std::monostate;
    using brush_expr_value = std::variant<nil_value, double, strokes, random_func>;
    using variables_map = std::unordered_map<std::string, brush_expr_value>;

    struct brush_context {
        polygon poly;
        cbrng_seed seed;
        variables_map variables;
        std::optional<ch::strokes> strokes;
    };

    class brush_expr;
    using brush_expr_ptr = std::shared_ptr<brush_expr>;

    class brush_expr {
    public:
        brush_expr() {}
        brush_expr(std::span<const brush_expr_ptr> children);
        virtual brush_expr_value eval(brush_context& ctxt) = 0;
    protected:
        std::vector<brush_expr_ptr> children_;
    };
    ranges::any_view<ch::point> transform(ranges::any_view<ch::point> poly, const ch::matrix& mat);
    stroke transform(stroke s, const ch::matrix& mat);
    strokes transform(strokes s, const ch::matrix& mat);
    void debug_brushes();
}