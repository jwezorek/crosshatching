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
#include <stdexcept>

namespace ch {

    struct stroke {
        ranges::any_view<ch::point> polyline;
        double pen_thickness;
    };
    using strokes = ranges::any_view<stroke>;

    class brush_expr;
    using brush_expr_ptr = std::shared_ptr<brush_expr>;
    
    using nil_value = std::monostate;
    using brush_expr_value = std::variant<nil_value, double, strokes, random_func, brush_expr_ptr>;
    using variables_map = std::unordered_map<std::string, brush_expr_value>;

    struct brush_context {
        polygon poly;
        variables_map variables;
        std::optional<ch::strokes> strokes;
    };

    class brush_expr {
    public:
        brush_expr() {}
        void set_children(std::span<const brush_expr_ptr> children);
        virtual brush_expr_value eval(brush_context& ctxt) = 0;
    protected:
        std::vector<brush_expr_ptr> children_;
    };

    std::variant<brush_expr_ptr, std::runtime_error> parse(const std::string& str);

    ranges::any_view<ch::point> transform(ranges::any_view<ch::point> poly, const ch::matrix& mat);
    stroke transform(stroke s, const ch::matrix& mat);
    strokes transform(strokes s, const ch::matrix& mat);
    void debug_brushes();
}