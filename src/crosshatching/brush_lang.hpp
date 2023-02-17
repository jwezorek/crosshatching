#pragma once

#include <map>
#include "geometry.hpp"
#include "util.hpp"
#include "strokes.hpp"
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

/*------------------------------------------------------------------------------------------------*/

namespace ch {

    class brush_expr;
    using brush_expr_ptr = std::shared_ptr<brush_expr>;
    
    using nil_value = std::monostate;
    using brush_expr_value = 
        std::variant<nil_value, double, drawing_comps_ptr, random_func, brush_expr_ptr>;
    using variables_map = std::unordered_map<std::string, brush_expr_value>;

    struct brush_context {
        polygon poly;
        variables_map variables;
        drawing_comps_ptr strokes;

        brush_context(const polygon& poly, double param, double value = 0.5, double flow = 0.0);
        brush_context clone() const;
        
        template<typename T>
        T get(const std::string& var) const {
            return std::get<T>(variables.at(var));
        }
    };

    class brush_expr {
    public:
        brush_expr() {}
        void set_children(std::span<const brush_expr_ptr> children);
        virtual brush_expr_value eval(brush_context& ctxt) = 0;
        virtual std::string short_string() const = 0;
        const std::vector<brush_expr_ptr>& children() const;
        void replace_child(const brush_expr_ptr& old_child, const brush_expr_ptr& new_child);
    protected:
        std::vector<brush_expr_ptr> children_;
    };

    std::variant<brush_expr_ptr, std::runtime_error> parse(const std::string& str);
    drawing_comps_ptr brush_expr_to_strokes(
        const brush_expr_ptr& expr, const brush_context& br_ctxt);
    std::string pretty_print(const brush_expr& expr);
    void debug_brushes();
}