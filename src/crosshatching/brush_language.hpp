#pragma once

#include "brush.hpp"
#include <vector>
#include <memory>
#include <string>
#include <optional>
#include <variant>
#include <iterator>

namespace ch {

    enum class symbol {
        true_,
        false_,
        pipe,
        linear_brush,
        norm_rnd,
        lerp,
        ramp,
        rotate,
        disintegrate,
        jiggle,
        merge
    };

    class brush_expr_base {
    public:
        virtual brush_pipeline_item eval() const = 0;
        virtual std::optional<symbol> sym_type() const = 0;
        virtual std::string to_string() const = 0;
        virtual std::optional<double> to_number() const = 0;
    };
    using brush_expr_ptr = std::shared_ptr<brush_expr_base>;

    class brush_expr : public brush_expr_base {
    public:
        brush_expr( ch::symbol op);

        template<typename Iter>
        brush_expr(ch::symbol op, Iter begin, Iter end) : brush_expr(op) {
            std::copy(begin, end, std::back_inserter(children_));
        }

        brush_pipeline_item eval() const override;
        std::optional<symbol> sym_type() const override;
        std::string to_string() const override;
        std::optional<double>  to_number() const override;
    
    private:
        ch::symbol op_;
        std::vector<brush_expr_ptr> children_;
    };

   
    class symbol_expr : public brush_expr_base {
    public:
        symbol_expr(symbol sym);
        brush_pipeline_item eval() const override;
        virtual std::optional<symbol> sym_type() const override;
        std::string to_string() const override;
        std::optional<double> to_number() const override;

    private:
        symbol sym_;
    };

    class num_expr : public brush_expr_base {
    public:
        num_expr(double val);
        brush_pipeline_item eval() const override;
        virtual std::optional<symbol> sym_type() const override;
        std::string to_string() const override;
        std::optional<double> to_number() const override;

    private:
        double val_;
    };

    std::variant<ch::brush_fn, std::string> parse_brush_language(const std::string& input);
};