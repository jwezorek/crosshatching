#pragma once

#include "brush.hpp"
#include <vector>
#include <memory>
#include <string>
#include <optional>
#include <iterator>

namespace ch {

    enum class symbol {
        true_,
        false_,
        pipe,
        norm_rnd,
        lerp,
        ramp,
        rotate,
        disintegrate
    };

    class brush_expr_base {
    public:
        virtual brush_pipeline_item eval() const = 0;
        virtual std::optional<symbol> sym_type() const = 0;
        virtual std::string to_string() const = 0;
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
    private:
        symbol sym_;
    };

    class num_expr : public brush_expr_base {
    public:
        num_expr(double val);
        brush_pipeline_item eval() const override;
        virtual std::optional<symbol> sym_type() const override;
        std::string to_string() const override;
    private:
        double val_;
    };

    brush_expr_ptr parse_brush_language(const std::string& input);
};