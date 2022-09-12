#pragma once

#include "old_brush.hpp"
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
        scatter_brush,
        norm_rnd,
        lerp,
        ramp,
        rotate,
        disintegrate,
        jiggle,
        merge
    };

    class brush_expression_base;
    using brush_expression_ptr = std::shared_ptr<brush_expression_base>;

    class brush_expression_base {
    public:
        virtual brush_pipeline_item eval() const = 0;
        virtual std::optional<symbol> sym_type() const = 0;
        virtual std::optional<double> to_number() const = 0;
        virtual std::string to_short_string() const = 0;
        virtual std::string to_string() const = 0;
        virtual std::string to_formatted_string(int n = 0) const;
        virtual const std::vector<brush_expression_ptr>* children() const;
        virtual bool is_expression() const;
    };

    class brush_expression : public brush_expression_base {
    public:
        brush_expression( ch::symbol op);

        template<typename Iter>
        brush_expression(ch::symbol op, Iter begin, Iter end) : brush_expression(op) {
            std::copy(begin, end, std::back_inserter(children_));
        }

        brush_pipeline_item eval() const override;
        std::optional<symbol> sym_type() const override;
        std::string to_short_string() const override;
        std::string to_string() const override;
        std::string to_formatted_string(int n) const override;
        std::optional<double>  to_number() const override;
        const std::vector<brush_expression_ptr>* children() const override;

        bool is_expression() const override;
        bool is_one_liner() const;
        void replace_child(brush_expression_ptr old_child, brush_expression_ptr new_child);
    
    private:
        ch::symbol op_;
        std::vector<brush_expression_ptr> children_;
    };

   
    class symbol_expr : public brush_expression_base {
    public:
        symbol_expr(symbol sym);
        brush_pipeline_item eval() const override;
        virtual std::optional<symbol> sym_type() const override;
        std::string to_short_string() const override;
        std::string to_string() const override;
        std::optional<double> to_number() const override;

    private:
        symbol sym_;
    };

    class num_expr : public brush_expression_base {
    public:
        num_expr(double val);
        brush_pipeline_item eval() const override;
        virtual std::optional<symbol> sym_type() const override;
        std::string to_short_string() const override;
        std::string to_string() const override;
        std::optional<double> to_number() const override;

    private:
        double val_;
    };

    std::variant<ch::brush_fn, std::string> brush_language_to_func(const std::string& input);
    std::variant<ch::brush_expression_ptr, std::string> brush_language_to_expr(const std::string& input);
};