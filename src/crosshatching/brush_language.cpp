#include "brush_language.hpp"
#include "peglib.h"
#include <assert.h>
#include <iostream>
#include <unordered_map>
#include <sstream>

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    using expr_args = std::vector<ch::brush_expr_ptr>;
    using brush_expr_evaluator = std::function< ch::brush_pipeline_item(const expr_args&)>;

    bool parser_initialized = false;

    peg::parser parser(R"(
        ROOT <- Expr
        Symbol <- < [_a-z]+ >
        Number <- < [.0-9]+ >
        Op <- 'pipe' | 'norm_rnd' | 'lerp' | 'ramp' | 'rot' | 'dis' | 'jiggle' | 'merge' | 'lin_brush' | 'scatter'
        Expr <- '(' Op (Expr / Symbol / Number)+ ')'
        %whitespace <- [ \t\r\n]*
      )"
    );

    ch::symbol string_to_sym(const std::string& str) {
        static std::unordered_map<std::string, ch::symbol> tbl = {
            {"true",     ch::symbol::true_},
            {"false",    ch::symbol::false_},
            {"pipe",     ch::symbol::pipe},
            {"norm_rnd", ch::symbol::norm_rnd},
            {"lerp",     ch::symbol::lerp},
            {"ramp",     ch::symbol::ramp},
            {"rot",      ch::symbol::rotate},
            {"dis",      ch::symbol::disintegrate},
            {"merge",    ch::symbol::merge},
            {"jiggle",   ch::symbol::jiggle},
            {"lin_brush",ch::symbol::linear_brush},
            {"scatter", ch::symbol::scatter_brush}
        };
        return tbl.at(str);
    }

    std::string sym_to_string(ch::symbol sym) {
        static std::unordered_map<ch::symbol, std::string> tbl = {
            {ch::symbol::true_,        "true"},
            {ch::symbol::false_,       "false"},
            {ch::symbol::pipe,         "pipe"},
            {ch::symbol::norm_rnd,     "norm_rnd"},
            {ch::symbol::lerp,         "lerp"},
            {ch::symbol::ramp,         "ramp"},
            {ch::symbol::rotate,       "rot"},
            {ch::symbol::disintegrate, "dis"},
            {ch::symbol::merge,        "merge"},
            {ch::symbol::jiggle,       "jiggle"},
            {ch::symbol::linear_brush, "lin_brush"},
            {ch::symbol::scatter_brush, "scatter"}
        };
        return tbl.at(sym);
    }

    void init_parser() {
        ::parser["Symbol"] = [](const peg::SemanticValues& vs)->ch::brush_expr_ptr {
            auto val = vs.token_to_string();
            return std::make_shared<ch::symbol_expr>(string_to_sym(val));
        };
        
        ::parser["Number"] = [](const peg::SemanticValues& vs)->ch::brush_expr_ptr {
            auto val = vs.token_to_number<double>();
            return std::make_shared<ch::num_expr>( val );
        };
        
        ::parser["Op"] = [](const peg::SemanticValues& vs)->ch::brush_expr_ptr {
            auto val = vs.token_to_string();
            return std::make_shared<ch::symbol_expr>(string_to_sym(val));
        };

        ::parser["Expr"] = [](const peg::SemanticValues& vs)->ch::brush_expr_ptr {
            std::vector<ch::brush_expr_ptr> args(vs.size());
            std::transform(vs.begin(), vs.end(), args.begin(),
                [](auto vsi)->ch::brush_expr_ptr {
                    return std::any_cast<ch::brush_expr_ptr>(vsi);
                }
            );
            auto op = args.front()->sym_type();
            if (!op) {
                throw std::runtime_error("bad expression: no operation");
            }
            return std::make_shared<ch::brush_expr>(*op, std::next(args.begin()), args.end());
        };

        parser_initialized = true;
    }

    ch::param_adapter_fn expr_to_parameter_adapter(ch::brush_expr_ptr expr) {
        auto numeric_value = expr->to_number();
        if (numeric_value) {
            return ch::make_constant_fn(*numeric_value);
        } else {
            return std::get< ch::param_adapter_fn>(expr->eval());
        }
    }

    bool is_symbol(ch::symbol sym, ch::brush_expr_ptr expr) {
        auto sym_val = expr->sym_type();
        return sym_val ? *sym_val == sym : false;
    }

    ch::brush_pipeline_item eval_pipe_expr(const expr_args& args) {
        return ch::make_run_pipeline_fn(
            args |
            rv::transform(
                [](ch::brush_expr_ptr expr) {
                    return expr->eval();
                }
            ) |
            r::to_vector
        );
    }

    ch::brush_pipeline_item eval_norm_rnd_expr(const expr_args& args) {
        return ch::make_normal_dist_fn(
            expr_to_parameter_adapter(args[0]),
            expr_to_parameter_adapter(args[1])
        );
    }

    ch::brush_pipeline_item eval_lerp_expr(const expr_args& args) {
        return ch::make_lerp_fn(
            expr_to_parameter_adapter(args[0]),
            expr_to_parameter_adapter(args[1])
        );
    }

    ch::brush_pipeline_item eval_ramp_expr(const expr_args& args) {
        return ch::make_ramp_fn(
            expr_to_parameter_adapter(args[0]),
            is_symbol(ch::symbol::true_, args[1]),
            is_symbol(ch::symbol::true_, args[2])
        );
    }

    ch::brush_pipeline_item eval_rotate_expr(const expr_args& args) {
        return ch::make_one_param_brush_adaptor(ch::rotate_in_degrees, expr_to_parameter_adapter(args[0]));
    }

    ch::brush_pipeline_item eval_disintegrate_expr(const expr_args& args) {
        return ch::make_one_param_brush_adaptor(ch::disintegrate, expr_to_parameter_adapter(args[0]));
    }

    ch::brush_pipeline_item eval_merge_expr(const expr_args& args) {
        return ch::make_merge_fn(
            args |
            rv::transform(
                [](ch::brush_expr_ptr expr) {
                    return std::get<ch::brush_fn>(expr->eval());
                }
            ) |
            r::to_vector
        );
    }

    ch::brush_pipeline_item eval_jiggle_expr(const expr_args& args) {
        return make_random_brush_adaptor(ch::jiggle, expr_to_parameter_adapter(args[0]));
    }

    ch::brush_pipeline_item eval_linear_brush_expr(const expr_args& args) {
        return ch::make_linear_hatching_brush_fn(
            expr_to_parameter_adapter(args[0]),
            expr_to_parameter_adapter(args[1]),
            expr_to_parameter_adapter(args[2]),
            ch::make_default_hatching_unit()
        );
    }

    ch::brush_pipeline_item eval_scatter_brush_expr(const expr_args& args) {
        return ch::make_scatter_hatching_brush_fn(
            expr_to_parameter_adapter(args[0])
        );
    }
}

/*----------------------------------------------------------------------------------------------------*/

ch::brush_expr::brush_expr(ch::symbol op) :
    op_(op)
{}

ch::brush_pipeline_item ch::brush_expr::eval() const {
    static std::unordered_map<ch::symbol, brush_expr_evaluator> tbl = {
        {ch::symbol::pipe,          eval_pipe_expr},
        {ch::symbol::norm_rnd,      eval_norm_rnd_expr},
        {ch::symbol::lerp,          eval_lerp_expr},
        {ch::symbol::ramp,          eval_ramp_expr},
        {ch::symbol::rotate,        eval_rotate_expr},
        {ch::symbol::disintegrate,  eval_disintegrate_expr},
        {ch::symbol::merge,         eval_merge_expr},
        {ch::symbol::jiggle,        eval_jiggle_expr},
        {ch::symbol::linear_brush,  eval_linear_brush_expr},
        {ch::symbol::scatter_brush,  eval_scatter_brush_expr}
    };
    const auto& evaluator = tbl.at(op_);
    return evaluator(children_);
}

std::optional<ch::symbol> ch::brush_expr::sym_type() const {
    return {};
}

std::string ch::brush_expr::to_short_string() const
{
    return sym_to_string(op_);
}

std::string ch::brush_expr::to_string() const {
    std::stringstream ss;
    
    ss << "(" << sym_to_string(op_) << " ";
    for (auto iter = children_.begin(); iter != children_.end(); ++iter) {
        ss << (*iter)->to_string();
        if (iter != std::prev(children_.end())) {
            ss << " ";
        }
    }
    ss << ")";

    return ss.str();
}

std::optional<double> ch::brush_expr::to_number() const {
    return {};
}

const std::vector<ch::brush_expr_ptr>* ch::brush_expr::children() const  {
    return &children_;
}

bool ch::brush_expr::is_expression() const {
    return true;
}

/*----------------------------------------------------------------------------------------------------*/

ch::symbol_expr::symbol_expr(ch::symbol sym) :
    sym_(sym)
{}

ch::brush_pipeline_item ch::symbol_expr::eval() const {
    return {};
}

std::optional<ch::symbol> ch::symbol_expr::sym_type() const {
    return sym_;
}

std::string ch::symbol_expr::to_short_string() const
{
    return to_string();
}

std::string ch::symbol_expr::to_string() const {
    return sym_to_string(sym_);
}

std::optional<double> ch::symbol_expr::to_number() const {
    return {};
}

/*----------------------------------------------------------------------------------------------------*/

ch::num_expr::num_expr(double val) :
    val_(val)
{}

ch::brush_pipeline_item ch::num_expr::eval() const {
    return {};
}

std::optional<ch::symbol> ch::num_expr::sym_type() const {
    return {};
}

std::string ch::num_expr::to_short_string() const
{
    return to_string();
}

std::string ch::num_expr::to_string() const  {
    return std::to_string(val_);
}

std::optional<double> ch::num_expr::to_number() const {
    return val_;
}

/*----------------------------------------------------------------------------------------------------*/

std::variant<ch::brush_fn, std::string> ch::brush_language_to_func(const std::string& input) {
    auto result = brush_language_to_expr(input);
    if (std::holds_alternative<std::string>(result)) {
        return { std::get<std::string>(result) };
    }
    auto item = std::get<ch::brush_expr_ptr>(result)->eval();
    return std::get<ch::brush_fn>(item);
}

std::variant<ch::brush_expr_ptr, std::string> ch::brush_language_to_expr(const std::string& input) {
    try {
        if (!parser_initialized) {
            init_parser();
        }
        ch::brush_expr_ptr expr;
        bool success = ::parser.parse(input, expr);
        if (!success) {
            return std::string("error parsing brush");
        }
        return expr;
    } catch (std::runtime_error e) {
        return std::string(e.what());
    } catch (...) {
        return std::string("unknown error");
    };
}