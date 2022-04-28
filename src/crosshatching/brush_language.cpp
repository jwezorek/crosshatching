#include "brush_language.hpp"
#include "peglib.h"
#include <assert.h>
#include <iostream>
#include <unordered_map>
#include <sstream>

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    bool parser_initialized = false;

    peg::parser parser(R"(
        ROOT <- Expr
        Symbol <- < [_a-z]+ >
        Number <- < [.0-9]+ >
        Expr <- '(' (Expr / Symbol / Number)* ')'
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
            {"dis",      ch::symbol::disintegrate}
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
            {ch::symbol::disintegrate, "dis"}
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
}

/*----------------------------------------------------------------------------------------------------*/

ch::brush_expr::brush_expr(ch::symbol op) :
    op_(op)
{}

ch::brush_pipeline_item ch::brush_expr::eval() const {
    return {};
}

std::optional<ch::symbol> ch::brush_expr::sym_type() const {
    return {};
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

std::string ch::symbol_expr::to_string() const {
    return sym_to_string(sym_);
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

std::string ch::num_expr::to_string() const  {
    return std::to_string(val_);
}

/*----------------------------------------------------------------------------------------------------*/

ch::brush_expr_ptr ch::parse_brush_language(const std::string& input) {
    if (!parser_initialized) {
        init_parser();
    }
    ch::brush_expr_ptr expr;
    bool success = ::parser.parse(input, expr);
    auto test = expr->to_string();
    return expr;
}