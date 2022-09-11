#include "brush_lang.hpp"
#include "point_set.hpp"
#include "peglib.h"
#include "qdebug.h"
#include <tuple>
#include <map>
#include <stdexcept>
#include <span>
#include <memory>
#include <optional>
#include <sstream>
#include <fstream>
#include <numeric>

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    double x_on_line_segement_at_y(double y, const ch::line_segment& line) {
        auto [u, v] = line;
        return (u.y * v.x - u.x * v.y + u.x * y - v.x * y) / (u.y - v.y);
    }

    class horz_hull_slicer {
    public:
        horz_hull_slicer(std::span<ch::point> hull) {
            populate_left_and_right_sides(hull);
        }

        double top() const {
            return left_side_.back().y;
        }

        double bottom() const {
            return left_side_.front().y;
        }

        std::tuple<double, double> vertical_extent() const {
            return { bottom(), top() };
        }

        std::optional<std::tuple<double, double>> slice(double y) {
            if (y < bottom() || y > top()) {
                return {};
            }
            auto left = intersecting_line_segment(left_side_, y);
            auto right = intersecting_line_segment(right_side_, y);
            return {{ x_on_line_segement_at_y(y, left), x_on_line_segement_at_y(y, right) }};
        }

    private:
        std::vector<ch::point> left_side_;
        std::vector<ch::point> right_side_;

        ch::line_segment intersecting_line_segment(const std::vector<ch::point>& side, double y) {
            auto iter = std::lower_bound(side.begin(), side.end(), y,
                [](const ch::point& lhs, double rhs) {
                    return lhs.y < rhs;
                }
            );
            return { *iter, *std::prev(iter) };
        }

        void populate_left_and_right_sides(std::span<ch::point> hull) {
            auto edges = hull | rv::cycle | rv::sliding(2) |
                rv::transform(
                    [](auto rng)->std::tuple< ch::point, ch::point> {
                        return { rng[0], rng[1] };
                    }
                );
            ch::point_map<ch::point> left_side_set;
            ch::point_map<ch::point> right_side_set;
            for (const auto& [pt1, pt2] : edges | rv::take(hull.size())) {
                if (pt1.y < pt2.y) {
                    left_side_set[pt1] = pt1;
                    left_side_set[pt2] = pt2;
                } else if (pt1.y > pt2.y) {
                    right_side_set[pt1] = pt1;
                    right_side_set[pt2] = pt2;
                }
            }
            left_side_ = left_side_set | 
                rv::transform([](const auto& p) {return p.second; }) | r::to_vector;
            right_side_ = right_side_set |
                rv::transform([](const auto& p) {return p.second; }) | r::to_vector;

            auto compare_y = [](const ch::point& lhs, const ch::point& rhs) {
                return lhs.y < rhs.y;
            };
            r::sort(left_side_, compare_y);
            r::sort(right_side_, compare_y);
        }
    };

    auto running_sum(uint32_t seed, const ch::random_func& gen_fn,
            double init_sum, uint32_t initial_counter) {
        return 
            rv::concat(
                rv::single(0),
                rv::iota(initial_counter) |
                    rv::transform(
                        [seed, gen_fn](int counter)->double {
                            return  gen_fn(ch::cbrng_state(seed, counter));
                        }
                    )
            )|
            rv::partial_sum() |
            rv::transform(
                [init_sum](double sum) {
                    return sum + init_sum;
                }
            );
    }

    struct run_and_space {
        double run;
        double space;
        double x;
    };

    auto row_of_strokes(double x1, double x2, double y, 
            uint32_t seed, uint32_t rnd_key_1,
            const ch::random_func& run_length, const ch::random_func& space_length) {
        x1 -= run_length(ch::cbrng_state(seed, rnd_key_1, 1));
        x2 += run_length(ch::cbrng_state(seed, rnd_key_1, 2));
        return rv::iota(1) |
            rv::transform(
                [=](int i)->run_and_space {
                    return run_and_space{
                        run_length(ch::cbrng_state(seed, rnd_key_1, 2 * i + 1)),
                        space_length(ch::cbrng_state(seed, rnd_key_1, 2 * i + 2)),
                        x1
                    };
                }
            ) |
            rv::partial_sum(
                [](run_and_space a, run_and_space b)->run_and_space {
                    return { b.run, b.space, a.x + a.run + a.space };
                }
            ) |
            rv::transform(
                    [y](run_and_space rs) {
                        auto x1 = rs.x;
                        auto x2 = rs.x + rs.run - 1;
                        return rv::concat(
                            rv::single(ch::point(x1, y)),
                            rv::single(ch::point(x2, y))
                        );
                    }
            ) |
            rv::take_while(
                [x2](auto line_segment) {
                    return line_segment[1].x < x2;
                }
            );
    }

    template<typename R>
    auto horz_strokes_from_y_positions(uint32_t seed, R row_y_positions, 
            std::shared_ptr<horz_hull_slicer> slicer, ch::random_func run_length, 
            ch::random_func space_length) {
        return row_y_positions | // given a view of y-positions
            rv::enumerate | // tag each with an integer index
            rv::transform(  
                    // turn them into possibly empty horizontal slices of the convex polygon at 
                    // the given y and also pass along the integer index which will be used as a 
                    // key in a counter-based RNG.
                [slicer, run_length, space_length, seed](const std::tuple<int, double>& pair) ->
                        std::optional<std::tuple<int, double,double,double>> {
                    auto [key, y] = pair;
                    auto slice = slicer->slice(y);
                    if (!slice) {
                        return {};
                    } else {
                        auto [x1, x2] = *slice;
                        return { {key, x1, x2, y } };
                    }
                }
            ) | 
            rv::remove_if( 
                    // remove any empty slices that may have occured above or below the polygon
                [](const auto& maybe_row) {
                    return !maybe_row.has_value();
                }
            ) | 
            rv::transform(  
                    // generate a range view of horizontal crosshatching in the slice
                [seed, run_length, space_length](const auto& maybe_row) {

                    auto [key, x1, x2, y] = *maybe_row;
                    return row_of_strokes(x1, x2, y, seed, key, run_length, space_length);
                }
            ) |
            rv::join; // flatten all the row views above into one view.
    }

    ch::strokes linear_strokes(double rotation, double thickness,
            const std::vector<ch::point>& region, ch::random_func run_length, 
            ch::random_func space_length, ch::random_func vert_space) {
        auto seed = ch::random_seed();
        auto rot_matrix = ch::rotation_matrix(rotation);
        auto hull = ch::convex_hull(region);
        hull = hull |
            rv::transform(
                    [rot_matrix](const auto& pt) {return ch::transform(pt, rot_matrix); }
                ) |
            r::to_vector;
        auto slicer = std::make_shared<horz_hull_slicer>(hull);
        auto [bottom, top] = slicer->vertical_extent();
        bottom -= vert_space(ch::cbrng_state(seed, 0));
        top += vert_space(ch::cbrng_state(seed, 1));
        auto rows = running_sum(seed, vert_space, bottom, 2) |
            rv::take_while( [top](double y) {return y < top;} );

        auto inverse_rot_matrix = ch::rotation_matrix(-rotation);
        return horz_strokes_from_y_positions(seed, rows, slicer, run_length, space_length) |
            rv::transform(
                [inverse_rot_matrix, thickness = thickness](auto polyline)->ch::stroke {
                    return {
                        polyline | rv::transform(
                            [inverse_rot_matrix](const auto& pt) {
                                return ch::transform(pt, inverse_rot_matrix); 
                            }
                        ),
                        thickness
                    };
                }
            );
    }

    struct drawing2 {
        ch::strokes content;
        ch::dimensions<int> size;
    };

    void write_to_svg(const std::string& filename, drawing2 d,
            std::function<void(double)> update_progress) {
        std::ofstream outfile(filename);

        outfile << ch::svg_header(static_cast<int>(d.size.wd), static_cast<int>(d.size.hgt));

        auto strokes = d.content | r::to_vector;
        auto n = static_cast<int>(strokes.size());
        
        for (const auto& [index, s] : rv::enumerate(strokes)) {
            auto polyline = s.polyline | r::to_vector;
            outfile << ch::polyline_to_svg(polyline, s.pen_thickness) << std::endl;
            if (update_progress) {
                update_progress(index / n);
            }
        }

        outfile << "</svg>" << std::endl;
        outfile.close();
    }

    constexpr auto k_brush_param_var = "<<t>>";
    constexpr auto k_rotation_var = "<<rotation>>";
    constexpr auto k_pen_thickness_var = "<<pen_thickness>>";

    template<typename T>
    T variable_value(const ch::variables_map& vars,  const std::string& key, T default_val) {
        auto iter = vars.find(key);
        if (iter == vars.end()) {
            return default_val;
        }
        return std::get<T>(iter->second);
    }

    class pipe_expr : public ch::brush_expr {

        class expr_val_visitor {
            ch::brush_context ctxt_;
        public:
            expr_val_visitor(ch::brush_context& ctxt) :
                ctxt_(ctxt)
            {}

            void operator()(ch::nil_value nil) {
            }

            void operator()(double value) {
                ctxt_.variables[k_brush_param_var] = value;
            }

            void operator()(ch::strokes strokes) {
                if (ctxt_.strokes) {
                    ctxt_.strokes = rv::concat(*ctxt_.strokes, strokes);
                } else {
                    ctxt_.strokes = strokes;
                }
            }

            void operator()(ch::random_func rnd){
                throw std::runtime_error("random number generator at pipeline scope.");
            }
        };

    public:

        ch::brush_expr_value eval(ch::brush_context& ctxt) {

            ch::brush_context local_ctxt = ctxt;
            expr_val_visitor visitor(local_ctxt);

            for (const auto& child : children_) {
                auto value = child->eval(local_ctxt);
                std::visit(visitor, value);
            }
            if (ctxt.strokes) {
                throw std::runtime_error("pipeline did not emit crosshatching.");
            }
            return *ctxt.strokes;
        }
    };

    class number_expr : public ch::brush_expr {
        double value_;
    public:
        number_expr(double num) : value_(num)
        {}

        ch::brush_expr_value eval(ch::brush_context& ctxt) {
            return value_;
        }
    };

    class variable_expr : public ch::brush_expr {
        std::string var_;
    public:
        variable_expr(const std::string& var) : var_(var)
        {}

        std::string var() const {
            return var_;
        }

        ch::brush_expr_value eval(ch::brush_context& ctxt) {
            auto iter = ctxt.variables.find(var_);
            if (iter == ctxt.variables.end()) {
                std::string err = "Undefined variable: " + var_;
                throw std::runtime_error(err);
            }
            return iter->second;
        }
    };

    class define_expr : public ch::brush_expr {
    
    public:

        ch::brush_expr_value eval(ch::brush_context& ctxt) {
            try {
                std::string var = std::dynamic_pointer_cast<variable_expr>(children_[0])->var();
                ctxt.variables[var] = children_[1]->eval(ctxt);
            } catch (...) {
                throw std::runtime_error("bad define expression");
            }
            return {};
        }
    };

    std::string error_msg_from_tag(const std::string& msg, const std::string& err_tag) {
        if (err_tag.empty()) {
            return msg;
        }
        return err_tag + " : " + msg;
    }

    template<typename... Ts>
    std::tuple<Ts...> evaluate_to_tuple(ch::brush_context& ctxt, 
            const std::vector<ch::brush_expr_ptr>& exprs, const std::string& err_tag = "") {
        if (sizeof...(Ts) != exprs.size()) {
            throw std::runtime_error( error_msg_from_tag("wrong number of args", err_tag) );
        }
        int i = 0;
        try {
            return  std::tuple<Ts...>{
                std::get<Ts>(exprs[i++]->eval(ctxt)) ...
            };
        } catch (...) {
            throw std::runtime_error( error_msg_from_tag("bad argument(s)", err_tag) );
        }
    }

    template<typename T>
    std::vector<T> evaluate_to_vector(ch::brush_context& ctxt,
            const std::vector<ch::brush_expr_ptr>& exprs, const std::string& err_tag = "") {
        try {
            return exprs |
                rv::transform(
                    [&ctxt](ch::brush_expr_ptr ptr) {
                        return std::get<T>(ptr->eval(ctxt));
                    }
                ) | r::to_vector;
        } catch (...) {
            throw std::runtime_error(error_msg_from_tag("bad argument(s)", err_tag));
        }
    }

    class ramp_expr : public ch::brush_expr {
    public:

        ch::brush_expr_value eval(ch::brush_context& ctxt) {
            if (children_.size() == 4) {
                auto [k, right, up, t] = evaluate_to_tuple<double, double, double, double>(
                    ctxt, children_,"ramp"
                );
                return ch::ramp(t, k, right != 0.0, up != 0.0);
            } 
            auto [k, right, up] = evaluate_to_tuple<double, double, double>(ctxt, children_);
            auto t = variable_value(ctxt.variables, k_brush_param_var, 0.0);
            return ch::ramp(t, k, right != 0.0, up != 0.0);
        }
    };

    class norm_rnd_expr : public ch::brush_expr {
    public:

        ch::brush_expr_value eval(ch::brush_context& ctxt) {
            auto [mean, stddev] = evaluate_to_tuple<double, double>(
                ctxt, children_, "norm_rnd"
            );
            return [=](const ch::cbrng_state& rnd)->double {
                return ch::normal_random(rnd, mean, stddev);
            };
        }
    };

    class lerp_expr : public  ch::brush_expr {
    public:

        ch::brush_expr_value eval(ch::brush_context& ctxt) {
            if (children_.size() == 3) {
                auto [from, to, t] = evaluate_to_tuple<double, double, double>(
                    ctxt, children_, "lerp"
                );
                return std::lerp(from, to, t);
            }
            auto [from, to] = evaluate_to_tuple<double, double>(
                ctxt, children_, "lerp"
            );
            auto t = variable_value(ctxt.variables, k_brush_param_var, 0.0);
            return std::lerp(from, to, t);
        }
    };

    class linear_brush_expr : public ch::brush_expr {
    public:
        ch::brush_expr_value eval(ch::brush_context& ctxt) {
            auto [run_length, space_length, vert_space] =
                evaluate_to_tuple<ch::random_func, ch::random_func, ch::random_func>(
                    ctxt, children_, "strokes"
                );
            auto pen_thickness = variable_value(ctxt.variables, k_pen_thickness_var, 1.0);
            auto rotation = variable_value(ctxt.variables, k_brush_param_var, 0.0);
            return linear_strokes(
                rotation, 
                pen_thickness,
                ctxt.poly.outer(), 
                run_length, 
                space_length, 
                vert_space
            );
        }
    };

    class add_expr : public ch::brush_expr {
    public:
        ch::brush_expr_value eval(ch::brush_context& ctxt) {
            auto args = evaluate_to_vector<double>(ctxt, children_, "add");
            return r::accumulate(args, 0.0);
        }
    };

    class multiply_expr : public ch::brush_expr {
    public:
        ch::brush_expr_value eval(ch::brush_context& ctxt) {
            auto args = evaluate_to_vector<double>(ctxt, children_, "multiply");
            return r::accumulate(args, 1.0, std::multiplies<double>());
        }
    };

    class subtract_expr : public ch::brush_expr {
    public:
        ch::brush_expr_value eval(ch::brush_context& ctxt) {
            auto args = evaluate_to_vector<double>(ctxt, children_, "subtract");
            return args.front() - r::accumulate(args | rv::drop(1), 0.0);
        }
    };

    class divide_expr : public ch::brush_expr {
    public:
        ch::brush_expr_value eval(ch::brush_context& ctxt) {
            auto args = evaluate_to_vector<double>(ctxt, children_, "divide");
            return args.front() / r::accumulate(args | rv::drop(1), 1.0, std::multiplies<double>());
        }
    };

    class disintegrate_expr : public ch::brush_expr {
    public:
        ch::brush_expr_value eval(ch::brush_context& ctxt) {
            auto result = children_.front()->eval(ctxt);
            if (!std::holds_alternative<double>(result)) {
                throw std::runtime_error("disintegrate : non-numeric argument");
            }
            auto probability = std::clamp(std::get<double>(result), 0.0, 1.0);
            auto seed = ch::random_seed();
            auto disintegrated =
                rv::enumerate(*ctxt.strokes) |
                rv::filter(
                    [probability, seed](const auto& tup)->bool {
                        auto [index, stroke] = tup;
                        return ch::uniform_rnd(ch::cbrng_state(seed, index), 0.0, 1.0) <= probability;
                    }
                ) |
                rv::transform(
                    [](const auto& tup) {
                        return std::get<1>(tup);
                    }
                );

            return ch::strokes{ disintegrated };
        }
    };
    
    class jiggle_expr : public ch::brush_expr {
    public:
        ch::brush_expr_value eval(ch::brush_context& ctxt) {
            auto result = children_.front()->eval(ctxt);
            if (!std::holds_alternative<ch::random_func>(result)) {
                throw std::runtime_error("jiggle : invalid argument");
            }
            auto seed = ch::random_seed();

            return ch::strokes{
                rv::enumerate(*ctxt.strokes) |
                rv::transform(
                    [seed, rand = std::get<ch::random_func>(result)](const auto& tup){
                        auto [index, stroke] = tup;
                        auto theta = ch::degrees_to_radians(
                            rand(ch::cbrng_state(seed, index))
                        );
                        auto rot_matrix = ch::rotation_matrix(theta);
                        return ch::transform(stroke, rot_matrix);
                    }
                )
            };
        }
    };

    enum class operation {
        define,
        brush,
        horz_strokes,
        norm_rnd,
        lerp,
        ramp,
        rotate,
        disintegrate,
        jiggle,
        pen_thickness,
        add,
        subtract,
        multiply,
        divide,
    };

    using expr_creation_fn = std::function<ch::brush_expr_ptr(std::span<const ch::brush_expr_ptr>)>;
    struct op_info {
        operation op;
        std::string str;
        expr_creation_fn create;
    };

    template<typename T>
    expr_creation_fn make_op_create_fn() {
        return [](std::span<const ch::brush_expr_ptr> args)->ch::brush_expr_ptr {
            auto expr = std::make_shared<T>();
            expr->set_children(args);
            return expr;
        };
    }

    expr_creation_fn make_rot_create_fn(){
        return [](std::span<const ch::brush_expr_ptr> args)->ch::brush_expr_ptr {
            std::vector<ch::brush_expr_ptr> new_args = {
                std::make_shared<variable_expr>(k_rotation_var),
                args.front()
            };
            auto def_expr = std::make_shared<define_expr>();
            def_expr->set_children(new_args);
            return def_expr;
        };
    }

    std::vector<op_info> g_ops{
        {operation::define,       "def",          make_op_create_fn<define_expr>()},
        {operation::brush,        "brush",        make_op_create_fn<pipe_expr>()},
        {operation::norm_rnd,     "norm_rnd",     make_op_create_fn<norm_rnd_expr>()},
        {operation::lerp,         "lerp",         make_op_create_fn<lerp_expr>()},
        {operation::ramp,         "ramp",         make_op_create_fn<ramp_expr>()},
        {operation::rotate,       "rot",          make_rot_create_fn()},
        {operation::disintegrate, "dis",          make_op_create_fn<disintegrate_expr>()},
        {operation::jiggle,       "jiggle",       make_op_create_fn<jiggle_expr>()},
        {operation::horz_strokes, "horz_strokes", make_op_create_fn<linear_brush_expr>()},
        {operation::add     ,     "+",            make_op_create_fn<add_expr>()},
        {operation::subtract,     "-",            make_op_create_fn<subtract_expr>()},
        {operation::multiply,     "*",            make_op_create_fn<multiply_expr>()},
        {operation::divide  ,     "/",            make_op_create_fn<divide_expr>()}
    };

    operation string_to_op(const std::string& str) {
        static std::unordered_map<std::string, operation> tbl;
        if (tbl.empty()) {
            for (const auto& opi : g_ops) {
                tbl[opi.str] = opi.op;
            }
        }
        return tbl.at(str);
    }

    std::string op_to_string(operation sym) {
        static std::unordered_map<operation, std::string> tbl;
        if (tbl.empty()) {
            for (const auto& opi : g_ops) {
                tbl[opi.op] = opi.str;
            }
        }
        return tbl.at(sym);
    }

    ch::brush_expr_ptr create_expr(operation sym, std::span<const ch::brush_expr_ptr> args) {
        static std::unordered_map<operation, expr_creation_fn> tbl;
        if (tbl.empty()) {
            for (const auto& opi : g_ops) {
                tbl[opi.op] = opi.create;
            }
        }
        auto create = tbl.at(sym);
        return create(args);
    }

    auto operations() {
        return g_ops |
            rv::transform(
                [](const op_info& opi)->operation {
                    return opi.op;
                }
        );
    }

    int num_ops() {
        return static_cast<int>(g_ops.size());
    }

    void insert_semantic_actions(peg::parser& parser) {
        parser["Op"] = [](const peg::SemanticValues& vs)->operation {
            return string_to_op(vs.token_to_string());
        };

        parser["Symbol"] = [](const peg::SemanticValues& vs)->ch::brush_expr_ptr {
            auto symbol = vs.token_to_string();
            if (symbol == "true") {
                return std::make_shared<number_expr>(1);
            } else if (symbol == "false") {
                return std::make_shared<number_expr>(0);
            }
            return std::make_shared<variable_expr>(symbol);
        };
        parser["Number"] = [](const peg::SemanticValues& vs)->ch::brush_expr_ptr {
            auto val = vs.token_to_number<double>();
            return std::make_shared<number_expr>(val);
        };

        parser["Expr"] = [](const peg::SemanticValues& vs)->ch::brush_expr_ptr {
            std::vector<ch::brush_expr_ptr> args(vs.size());
            std::transform(std::next(vs.begin()), vs.end(), args.begin(),
                [](auto vsi)->ch::brush_expr_ptr {
                    return std::any_cast<ch::brush_expr_ptr>(vsi);
                }
            );
            auto op = std::any_cast<operation>(vs.front());
            return  create_expr(op, args);
        };
    }

    peg::parser make_parser() {
        auto make_op_rule = []() {
            auto quoted_strings =
                g_ops |
                rv::transform(
                     [](const op_info& opi)->std::string {
                         return " '" + opi.str + "' ";
                     }
                ) | r::to_vector;
            std::string rhs = quoted_strings | rv::join('|') | r::to_<std::string>();
            return "Op <-" + rhs;
        };

        std::stringstream ss;
        ss << "ROOT <- Expr\n";
        ss << make_op_rule() << "\n";
        ss << "Symbol <- < [_a-z]+ >\n";
        ss << "Number <- < [.0-9]+ >\n";
        ss << "Expr <- '(' Op (Expr / Symbol / Number)+ ')'\n";

        auto parser = peg::parser(ss.str());
        insert_semantic_actions(parser);

        return parser;
    }
}

ranges::any_view<ch::point> ch::transform(r::any_view<ch::point> poly, const ch::matrix& mat) {
    return poly |
        rv::transform(
            [mat](const ch::point& pt) {
                return ch::transform(pt, mat);
            }
    );
}

ch::stroke ch::transform(ch::stroke s, const ch::matrix& mat) {
    return ch::stroke{
        ch::transform(s.polyline, mat),
        s.pen_thickness
    };
}

ch::strokes ch::transform(ch::strokes strokes, const ch::matrix& mat) {
    return strokes |
        rv::transform(
            [mat](auto stroke) {
                return ch::transform(stroke, mat);
            }
        );
}

void ch::brush_expr::set_children(std::span<const ch::brush_expr_ptr> children) {
    children_ = children | r::to_vector;
}

std::variant<ch::brush_expr_ptr, std::runtime_error> ch::parse(const std::string& str) {
    return std::runtime_error("TODO");
}

void ch::debug_brushes() {

    auto parser = make_parser();

    std::vector<ch::point> pts{
        {100, 200}, {150,200},
        {50,150}, {200,150},
        {50,100},{200,100},
        {100,50},{150,50}
    };
    auto strokes = linear_strokes(
        ch::degrees_to_radians(45),
        1,
        pts,
        normal_rnd_func(50, 15),
        normal_rnd_func(2, 0.25),
        normal_rnd_func(3, 0.5)
    );

    auto hull = convex_hull(pts) | r::to_vector;
    strokes = rv::concat(rv::single(ch::stroke{hull, 2 }), strokes);
    
    write_to_svg(
        "C:\\test\\foo.svg",
        drawing2{ strokes, {500,500} },
        std::function<void(double)>{}
    );
}