#include "brushes.hpp"
#include "point_set.hpp"
#include "qdebug.h"
#include <tuple>
#include <map>
#include <stdexcept>
#include <span>
#include <memory>
#include <optional>
#include <sstream>
#include <fstream>

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

    auto running_sum(const ch::cbrng_seed& seed, const ch::random_func& gen_fn, double init_sum, uint32_t initial_counter) {
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

    auto row_of_strokes(double x1, double x2, double y, int rnd_key_1,
            ch::random_func run_length, ch::random_func space_length) {
        qDebug() << x1 << " " << x2 << " " << y << " " << rnd_key_1 << "\n";
        return 0;
    }

    struct run_and_space {
        double run;
        double space;
        double x;
    };

    auto row_of_strokes(double x1, double x2, double y, const ch::cbrng_seed& seed, uint32_t rnd_key_1, const ch::random_func& run_length, const ch::random_func& space_length) {
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
    auto horz_strokes_from_y_positions(const ch::cbrng_seed& seed, R row_y_positions, 
            std::shared_ptr<horz_hull_slicer> slicer, ch::random_func run_length, ch::random_func space_length) {
        return row_y_positions | // given a view of y-positions
            rv::enumerate | // tag each with an integer index
            rv::transform(  
                    // turn them into possibly empty horizontal slices of the convex polygon at the given y
                    // and also pass along the integer index which will be used as a key in a counter-based RNG.
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
            rv::remove_if( // remove any empty slices that may have occured above or below the polygon
                [](const auto& maybe_row) {
                    return !maybe_row.has_value();
                }
            ) | 
            rv::transform(  // generate a range view of horizontal crosshatching in the slice
                [seed, run_length, space_length](const auto& maybe_row) {

                    auto [key, x1, x2, y] = *maybe_row;
                    return row_of_strokes(x1, x2, y, seed, key, run_length, space_length);
                }
            ) |
            rv::join; // flatten all the row views above into one view.
    }

    ch::strokes linear_strokes(const ch::cbrng_seed& seed, double rotation, double thickness,
            const std::vector<ch::point>& region, ch::random_func run_length, 
            ch::random_func space_length, ch::random_func vert_space) {
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

    void write_to_svg(const std::string& filename, ch::drawing2 d,
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
}

void ch::debug_brushes() {

    std::vector<ch::point> pts{
        {100, 200}, {150,200},
        {50,150}, {200,150},
        {50,100},{200,100},
        {100,50},{150,50}
    };
    auto strokes = linear_strokes(
        cbrng_seed(17, 0, 0),
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
        ch::drawing2{ strokes, {500,500} },
        std::function<void(double)>{}
    );
}