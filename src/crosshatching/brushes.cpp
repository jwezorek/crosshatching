#include "brushes.hpp"
#include "point_set.hpp"
#include <tuple>
#include <map>
#include <stdexcept>

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
            return left_side_.front().y;
        }

        double bottom() const {
            return left_side_.back().y;
        }

        std::tuple<double, double> slice(double y) {
            auto left = intersecting_line_segment(left_side_, y);
            auto right = intersecting_line_segment(right_side_, y);
            return { x_on_line_segement_at_y(y, left), x_on_line_segement_at_y(y, right) };
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

    ch::stroked_region linear_strokes(ch::rand_number_state& rng, const ch::drawing_context& ctxt,
            const ch::ring& region) {
        
    }

}

ch::brush_func ch::make_linear_strokes_fn( rnd_fn run_length,  rnd_fn space_length, 
        rnd_fn vert_space) {
    return {};
}

void ch::debug_brushes() {
    std::vector<ch::point> pts{
        {10, 20}, {15,20},
        {5,15}, {20,15},
        {5,10},{20,10},
        {10,5},{15,5}
    };
    auto hull = ch::convex_hull(pts);
    horz_hull_slicer slicer(hull);
    auto [x1, x2] = slicer.slice(7.5);
    int aaa;
    aaa = 5;
}