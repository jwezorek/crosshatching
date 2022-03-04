#pragma once

#include "geometry.hpp"
#include <opencv2/core.hpp>
#include <range/v3/all.hpp>
#include <tuple>
#include <functional>

namespace ch {
    using poly_range = ranges::any_view<point>;
    using hatching_range = ranges::any_view<polyline>;
    using unit_of_hatching_fn = std::function<hatching_range(double, double, double, double)>;

    unit_of_hatching_fn make_shading_stroke(double sz_pcnt, double variance, bool centered);
    unit_of_hatching_fn make_brick_stroke();
    hatching_range one_horz_stroke(double x1, double x2, double y, double hgt);

    hatching_range linear_crosshatching(int dim, const ch::rnd_fn& run_length, const ch::rnd_fn& space_length, const ch::rnd_fn& vert_space,
        const unit_of_hatching_fn& h_fn = one_horz_stroke);
    hatching_range rotated_crosshatching(double theta, int dim, const ch::rnd_fn& run_length, const ch::rnd_fn& run_space_length, const ch::rnd_fn& line_space,
        const unit_of_hatching_fn& h_fn = one_horz_stroke);
    hatching_range apply_jitter(hatching_range rng, const ch::rnd_fn& run_length, const ch::rnd_fn& jitter);

    struct fragment_param {
        double min_length;
        ch::rnd_fn run_length;
    };

    hatching_range disintegrate(hatching_range rng, double amount, const std::optional<fragment_param>& frag = {});
    cv::Mat paint_cross_hatching(int thickness, int dim, hatching_range rng);
    double gray_level(int thickness, hatching_range rng);
    void to_svg(const std::string& filename, int dim, int thickness, hatching_range rng);

    void debug();
}