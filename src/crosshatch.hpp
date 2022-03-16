#pragma once

#include "geometry.hpp"
#include "util.hpp"
#include <opencv2/core.hpp>
#include <range/v3/all.hpp>
#include <tuple>
#include <functional>

namespace ch {

    constexpr int k_swatch_sz = 256;

    using hatching_range = ranges::any_view<polyline>;
    using unit_of_hatching_fn = std::function<hatching_range(double, double, double, double)>;

    unit_of_hatching_fn make_shading_stroke(double sz_pcnt, double variance, bool centered);
    unit_of_hatching_fn make_brick_stroke();
    hatching_range one_horz_stroke(double x1, double x2, double y, double hgt);

    hatching_range linear_crosshatching(const ch::rnd_fn& run_length, const ch::rnd_fn& space_length, const ch::rnd_fn& vert_space,
        const unit_of_hatching_fn& h_fn = one_horz_stroke, int swatch_sz = k_swatch_sz);
    hatching_range fragment(hatching_range rng, const ch::rnd_fn& frag);
    hatching_range jitter(hatching_range rng, const ch::rnd_fn& jitter);
    hatching_range jiggle(hatching_range rng, const ch::rnd_fn& jiggle);
    hatching_range rotate(hatching_range rng, double theta);
    hatching_range disintegrate(hatching_range rng, double amount);

    cv::Mat paint_cross_hatching(int thickness, hatching_range rng, int swatch_sz = k_swatch_sz);
    double gray_level(int thickness, hatching_range rng);
    void to_svg(const std::string& filename, int thickness, hatching_range rng, int swatch_sz = k_swatch_sz);
   
}