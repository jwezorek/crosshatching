#pragma once

#include "geometry.hpp"
#include "util.hpp"
#include <opencv2/core.hpp>
#include <range/v3/all.hpp>
#include <tuple>
#include <functional>

namespace ch {

    struct dimensions {
        double wd;
        double hgt;
             
        dimensions(double d = 0.0);
        dimensions(double w, double h);
    };

    using crosshatching_range = ranges::any_view<polyline>;
    struct crosshatching_swatch {
        crosshatching_range content;
        dimensions sz;
        double stroke_wd;
    };

    using unit_of_hatching_fn = std::function<crosshatching_range(double, double, double, double)>;

    unit_of_hatching_fn make_shading_stroke(double sz_pcnt, double variance, bool centered);
    unit_of_hatching_fn make_brick_stroke();
    crosshatching_range one_horz_stroke(double x1, double x2, double y, double hgt);

    crosshatching_swatch linear_crosshatching( rnd_fn run_length,  rnd_fn space_length, rnd_fn vert_space,
        unit_of_hatching_fn h_fn = one_horz_stroke, dimensions sz = { 512,512 }, double stroke_wd = 1);
    crosshatching_swatch fragment(crosshatching_swatch rng, ch::rnd_fn frag);
    crosshatching_swatch jitter(crosshatching_swatch rng, ch::rnd_fn jitter);
    crosshatching_swatch jiggle(crosshatching_swatch rng, ch::rnd_fn jiggle);
    crosshatching_swatch rotate(crosshatching_swatch rng, double theta);
    crosshatching_swatch disintegrate(crosshatching_swatch rng, double amount);

    cv::Mat paint_cross_hatching(crosshatching_swatch swatch);
    double gray_level(crosshatching_swatch rng);
    void to_svg(const std::string& filename, crosshatching_swatch swatch);
   
}