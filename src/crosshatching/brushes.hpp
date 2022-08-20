#pragma once

#include <map>
#include "geometry.hpp"
#include "util.hpp"
#include <opencv2/core.hpp>
#include <range/v3/all.hpp>
#include <tuple>
#include <functional>
#include <variant>
#include <vector>
#include <memory>

namespace ch {

    struct stroke {
        ranges::any_view<ch::point> polyline;
        double pen_thickness;
    };
    using strokes = ranges::any_view<stroke>;

    struct stroked_region {
        strokes content;
        ring region;
    };

    struct drawing_context {
        double t;
        double pen_direction;
        double pen_thickness;
    };

    using brush_func = std::function<stroked_region(rand_number_state&, 
            const drawing_context&, const ring&)>;
    brush_func make_linear_strokes_fn(rnd_fn run_length, rnd_fn space_length, rnd_fn vert_space);

    void debug_brushes();
}