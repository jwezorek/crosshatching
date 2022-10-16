#pragma once

#include <range/v3/all.hpp>
#include "qpainter.h"
#include "geometry.hpp"

namespace ch {

    using stroke_range = ranges::any_view<ch::point>;
    using stroke_ranges = ranges::any_view<stroke_range>;
    struct stroke_cluster {
        stroke_ranges strokes;
        double thickness;
    };
    using strokes = ranges::any_view<stroke_cluster>;

    struct drawn_stroke_cluster {
        ch::polylines strokes;
        double thickness;
    };
    using drawn_strokes = std::vector<drawn_stroke_cluster>;

    stroke_cluster transform(stroke_cluster s, const ch::matrix& mat);
    strokes transform(strokes s, const ch::matrix& mat);
    drawn_strokes to_drawn_strokes(strokes strks);
    ch::polylines to_polylines(ch::stroke_ranges ranges);
    void paint_strokes(QPainter& g, const drawn_strokes& str);
}