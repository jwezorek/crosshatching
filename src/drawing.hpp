#pragma once

#include "imgproc.hpp"
#include "brush.hpp"

namespace ch {
    struct drawing {
        std::vector<polyline> strokes;
        dimensions sz;
    };

    std::vector<polyline> crosshatched_poly_with_holes(const polygon_with_holes& poly, uchar color, brush& br);
    void to_svg(const std::string& filename, const drawing& d);

}