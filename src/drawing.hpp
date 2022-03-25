#pragma once

#include "imgproc.hpp"
#include "brush.hpp"
#include <string>

namespace ch {

    struct segmentation_params {
        int sigmaS;
        float sigmaR;
        int minSize;
    };

    struct drawing {
        std::vector<polyline> strokes;
        dimensions sz;
        double stroke_wd;
    };

    drawing generate_crosshatched_drawing(const std::string& image_file, segmentation_params params, double scale, brush& br);
    std::vector<polyline> crosshatched_poly_with_holes(const polygon_with_holes& poly, double color, brush& br);
    void to_svg(const std::string& filename, const drawing& d);

}