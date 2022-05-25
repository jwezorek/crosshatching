#pragma once
#include "brush.hpp"
#include "geometry.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#include <string>

namespace ch {

    struct drawing {
        std::vector<polyline> strokes;
        dimensions sz;
        double stroke_wd;
    };

    drawing generate_crosshatched_drawing(cv::Mat img, cv::Mat label_img, double scale, const std::vector<std::tuple<ch::brush_fn, double>>& brushes);
    drawing generate_crosshatched_drawing(cv::Mat img, double scale, const std::vector<std::tuple<ch::brush_fn, double>>& brushes);
    void to_svg(const std::string& filename, const drawing& d);
    void debug(cv::Mat img, cv::Mat labels);

}