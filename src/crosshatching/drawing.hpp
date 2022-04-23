#pragma once
#include "brush.hpp"
#include "geometry.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
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

    struct polygon_with_holes {
        polyline border;
        std::vector<polyline> holes;
    };

    struct gray_level {
        double value;
        std::vector<polygon_with_holes> blobs;
    };

    std::vector<gray_level> extract_gray_levels(const cv::Mat& gray_scale_img, bool hierarchical);
    std::vector<gray_level> scale(const std::vector<gray_level>& planes, double scale);
    void write_to_svg(const std::string& filename, const std::vector<ch::gray_level>& levels, int wd, int hgt, double scale);

    drawing generate_crosshatched_drawing(const std::string& image_file, segmentation_params params, double scale, brush& br);
    std::vector<polyline> crosshatched_poly_with_holes(const polygon_with_holes& poly, double color, brush& br);
    std::vector<polyline> crosshatched_poly_with_holes(const polygon_with_holes& poly, double color, hierarchical_brush& br);
    void to_svg(const std::string& filename, const drawing& d);

    drawing generate_hierarchical_drawing(cv::Mat image, double scale, const std::vector<brush_fn>& brush_fns, const std::vector<double>& gray_intervals,
        int line_thickness = 1, double epsilon = k_epsilon, dimensions swatch_sz = { k_swatch_sz });
}