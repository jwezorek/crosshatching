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

    struct gray_level_plane {
        uchar gray;
        std::vector<polygon_with_holes> blobs;
    };

    struct gray_levels {
        double value;
        std::vector<polygon_with_holes> blobs;
    };

    std::vector<gray_level_plane> extract_gray_level_planes(const cv::Mat& gray_scale_img);
    std::vector<gray_levels> extract_gray_levels(const cv::Mat& gray_scale_img);
    std::vector<gray_level_plane> scale(const std::vector<gray_level_plane>& planes, double scale);
    void write_to_svg(const std::string& filename, const std::vector<gray_level_plane>& levels, int wd, int hgt, double scale);
    void debug();

    drawing generate_crosshatched_drawing(const std::string& image_file, segmentation_params params, double scale, brush& br);
    std::vector<polyline> crosshatched_poly_with_holes(const polygon_with_holes& poly, double color, brush& br);
    void to_svg(const std::string& filename, const drawing& d);

}