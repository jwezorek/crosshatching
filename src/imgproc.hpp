#pragma once
#include "geometry.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#include <string>

namespace ch {

    struct polygon_with_holes {
        polyline border;
        std::vector<polyline> holes;
    };

    struct gray_level_plane {
        uchar gray;
        std::vector<polygon_with_holes> blobs;
    };

    cv::Mat do_segmentation(const cv::Mat& input, int sigmaS, float sigmaR, int minSize);
    std::vector<gray_level_plane> extract_gray_level_planes(const cv::Mat& gray_scale_img);
    void write_to_svg(const std::string& filename, const std::vector<gray_level_plane>& levels, int wd, int hgt, double scale);

    void debug();
}