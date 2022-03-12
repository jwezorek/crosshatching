#pragma once

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

namespace ch {

    cv::Mat do_segmentation(const cv::Mat& input, int sigmaS, float sigmaR, int minSize);
    void debug();
}