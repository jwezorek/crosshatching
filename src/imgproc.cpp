#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include "imgproc.hpp"
#include "MeanShiftSegmentation.hpp"

cv::Mat ch::do_segmentation(const cv::Mat& input, float sigmaS, float sigmaR, int minSize) {
    cv::Mat output;
    cv::Mat labels;
    auto mss = createMeanShiftSegmentation(sigmaS, sigmaR, minSize, true);
    mss->processImage(input, output, labels);

    cv::Mat gray;
    cv::cvtColor(output, gray, cv::COLOR_BGR2GRAY);

    return gray;
}