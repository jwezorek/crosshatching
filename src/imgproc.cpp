#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include "imgproc.hpp"
#include "meanshift.hpp"
#include <range/v3/all.hpp>
#include <set>
#include <map>
#include <stdexcept>
#include <array>
#include <iostream>

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    std::vector<uchar> get_gray_levels(const cv::Mat& input) {
        if (input.channels() != 1) {
            throw std::runtime_error("called get_gray_levels on color image");
        }
        std::array<bool,256> grays = {};
        input.forEach<uchar>( [&grays](uchar gray, const int* pos) { grays[gray] = true;  } );
        return rv::iota(0) | 
            rv::take(256) | 
            rv::filter([&grays](int g) {return grays[g]; }) | 
            r::to<std::vector<uchar>>();
    }

    std::map<uchar, cv::Mat> gray_planes(const cv::Mat& input) {
        auto levels = get_gray_levels(input);
        return
            levels |
            rv::transform(
                [&input](auto val)->std::map<uchar, cv::Mat>::value_type {
                    cv::Mat output;
                    cv::inRange(input, val, val, output);
                    return { val, output };
                }
            ) |
            r::to<std::map<uchar, cv::Mat>>();
    }

}

cv::Mat ch::do_segmentation(const cv::Mat& input, int sigmaS, float sigmaR, int minSize) {
    cv::Mat output;
    cv::Mat labels;
    auto mss = createMeanShiftSegmentation(sigmaS, sigmaR, minSize, true);
    mss->processImage(input, output, labels);

    cv::Mat gray;
    cv::cvtColor(output, gray, cv::COLOR_BGR2GRAY);

    return gray;
}

void ch::debug() {
    cv::Mat img = cv::imread("C:\\test\\dunjon.png");

    cv::Mat mat = ch::do_segmentation(img, 8.0f, 3.0f, 12);
    auto levels = gray_planes(mat);

    for (const auto [lev, img] : levels) {
        cv::imshow("level", img);
        cv::waitKey(0);
    }

    std::cout << levels.size() << "\n";
}