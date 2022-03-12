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

    struct contours_trees {
        std::vector<std::vector<cv::Point>> contours;
        std::vector<cv::Vec4i> hierarchy;
    };

    std::map<uchar, contours_trees> gray_level_contours(const cv::Mat& input) {
        auto planes = gray_planes(input);
        return planes |
            rv::transform(
                [](const auto& pair)->std::map<uchar, contours_trees>::value_type {
                    const auto& [gray, mat] = pair;
                    contours_trees output;
                    cv::findContours(mat, output.contours, output.hierarchy, cv::RETR_CCOMP, cv::CHAIN_APPROX_SIMPLE);
                    return { gray, std::move(output) };
                }
            ) |
            r::to<std::map<uchar, contours_trees>>();
    }

    cv::Mat paint_contours(int wd, int hgt, const std::map<uchar, contours_trees>& cont) {
        auto img = cv::Mat(hgt, wd, CV_8UC1, 255);
        for (const auto& [gray, c] : cont) {
            int n = static_cast<int>(c.contours.size());
            for ( int i = 0; i < n; ++i) {
                cv::drawContours(img, c.contours, i, gray, cv::FILLED, 8, c.hierarchy);
            }
        }
        return img;
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
    auto contours = gray_level_contours(mat);
    auto output = paint_contours(mat.cols, mat.rows, contours);

    cv::imshow("contours", output);
}