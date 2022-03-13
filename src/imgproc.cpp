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

    cv::Mat dilate_blob(const cv::Mat& blob) {
        static cv::Mat strel = (
            cv::Mat_<uchar>(3, 3) <<
                0, 0, 0,
                0, 1, 1,
                0, 1, 1
        );
        cv::Mat output;
        cv::dilate(blob, output, strel, { 1,1 });
        return output;
    }

    std::map<uchar, contours_trees> gray_level_contours(const cv::Mat& input) {
        auto planes = gray_planes(input);
        return planes |
            rv::transform(
                [](const auto& pair)->std::map<uchar, contours_trees>::value_type {
                    const auto& [gray, mat] = pair;
                    auto dilated_mat = dilate_blob(mat);
                    contours_trees output;
                    cv::findContours(dilated_mat, output.contours, output.hierarchy, cv::RETR_CCOMP, cv::CHAIN_APPROX_SIMPLE);
                    return { gray, std::move(output) };
                }
            ) |
            r::to<std::map<uchar, contours_trees>>();
    }

    std::vector<std::vector<cv::Point>> apply_douglas_peucker(const std::vector<std::vector<cv::Point>>& polys, float param) {
        return polys |
            rv::transform(
                [param](const std::vector<cv::Point>& poly)->std::vector<cv::Point> {
                    std::vector<cv::Point> output;
                    cv::approxPolyDP(poly, output, param, true);
                    return output;
                }
            ) |
            r::to<std::vector<std::vector<cv::Point>>>();
    }

    void apply_douglas_peucker(std::map<uchar, contours_trees>& ct, float param) {
        for (auto& [k, v] : ct) {
            v.contours = apply_douglas_peucker(v.contours, param);
        }
    }

    void scale(std::map<uchar, contours_trees>& cont_trees, double s) {
        for (auto& [k, cont_tree] : cont_trees) {
            for (auto& polyline : cont_tree.contours) {
                for (auto& pt : polyline) {
                    pt = { 
                        static_cast<int>(std::round(s * pt.x)), 
                        static_cast<int>(std::round(s * pt.y))
                    };
                }
            }
        }
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

    cv::Mat mat = ch::do_segmentation(img, 8, 3.0f, 12);
    auto contours = gray_level_contours(mat);
    apply_douglas_peucker(contours, 0.5);
    scale(contours, 3.0);
    auto output = paint_contours(3*mat.cols, 3*mat.rows, contours);

    cv::imshow("contours", output);
}