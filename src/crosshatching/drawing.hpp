#pragma once
#include "util.hpp"
#include "brush.hpp"
#include "geometry.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#include <string>
#include <functional>

namespace ch {

    struct drawing {
        std::vector<polyline> strokes;
        dimensions<double> sz;
        double stroke_wd;
    };

    struct parameters {
        double scale;
        int stroke_width; 
        double epsilon;
        int swatch_sz;
        bool use_true_black;

        parameters(double sc = 5.0, int sw = 1, double eps = k_epsilon, 
                int sz = k_swatch_sz, bool black = false) :
            scale(sc),
            stroke_width(sw),
            epsilon(eps),
            swatch_sz(sz),
            use_true_black(black)
        {}
    };

    struct crosshatching_job {
        std::string title;
        std::vector<std::tuple<ch::brush_fn, cv::Mat>> layers;
        parameters params;
    };

    struct callbacks {
        std::function<void(double)> update_progress_cb;
        std::function<void(const std::string&)> update_status_cb;
        std::function<void(const std::string&)> log_message_cb;
    };

    std::vector<std::tuple<ch::brush_fn, cv::Mat>> generate_ink_layer_images(
        cv::Mat img, cv::Mat label_img,
        const std::vector<std::tuple<ch::brush_fn, double>>& brush_intervals
    );
    drawing generate_crosshatched_drawing(const crosshatching_job& job, const callbacks& cbs = {});
    void write_to_svg(const std::string& filename, const drawing& d, 
        std::function<void(double)> update_progress_cb = {});
    drawing scale(const drawing& d, double scale);
    cv::Mat paint_drawing(const drawing& d, std::function<void(double)> update_progress_cb = {});

    namespace detail {
        std::vector<std::tuple<uchar,polygon>> to_blobs_from_1channel_image(const cv::Mat& input);
        std::vector<std::tuple<color, polygon>> to_blobs_from_3channel_image(const cv::Mat& input);
    }

    template<typename T>
    std::vector<std::tuple<T, polygon>> to_blobs(const cv::Mat& img) {
        if constexpr (std::is_same<T, uchar>::value) {
            return detail::to_blobs_from_1channel_image(img);
        } else {
            return detail::to_blobs_from_3channel_image(img);
        }
    }

    void debug(const cv::Mat& mat);
}