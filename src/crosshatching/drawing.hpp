#pragma once
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
        dimensions sz;
        double stroke_wd;
    };

    struct crosshatching_params {
        double scale;
        int stroke_width; 
        double epsilon;
        double swatch_sz;

        crosshatching_params(double sc = 1.0, int sw = 1, double eps = k_epsilon, double sz = k_swatch_sz) :
            scale(sc),
            stroke_width(sw),
            epsilon(eps),
            swatch_sz(sz)
        {}
    };

    struct crosshatching_job {
        cv::Mat img; 
        cv::Mat label_img;
        std::vector<std::tuple<ch::brush_fn, double>> layers;
        crosshatching_params params;
    };

    struct prog_callbacks {
        std::function<void(double)> update_progress_cb;
        std::function<void(const std::string&)> update_status_cb;
        std::function<void(const std::string&)> log_message_cb;
    };

    drawing generate_crosshatched_drawing(const crosshatching_job& job, const prog_callbacks& cbs = {});
    drawing generate_crosshatched_drawing(cv::Mat img, cv::Mat label_img, const std::vector<std::tuple<ch::brush_fn, double>>& brushes, const crosshatching_params& params = {}, const prog_callbacks& cbs = {});
    drawing generate_crosshatched_drawing(cv::Mat img, const std::vector<std::tuple<ch::brush_fn, double>>& brushes, const crosshatching_params& params = {}, const prog_callbacks& cbs = {});
    void to_svg(const std::string& filename, const drawing& d);
    void debug(cv::Mat img, cv::Mat labels);

}