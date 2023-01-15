#pragma once
#include "util.hpp"
#include "strokes.hpp"
#include "brush.hpp"
#include "geometry.hpp"
#include "ink_layers.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#include <string>
#include <functional>

/*------------------------------------------------------------------------------------------------*/

namespace ch {

    struct drawing {
        std::vector<polygon> black_regions_;
        drawing_comps content;
        dimensions<double> sz;
        size_t stroke_count() const;
    };

    struct parameters {
        double scale;
        double epsilon;
        int swatch_sz;
        int num_samples;
        bool use_true_black;

        parameters(double sc = 10.0, double eps = k_epsilon,
                int sz = k_swatch_sz, int ns = k_default_num_samples, 
                bool black = false) :
            scale(sc),
            epsilon(eps),
            swatch_sz(sz),
            num_samples(ns),
            use_true_black(black)
        {}
    };

    struct crosshatching_job {
        std::string title;
        ink_layers layers;
        parameters params;
    };

    struct callbacks {
        std::function<void(double)> update_progress_cb;
        std::function<void(const std::string&)> update_status_cb;
        std::function<void(const std::string&)> log_message_cb;
    };

    drawing generate_crosshatched_drawing(const crosshatching_job& job, const callbacks& cbs = {});
    void write_to_svg(const std::string& filename, const drawing& d, 
        std::function<void(double)> update_progress_cb = {});
    //drawing scale(const drawing& d, double scale);
    cv::Mat paint_drawing(const drawing& d, std::function<void(double)> update_progress_cb = {});

    void debug_drawing();
}