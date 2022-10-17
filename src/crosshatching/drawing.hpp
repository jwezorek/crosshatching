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
        stroke_groups content;
        dimensions<double> sz;
        size_t stroke_count() const;

        inline auto strokes() const {
            namespace r = ranges;
            namespace rv = ranges::views;
            return content |
                rv::transform(
                    [](const stroke_group& sg) {
                        return sg.strokes |
                            rv::transform(
                                [thickness = sg.thickness]
                                (const polyline& poly)->std::tuple<polyline, double> {
                                    return { poly, thickness };
                                }
                        );
                    }
                ) | rv::join;
        }

    };

    struct parameters {
        double scale;
        double epsilon;
        int swatch_sz;
        bool use_true_black;

        parameters(double sc = 5.0, double eps = k_epsilon,
                int sz = k_swatch_sz, bool black = false) :
            scale(sc),
            epsilon(eps),
            swatch_sz(sz),
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