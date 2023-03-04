#pragma once

#include "brush_lang.hpp"
#include <opencv2/core.hpp>
#include <range/v3/all.hpp>
#include <map>
#include <unordered_map>
#include <memory>
#include <deque>

/*------------------------------------------------------------------------------------------------*/

namespace ch {

    using bkgd_swatches = std::vector<cv::Mat>;

    constexpr int k_swatch_sz = 600;
    constexpr double k_epsilon = 0.001;
    constexpr auto k_default_num_samples = 4;

    class brush {
    private:
        using gray_map_iter = std::map<double, double>::iterator;
        std::map<double, double> gray_to_param_;
        std::unordered_map<double, double> param_to_gray_;
        brush_expr_ptr brush_expr_;
        int line_thickness_;
        dimensions<double> swatch_sz_;
        double eps_;
        bkgd_swatches bkgds_;
        int num_samples_;

        bool is_uinitiailized() const;
        double get_or_sample_param(double param);
        double build_between(double v, gray_map_iter left, gray_map_iter right,
            int depth = 0);
        double build_to_gray_level(double gray_level);

    public:
        brush();
        brush(brush_expr_ptr expr, double epsilon = k_epsilon, 
            int num_samples = k_default_num_samples,
            dimensions<float> swatch_sz = { static_cast<float>(k_swatch_sz) }, 
            const bkgd_swatches& bkgds = {});

        void build_n(int n);
        double gray_value_to_param(double gray_val);
        bkgd_swatches render_swatches(double gray_level, int n);
        drawing_comps_ptr draw_strokes(const brush_context& ctxt, bool clip_to_poly = true);
        double min_gray_level() const;
        double max_gray_level() const;
        int num_samples() const;
        int swatch_dim() const;
    };

    cv::Mat test_brush(brush_expr_ptr expr);

    using brush_ptr = std::shared_ptr<brush>;
}