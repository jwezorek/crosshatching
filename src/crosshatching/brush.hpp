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

    constexpr int k_swatch_sz = 800;
    constexpr double k_epsilon = 0.01;
    constexpr auto k_num_samples = 8;

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

        bool is_uinitiailized() const;
        double get_or_sample_param(double param);
        double build_between(double v, gray_map_iter left, gray_map_iter right,
            std::deque<double>& history);
        double build_to_gray_level(double gray_level);

    public:
        brush();
        brush(brush_expr_ptr expr, double epsilon = k_epsilon,
            dimensions<double> swatch_sz = { static_cast<double>(k_swatch_sz) }, 
            const bkgd_swatches& bkgds = {});

        void build_n(int n);
        double gray_value_to_param(double gray_val);
        bkgd_swatches render_swatches(double gray_level, int n = k_num_samples);
        strokes_ptr draw_strokes(const polygon& poly, double gray_level, bool clip_to_poly = true);
        double min_gray_level() const;
        double max_gray_level() const;
        static int num_samples();
    };

    using brush_ptr = std::shared_ptr<brush>;
}