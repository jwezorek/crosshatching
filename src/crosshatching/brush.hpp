#pragma once

#include <map>
#include <unordered_map>
#include <opencv2/core.hpp>
#include <range/v3/all.hpp>
#include "brush_lang.hpp"

/*------------------------------------------------------------------------------------------------------*/


namespace ch {

    using bkgd_swatches = std::vector<cv::Mat>;

    constexpr int k_swatch_sz = 512;
    constexpr double k_epsilon = 0.001;

    struct swatch {
        strokes content;
        dimensions<double> sz;
    };

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
        double build_between(double v, gray_map_iter left, gray_map_iter right);
        double build_to_gray_level(double gray_level);

    public:
        brush();
        brush(brush_expr_ptr expr, double epsilon = k_epsilon,
            dimensions<double> swatch_sz = { static_cast<double>(k_swatch_sz) }, 
            const bkgd_swatches& bkgds = {});

        void build_n(int n);
        double gray_value_to_param(double gray_val);
        swatch get_hatching(double gray_val, dimensions<double> sz);
        bkgd_swatches render_swatches(double gray_level);
        cv::Mat swatch(double gray_level);
        double min_gray_level() const;
        double max_gray_level() const;
        static int num_samples();
    };
}