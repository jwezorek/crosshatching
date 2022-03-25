#pragma once

#include <string>
#include <functional>
#include <range/v3/all.hpp>
#include <opencv2/core.hpp>

namespace ch {

    using point = cv::Point2d;
    using polyline = std::vector<point>;

    std::string svg_header(int wd, int hgt, bool bkgd_rect = false);
    std::string polyline_to_svg(const ch::polyline& poly, int thickness);
    std::string gray_to_svg_color(unsigned char gray);

    using rnd_fn = std::function<double()>;
    double normal_rnd(double mean, double stddev);
    double uniform_rnd(double lower_bound, double upper_bound);
    rnd_fn normal_rnd_fn(double mean, double stddev);
    ch::rnd_fn const_rnd_fn(double val);

    double ramp(double t, double k, bool right, bool up);

    template<typename R>
    auto rotate_view(R rng, int pivot) {
        return ranges::views::concat(
            rng | ranges::views::drop(pivot), 
            rng | ranges::views::take(pivot)
        );
    }

}