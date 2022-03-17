#pragma once

#include <string>
#include <functional>

namespace ch {

    std::string svg_header(int wd, int hgt, bool bkgd_rect = false);
    std::string gray_to_svg_color(unsigned char gray);

    using rnd_fn = std::function<double()>;
    double normal_rnd(double mean, double stddev);
    double uniform_rnd(double lower_bound, double upper_bound);
    rnd_fn normal_rnd_fn(double mean, double stddev);
    ch::rnd_fn const_rnd_fn(double val);

}