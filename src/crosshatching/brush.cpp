#include "brush.hpp"

bool ch::brush::is_uinitiailized() const {
    return false;
}

double ch::brush::get_or_sample_param(double param) {
    return 0;
}
double ch::brush::build_between(double v, gray_map_iter left, gray_map_iter right) {
    return 0;
}

double ch::brush::build_to_gray_level(double gray_level) {
    return 0;
}

ch::brush::brush() {
};

ch::brush::brush(ch::brush_expr_ptr expr, double epsilon,
        ch::dimensions<double> swatch_sz, const ch::bkgd_swatches& bkgds) {

}

void ch::brush::build() {

}

void ch::brush::build_n(int n) {

}

double ch::brush::gray_value_to_param(double gray_val) {
    return 0;
}

ch::crosshatching_swatch ch::brush::get_hatching(double gray_val, ch::dimensions<double> sz) {
    return {};
}

ch::bkgd_swatches ch::brush::render_swatches(double gray_level) {
    return {};
}

cv::Mat ch::brush::swatch(double gray_level) {
    return {};
}

double ch::brush::min_gray_level() const {
    return 0;
}

double ch::brush::max_gray_level() const {
    return 0;
}

int ch::brush::num_samples() {
    return 0;
}