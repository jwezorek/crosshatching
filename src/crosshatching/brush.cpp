#include "brush.hpp"
#include "geometry.hpp"
#include <opencv2/imgproc.hpp>
#include <future>
#include <numeric>

/*------------------------------------------------------------------------------------------------*/

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    constexpr auto k_num_samples = 4;

    ch::swatch create_swatch(ch::brush_expr_ptr expr, const ch::dimensions<double>& sz, double t) {
        double half_wd = sz.wd / 2.0;
        double half_hgt = sz.hgt / 2.0;
        std::array<ch::point, 4> rect_verts = {
            ch::point{-half_wd,-half_hgt}, ch::point{half_wd,-half_hgt},
            ch::point{half_wd,half_hgt}, ch::point{-half_wd,half_hgt},
        };
        auto strokes = ch::brush_expr_to_strokes(expr, ch::make_polygon(rect_verts), t);
        return {
            ch::transform(strokes, ch::translation_matrix(half_wd, half_hgt)),
            sz
        };
    }

    void paint_stroke(cv::Mat mat, ch::stroke stroke) {
        auto poly = stroke.polyline | r::to_vector;
        std::vector<cv::Point> int_pts(poly.size());
        std::transform(poly.begin(), poly.end(), int_pts.begin(),
            [](const auto& p)->cv::Point {
                return cv::Point(
                    static_cast<int>(std::round(p.x)),
                    static_cast<int>(std::round(p.y))
                );
            }
        );
        auto npts = int_pts.size();
        cv::polylines(mat, int_pts, false, 0, stroke.pen_thickness, 8, 0);
    }

    cv::Mat paint_swatch(ch::swatch swatch, cv::Mat bkgd) {
        cv::Mat mat = (bkgd.empty()) ?
            cv::Mat(static_cast<int>(swatch.sz.hgt), static_cast<int>(swatch.sz.wd), CV_8U, 255) :
            bkgd.clone();
        auto qimg = ch::mat_to_qimage(mat, false);
        QPainter g(&qimg);
        g.setRenderHint(QPainter::Antialiasing, true);
        ch::paint_strokes(g, swatch.content);
        return mat;
    }

    double measure_gray_level(ch::swatch swatch, cv::Mat bkgd) {
        auto mat = paint_swatch(swatch, bkgd);
        double mean_val = cv::mean(mat).val[0];
        return (255.0 - mean_val) / 255.0;
    }

    double sample(ch::dimensions<double> sz, ch::brush_expr_ptr expr, 
            double t, int n, const std::vector<cv::Mat>& bkgds) {

        auto compute_gray_level = [sz](ch::brush_expr_ptr expr, double t, cv::Mat bkgd) -> double {
            return measure_gray_level(create_swatch(expr, sz, t), bkgd);
        };

        std::vector<std::future<double>> samples = rv::iota(0, n) |
            rv::transform(
                [=, &bkgds](int i)->std::future<double> {
                    cv::Mat bkgd = (bkgds.empty()) ? cv::Mat() : bkgds.at(i);
                    return std::async(std::launch::async, compute_gray_level, expr, t, bkgd);
                }
            ) | r::to_vector;

        double sum = std::accumulate(samples.begin(), samples.end(), 0.0,
            [](double s, std::future<double>& fut) {
                return s + fut.get();
            }
        );

        auto gray = sum / n;
        return gray;
    }

}

bool ch::brush::is_uinitiailized() const {
    return gray_to_param_.empty();
}

double ch::brush::get_or_sample_param(double param) {
    auto iter = param_to_gray_.find(param);
    if (iter != param_to_gray_.end()) {
        return iter->second;
    }
    auto gray = sample(swatch_sz_, brush_expr_, param, k_num_samples, bkgds_);
    gray_to_param_[gray] = param;
    param_to_gray_[param] = gray;

    return gray;
}

double ch::brush::build_between(double gray, gray_map_iter left, gray_map_iter right) {
    double mid_param = (left->second + right->second) / 2.0;
    auto mid_gray_value = get_or_sample_param(mid_param);
    if (std::abs(mid_gray_value - gray) < eps_) {
        return mid_param;
    }
    auto mid_iter = gray_to_param_.find(mid_gray_value);
    if (gray < mid_gray_value) {
        return build_between(gray, left, mid_iter);
    } else {
        return build_between(gray, mid_iter, right);
    }
}

double ch::brush::build_to_gray_level(double gray_level) {
    auto right = gray_to_param_.lower_bound(gray_level);
    if (right == gray_to_param_.end()) {
        return gray_to_param_.rbegin()->second;
    }

    auto left = (right != gray_to_param_.begin()) ? std::prev(right) : right;
    if (std::abs(gray_level - left->first) < eps_) {
        return left->second;
    }
    if (std::abs(gray_level - right->first) < eps_) {
        return right->second;
    }
    return build_between(gray_level, left, right);
}

ch::brush::brush() {
};

ch::brush::brush(ch::brush_expr_ptr expr, double epsilon,
        ch::dimensions<double> swatch_sz, const ch::bkgd_swatches& bkgds) :
    brush_expr_(expr),
    swatch_sz_{ swatch_sz },
    eps_(epsilon),
    bkgds_(bkgds) {

}

void ch::brush::build_n(int n) {
    if (!is_uinitiailized()) {
        throw std::runtime_error("brush is already built");
    }
    double delta = 1.0 / (static_cast<double>(n) - 1.0);
    for (int i = 0; i < n; ++i) {
        get_or_sample_param(delta * i);
    }
}

double ch::brush::gray_value_to_param(double gray_val) {
    if (is_uinitiailized()) {
        throw std::runtime_error("brush is not built");
    }
    return build_to_gray_level(gray_val);
}

ch::swatch ch::brush::get_hatching(double gray_level, ch::dimensions<double> sz) {
    if (is_uinitiailized()) {
        throw std::runtime_error("brush is not built");
    }
    if (gray_level < min_gray_level()) {
        return get_hatching(min_gray_level(), sz);
    } else if (gray_level > max_gray_level()) {
        return get_hatching(max_gray_level(), sz);
    }
    auto param = build_to_gray_level(gray_level);
    return create_swatch(brush_expr_, sz, param);
}

ch::bkgd_swatches ch::brush::render_swatches(double gray_value) {
    auto bkgds = bkgds_;
    if (bkgds.empty()) {
        bkgds.resize(k_num_samples);
    }
    return rv::iota(0, k_num_samples) |
        rv::transform(
            [&](int i)->cv::Mat {
                auto bkgd = bkgds.at(i);
                auto swatch = this->get_hatching(gray_value, swatch_sz_);
                return paint_swatch(swatch, bkgd);
            }
        ) |
        r::to_vector;
}

cv::Mat ch::brush::swatch(double gray_value) {
    cv::Mat bkgd = bkgds_.empty() ? cv::Mat() : bkgds_.at(0);
    auto swatch = this->get_hatching(gray_value, swatch_sz_);
    return paint_swatch(swatch, bkgd);
}

double ch::brush::min_gray_level() const {
    if (is_uinitiailized()) {
        throw std::runtime_error("brush is not built");
    }
    return gray_to_param_.begin()->first;
}

double ch::brush::max_gray_level() const {
    if (is_uinitiailized()) {
        throw std::runtime_error("brush is not built");
    }
    return std::prev(gray_to_param_.end())->first;
}

int ch::brush::num_samples() {
    return k_num_samples;
}