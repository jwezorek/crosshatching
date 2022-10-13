#include "brush.hpp"
#include "brush_lang.hpp"
#include "geometry.hpp"
#include <opencv2/imgproc.hpp>
#include <qdebug.h>
#include <future>
#include <numeric>
#include <sstream>
/*------------------------------------------------------------------------------------------------*/

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    ch::polygon make_rectangle(const ch::dimensions<double>& sz) {
        return ch::make_polygon(
            { { {0.0,0.0}, {sz.wd, 0.0}, {sz.wd,sz.hgt}, {0.0, sz.hgt} } }
        );
    }

    ch::drawn_strokes stroke_rectangle(ch::brush_expr_ptr brush_expr, const ch::dimensions<double>& sz, double t) {

        auto rect = make_rectangle(sz);
        return ch::to_drawn_strokes(
            ch::brush_expr_to_strokes(brush_expr, rect, t)
        );
    }

    cv::Mat drawn_strokes_to_mat(ch::drawn_strokes strokes, cv::Mat mat) {
        auto qimg = ch::mat_to_qimage(mat, false);
        QPainter g(&qimg);
        g.setRenderHint(QPainter::Antialiasing, true);
        ch::paint_strokes(g, strokes);
        return mat;
    }

    cv::Mat render_brush_expr_to_mat(ch::brush_expr_ptr brush_expr, cv::Mat mat, double t) {
        return drawn_strokes_to_mat(
            stroke_rectangle(  brush_expr, ch::mat_dimensions(mat), t), 
            mat
        );
    }

    double measure_gray_level(cv::Mat swatch, cv::Mat bkgd) {
        double mean_val = cv::mean(swatch).val[0];
        return (255.0 - mean_val) / 255.0;
    }

    template<typename T>
    cv::Mat blank_mat(const ch::dimensions<T>& sz) {
        return cv::Mat(static_cast<int>(sz.hgt), static_cast<int>(sz.wd), CV_8U, 255);
    }

    double sample(ch::dimensions<double> sz, ch::brush_expr_ptr expr, 
            double t, int n, const std::vector<cv::Mat>& bkgds) {

        auto compute_gray_level = [sz](ch::brush_expr_ptr expr, double t, cv::Mat bkgd) -> double {
            auto swatch = render_brush_expr_to_mat(expr, bkgd, t);
            return measure_gray_level(swatch, bkgd);
        };

        std::vector<std::future<double>> samples = rv::iota(0, n) |
            rv::transform(
                [=, &bkgds](int i)->std::future<double> {
                    cv::Mat bkgd = (bkgds.empty()) ? 
                        blank_mat(sz) : 
                        bkgds.at(i);
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

    ch::drawn_strokes clip_drawn_strokes(const ch::polygon& poly, const ch::drawn_strokes& strokes) {
        return strokes |
            rv::transform(
                [&poly](const ch::drawn_stroke_cluster& cluster)->ch::drawn_stroke_cluster {
                    return {
                        ch::clip_polylines_to_poly(cluster.strokes, poly),
                        cluster.thickness
                    };
                }
            ) | r::to_vector;
    }

    ch::drawn_stroke_cluster transform(const ch::drawn_stroke_cluster& dsc, const ch::matrix& mat) {
        auto polyline_range = dsc.strokes |
            rv::transform(
                [&mat](const ch::polyline& poly)->ch::polyline {
                    return rv::all(poly) | 
                        rv::transform(
                            [&mat](ch::point pt) {
                                return ch::transform(pt, mat);
                            }
                        ) |
                        r::to< ch::polyline>;
                }
            );

        ch::polylines output;
        output.resize(dsc.strokes.size());
        for (auto&& [i, poly] : rv::enumerate(polyline_range)) {
            output[i] = std::move(poly);
        }

        return {
            std::move(output),
            dsc.thickness
        };
    }

    ch::drawn_strokes transform(const ch::drawn_strokes& strokes, const ch::matrix& mat) {
        return strokes |
            rv::transform(
                [&mat](const ch::drawn_stroke_cluster& dsc) {
                    return ::transform(dsc, mat);
                }
            ) | r::to_vector;
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

ch::drawn_strokes ch::brush::draw_strokes(const polygon& poly, double gray_level, bool clip_to_poly) {
    
    if (is_uinitiailized()) {
        throw std::runtime_error("brush is not built");
    }
    
    if (gray_level < min_gray_level()) {
        return draw_strokes(poly, min_gray_level(), clip_to_poly);
    } else if (gray_level > max_gray_level()) {
        return draw_strokes(poly, max_gray_level(), clip_to_poly);
    }

    auto param = build_to_gray_level(gray_level);
    auto centroid = mean_point(poly.outer());
    auto canonicalized = ch::transform(poly, ch::translation_matrix(-centroid));

    auto strokes = to_drawn_strokes(
        brush_expr_to_strokes(brush_expr_, canonicalized, param)
    );

    if (clip_to_poly) {
        strokes = clip_drawn_strokes(poly, strokes);
    }

    return ::transform(strokes, ch::translation_matrix(centroid));
}

ch::bkgd_swatches ch::brush::render_swatches(double gray_value, int n) {
    auto bkgds = bkgds_;
    if (bkgds.empty()) {
        bkgds.resize(n, blank_mat(swatch_sz_));
    }

    return rv::iota(0, n) |
        rv::transform(
            [&](int i)->cv::Mat {
                auto bkgd = bkgds.at(i);
                auto swatch = draw_strokes(make_rectangle(swatch_sz_), gray_value, false);
                return drawn_strokes_to_mat(swatch, bkgd);
            }
        ) |
        r::to_vector;
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

