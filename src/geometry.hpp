#pragma once

#include <opencv2/core.hpp>
#include <Eigen/Dense>
#include <vector>
#include <range/v3/all.hpp>

namespace ch {

    using point = cv::Point2d;
    using matrix = Eigen::Matrix<double, 3, 3>;
    using vec = Eigen::Matrix<double, 3, 1>;
    using polyline = std::vector<point>;
    using rnd_fn = std::function<double()>;

    matrix rotation_matrix(double theta);
    matrix rotation_matrix(double cos_theta, double sin_theta);
    matrix translation_matrix(double x, double y);
    matrix scale_matrix(double x_scale, double y_scale);

    ranges::any_view<polyline> transform(ranges::any_view<polyline> polys, const matrix& mat);
    std::vector<polyline> transform(const std::vector<polyline>& poly, const matrix& mat);
    polyline transform(const polyline& poly, const matrix& mat);
    point transform(const point& pt, const matrix& mat);
    void paint_polyline(cv::Mat& mat, const polyline& p, int thickness, int color, point offset = { 0,0 });

    double euclidean_distance(const point& pt1, const point& pt2);
    double normal_rnd(double mean, double stddev);
    double uniform_rnd(double lower_bound, double upper_bound);
    rnd_fn normal_rnd_fn(double mean, double stddev);
    ch::rnd_fn const_rnd_fn(double val);
}