#pragma once

#include <opencv2/core.hpp>
#include <Eigen/Dense>
#include <vector>
#include <range/v3/all.hpp>
#include <sstream>

namespace ch {

    using point = cv::Point2d;
    using matrix = Eigen::Matrix<double, 3, 3>;
    using vec = Eigen::Matrix<double, 3, 1>;
    using polyline = std::vector<point>;

    matrix rotation_matrix(double theta);
    matrix rotation_matrix(double cos_theta, double sin_theta);
    matrix translation_matrix(double x, double y);
    matrix scale_matrix(double x_scale, double y_scale);

    ranges::any_view<polyline> transform(ranges::any_view<polyline> polys, const matrix& mat);
    std::vector<polyline> transform(const std::vector<polyline>& poly, const matrix& mat);
    point mean_point(const polyline& poly);
    polyline transform(const polyline& poly, const matrix& mat);
    point transform(const point& pt, const matrix& mat);
    void paint_polyline(cv::Mat& mat, const polyline& p, double thickness, int color, point offset = { 0,0 });

    double euclidean_distance(const point& pt1, const point& pt2);

    template<typename P>
    std::string poly_to_string(const std::vector<P>& polyline) {
        std::stringstream ss;
        ss << "[ ";
        for (const auto& pt : polyline) {
            ss << pt.x << "," << pt.y << " ";
        }
        ss << "]";
        return ss.str();
    }

    template<typename T>
    cv::Point_<T> normalize_offset(const cv::Point_<T>& pt) {
        auto offset = pt;
        offset /= std::max(std::abs(pt.x), std::abs(pt.y));
        return offset;
    }

    struct point_hasher
    {
        std::size_t operator()(const cv::Point& p) const;
    };

}