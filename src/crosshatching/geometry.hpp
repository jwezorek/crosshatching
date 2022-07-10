#pragma once

#include <opencv2/core.hpp>
#include <Eigen/Dense>
#include <vector>
#include <range/v3/all.hpp>
#include <sstream>
#include <optional>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/multi/geometries/multi_linestring.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/ring.hpp>
#include <boost/geometry/geometries/register/linestring.hpp>
#include <boost/geometry/geometries/register/multi_linestring.hpp>

namespace ch { using point = cv::Point2d; }
BOOST_GEOMETRY_REGISTER_POINT_2D(ch::point, double, boost::geometry::cs::cartesian, x, y);

namespace ch {

    using polygon = boost::geometry::model::polygon<point, true, false>;
    using ring = boost::geometry::model::ring<point, true, false>;
    using polyline = boost::geometry::model::linestring<point>;
    using polylines = boost::geometry::model::multi_linestring<polyline>;
    using line_segment = std::tuple<point, point>;
    using rectangle = std::tuple<double, double, double, double>;
    using matrix = Eigen::Matrix<double, 3, 3>;

    polyline make_polyline(size_t sz);
    ring make_ring(size_t);
    polygon make_polygon(const ring& outer, const std::vector<ring>& inners);
    matrix rotation_matrix(double theta);
    matrix rotation_matrix(double cos_theta, double sin_theta);
    matrix translation_matrix(double x, double y);
    matrix translation_matrix(const point& pt);
    matrix scale_matrix(double x_scale, double y_scale);
    ranges::any_view<polyline> transform(ranges::any_view<polyline> polys, const matrix& mat);
    std::vector<polyline> transform(const std::vector<polyline>& poly, const matrix& mat);
    point mean_point(const polyline& poly);
    polyline transform(const polyline& poly, const matrix& mat);
    ring transform(const ring& r, const matrix& mat);
    point transform(const point& pt, const matrix& mat);
    polygon transform(const polygon& poly, const matrix& mat);
    void paint_polyline(cv::Mat& mat, const polyline& p, double thickness, int color, point offset = { 0,0 });
    void paint_polyline_aa(cv::Mat& mat, const polyline& p, double thickness, int color, point offset = { 0,0 });
    double euclidean_distance(const point& pt1, const point& pt2);
    polylines clip_lines_to_poly(const polylines& strokes, const polygon& poly);
    ring scale(const ring& r, double scale);
    polyline scale(const polyline& poly, double scale);
    polygon scale(const polygon& poly, double scale);
    ch::ring simplify_rectilinear_ring(const ring& poly);
    rectangle bounding_rectangle(const polyline& poly);
    rectangle bounding_rectangle(const ring& poly);
    cv::Rect union_rect_and_pt(const cv::Rect& r, cv::Point2i pt);
    std::optional<line_segment> linesegment_rectangle_intersection(const line_segment& line_seg, const rectangle& rect);

}