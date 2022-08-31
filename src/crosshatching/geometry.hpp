#pragma once

#include <vector>
#include <sstream>
#include <optional>
#include <unordered_map>
#include <span>
#include <unordered_set>
#include <opencv2/core.hpp>
#include <Eigen/Dense>
#include <range/v3/all.hpp>
#include <boost/functional/hash.hpp>
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
    using int_point = cv::Point;

    template <typename T>
    struct dimensions {
        T wd;
        T hgt;

        dimensions(T d = {}) : wd(d), hgt(d) {}
        dimensions(T w, T h) : wd(w), hgt(h) {}

        template<typename U>
        dimensions(U w, U h) : 
            wd( static_cast<T>(w)), hgt(static_cast<T>(h))
        {}

        template<typename U>
        dimensions(const dimensions<U>& d) :
            wd(static_cast<T>(d.wd)), hgt(static_cast<T>(d.hgt))
        {}

        T area() const {
            return wd * hgt;
        }
    };

    template <typename T>
    dimensions<T> operator*(T left, const dimensions<T>& right) {
        return { left * right.wd, left * right.hgt };
    }

    polyline make_polyline(size_t sz);
    ring make_ring(size_t);
    polygon make_polygon(const ring& outer, const std::vector<ring>& inners);
    matrix rotation_matrix(double theta);
    matrix rotation_matrix(double cos_theta, double sin_theta);
    matrix translation_matrix(double x, double y);
    matrix translation_matrix(const point& pt);
    matrix scale_matrix(double x_scale, double y_scale);
    ranges::any_view<polyline> transform_view(ranges::any_view<polyline> polys, const matrix& mat);
    std::vector<polyline> transform( std::span<const polyline> polys, const matrix& mat);
    point mean_point(std::span<const point> points);
    polyline transform(std::span<const point>, const matrix& mat);
    ring transform(const ring& r, const matrix& mat);
    point transform(const point& pt, const matrix& mat);
    polygon transform(const polygon& poly, const matrix& mat);
    void paint_polyline(cv::Mat& mat, const polyline& p, double thickness, 
        int color, point offset = { 0,0 });
    void paint_polyline_aa(cv::Mat& mat, const polyline& p, double thickness, 
        int color, point offset = { 0,0 });
    double euclidean_distance(const point& pt1, const point& pt2);
    polylines clip_lines_to_poly(const polylines& strokes, const polygon& poly);
    ring scale(const ring& r, double scale);
    polyline scale(const polyline& poly, double scale);
    polygon scale(const polygon& poly, double scale);
    ch::ring simplify_rectilinear_ring(const ring& poly);
    rectangle bounding_rectangle(const polyline& poly);
    rectangle bounding_rectangle(const ring& poly);
    rectangle bounding_rectangle(const polygon& poly);
    rectangle bounding_rectangle(std::span<const polygon> polys);
    rectangle union_rect(const rectangle& r1, const rectangle& r2);
    cv::Rect union_rect_and_pt(const cv::Rect& r, cv::Point2i pt);
    std::optional<line_segment> linesegment_rectangle_intersection(
        const line_segment& line_seg, const rectangle& rect);
    ch::point southeast_most_point(std::span<const ch::point> points);
    std::vector<ch::polygon> simplify_polygons(
        std::span<const ch::polygon> dissection, double param);
    size_t vert_count(const ch::polygon& poly);
    std::vector<polygon> buffer(const ch::polygon& poly, double amt);
    std::vector<polygon> buffer(std::span<const ch::polygon> polys, double amt);
    std::vector<ch::point> convex_hull(std::span<const ch::point> points);

    template<typename... Args>
    auto first(const std::tuple<Args...>& tup)->decltype(auto) {
        return std::get<0>(tup);
    }

    template<typename... Args>
    auto second(const std::tuple<Args...>& tup)->decltype(auto) {
        return std::get<1>(tup);
    }

    template<typename T>
    std::vector<std::tuple<T, polygon>> simplify_colored_polygons(
            std::span< std::tuple<T, polygon>> blobs, double param) {
        namespace r = ranges;
        namespace rv = ranges::views;

        auto polys = blobs |
            rv::transform( second<T,polygon> ) |
            r::to_vector;

        polys = simplify_polygons(polys, param);

        return rv::zip(
                blobs | rv::transform( first<T, polygon> ),
                polys
            ) | rv::transform(
                 [](const auto& pair)->std::tuple<T, polygon> {
                     return { pair.first, pair.second };
                 }
            ) | r::to_vector;
    }

    void debug_geom(cv::Mat mat);
}