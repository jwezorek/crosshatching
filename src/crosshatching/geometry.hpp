#pragma once

#include <vector>
#include <sstream>
#include <optional>
#include <unordered_map>
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
    rectangle bounding_rectangle(const polygon& poly);
    rectangle bounding_rectangle(const std::vector<polygon>& polys);
    rectangle union_rect(const rectangle& r1, const rectangle& r2);
    cv::Rect union_rect_and_pt(const cv::Rect& r, cv::Point2i pt);
    std::optional<line_segment> linesegment_rectangle_intersection(const line_segment& line_seg, const rectangle& rect);
    std::vector<ch::polygon> simplify_rectangle_dissection(const std::vector<ch::polygon>& dissection, 
        const dimensions<double>& rect, double param);

    namespace detail {
        struct point_hasher
        {
            size_t operator()(const int_point& p) const
            {
                std::size_t seed = 0;
                boost::hash_combine(seed, p.x);
                boost::hash_combine(seed, p.y);

                return seed;
            }
        };
    }
    template<typename V>
    using point_map = std::unordered_map<int_point, V, detail::point_hasher>;
    using point_set = std::unordered_set<int_point, detail::point_hasher>;

    void debug_geom(cv::Mat mat, const std::vector<polygon>& polygons);
}