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
#include <boost/geometry/algorithms/envelope.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/multi/geometries/multi_linestring.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/ring.hpp>
#include <boost/geometry/geometries/register/linestring.hpp>
#include <boost/geometry/geometries/register/multi_linestring.hpp>

/*------------------------------------------------------------------------------------------------*/

namespace ch { using point = cv::Point2f; }
BOOST_GEOMETRY_REGISTER_POINT_2D(ch::point, double, boost::geometry::cs::cartesian, x, y);

namespace ch {

    using polygon = boost::geometry::model::polygon<point, true, false>;
    using polygons = boost::geometry::model::multi_polygon<polygon>;
    using ring = boost::geometry::model::ring<point, true, false>;
    using polyline = boost::geometry::model::linestring<point>;
    using polylines = boost::geometry::model::multi_linestring<polyline>;
    using line_segment = boost::geometry::model::segment<point>;
    using box = boost::geometry::model::box<point>;
    using rectangle = std::tuple<float, float, float, float>;
    using matrix = Eigen::Matrix<float, 3, 3>;
    using int_point = cv::Point;

    template<typename T, typename F>
    class polygon_tree {

    private:

        
        using value_type = std::pair<box, T>;
        using poly_rtree = boost::geometry::index::rtree<value_type, boost::geometry::index::rstar<16>>;

        poly_rtree impl_;
    public:
        void insert(T val) {
            box b;
            F get_poly;
            boost::geometry::envelope(get_poly(val), b);
            impl_.insert( value_type{ b, val } );
        }
        template<typename U>
        std::vector<T> query(const U& inp) const {
            namespace r = ranges;
            namespace rv = ranges::views;
            namespace bg = boost::geometry;
            namespace bgi = boost::geometry::index;
            std::vector<value_type> result;
            impl_.query(bgi::intersects(inp), std::back_inserter(result));
            return result | 
                rv::remove_if(
                    [&inp](const auto& p) {
                        F get_poly;
                        return !bg::intersects(inp, get_poly(p.second));
                    }
                ) | rv::transform(
                    [](const auto& p) {return p.second; }
                ) | r::to_vector;
        }
        std::optional<T> find(const point& pt) const {
            namespace r = ranges;
            namespace rv = ranges::views;
            namespace bg = boost::geometry;
            namespace bgi = boost::geometry::index;
            F get_poly;
            std::vector<value_type> result;
            impl_.query(bgi::intersects(pt), std::back_inserter(result));
            for (const auto& [k, v] : result) {
                if (bg::within(pt, get_poly(v))) {
                    return v;
                }
            }
            return {};
        }
    };

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
    dimensions<T> operator*(int left, const dimensions<T>& right) {
        return { left * right.wd, left * right.hgt };
    }

    template <typename T>
    dimensions<T> operator*(double left, const dimensions<T>& right) {
        return { static_cast<T>(left * right.wd), static_cast<T>(left * right.hgt) };
    }

    polyline make_polyline(size_t sz);
    polyline make_polyline(std::span<const point> pts);
    ring make_ring(size_t);
    polygon make_polygon(const ring& outer, const std::vector<ring>& inners);
    polygon make_polygon(std::span<const point> verts);
    matrix rotation_matrix(double theta);
    matrix rotation_matrix(double cos_theta, double sin_theta);
    matrix translation_matrix(double x, double y);
    matrix translation_matrix(const point& pt);
    matrix scale_matrix(double x_scale, double y_scale);

    point mean_point(std::span<const point> points);
    polyline transform(std::span<const point>, const matrix& mat);
    ring transform(const ring& r, const matrix& mat);
    point transform(const point& pt, const matrix& mat);
    polygon transform(const polygon& poly, const matrix& mat);
    polylines transform(const polylines& pl, const matrix& mat);
    void paint_polyline(cv::Mat& mat, const polyline& p, double thickness, 
        int color, point offset = { 0,0 });
    void paint_polyline_aa(cv::Mat& mat, const polyline& p, double thickness, 
        int color, point offset = { 0,0 });
    double euclidean_distance(const point& pt1, const point& pt2);
    polylines clip_polylines_to_poly(const polylines& strokes, const polygon& poly);
    polylines clip_polyline_to_poly(const polyline& stroke, const polygon& poly);
    ring scale(const ring& r, double scale);
    polyline scale(const polyline& poly, double scale);
    polygon scale(const polygon& poly, double scale);
    std::vector<polygon> scale(const std::vector<polygon>& polys, double scale);
    ch::ring simplify_rectilinear_ring(const ring& poly);
    rectangle bounding_rectangle(std::span<const ch::point> poly);
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
    bool is_point_in_rect(const point& pt, const rectangle& rect);
    bool is_point_in_polygon(const point& pt, const polygon& poly);
    point unit_vector(float theta);
    point normalize(const point& pt);
    point representative_point(const polygon& poly);
    bool is_degenerate_ring(const ch::ring& r);

    struct triangle {
        ch::point a;
        ch::point b;
        ch::point c;

        double area() const;
        point centroid() const;
    };

    std::vector<triangle> triangulate(const polygon& poly);

    using edge = std::tuple<point, point>;

    inline auto all_edges(const ring& r) {
        namespace rv = ranges::views;
        return rv::concat(r, rv::single(r.front())) | rv::sliding(2);
    }

    inline auto all_edges(const polygon& poly) {
        namespace r = ranges;
        namespace rv = ranges::views;

        return rv::concat(
            all_edges(poly.outer()),
            poly.inners() |
                rv::transform(
                    [](const auto& r) {
                        return all_edges(r);
                    }
                ) |
                rv::join
        ) | rv::transform(
            [](auto rng)->edge {
                return { rng[0], rng[1] };
            }
        );
    }

    inline auto all_vertices(const polygon& poly) {
        namespace r = ranges;
        namespace rv = ranges::views;

        return rv::concat(
            rv::all(poly.outer()),
            poly.inners() | 
                rv::transform(
                    [](const auto& r) {
                        return rv::all(r);
                    }
                ) | rv::join
        ); 
    }

    inline auto points_on_linesegment(const point& u, const point& v, float delta) {
        namespace r = ranges;
        namespace rv = ranges::views;

        auto length = euclidean_distance(u, v);
        auto n = static_cast<int>(length / delta);
        ch::point vec = delta * ((v - u) / length);
        return rv::iota(0, n) |
            rv::transform(
                [=](int i)->ch::point {
                    return u + i * vec;
                }
        );
    }

    void debug_geom();
}