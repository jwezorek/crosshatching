#include "geometry.hpp"
#include "segmentation.hpp"
#include "util.hpp"
#include "point_set.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <range/v3/all.hpp>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <array>
#include <numeric>
#include <qdebug.h>
#include "drawing.hpp"

namespace r = ranges;
namespace rv = ranges::views;
namespace bg = boost::geometry;

namespace {

	using vec = Eigen::Matrix<double, 3, 1>;

	template<typename T>
	T scale_geometry(const T& poly, double scale)
	{
		return poly |
			rv::transform([scale](const auto& pt) { return scale * pt; }) |
			r::to_<T>();
	}

	template<typename T>
	ch::rectangle bounding_rect_of_geometry(const T& poly) {
		auto floats = poly | 
			rv::transform([](const auto& pt) {return cv::Point2f(pt); }) |
			r::to_vector;
		cv::Rect2d rect = cv::boundingRect(floats);
		return {
			static_cast<double>(rect.x),
			static_cast<double>(rect.y),
			static_cast<double>(rect.x + rect.width),
			static_cast<double>(rect.y + rect.height)
		};
	}

	bool are_parallel(ch::point p1, ch::point p2, ch::point p3) {
		if (p1.x == p2.x && p2.x == p3.x) {
			return true;
		}
		if (p1.y == p2.y && p2.y == p3.y) {
			return true;
		}
		return false;
	}

	ch::ring elide_adjacent_parallel_edges(const ch::ring& poly) {
		auto first = poly.front();
		auto last = poly.back();
		return rv::concat(rv::concat(rv::single(last), poly), rv::single(first)) |
			rv::sliding(3) |
			rv::remove_if(
				[](auto rng3)->bool {
					return are_parallel(rng3[0], rng3[1], rng3[2]);
				}
			) |
			rv::transform(
				[](auto rng3)->ch::point {
					return rng3[1];
				}
			) | r::to_<ch::ring>();
	}

	namespace cohen_sutherland {

		typedef int out_code;

		const int INSIDE = 0; // 0000
		const int LEFT = 1;   // 0001
		const int RIGHT = 2;  // 0010
		const int BOTTOM = 4; // 0100
		const int TOP = 8;    // 1000

		// Compute the bit code for a point (x, y) using the clip
		// bounded diagonally by (xmin, ymin), and (xmax, ymax)

		// ASSUME THAT xmax, xmin, ymax and ymin are global constants.

		out_code compute_out_code(double x, double y, double xmin, double ymin, double xmax, double ymax)
		{
			out_code code;

			code = INSIDE;          // initialised as being inside of [[clip window]]

			if (x < xmin)           // to the left of clip window
				code |= LEFT;
			else if (x > xmax)      // to the right of clip window
				code |= RIGHT;
			if (y < ymin)           // below the clip window
				code |= BOTTOM;
			else if (y > ymax)      // above the clip window
				code |= TOP;

			return code;
		}

		// Cohen-Sutherland clipping algorithm clips a line from
		// P0 = (x0, y0) to P1 = (x1, y1) against a rectangle with 
		// diagonal from (xmin, ymin) to (xmax, ymax).
		std::optional<ch::line_segment>  clip(double x0, double y0, double x1, double y1,
			double xmin, double ymin, double xmax, double ymax)
		{
			// compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
			out_code outcode0 = compute_out_code(x0, y0, xmin, ymin, xmax, ymax);
			out_code outcode1 = compute_out_code(x1, y1, xmin, ymin, xmax, ymax);
			bool accept = false;

			while (true) {
				if (!(outcode0 | outcode1)) {
					// bitwise OR is 0: both points inside window; trivially accept and exit loop
					using p = ch::point;
					return ch::line_segment(p{ x0, y0 }, p{ x1, y1 });
				} else if (outcode0 & outcode1) {
					// bitwise AND is not 0: both points share an outside zone (LEFT, RIGHT, TOP,
					// or BOTTOM), so both must be outside window; exit loop (accept is false)
					return {};
				} else {
					// failed both tests, so calculate the line segment to clip
					// from an outside point to an intersection with clip edge
					double x, y;

					// At least one endpoint is outside the clip rectangle; pick it.
					out_code outcodeOut = outcode1 > outcode0 ? outcode1 : outcode0;

					// Now find the intersection point;
					// use formulas:
					//   slope = (y1 - y0) / (x1 - x0)
					//   x = x0 + (1 / slope) * (ym - y0), where ym is ymin or ymax
					//   y = y0 + slope * (xm - x0), where xm is xmin or xmax
					// No need to worry about divide-by-zero because, in each case, the
					// outcode bit being tested guarantees the denominator is non-zero
					if (outcodeOut & TOP) {           // point is above the clip window
						x = x0 + (x1 - x0) * (ymax - y0) / (y1 - y0);
						y = ymax;
					} else if (outcodeOut & BOTTOM) { // point is below the clip window
						x = x0 + (x1 - x0) * (ymin - y0) / (y1 - y0);
						y = ymin;
					} else if (outcodeOut & RIGHT) {  // point is to the right of clip window
						y = y0 + (y1 - y0) * (xmax - x0) / (x1 - x0);
						x = xmax;
					} else if (outcodeOut & LEFT) {   // point is to the left of clip window
						y = y0 + (y1 - y0) * (xmin - x0) / (x1 - x0);
						x = xmin;
					}

					// Now we move outside point to intersection point to clip
					// and get ready for next pass.
					if (outcodeOut == outcode0) {
						x0 = x;
						y0 = y;
						outcode0 = compute_out_code(x0, y0, xmin, ymin, xmax, ymax);
					} else {
						x1 = x;
						y1 = y;
						outcode1 = compute_out_code(x1, y1, xmin, ymin, xmax, ymax);
					}
				}
			}
		}
	}

	int update_vertex_usage_count(ch::point_map<int>& vertex_count, const ch::point& pt,
			const ch::dimensions<int>& dims) {
		bool not_in_table = vertex_count.find(pt) == vertex_count.end();
		if (not_in_table && (pt.x == 0 || pt.x == dims.wd)) {
			++vertex_count[pt];
		}
		if (not_in_table && (pt.y == 0 || pt.y == dims.hgt)) {
			++vertex_count[pt];
		}
		return ++vertex_count[pt];
	}

	r::any_view<ch::point> update_vertex_usage_counts(const ch::ring& ring,
			ch::point_map<int>& vertex_count, const ch::dimensions<int>& dims) {
		return ring | rv::remove_if(
				[&](const ch::point& pt)->bool {
					int count = update_vertex_usage_count(vertex_count, pt, dims);
					return count != 3;
				}
			);
	}

	r::any_view<ch::point> update_vertex_usage_counts(const ch::polygon& poly,
			ch::point_map<int>& vertex_count, const ch::dimensions<int>& dims) {
		return rv::concat(
			rv::single(update_vertex_usage_counts(poly.outer(), vertex_count, dims)),
			poly.inners() |
			rv::transform(
				[&](const ch::ring& r) {
					return update_vertex_usage_counts(r, vertex_count, dims);
				}
			)
		) | rv::join;
	}
	
	bool is_critical_point_free(const ch::point_set& set, std::span<const ch::point> pts) {
		for (ch::point pt : pts) {
			if (set.find(pt) != set.end()) {
				return false;
			}
		}
		return true;
	}

	void insert_island_critical_points(ch::point_set& set, const ch::ring& r) {
		if (is_critical_point_free(set, r)) {
			set.insert(r.front());
		}
	}

	void insert_island_critical_points(
			ch::point_set& set, std::span<const ch::polygon> polygons) {
		for (const auto& poly : polygons) {
			insert_island_critical_points(set, poly.outer());
			for (const auto& hole : poly.inners()) {
				insert_island_critical_points(set, hole);
			}
		}
	}
	
	ch::point_set generate_critical_points_set( 
			std::span<const ch::polygon> polygons, const ch::dimensions<double>& rect) 
	{
		ch::dimensions<int> dims(rect);
		ch::point_map<int> vertex_usage_counts;
		vertex_usage_counts.reserve(polygons.size() * 10);

		auto set = polygons |
			rv::transform(
				[&vertex_usage_counts, &dims](const auto& poly) {
					return update_vertex_usage_counts(poly, vertex_usage_counts, dims);
				}
			) | 
			rv::join |
			r::to_<ch::point_set>();

		insert_island_critical_points(set, polygons);

		return set;
	}

}

ch::polyline ch::make_polyline(size_t sz)
{
	polyline poly;
	poly.resize(sz);
	return poly;
}

ch::ring ch::make_ring(size_t sz) {
	ring r;
	r.resize(sz);
	return r;
}

ch::polygon ch::make_polygon(const ch::ring& outer, const std::vector<ch::ring>& inners) {
	ch::polygon poly;
	poly.outer() = outer;
	poly.inners() = inners;
	return poly;
}

ch::matrix ch::rotation_matrix(double theta)
{
	return rotation_matrix(std::cos(theta), std::sin(theta));
}

ch::matrix ch::rotation_matrix(double cos_theta, double sin_theta) {
	ch::matrix rotation;
	rotation <<
		cos_theta, -sin_theta, 0,
		sin_theta, cos_theta, 0,
		0, 0, 1;
	return rotation;
}

ch::matrix ch::translation_matrix(double x, double y) {
	ch::matrix translation;
	translation <<
		1, 0, x,
		0, 1, y,
		0, 0, 1;
	return translation;
}

ch::matrix ch::translation_matrix(const point& pt)
{
	return translation_matrix(pt.x, pt.y);
}

ch::matrix ch::scale_matrix(double x_scale, double y_scale) {
	ch::matrix scale;
	scale <<
		x_scale, 0, 0,
		0, y_scale, 0,
		0, 0, 1;
	return scale;
}

ranges::any_view<ch::polyline> ch::transform_view(ranges::any_view<polyline> polys, const matrix& mat)
{
	return polys | rv::transform([=](const auto& poly) { return ch::transform(poly, mat); });
}

std::vector<ch::polyline> ch::transform( std::span<const polyline> polys, const matrix& mat)
{
	std::vector<ch::polyline> output(polys.size());
	std::transform(polys.begin(), polys.end(), output.begin(),
		[&mat](const auto& poly) {
			return ch::transform(poly, mat);
		}
	);
	return output;
}

ch::point ch::mean_point(std::span<const point> points)
{
	auto x = 0.0;
	auto y = 0.0;
	for (const auto& pt : points) {
		x += pt.x;
		y += pt.y;
	}
	return {
		x / points.size(),
		y / points.size()
	};
}

ch::polyline ch::transform(std::span<const point> poly, const matrix& mat)
{
	ch::polyline output = make_polyline(poly.size());
	std::transform(poly.begin(), poly.end(), output.begin(),
		[&mat](const auto& p) {return transform(p, mat); }
	);
	return output;
}

ch::point ch::transform(const point& pt, const matrix& mat)
{
	vec v;
	v << pt.x, pt.y, 1.0;
	v = mat * v;
	return { v[0], v[1] }; 
}

ch::ring ch::transform(const ring& r, const matrix& mat) {
	ch::ring output = make_ring(r.size());
	std::transform(r.begin(), r.end(), output.begin(),
		[&mat](const auto& p) {return transform(p, mat); }
	);
	return output;
}

ch::polygon ch::transform(const ch::polygon& p, const ch::matrix& mat) {
	return make_polygon(
		transform(p.outer(), mat),
		p.inners() | rv::transform([&mat](const auto& hole) {return ch::transform(hole, mat); }) | r::to_vector
	);
}

void ch::paint_polyline(cv::Mat& mat, const polyline& poly, double thickness, int color, point offset)
{
	std::vector<cv::Point> int_pts(poly.size());
	std::transform(poly.begin(), poly.end(), int_pts.begin(),
		[offset](const auto& p) {
			return cv::Point(
				static_cast<int>(std::round(p.x + offset.x)),
				static_cast<int>(std::round(p.y + offset.y))
			); 
		}
	);
	auto npts = int_pts.size();
	cv::polylines(mat, int_pts, false, color, thickness, 8, 0);
}

void ch::paint_polyline_aa(cv::Mat& mat, const polyline& poly, double thickness, int color, point offset)
{
	std::vector<cv::Point> int_pts(poly.size());
	std::transform(poly.begin(), poly.end(), int_pts.begin(),
		[offset](const auto& p) {
			return cv::Point(
				static_cast<int>(std::round(p.x + offset.x)),
				static_cast<int>(std::round(p.y + offset.y))
			);
		}
	);
	auto npts = int_pts.size();
	cv::polylines(mat, int_pts, false, color, thickness, cv::LINE_AA, 0);
}

double ch::euclidean_distance(const point& pt1, const point& pt2)
{
	auto x_diff = pt2.x - pt1.x;
	auto y_diff = pt2.y - pt1.y;
	return std::sqrt(x_diff * x_diff + y_diff * y_diff);
}

ch::polylines ch::clip_lines_to_poly(const ch::polylines& strokes, const ch::polygon& poly) {
	polylines result;
	bg::intersection(poly, strokes, result);
	return result;
}

ch::ring ch::scale(const ch::ring& r, double scale)
{
	return scale_geometry<ring>(r, scale);
}

ch::polyline ch::scale(const polyline& poly, double scale)
{
	return scale_geometry<polyline>(poly, scale);
}

ch::polygon ch::scale(const polygon& poly, double val) {
	return make_polygon(
		scale(poly.outer(), val),
		poly.inners() |
		rv::transform([val](const auto& r) {return scale(r, val); }) |
		r::to_vector
	);
}

ch::ring ch::simplify_rectilinear_ring(const ch::ring& poly) {
	return elide_adjacent_parallel_edges(poly);
}

std::optional<ch::line_segment> ch::linesegment_rectangle_intersection(const ch::line_segment& line_seg, const ch::rectangle& rect) {
	auto [p1, p2] = line_seg;
	auto [x1, y1, x2, y2] = rect;
	return cohen_sutherland::clip(
		p1.x, p1.y, p2.x, p2.y,
		x1, y1, x2, y2
	);
}

ch::rectangle ch::bounding_rectangle(const ch::polyline& poly) {
	return bounding_rect_of_geometry(poly);
}

ch::rectangle ch::bounding_rectangle(const ch::ring& poly) {
	return bounding_rect_of_geometry(poly);
}

ch::rectangle ch::bounding_rectangle(const ch::polygon& poly) {
	return bounding_rectangle(poly.outer());
}

ch::rectangle ch::bounding_rectangle(std::span<const polygon> polys) {
	return  std::accumulate(polys.begin(), polys.end(), bounding_rectangle(polys.front()),
		[](const rectangle& r, const polygon& poly)->rectangle {
			return union_rect(bounding_rectangle(poly), r);
		}
	);
}

ch::rectangle ch::union_rect(const ch::rectangle & r1, const ch::rectangle& r2) {

	std::array<double, 4> x_vals = { std::get<0>(r1), std::get<2>(r1),
		std::get<0>(r2), std::get<2>(r2) };
	std::array<double, 4> y_vals = { std::get<1>(r1), std::get<3>(r1),
		std::get<1>(r2), std::get<3>(r2) };

	auto x_pair = std::minmax_element(x_vals.begin(), x_vals.end());
	auto y_pair = std::minmax_element(y_vals.begin(), y_vals.end());

	return { *x_pair.first, *y_pair.first, *x_pair.second, *y_pair.second};
}

cv::Rect ch::union_rect_and_pt(const cv::Rect& r, cv::Point2i pt) {
	int x1 = r.x;
	int y1 = r.y;
	int x2 = r.x + r.width - 1;
	int y2 = r.y + r.height - 1;

	x1 = std::min(x1, pt.x);
	y1 = std::min(y1, pt.y);
	x2 = std::max(x2, pt.x);
	y2 = std::max(y2, pt.y);

	return {
		x1,
		y1,
		x2 - x1 + 1,
		y2 - y1 + 1
	};
}

using edge = std::tuple<ch::ring::const_iterator, ch::ring::const_iterator>;
/*
std::vector<edge> get_shared_edges(
		const ch::ring& verts, const ch::path_table& tbl) {

	auto iter_list = rv::iota(0, static_cast<int>(verts.size())) |
		rv::transform(
			[&verts](int index)->ch::ring::const_iterator {
				return verts.cbegin() + index;
			}
		) | r::to_vector;

	auto critical_pts = iter_list | rv::remove_if(
			[&tbl](const auto& i) {
				return tbl.find(*i) == tbl.end();
			}
		) | r::to_vector;

	return critical_pts |
		rv::cycle |
		rv::sliding(2) |
		rv::take(critical_pts.size()) |
		rv::transform(
			[](auto pair)->edge {
				return { pair[0] , pair[1] };
			}
		) | r::to_vector;
}

std::vector<std::vector<ch::point>> get_paths( const ch::ring& ring, 
		const path_table& tbl) {

	return {};
}

ch::ring simplify_ring(const ch::ring& ring, path_table& paths, double param) {
	return {};
}

ch::polygon simplify_polygon(const ch::polygon& poly, path_table& paths, double param) {
	return {};
}

std::vector<ch::polygon> ch::simplify_rectangle_dissection(const std::vector<ch::polygon>& dissection,
		const ch::dimensions<double>& rect, double param) {

	auto path_tbl = generate_critical_points_table(dissection, rect);

	return dissection;
}

void ch::debug_geometry() {
	ch::ring r = { {34,15},{7,7},{16,1},{18,12}, {100,100}, {3,3}, {14,14} , {34,21} };
	path_table tbl;
	tbl[{7, 7}] = {};
	tbl[{100, 100}] = {};
	tbl[{14, 14}] = {};
	auto foo = get_shared_edges(r, tbl);
	for (auto [p1, p2] : foo) {
		qDebug() << "(" << p1->x << "," << p1->y << ") (" << p2->x << "," << p2->y << ")";
	}
}
*/

std::vector<ch::polygon> ch::simplify_rectangle_dissection(const std::vector<ch::polygon>& dissection,
		const ch::dimensions<double>& rect, double param) {

	return {};
}

void ch::debug_geometry() {
}


void ch::debug_geom(cv::Mat mat) {

	auto blobs = ch::to_blobs<ch::color>(mat);
	auto polygons = blobs |
		rv::transform([](const auto& tup) {return std::get<1>(tup); }) |
		r::to_vector;
	auto dims = ch::mat_dimensions(mat);
	auto critical_points = generate_critical_points_set(polygons, dims);
	cv::Mat debug_img(dims.hgt+1, dims.wd+1, CV_8UC3, cv::Scalar(255, 255, 255));
	mat.copyTo(debug_img(cv::Rect(0, 0, dims.wd, dims.hgt)));
	for (const auto& pt : critical_points | 
		rv::transform([](const auto& p) {return p.to_point(); })) {
		debug_img.at<cv::Vec3b>(int_point(pt)) = cv::Vec3b(0, 0, 255);
	}
	cv::imwrite("C:\\test\\test_blobs_cp.png", debug_img);
}
