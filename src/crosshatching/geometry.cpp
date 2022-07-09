#include "geometry.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <boost/functional/hash.hpp>
#include <range/v3/all.hpp>
#include <cmath>
#include <algorithm>

namespace r = ranges;
namespace rv = ranges::views;

namespace {

	bool are_parallel(ch::point p1, ch::point p2, ch::point p3) {
		if (p1.x == p2.x && p2.x == p3.x) {
			return true;
		}
		if (p1.y == p2.y && p2.y == p3.y) {
			return true;
		}
		return false;
	}

	ch::polyline elide_adjacent_parallel_edges(const ch::polyline& poly) {
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
			) |
					r::to_vector;
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

ch::matrix ch::translation_matrix(const cv::Point2d& pt)
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

ranges::any_view<ch::polyline> ch::transform(ranges::any_view<polyline> polys, const matrix& mat)
{
	return polys | rv::transform([=](const auto& poly) { return ch::transform(poly, mat); });
}

std::vector<ch::polyline> ch::transform(const std::vector<polyline>& polys, const matrix& mat)
{
	std::vector<ch::polyline> output(polys.size());
	std::transform(polys.begin(), polys.end(), output.begin(),
		[&mat](const auto& poly) {
			return ch::transform(poly, mat);
		}
	);
	return output;
}

ch::point ch::mean_point(const polyline& poly)
{
	auto x = 0.0;
	auto y = 0.0;
	for (const auto& pt : poly) {
		x += pt.x;
		y += pt.y;
	}
	return {
		x / poly.size(),
		y / poly.size()
	};
}

ch::polyline ch::transform(const polyline& poly, const matrix& mat)
{
	ch::polyline output(poly.size());
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

std::size_t ch::point_hasher::operator()(const cv::Point& p) const
{
	std::size_t seed = 0;
	boost::hash_combine(seed, p.x);
	boost::hash_combine(seed, p.y);

	return seed;
}

ch::polyline ch::scale(const polyline& poly, double scale)
{
	return poly |
		rv::transform([scale](const auto& pt) { return scale * pt; }) |
		r::to_vector;
}

ch::polyline ch::simplify_rectilinear_polygon(const ch::polyline& poly) {
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
	auto floats = poly | rv::transform([](const cv::Point2d& pt) {return cv::Point2f(pt); }) | r::to_vector;
	cv::Rect2d rect = cv::boundingRect(floats);
	return {
		static_cast<double>(rect.x),
		static_cast<double>(rect.y),
		static_cast<double>(rect.x + rect.width),
		static_cast<double>(rect.y + rect.height)
	};
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