#include "geometry.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <cmath>
#include <algorithm>

namespace r = ranges;
namespace rv = ranges::views;

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

void ch::paint_polyline(cv::Mat& mat, const polyline& poly, int thickness, int color, point offset)
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

double ch::euclidean_distance(const point& pt1, const point& pt2)
{
	auto x_diff = pt2.x - pt1.x;
	auto y_diff = pt2.y - pt1.y;
	return std::sqrt(x_diff * x_diff + y_diff * y_diff);
}

