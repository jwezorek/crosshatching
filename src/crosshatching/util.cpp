#include "util.hpp"
#include <opencv2/imgproc.hpp>
#include <sstream>
#include <random>
#include <iomanip>

namespace r = ranges;
namespace rv = ranges::views;

namespace {
	std::random_device rd;
	std::mt19937 random(rd());

	double ramp_up_right(double t, double k) {
		return (t <= k) ? 0.0 : (t - k) / (1.0 - k);
	}
}

std::string ch::svg_header(int wd, int hgt, bool bkgd_rect)
{
    std::stringstream ss;

    ss << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
    ss << "<svg width=\"" + std::to_string(wd) + "px\" height=\"" + std::to_string(hgt) + "px\"  xmlns = \"http://www.w3.org/2000/svg\" version = \"1.1\">\n";
	if (bkgd_rect) {
		ss << "<rect width=\"" + std::to_string(wd) + "\" height=\"" + std::to_string(hgt) + "\"  fill=\"white\"/>\n";
	}
    return ss.str();
}

std::string uchar_to_hex(unsigned char uc) {
	std::stringstream ss;
	ss << std::setw(2) << std::setfill('0') << std::hex << static_cast<int>(uc);
	return ss.str();
}

std::string ch::gray_to_svg_color(unsigned char gray)
{
	auto hex_val = uchar_to_hex(gray);
	std::stringstream ss;
	ss << "#" << hex_val << hex_val << hex_val;
	return ss.str();
}

ch::polyline ch::scale(const polyline& poly, double scale)
{
	return poly |
		rv::transform([scale](const auto& pt) { return scale * pt; }) |
		r::to_vector;
}

double ch::normal_rnd(double mean, double stddev)
{
	std::normal_distribution<double> nd(mean, stddev);
	return nd(random);
}

double ch::uniform_rnd(double lower_bound, double upper_bound)
{
	std::uniform_real_distribution<> ud(lower_bound, upper_bound);
	return ud(random);
}

ch::rnd_fn ch::normal_rnd_fn(double mean, double stddev)
{
	return [mean, stddev]() {return normal_rnd(mean, stddev); };
}

ch::rnd_fn ch::const_rnd_fn(double val)
{
	return [val]() {return val; };
}

double ch::ramp(double t, double k, bool right, bool up) {
	if (right) {
		return up ? ramp_up_right(t, k) : 1.0 - ramp_up_right(t, k);
	} else {
		t = 1.0 - t;
		k = 1.0 - k;
		return up ? 1.0 - ramp_up_right(t, k) : ramp_up_right(t, k);
	}
}

double ch::sigmoidal_contrast(double u, double contrast, double thresh)
{
	return (1.0 / (1.0 + std::exp(contrast * (thresh - u))) - 1.0 / (1.0 + std::exp(contrast * thresh))) /
		(1.0 / (1.0 + std::exp(contrast * (thresh - 1.0))) - 1.0 / (1.0 + std::exp(contrast * thresh)));
}

uchar ch::sigmoidal_contrast(uchar u, double contrast, double thresh) {
	return static_cast<uchar>(
		255.0 * sigmoidal_contrast(
			static_cast<double>(u) / 255.0, contrast, thresh
		) 
	);
}

cv::Mat ch::apply_contrast(cv::Mat img, double beta, double sigma) {
	return img;
}

std::string ch::polyline_to_svg(const ch::polyline& poly, double thickness) {
	std::stringstream ss;
	ss << "<polyline points=\"";
	for (const auto& pt : poly) {

		ss << " " << pt.x << "," << pt.y;
	}
	ss << "\" style=\"fill:none;stroke:black;stroke-width:" << thickness << "px\" />";
	return ss.str();
}

std::string ch::to_string(double val, int precision) {
	std::stringstream ss;
	ss << std::fixed << std::setprecision(2) << val;
	return ss.str();
}

cv::Mat ch::scale(cv::Mat mat, double scale) {
	cv::Mat scaled;
	int scaled_wd = static_cast<int>(scale * mat.cols);
	int scaled_hgt = static_cast<int>(scale * mat.rows);
	cv::resize(mat, scaled, cv::Size(scaled_wd, scaled_hgt), 0, 0, cv::INTER_CUBIC);
	return scaled;
}