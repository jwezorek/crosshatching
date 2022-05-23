#include "util.hpp"
#include "meanshift.hpp"
#include <opencv2/imgproc.hpp>
#include <opencv2/ximgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <sstream>
#include <random>
#include <iomanip>
#include <numbers>

namespace r = ranges;
namespace rv = ranges::views;

namespace {
	std::random_device rd;
	std::mt19937 random(rd());

	double ramp_up_right(double t, double k) {
		return (t <= k) ? 0.0 : (t - k) / (1.0 - k);
	}

	double sigmoidal_contrast(double u, double contrast, double thresh)
	{
		
	}

	uchar sigmoidal_contrast(uchar u, double contrast, double thresh) {
		auto sigmoidal_contrast = [contrast, thresh](double u) {
			return (1.0 / (1.0 + std::exp(contrast * (thresh - u))) - 1.0 / (1.0 + std::exp(contrast * thresh))) /
				(1.0 / (1.0 + std::exp(contrast * (thresh - 1.0))) - 1.0 / (1.0 + std::exp(contrast * thresh)));
		};
		return static_cast<uchar>(
			255.0 * sigmoidal_contrast( static_cast<double>(u) / 255.0 )
		);
	}

	std::vector<uchar> make_contrast_lookup_table(double beta, double sigma) {
		return rv::iota(0,256) |
			rv::transform(
				[beta, sigma](int u) {
					return sigmoidal_contrast(static_cast<uchar>(u), beta, sigma);
				}
			) | 
			r::to_vector;
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

	ch::polyline elide_adjacent_parallel_edges(const ch::polyline& poly) {
		auto first = poly.front();
		auto last = poly.back();
		return rv::concat(rv::concat(rv::single(last), poly), rv::single(first) ) |
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

int ch::uniform_rnd_int(int low, int high) {
	std::uniform_int_distribution<> ud(low, high);
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


cv::Mat ch::apply_contrast(cv::Mat img, double beta, double sigma) {
	using hsv_pixel = cv::Point3_<uchar>;
	auto lookup_table = make_contrast_lookup_table(beta, sigma);
	cv::Mat hsv;
	cv::cvtColor(img, hsv, cv::COLOR_BGR2HSV);

	hsv.forEach<hsv_pixel>(
		[&lookup_table](hsv_pixel& pix, const int* position) {
			pix.z = lookup_table[pix.z];
		}
	);

	cv::Mat output;
	cv::cvtColor(hsv, output, cv::COLOR_HSV2BGR);
	return output;
}

std::string ch::polyline_to_svg(const ch::polyline& poly, double thickness, bool closed) {
	std::stringstream ss;
	std::string cmd = (closed) ? "polygon" : "polyline";
	ss << "<" << cmd << " points = \"";
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

cv::Mat ch::convert_to_3channel_grayscale(cv::Mat img) {
	cv::Mat gray;
	cv::cvtColor(img, gray, cv::COLOR_BGR2GRAY);
	cv::Mat output;
	cv::cvtColor(gray, output, cv::COLOR_GRAY2BGR);
	return output;
}

/* ==============================================
*   Coherence-Enhancing Shock Filters
*  Author:WinCoder@qq.com
*  inspired by
*  Joachim Weickert "Coherence-Enhancing Shock Filters"
*  http://www.mia.uni-saarland.de/Publications/weickert-dagm03.pdf
*
*   Paras:
*   @img        : input image ranging value from 0 to 255.
*   @sigma      : sobel kernel size.
*   @str_sigma  : neighborhood size,see detail in reference[2]
*   @belnd      : blending coefficient.default value 0.5.
*   @iter       : number of iteration.
*
*   Example:
*   Mat dst = CoherenceFilter(I,11,11,0.5,4);
*   imshow("shock filter",dst);
*/

cv::Mat ch::coherence_filter(cv::Mat img, int sigma, int str_sigma, float blend, int iter)
{
	cv::Mat I = img.clone();
	int height = I.rows;
	int width = I.cols;

	for (int i = 0; i < iter; i++)
	{
		cv::Mat gray;
		cv::cvtColor(I, gray, cv::COLOR_BGR2GRAY);
		cv::Mat eigen;
		cv::cornerEigenValsAndVecs(gray, eigen, str_sigma, 3);

		std::vector<cv::Mat> vec;
		cv::split(eigen, vec);

		cv::Mat x, y;
		x = vec[2]; // c
		y = vec[3]; // s

		cv::Mat gxx, gxy, gyy;
		cv::Sobel(gray, gxx, CV_32F, 2, 0, sigma);
		cv::Sobel(gray, gxy, CV_32F, 1, 1, sigma);
		cv::Sobel(gray, gyy, CV_32F, 0, 2, sigma);

		cv::Mat ero;
		cv::Mat dil;
		cv::erode(I, ero, cv::Mat());
		cv::dilate(I, dil, cv::Mat()); 

		cv::Mat img1 = ero;
		for (int nY = 0; nY < height; nY++)
		{
			for (int nX = 0; nX < width; nX++)
			{
				if (x.at<float>(nY, nX) * x.at<float>(nY, nX) * gxx.at<float>(nY, nX)
					+ 2 * x.at<float>(nY, nX) * y.at<float>(nY, nX) * gxy.at<float>(nY, nX)
					+ y.at<float>(nY, nX) * y.at<float>(nY, nX) * gyy.at<float>(nY, nX) < 0)
				{
					img1.at<cv::Vec3b>(nY, nX) = dil.at<cv::Vec3b>(nY, nX); 
				}
			}
		}
		I = I * (1.0 - blend) + img1 * blend;
	}
	return I;
}

cv::Mat ch::anisotropic_diffusion(cv::Mat img, double alpha, double k, int iters) {
	cv::Mat output;
	cv::ximgproc::anisotropicDiffusion(img, output, alpha, k, iters);
	return output;
}

std::tuple<cv::Mat, cv::Mat> ch::meanshift_segmentation(const cv::Mat& input, int sigmaS, float sigmaR, int minSize) {
	cv::Mat output;
	cv::Mat labels;
	auto mss = createMeanShiftSegmentation(sigmaS, sigmaR, minSize, 4, true);
	mss->processImage(input, output, labels);
	return { output, labels };
}

cv::Mat ch::convert_to_1channel_gray(const cv::Mat& color) {
	cv::Mat gray;
	cv::cvtColor(color, gray, cv::COLOR_BGR2GRAY);
	return gray;
}

double ch::degrees_to_radians(double degrees)
{
	return (std::numbers::pi * degrees) / 180.0;
}

ch::polyline ch::simplify_rectilinear_polygon(const ch::polyline& poly) {
	return elide_adjacent_parallel_edges(poly);
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

	// Cohen–Sutherland clipping algorithm clips a line from
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

void ch::write_label_map_visualization(cv::Mat img, const std::string& output_file) {

	using pixel = cv::Point3_<uint8_t>;
	static auto random_color = []()->pixel {
		return pixel(ch::uniform_rnd_int(0, 255), ch::uniform_rnd_int(0, 255), ch::uniform_rnd_int(0, 255));
	};

	cv::Mat mat(img.rows, img.cols, CV_8UC3);
	std::unordered_map<int, pixel> tbl;
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			int label = img.at<int>(i, j);
			pixel color;
			auto iter = tbl.find(label);
			if (iter != tbl.end()) {
				color = iter->second;
			} else {
				color = random_color();
				tbl[label] = color;
			}
			mat.at<pixel>(i, j) = color;
		}
	}
	cv::imwrite(output_file, mat);
}

std::vector<uchar> ch::unique_gray_values(const cv::Mat& input) {
	if (input.channels() != 1) {
		throw std::runtime_error("called get_gray_levels on color image");
	}
	std::array<bool, 256> grays = {};
	input.forEach<uchar>([&grays](uchar gray, const int* pos) { grays[gray] = true;  });
	return rv::iota(0) |
		rv::take(256) |
		rv::filter([&grays](int g) {return grays[g]; }) |
		r::to<std::vector<uchar>>() |
		r::action::reverse;
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

int ch::max_val_in_mat(cv::Mat mat) {
	double min_val;
	double max_val;
	cv::Point minLoc;
	cv::Point maxLoc;

	cv::minMaxLoc(mat, &min_val, &max_val, &minLoc, &maxLoc);
	return static_cast<int>(max_val);
}