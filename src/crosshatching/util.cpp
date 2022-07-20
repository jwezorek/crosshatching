#include "util.hpp"
#include "meanshift.hpp"
#include <opencv2/imgproc.hpp>
#include <opencv2/ximgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <sstream>
#include <random>
#include <iomanip>
#include <numbers>
#include <fstream>
#include <unordered_set>
#include <qdebug.h>

namespace r = ranges;
namespace rv = ranges::views;

namespace {
	std::random_device rd;
	std::mt19937 random(rd());

	struct color_hasher {
		size_t operator()(const ch::color& c) const
		{
			std::size_t seed = 0;
			boost::hash_combine(seed, c[0]);
			boost::hash_combine(seed, c[1]);
			boost::hash_combine(seed, c[2]);

			return seed;
		}
	};

	using color_set = std::unordered_set<ch::color, color_hasher>;

	double ramp_up_right(double t, double k) {
		return (t <= k) ? 0.0 : (t - k) / (1.0 - k);
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

	std::string loop_to_path_commands(const ch::ring& poly, double scale) {
		std::stringstream ss;
		ss << "M " << scale * poly[0].x << "," << scale * poly[0].y << " L";
		for (const auto& pt : rv::tail(poly)) {
			ss << " " << scale * pt.x << "," << scale * pt.y;
		}
		ss << " Z";
		return ss.str();
	}

	std::string svg_path_commands(const ch::polygon& poly, double scale) {
		std::stringstream ss;
		ss << loop_to_path_commands(poly.outer(), scale);
		for (const auto& hole : poly.inners()) {
			ss << " " << loop_to_path_commands(hole, scale);
		}
		return ss.str();
	}

	template<typename T>
	std::string polygon_to_svg(const ch::polygon& poly, T color, double scale) {
		std::stringstream ss;
		ss << "<path fill-rule=\"evenodd\" stroke=\"none\" fill=\"";
		ss << ch::to_svg_color(color) << "\" d=\"";
		ss << svg_path_commands(poly, scale);
		ss << "\" />";
		return ss.str();
	}

	template<typename T>
	void polygons_to_svg(const std::string& output_file,
			std::span<std::tuple<T, ch::polygon>> polys,
			double scale) {

		auto [x1, y1, wd, hgt] = ch::bounding_rectangle(
			polys |
			rv::transform([](const auto& tup) {return std::get<1>(tup); }) |
			r::to_vector
		);

		std::ofstream outfile(output_file);
		outfile << ch::svg_header(static_cast<int>(wd), static_cast<int>(hgt));
		for (const auto& [color, poly] : polys) {
			outfile << polygon_to_svg(poly, color, scale) << std::endl;
		}
		outfile << "</svg>" << std::endl;
		outfile.close();
	}
}

std::string uchar_to_hex(unsigned char uc) {
	std::stringstream ss;
	ss << std::setw(2) << std::setfill('0') << std::hex << static_cast<int>(uc);
	return ss.str();
}

std::string ch::to_svg_color(uchar gray)
{
	auto hex_val = uchar_to_hex(gray);
	std::stringstream ss;
	ss << "#" << hex_val << hex_val << hex_val;
	return ss.str();
}

std::string ch::to_svg_color(const color& c) {
	std::stringstream ss;
	ss << "#" << uchar_to_hex(c[2]) << uchar_to_hex(c[1]) << uchar_to_hex(c[0]);
	return ss.str();
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

void ch::detail::polygons_to_svg_aux(const std::string& output_file,
		std::span<std::tuple<color, polygon>> polys,
		double scale) {
	::polygons_to_svg<color>(output_file, polys, scale);
}

void ch::detail::polygons_to_svg_aux(const std::string& output_file,
		std::span<std::tuple<uchar, polygon>> polys,
		double scale) {
	::polygons_to_svg<uchar>(output_file, polys, scale);
}

std::string ch::svg_header(int wd, int hgt, bool bkgd_rect)
{
	std::stringstream ss;

	ss << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
	ss << "<svg width=\"" + std::to_string(wd) + "px\" height=\"" + 
		std::to_string(hgt) + "px\"  xmlns = \"http://www.w3.org/2000/svg\" version = \"1.1\">\n";
	if (bkgd_rect) {
		ss << "<rect width=\"" + std::to_string(wd) + "\" height=\"" + 
			std::to_string(hgt) + "\"  fill=\"white\"/>\n";
	}
	return ss.str();
}

std::string ch::polyline_to_svg(std::span<const point> poly, double thickness, bool closed) {
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
	if (img.channels() == 3) {
		cv::cvtColor(img, gray, cv::COLOR_BGR2GRAY);
	} else {
		gray = img;
	}
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

std::tuple<cv::Mat, cv::Mat> ch::meanshift_segmentation(const cv::Mat& input, int sigmaS, 
		float sigmaR, int minSize) {
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

void ch::label_map_to_visualization_img(cv::Mat img, const std::string& output_file) {

	using pixel = cv::Point3_<uint8_t>;
	static auto random_color = []()->pixel {
		return pixel(
			ch::uniform_rnd_int(0, 255), 
			ch::uniform_rnd_int(0, 255), 
			ch::uniform_rnd_int(0, 255)
		);
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



ch::dimensions<int> ch::mat_dimensions(cv::Mat mat) {
	return { mat.cols, mat.rows };
}

std::vector<uchar> ch::detail::unique_1channel_values(const cv::Mat& input) {
	if (input.channels() != 1) {
		throw std::runtime_error("called unique_1channel_values on color image");
	}
	std::array<bool, 256> grays = {};
	input.forEach<uchar>([&grays](uchar gray, const int* pos) { grays[gray] = true;  });
	return rv::iota(0) |
		rv::take(256) |
		rv::filter([&grays](int g) {return grays[g]; }) |
		r::to<std::vector<uchar>>() |
		r::action::reverse;
}

std::vector<ch::color> ch::detail::unique_3channel_values(const cv::Mat& input) {
	if (input.channels() != 3) {
		throw std::runtime_error("called unique_3channel_values on color image");
	}
	color_set colors;
	for (int y = 0; y < input.rows; ++y) {
		for (int x = 0; x < input.cols; ++x) {
			auto pixel = input.at<ch::color>(y, x);
			colors.insert(pixel);
		}
	}
	return colors | r::to_vector;
}

int ch::max_val_in_mat(cv::Mat mat) {
	double min_val;
	double max_val;
	cv::Point minLoc;
	cv::Point maxLoc;

	cv::minMaxLoc(mat, &min_val, &max_val, &minLoc, &maxLoc);
	return static_cast<int>(max_val);
}


