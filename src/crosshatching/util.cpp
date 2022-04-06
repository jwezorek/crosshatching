#include "util.hpp"
#include "meanshift.hpp"
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

cv::Mat ch::do_segmentation(const cv::Mat& input, int sigmaS, float sigmaR, int minSize) {
	cv::Mat output;
	cv::Mat labels;
	auto mss = createMeanShiftSegmentation(sigmaS, sigmaR, minSize, true);
	mss->processImage(input, output, labels);
	return output;
}

cv::Mat ch::convert_to_gray(const cv::Mat& color) {
	cv::Mat gray;
	cv::cvtColor(color, gray, cv::COLOR_BGR2GRAY);
	return gray;
}