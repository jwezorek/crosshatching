#include "util.hpp"
#include "meanshift.hpp"
#include "random/counter_based_engine.hpp"
#include "random/philox_prf.hpp"
#include "perlin_noise.hpp"
#include "qpainterpath.h"
#include "qvector.h"
#include <opencv2/imgproc.hpp>
#include <opencv2/ximgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/photo.hpp>
#include <sstream>
#include <random>
#include <iomanip>
#include <numbers>
#include <fstream>
#include <unordered_set>
#include <qdebug.h>

/*------------------------------------------------------------------------------------------------*/

namespace r = ranges;
namespace rv = ranges::views;

namespace {
    std::random_device rd;
    std::mt19937 random(rd());

    uint32_t random_uint32() {
        std::uniform_int_distribution<uint32_t> distr;
        return distr(random);
    }

    double ramp_up_right(double t, double k) {
        return (t <= k) ? 0.0 : (t - k) / (1.0 - k);
    }

    uchar sigmoidal_contrast(uchar u, double contrast, double thresh) {
        using namespace std;
        auto sigmoidal_contrast = [k = contrast, mu = thresh](double u) {
            return (1.0 / (1.0 + exp(k * (mu - u))) - 1.0 / (1.0 + exp(k * mu))) /
                (1.0 / (1.0 + exp(k * (mu - 1.0))) - 1.0 / (1.0 + exp(k * mu)));
        };
        return static_cast<uchar>(
            255.0 * sigmoidal_contrast(static_cast<double>(u) / 255.0)
            );
    }

    std::vector<uchar> make_contrast_lookup_table(double beta, double sigma,
        double white_cutoff, double black_cutoff) {
        uchar white = static_cast<uchar>((1.0 - white_cutoff) * 255.0);
        uchar black = static_cast<uchar>((1.0 - black_cutoff) * 255.0);
        return rv::iota(0, 256) |
            rv::transform(
                [beta, sigma, white, black](int u)->uchar {
                    auto val = sigmoidal_contrast(static_cast<uchar>(u), beta, sigma);
                    if (val <= black) {
                        return 0;
                    } else if (val >= white) {
                        return 255;
                    }
                    return val;
                }
            ) |
            r::to_vector;
    }

    std::string loop_to_path_commands(std::span<const ch::point> poly, double scale) {
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
	void polygons_to_svg(const std::string& output_file,
			std::span<std::tuple<T, ch::polygon>> polys,
			double scale) {

		auto [x1, y1, wd, hgt] = ch::bounding_rectangle(
			polys |
			rv::transform([](const auto& tup) {return std::get<1>(tup); }) |
			r::to_vector
		);

		std::ofstream outfile(output_file);
		outfile << ch::svg_header(static_cast<int>(scale*wd), static_cast<int>(scale * hgt));
		for (const auto& [color, poly] : polys) {
			outfile << ch::polygon_to_svg(poly, ch::to_svg_color(color), scale) << std::endl;
		}
		outfile << "</svg>" << std::endl;
		outfile.close();
	}

	QPolygonF ring_to_qpolygon(const ch::ring& r) {
		auto qpoly = QPolygonF(
			r |
			rv::transform(
				[](const ch::point pt)->QPointF {
					return {
						static_cast<float>(pt.x),
						static_cast<float>(pt.y)
					};
				}
			) | r::to_<QVector>
		);
        if (!qpoly.isClosed()) {
            qpoly << qpoly.first();
        }
        return qpoly;
	}

	QColor cv_color_to_qcolor(ch::color c) {
		auto red = c[2];
		auto green = c[1];
		auto blue = c[0];
		return QColor(red, green, blue);
	}

    std::string uchar_to_hex(unsigned char uc) {
        std::stringstream ss;
        ss << std::setw(2) << std::setfill('0') << std::hex << static_cast<int>(uc);
        return ss.str();
    }

    std::vector<ch::color> some_colors = {
        ch::rgb(255,0,0),
        ch::rgb(255,97,3),
        ch::rgb(255,255,0),
        ch::rgb(0,255,0),
        ch::rgb(0,0,255),
        ch::rgb(138,43,226),
        ch::rgb(255, 192, 203),
        ch::rgb(255, 160, 122),
        ch::rgb(255, 250, 205),
        ch::rgb(152, 251, 152),
        ch::rgb(173, 216, 230),
        ch::rgb(230, 230, 250)
    };
}

std::string ch::polygon_to_svg(const ch::polygon& poly, const std::string& color, double scale) {
    std::stringstream ss;
    ss << "<path fill-rule=\"evenodd\" stroke=\"";
    ss << "none" << "\" fill=\"";
    ss << color << "\" d=\"";
    ss << svg_path_commands(poly, scale);
    ss << "\" />";
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


uint32_t ch::random_seed() {
	static std::random_device rd;    
	static std::mt19937_64 eng(rd()); 
	static std::uniform_int_distribution<uint32_t> distr;
	return distr(eng);
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

cv::Mat ch::apply_contrast(cv::Mat img, double beta, double sigma, double white_cutoff, double black_cutoff) {
	using hsv_pixel = cv::Point3_<uchar>;
	auto lookup_table = make_contrast_lookup_table(beta, sigma, white_cutoff, black_cutoff);
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

std::string ch::stippling_to_svg(std::span<const ch::point> points, double diameter) {
    std::stringstream ss;
    auto radius = diameter / 2.0;
    ss << "<g>";
    for (const auto& pt : points) {
        ss << "<circle cx=\"" << pt.x << "\" cy=\"" << pt.y << "\" r=\"" << radius << "\" />";
    }
    ss << "</g>";
    return ss.str();
}

std::string ch::to_string(double val, int precision) {
	std::stringstream ss;
	ss << std::fixed << std::setprecision(precision) << val;
	return ss.str();
}

std::string ch::to_string(const point& pt) {

	std::stringstream ss;
	ss << "( " << pt.x << " , " << pt.y << " )";
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

cv::Mat ch::stylize(cv::Mat img, double sigma_s, double sigma_r) {
    cv::Mat output;
    cv::stylization(img, output, sigma_s, sigma_r);
    return output;
}

cv::Mat ch::edge_preserving_smoothing(cv::Mat img, int flag, double sigma_s, double sigma_r) {
    cv::Mat output;
    cv::edgePreservingFilter(img, output, flag, sigma_s, sigma_r);
    return output;
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

cv::Mat ch::convert_to_1channel_gray(const cv::Mat& color, bool invert) {
	cv::Mat gray;
	cv::cvtColor(color, gray, cv::COLOR_BGR2GRAY);
	if (invert) {
		gray = cv::Scalar::all(255) - gray;
	}
	return gray;
}

double ch::degrees_to_radians(double degrees)
{
	return (std::numbers::pi * degrees) / 180.0;
}

void ch::label_map_to_visualization_img(cv::Mat img, const std::string& output_file) {
	auto uniform_rnd_int = [](int low, int high)->uint8_t {
			std::uniform_int_distribution<> ud(low, high);
			return static_cast<uint8_t>(ud(random));
		};
	using pixel = cv::Point3_<uint8_t>;
	static auto random_color = [uniform_rnd_int]()->pixel {
		return pixel(
			uniform_rnd_int(0, 255), 
			uniform_rnd_int(0, 255), 
			uniform_rnd_int(0, 255)
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

int ch::max_val_in_mat(cv::Mat mat) {
	double min_val;
	double max_val;
	cv::Point minLoc;
	cv::Point maxLoc;

	cv::minMaxLoc(mat, &min_val, &max_val, &minLoc, &maxLoc);
	return static_cast<int>(max_val);
}

double ch::normal_random(const cbrng_state& cbrng, double mean, double stddev) {
	std::counter_based_engine<std::philox4x32_prf, 1> philox{ cbrng.keys };
	std::normal_distribution<double> nd(mean, stddev);
	return nd(philox);
}

double ch::uniform_rnd(const cbrng_state& cbrng, double lower_bound, double upper_bound) {
	std::counter_based_engine<std::philox4x32_prf, 1> philox{ cbrng.keys };
	std::uniform_real_distribution<double> ud(lower_bound, upper_bound);
	return ud(philox);
}

int ch::uniform_rnd_int(const cbrng_state& cbrng, int low, int high) {
	std::counter_based_engine<std::philox4x32_prf, 1> philox{ cbrng.keys };
	std::uniform_int_distribution<> ud(low, high);
	return ud(philox);
}

ch::random_func ch::normal_rnd_func(double mean, double stddev) {
	return [=](const cbrng_state& state) {
		return normal_random(state, mean, stddev);
	};
}

ch::random_func ch::const_rnd_func(double val) {
	return [=](const cbrng_state& state) {
		return val;
	};
}

ch::cbrng_state::cbrng_state(uint32_t k1, uint32_t k2, uint32_t k3, uint32_t k4) :
	keys( {k1,k2,k3,k4})
{}

QImage ch::create_grayscale_qimage(int wd, int hgt) {
	QImage img(wd, hgt, QImage::Format::Format_Grayscale8);
	img.fill(Qt::white);
	return img;
}


QImage ch::create_compatible_qimage(int wd, int hgt) {
	QImage img(wd, hgt, QImage::Format::Format_BGR888);
	img.fill(Qt::white);
	return img;
}

cv::Mat ch::qimage_to_mat(QImage image, bool copy) {
	if (image.format() == QImage::Format::Format_BGR888) {
		cv::Mat mat(image.height(), image.width(), CV_8UC3, (cv::Scalar*)image.scanLine(0));
		return copy ? mat.clone() : mat;
	} else if (image.format() == QImage::Format::Format_Grayscale8) {
		cv::Mat mat(image.height(), image.width(), CV_8UC1, (cv::Scalar*)image.scanLine(0));
		return copy ? mat.clone() : mat;
	}
	return {};
}

QImage ch::mat_to_qimage(cv::Mat mat, bool copy) {
	QImage img;
	if (mat.channels() == 1) {
		img = QImage(mat.data, mat.cols, mat.rows, mat.step, QImage::Format_Grayscale8);
	} else if (mat.channels() == 3) {
		img = QImage(mat.data, mat.cols, mat.rows, mat.step, QImage::Format_BGR888);
	} else {
		return {};
	}
	return copy ? img.copy() : img;
}

QPen ch::create_pen(uchar color, double thickness) {
	return QPen(
		QBrush(QColor(color, color, color)),
		static_cast<float>(thickness)
	);
}

void ch::paint_polygon(QPainter& g, const polygon& poly, color col, bool filled, int thickness) {
	QPainterPath path;
	path.addPolygon(ring_to_qpolygon(poly.outer()));
	for (const auto& hole : poly.inners()) {
		path.addPolygon(ring_to_qpolygon(hole));
	}
	auto qcolor = cv_color_to_qcolor(col);
    g.setPen((thickness == 0) ? QPen(qcolor) : QPen(QBrush(qcolor), thickness));
	g.setBrush((filled) ? QBrush(qcolor) : QBrush(Qt::NoBrush));
	g.drawPath(path);
}

void ch::fill_polygon(QPainter& g, const polygon& poly, QBrush brush) {
    QPainterPath path;
    path.addPolygon(ring_to_qpolygon(poly.outer()));
    for (const auto& hole : poly.inners()) {
        path.addPolygon(ring_to_qpolygon(hole));
    }
    g.setPen(Qt::NoPen);
    g.setBrush(brush);
    g.drawPath(path);
}

cv::Mat ch::paint_polygons(const std::vector<std::tuple<color, polygon>>& polys,
		dimensions<int> sz) {
	cv::Mat mat(sz.hgt, sz.wd, CV_8UC3);
	mat.setTo(cv::Scalar(255, 255, 255));
	QImage img = mat_to_qimage(mat, false);
	QPainter g(&img);
	g.setRenderHint(QPainter::Antialiasing, true);

	for (const auto& [col, poly] : polys) {
		paint_polygon(g, poly, col);
	}

	return mat;
}

cv::Mat ch::paint_polygons(const std::vector<std::tuple<uchar, ch::polygon>>& gray_polys,
		ch::dimensions<int> sz, bool invert) {
	auto colored_polys = gray_polys |
		rv::transform(
			[invert](const auto& tup)->std::tuple<color, ch::polygon> {
				auto [gray, poly] = tup;
				gray = (invert) ? 255 - gray : gray;
				return { rgb(gray,gray,gray), poly };
			}
	) | r::to_vector;
	return paint_polygons(colored_polys, sz);
}

//  0.299 R + 0.587 G + 0.114 B
uchar ch::color_to_monochrome(ch::color col) {
	float red = static_cast<float>(col[2]);
	float green = static_cast<float>(col[1]);
	float blue = static_cast<float>(col[0]);
	float gray = 0.299 * red + 0.587 * green + 0.114 * blue;
	return static_cast<uchar>(gray);
}

cv::Mat ch::blank_monochrome_bitmap(int dim) {
    return  cv::Mat(dim, dim, CV_8U, 255);
}

double ch::measure_gray_level(cv::Mat swatch) {
    double mean_val = cv::mean(swatch).val[0];
    return (255.0 - mean_val) / 255.0;
}

/*
void ch::debug_polygons_to_svg(const std::string& output_file, double scale, 
        std::span<ch::polygon> polys) {
    auto colors = some_colors | rv::cycle | rv::take(polys.size());
    auto tups = rv::zip(colors, polys) |
        rv::transform(
            [](auto&& p)->std::tuple<ch::color, ch::polygon> {
                auto [c, poly] = p;
                return { c,poly, poly. };
            }
        ) | r::to_vector;
    ch::polygons_to_svg<color>(output_file, tups, scale);
}
*/

void ch::debug_polygons(const std::string& output_file, dimensions<int> sz,
		std::span<std::tuple<uchar, ch::polygon>> cpolys) {

	auto polys = cpolys | rv::transform([](const auto& tup) {return std::get<1>(tup); }) | r::to_vector;

	auto [x1, y1, x2, y2] = ch::bounding_rectangle(polys);
	int wd = sz.wd;
	int hgt = sz.hgt;

	auto n = std::min(cpolys.size(), some_colors.size());
	auto old_hgt = hgt;
	hgt += 35 * (n+2);
	auto colored_polys = rv::zip(some_colors, polys | rv::take(n)) |
		rv::transform(
			[](const auto& p) {return std::tuple(p.first, p.second); }
		) | r::to_vector;
	
	auto img = paint_polygons(colored_polys, dimensions<int>(wd, hgt));

	for (int i = 0; i < static_cast<int>(n); ++i) {
		cv::putText(img, //target image
			std::to_string(i).c_str(), //text
			cv::Point(50, old_hgt + 35 * (i+1)), //top-left position
			cv::FONT_HERSHEY_DUPLEX,
			1.0,
            some_colors[i], //font color
			2);
	}

	cv::imwrite(output_file, img);
}

ch::color ch::ink_shade_to_color(double ink_shade) {
    ink_shade = (1.0 - ink_shade) * 255.0;
    uchar ch = static_cast<uchar>(std::round(ink_shade));
    return rgb(ch, ch, ch);
}

ch::color ch::ink_shade_to_color(uchar ink_shade) {
    uchar channel = 255 - ink_shade;
    return rgb(channel, channel, channel);
}

ch::color ch::rgb(uchar r, uchar g, uchar b) {
	return { b,g,r };
}

cv::Mat ch::perlin_noise(const ch::dimensions<int>& sz, uint32_t seed, int octaves, double freq) {
    auto noise = cv::Mat(sz.hgt, sz.wd, CV_32FC1, 0.0f);

    siv::PerlinNoise perlin{ seed };
    auto dim = std::max(sz.wd, sz.hgt);
    double freq_per_pix = freq / dim;

    for (auto y = 0; y < sz.hgt; ++y) {
        for (auto x = 0; x < sz.wd; ++x) {
            auto value = perlin.octave2D_01(x * freq_per_pix, y * freq_per_pix, octaves);
            noise.at<float>(y,x) = value;
        }
    }

    return noise;
}

cv::Mat ch::perlin_flow_vector_field(const ch::dimensions<int>& sz, uint32_t seed1, uint32_t seed2,
        int octaves, double freq) {

    static auto to_flow_units = [](float v)->float {
        return 2.0f * (v - 0.5f);
    };

    auto field = cv::Mat( sz.hgt+1, sz.wd+1, CV_32FC2, cv::Vec2f(0,0) );
    siv::PerlinNoise perlin_x{ seed1 };
    siv::PerlinNoise perlin_y{ seed2 };
    auto dim = std::max(sz.wd, sz.hgt);
    double freq_per_pix = freq / dim;

    for (auto y = 0; y <= sz.hgt; ++y) {
        for (auto x = 0; x <= sz.wd; ++x) {
            auto x_value = perlin_x.octave2D_01(x * freq_per_pix, y * freq_per_pix, octaves);
            auto y_value = perlin_y.octave2D_01(x * freq_per_pix, y * freq_per_pix, octaves);
            field.at<cv::Vec2f>(y, x) = { to_flow_units(x_value), to_flow_units(y_value) };
        }
    }

    return field;
}

ch::dimensions<int> ch::vector_field_size(const cv::Mat& vector_field) {
    return {
        vector_field.cols - 1,
        vector_field.rows - 1
    };
}

cv::Mat ch::uniform_direction_vector_field(const ch::dimensions<int>& sz, double theta) {
    auto x_value = std::cosf(theta);
    auto y_value = std::sinf(theta);
    return cv::Mat(sz.hgt+1, sz.wd+1, CV_32FC2, cv::Vec2f(x_value, y_value));
}

cv::Mat ch::normalize_vector_field(const cv::Mat& input) {
    auto magnitude = [](const cv::Vec2f v)->float {
        return std::sqrtf(v[0] * v[0] + v[1] * v[1]);
    };
    auto largest_magnitude = 0.0f;
    for (auto y = 0; y < input.rows; ++y) {
        for (auto x = 0; x < input.cols; ++x) {
            auto mag = magnitude(input.at<cv::Vec2f>(y, x));
            if (mag > largest_magnitude) {
                largest_magnitude = mag;
            }
        }
    }
    auto normalized = input.clone();
    for (auto y = 0; y < input.rows; ++y) {
        for (auto x = 0; x < input.cols; ++x) {
            auto old = input.at<cv::Vec2f>(y, x);
            normalized.at<cv::Vec2f>(y, x) = old / largest_magnitude;
        }
    }
    return normalized;
}

cv::Mat ch::float_noise_to_grayscale(const cv::Mat mat) {
    cv::Mat scaled = 255.0f * mat;
    cv::Mat output;
    scaled.convertTo(output, CV_8UC1);
    return output;
}

float ch::interpolate_float_mat(const cv::Mat mat, const ch::point& pt) {
    ch::point pt1 = { std::trunc(pt.x), std::trunc(pt.y) };
    ch::point pt2 = pt1 + point{ 1.0f,1.0f };

    auto noise_fn = [&mat](float x, float y)->float {
        return mat.at<float>(static_cast<int>(x), static_cast<int>(y));
    };

    return bilinear_interpolation(
        noise_fn(pt1.x, pt1.y), noise_fn(pt1.x, pt2.y),
        noise_fn(pt2.x, pt1.y), noise_fn(pt2.x, pt2.y),
        pt1.x, pt2.x, pt1.y, pt2.y,
        pt.x, pt.y
    );
}

cv::Vec2f ch::interpolate_vector_field(const cv::Mat& field, const ch::point& pt) {
    ch::point pt1 = { std::trunc(pt.x), std::trunc(pt.y) };
    ch::point pt2 = pt1 + point{ 1.0f,1.0f };

    auto get_val = [&field](float x, float y, int i)->float {
        return field.at<cv::Vec2f>(static_cast<int>(x), static_cast<int>(y))[i];
    };

    return cv::Vec2f(
        bilinear_interpolation(
            get_val(pt1.x, pt1.y, 0), get_val(pt1.x, pt2.y, 0),
            get_val(pt2.x, pt1.y, 0), get_val(pt2.x, pt2.y, 0),
            pt1.x, pt2.x, pt1.y, pt2.y,
            pt.x, pt.y
        ),
        bilinear_interpolation(
            get_val(pt1.x, pt1.y, 1), get_val(pt1.x, pt2.y, 1),
            get_val(pt2.x, pt1.y, 1), get_val(pt2.x, pt2.y, 1),
            pt1.x, pt2.x, pt1.y, pt2.y,
            pt.x, pt.y
        )
    );
}

float ch::bilinear_interpolation(float q11, float q12, float q21, float q22,
        float x1, float x2, float y1, float y2,
        float x, float y) {
    float x2x1, y2y1, x2x, y2y, yy1, xx1;
    x2x1 = x2 - x1;
    y2y1 = y2 - y1;
    x2x = x2 - x;
    y2y = y2 - y;
    yy1 = y - y1;
    xx1 = x - x1;
    return 1.0 / (x2x1 * y2y1) * (
            q11 * x2x * y2y +
            q21 * xx1 * y2y +
            q12 * x2x * yy1 +
            q22 * xx1 * yy1
        );
}

void ch::debug_polys(std::string name, const std::vector<ch::polygon>& polys) {
    auto [x1, y1, wd, hgt] = ch::bounding_rectangle(polys);

    std::ofstream outfile(name);
    outfile << ch::svg_header(static_cast<int>( wd), static_cast<int>( hgt));
    for (const auto& poly : polys) {
        auto pt = representative_point(poly);

        outfile << "<path fill-rule=\"evenodd\" stroke=\"";
        outfile << "black" << "\" fill=\"";
        outfile << "none" << "\" d=\"";
        outfile << svg_path_commands(poly, 1.0);
        outfile << "\" />";
        outfile << "<circle fill=\"red\" cx=\"" << pt.x << "\" cy=\"" << pt.y << "\" r=\"" << 0.5 << "\" />";
    }
    outfile << "</svg>" << std::endl;
    outfile.close();
}


