#pragma once

#include <string>
#include <functional>
#include <tuple>
#include <optional>
#include <optional>
#include <span>
#include <array>
#include "qpainter.h"
#include "qimage.h"
#include <range/v3/all.hpp>
#include <opencv2/core.hpp>
#include "geometry.hpp"

/*------------------------------------------------------------------------------------------------*/

namespace ch {

    using color = cv::Vec3b;

    namespace detail {
        void polygons_to_svg_aux(const std::string& output_file,
            std::span<std::tuple<color, polygon>> polys,
            double scale);
        void polygons_to_svg_aux(const std::string& output_file,
            std::span<std::tuple<uchar, polygon>> polys,
            double scale);
    }

    // SVG
    std::string svg_header(int wd, int hgt, bool bkgd_rect = false);
    std::string polyline_to_svg(std::span<const point> poly, double thickness, 
            bool closed = false);
    std::string stippling_to_svg(std::span<const point> points, double thickness);
    std::string to_svg_color(uchar gray);
    std::string to_svg_color(const color& c);

    template<typename T>
    void polygons_to_svg(const std::string& output_file,
            std::span<std::tuple<T, polygon>> polys,
            double scale = 1.0) {
        detail::polygons_to_svg_aux(output_file, polys, scale);
    }

    void debug_polygons(const std::string& output_file, dimensions<int> sz,
        std::span<std::tuple<uchar, polygon>> polys);

    // counter-based RNG...
    struct cbrng_state {
        std::array<uint32_t, 4> keys;
        cbrng_state(uint32_t k1 = 0, uint32_t k2 = 0, uint32_t k3 = 0, uint32_t k4 = 0);
    };
    uint32_t random_seed();
    using random_func = std::function<double(const cbrng_state&)>;
    double normal_random(const cbrng_state& rnd, double mean, double stddev);
    double uniform_rnd(const cbrng_state& rnd, double lower_bound, double upper_bound);
    int uniform_rnd_int(const cbrng_state& rnd, int low, int high);
    random_func normal_rnd_func(double mean, double stddev);
    random_func const_rnd_func(double val);

    // image processing
    cv::Mat apply_contrast(cv::Mat img, double beta, double sigma);
    cv::Mat scale(cv::Mat mat, double scale);
    cv::Mat convert_to_3channel_grayscale(cv::Mat img);
    cv::Mat convert_to_1channel_gray(const cv::Mat& img, bool invert = false);
    cv::Mat coherence_filter(cv::Mat img, int sigma, int str_sigma, float blend, int iter);
    cv::Mat anisotropic_diffusion(cv::Mat img, double alpha, double k, int iters);
    std::tuple<cv::Mat,cv::Mat> meanshift_segmentation(
        const cv::Mat& input, int sigmaS, float sigmaR, int minSize);
    int max_val_in_mat(cv::Mat mat);
    void label_map_to_visualization_img(cv::Mat img, const std::string& output_file);
    dimensions<int> mat_dimensions(cv::Mat mat);
    uchar color_to_monochrome(color col);
    cv::Mat blank_monochrome_bitmap(int sz);

    // perlin noise

    cv::Mat perlin_noise(const ch::dimensions<int>& sz, uint32_t seed, int octaves, double freq);
    cv::Mat float_noise_to_grayscale(const cv::Mat mat);

    // painting with Qt...
    QImage create_grayscale_qimage(int wd, int hgt);
    QImage create_compatible_qimage(int wd, int hgt);
    cv::Mat qimage_to_mat(QImage img, bool copy = true);
    QImage mat_to_qimage(cv::Mat mat, bool copy = true);
    QPen create_pen(uchar color, double thickness);
    void paint_polygon(QPainter& g, const polygon& poly, color col);
    cv::Mat paint_polygons(const std::vector<std::tuple<color, polygon>>& polys,
        dimensions<int> sz);
    cv::Mat paint_polygons(const std::vector<std::tuple<uchar, polygon>>& polys,
        dimensions<int> sz, bool invert = false);

    // etc.
    std::string to_string(double val, int precision);
    std::string to_string(const point& pt);
    double degrees_to_radians(double degrees);
    double ramp(double t, double k, bool right, bool up);

    color rgb(uchar r, uchar g, uchar b);
}