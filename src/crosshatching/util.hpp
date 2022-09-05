#pragma once

#include <string>
#include <functional>
#include <range/v3/all.hpp>
#include <opencv2/core.hpp>
#include <tuple>
#include <optional>
#include <optional>
#include <span>
#include <array>
#include "geometry.hpp"

namespace ch {

    using color = cv::Vec3b;

    namespace detail {
        std::vector<uchar> unique_1channel_values(const cv::Mat& input);
        std::vector<color> unique_3channel_values(const cv::Mat& input);
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
    std::string to_svg_color(uchar gray);
    std::string to_svg_color(const color& c);

    template<typename T>
    void polygons_to_svg(const std::string& output_file,
            std::span<std::tuple<T, polygon>> polys,
            double scale = 1.0) {
        detail::polygons_to_svg_aux(output_file, polys, scale);
    }

    void debug_polygons(const std::string& output_file, std::span<polygon> polys);

    struct cbrng_state {
        std::array<uint32_t, 4> keys;
        cbrng_state(uint32_t k1 = 0, uint32_t k2 = 0, uint32_t k3 = 0, uint32_t k4 = 0);
    };

    using random_func = std::function<double(const cbrng_state&)>;
    double normal_random(const cbrng_state& rnd, double mean, double stddev);
    double uniform_rnd(const cbrng_state& rnd, double lower_bound, double upper_bound);
    int uniform_rnd_int(const cbrng_state& rnd, int low, int high);
    random_func normal_rnd_func(double mean, double stddev);
    random_func const_rnd_func(double val);

    using rnd_fn = std::function<double()>;
    uint32_t random_seed();
    double normal_rnd(double mean, double stddev);
    double uniform_rnd(double lower_bound, double upper_bound);
    int uniform_rnd_int(int low, int high);
    rnd_fn normal_rnd_fn(double mean, double stddev);
    ch::rnd_fn const_rnd_fn(double val);

    // image processing
    cv::Mat apply_contrast(cv::Mat img, double beta, double sigma);
    cv::Mat scale(cv::Mat mat, double scale);
    cv::Mat convert_to_3channel_grayscale(cv::Mat img);
    cv::Mat convert_to_1channel_gray(const cv::Mat& img);
    cv::Mat coherence_filter(cv::Mat img, int sigma, int str_sigma, float blend, int iter);
    cv::Mat anisotropic_diffusion(cv::Mat img, double alpha, double k, int iters);
    std::tuple<cv::Mat,cv::Mat> meanshift_segmentation(const cv::Mat& input, int sigmaS, float sigmaR, int minSize);
    int max_val_in_mat(cv::Mat mat);
    void label_map_to_visualization_img(cv::Mat img, const std::string& output_file);
    dimensions<int> mat_dimensions(cv::Mat mat);
    std::vector<ch::point> perform_douglas_peucker_simplification(
        const std::vector<ch::point>& pts, double param);

    template<typename T>
    std::vector<T> unique_values(const cv::Mat& img) {
        if constexpr (std::is_same<T, uchar>::value) {
            return detail::unique_1channel_values(img);
        } else {
            return detail::unique_3channel_values(img);
        }
    }

    // etc.
    std::string to_string(double val, int precision);
    double degrees_to_radians(double degrees);
    double ramp(double t, double k, bool right, bool up);
    std::vector<cv::Point_<float>> to_float_points(std::span<const ch::point> pts);
    std::vector<ch::point> from_float_points(std::span<const cv::Point_<float>> pts);
}