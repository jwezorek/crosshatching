#pragma once

#include <string>
#include <functional>
#include <range/v3/all.hpp>
#include <opencv2/core.hpp>
#include <tuple>
#include <optional>

namespace ch {

    using point = cv::Point2d;
    using polyline = std::vector<point>; 
    using line_segment = std::tuple< ch::point, ch::point>;
    using rectangle = std::tuple<double, double, double, double>;

    std::string svg_header(int wd, int hgt, bool bkgd_rect = false);
    std::string polyline_to_svg(const ch::polyline& poly, double thickness, bool closed = false);
    std::string gray_to_svg_color(unsigned char gray);

    polyline scale(const polyline& poly, double scale);

    using rnd_fn = std::function<double()>;
    double normal_rnd(double mean, double stddev);
    double uniform_rnd(double lower_bound, double upper_bound);
    int uniform_rnd_int(int low, int high);
    rnd_fn normal_rnd_fn(double mean, double stddev);
    ch::rnd_fn const_rnd_fn(double val);

    double ramp(double t, double k, bool right, bool up);

    template<typename R>
    auto rotate_view(R rng, int pivot) {
        return ranges::views::concat(
            rng | ranges::views::drop(pivot), 
            rng | ranges::views::take(pivot)
        );
    }

    cv::Mat apply_contrast(cv::Mat img, double beta, double sigma);
    std::string to_string(double val, int precision);
    cv::Mat scale(cv::Mat mat, double scale);
    cv::Mat convert_to_3channel_grayscale(cv::Mat img);
    cv::Mat convert_to_1channel_gray(const cv::Mat& color);
    cv::Mat coherence_filter(cv::Mat img, int sigma, int str_sigma, float blend, int iter);
    cv::Mat anisotropic_diffusion(cv::Mat img, double alpha, double k, int iters);
    std::tuple<cv::Mat,cv::Mat> meanshift_segmentation(const cv::Mat& input, int sigmaS, float sigmaR, int minSize);
    double degrees_to_radians(double degrees);
    ch::polyline simplify_rectilinear_polygon(const ch::polyline& poly);
    std::vector<uchar> unique_gray_values(const cv::Mat& input);
    std::optional<line_segment> linesegment_rectangle_intersection( const line_segment& line_seg, const rectangle& rect );
    rectangle bounding_rectangle(const polyline& poly);
    cv::Rect union_rect_and_pt(const cv::Rect& r, cv::Point2i pt);
    int max_val_in_mat(cv::Mat mat);

    void write_label_map_visualization(cv::Mat img, const std::string& output_file);
}