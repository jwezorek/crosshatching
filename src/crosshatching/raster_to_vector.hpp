#pragma once

#include "util.hpp"
#include "geometry.hpp"
#include <vector>
#include <tuple>

namespace ch {

    using colored_polygon = std::tuple<color, polygon>;
    using gray_polygon = std::tuple<uchar, polygon>;

    std::vector<gray_polygon> raster_to_vector_grayscale(cv::Mat mat, double param);
    std::vector<colored_polygon> raster_to_vector(cv::Mat mat, double param);

    std::vector<ch::point> perform_douglas_peucker_simplification(
        const std::vector<ch::point>& pts, double param);

    namespace detail {
        std::vector<uchar> unique_1channel_values(const cv::Mat& input);
        std::vector<color> unique_3channel_values(const cv::Mat& input);
    }

    template<typename T>
    std::vector<T> unique_values(const cv::Mat& img) {
        if constexpr (std::is_same<T, uchar>::value) {
            return detail::unique_1channel_values(img);
        } else {
            return detail::unique_3channel_values(img);
        }
    }

}