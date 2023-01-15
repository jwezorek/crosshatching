#pragma once

#include <vector>
#include <memory>
#include <range/v3/all.hpp>
#include "qpainter.h"
#include "geometry.hpp"

namespace ch {

    struct stroke_group {
        ch::polylines strokes;
        double thickness;
    };

    struct stippling_group {
        std::vector<ch::point> points;
        double thickness;
    };

    struct shaded_polygon {
        ch::polygon poly;
        double shade;
    };

    using drawing_component = std::variant<stroke_group, stippling_group, shaded_polygon>;
    using drawing_comps = std::vector<drawing_component>;
    using drawing_comps_ptr = std::shared_ptr<drawing_comps>;

    bool is_empty(const drawing_component& dc);
    drawing_comps transform(const drawing_comps& sg, const ch::matrix& mat);
    void paint_strokes(QPainter& g, const drawing_comps& str);
    void paint_strokes(QPainter& g, drawing_comps_ptr str);
    std::string to_svg(const drawing_comps& s);
    void append(drawing_comps_ptr dst, drawing_comps_ptr src);
    drawing_comps_ptr clip_strokes(const polygon& poly, drawing_comps_ptr strokes);
    drawing_comps_ptr transform(drawing_comps_ptr strokes, const ch::matrix& mat);
    cv::Mat strokes_to_mat(drawing_comps_ptr strokes, cv::Mat mat);
    cv::Mat strokes_to_mat(const drawing_comps& strokes, cv::Mat);

    namespace detail {
        struct to_polyline_tag {};
        struct to_strokes_tag {};
    }
}

constexpr ch::detail::to_polyline_tag to_polylines = {};
template<typename Range>
ch::polylines operator|(Range&& lines, ch::detail::to_polyline_tag) {
    ch::polylines polys;
    polys.resize(ranges::distance(lines));
    for (auto&& [i, line] : ranges::views::enumerate(lines)) {
        polys[i] = line;
    }
    return polys;
}

constexpr ch::detail::to_strokes_tag to_strokes = {};
template<typename Range>
ch::drawing_comps_ptr operator|(Range&& dcs, ch::detail::to_strokes_tag) {
    return std::make_shared<ch::drawing_comps>(
        dcs | ranges::views::remove_if(
            [](const auto& dc) {
                return ch::is_empty(dc);
            }
        ) | ranges::to_vector
   );
}

namespace ch {
    template<typename T>
    drawing_comps_ptr single_stroke_group(T&& sg) {
        auto dcs_ptr = std::make_shared<drawing_comps>();
        dcs_ptr->emplace_back(std::move(sg));
        return dcs_ptr;
    }
}