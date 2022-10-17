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
    using stroke_groups = std::vector<stroke_group>;
    using strokes_ptr = std::shared_ptr<stroke_groups>;

    stroke_groups transform(const stroke_groups& sg, const ch::matrix& mat);
    void paint_strokes(QPainter& g, const stroke_groups& str);
    void paint_strokes(QPainter& g, strokes_ptr str);
    std::string to_svg(const stroke_groups& s);
    void append(strokes_ptr dst, strokes_ptr src);
    strokes_ptr clip_strokes(const polygon& poly, strokes_ptr strokes);
    strokes_ptr transform(strokes_ptr strokes, const ch::matrix& mat);
    cv::Mat strokes_to_mat(strokes_ptr strokes, cv::Mat mat);
    

    namespace detail {
        struct to_polyline_tag {};
    }
}

template<typename Range>
ch::polylines operator|(Range&& lines, ch::detail::to_polyline_tag) {
    ch::polylines polys;
    polys.resize(ranges::distance(lines));
    for (auto&& [i, line] : ranges::views::enumerate(lines)) {
        polys[i] = line;
    }
    return polys;
}

constexpr ch::detail::to_polyline_tag to_polylines = {};

namespace ch {

    template<typename Rng>
    strokes_ptr single_stroke_group(Rng r, double thickness) {
        auto strokes = std::make_shared<stroke_groups>();
        strokes->emplace_back(r | to_polylines, thickness);
        return strokes;
    }

    template<typename Rng>
    strokes_ptr to_strokes(Rng r) {
        return std::make_shared<stroke_groups>(
            r | ranges::views::remove_if(
                [](const auto& sg) {
                    return sg.strokes.empty();
                }
            ) | ranges::to_vector
        );
    }

}