#pragma once

#include "geometry.hpp"
#include "raster_to_vector.hpp"
#include "brush.hpp"
#include <vector>
#include <span>

namespace ch {

    struct ink_layer_item {
        int id;
        uchar value;
        ch::polygon poly;
        ink_layer_item* parent;
    };

    struct ink_layer {
        ch::brush_expr_ptr brush;
        std::vector<ink_layer_item> content;
    };

    using ink_layers = std::vector<ink_layer>;

    ink_layers split_into_layers(const std::vector<gray_polygon>& polys,
        std::span<const ch::brush_expr_ptr> brushes,
        std::span<const uchar> gray_levels);

    std::vector<ch::gray_polygon> to_polygons(const ink_layer& il);

    void debug_layers(cv::Mat img);
}