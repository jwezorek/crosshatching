#pragma once

#include "geometry.hpp"
#include "raster_to_vector.hpp"
#include "brush.hpp"
#include <vector>
#include <span>

/*------------------------------------------------------------------------------------------------*/

namespace ch {

    using brush_token = uint64_t;

    struct ink_layer_item {
        int id;
        int layer_id;
        uchar value;
        ch::polygon poly;
        ink_layer_item* parent;

        brush_token token() const;
        brush_token parent_token() const;

    };

    struct ink_layer {
        ch::brush_expr_ptr brush;
        std::vector<ink_layer_item> content;
    };

    struct ink_layers {
        dimensions<double> sz;
        std::vector<ink_layer> content;

        ink_layers clone() const;
    };

    ink_layers scale(const ink_layers& il, double scale_factor);

    ink_layers split_into_layers(
        const dimensions<double>& sz,
        const std::vector<gray_polygon>& polys,
        std::span<const ch::brush_expr_ptr> brushes,
        std::span<const uchar> gray_levels);

    std::vector<ch::gray_polygon> to_polygons(const ink_layer& il);

    void debug_layers(cv::Mat img);
}