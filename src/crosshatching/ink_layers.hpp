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
        ch::brush_expr_ptr brush;
        ink_layer_item* parent;

        // hashable ID for this item's value, brush, and the
        // full chain of parent tokens
        brush_token token() const;

        // the above of this item's parent
        brush_token parent_token() const;

        // hashable id of this item's full chain of parent tokens and 
        // brush but not its value;
        brush_token brush_token() const;
    };

    using ink_layer = std::vector<ink_layer_item>;

    int ink_layer_index(const ink_layer& il);

    struct ink_layers {
        dimensions<double> sz;
        std::vector<ink_layer> content;

        ink_layers clone() const;
        bool empty() const;
        int count() const;
        void clear();
    };

    ink_layers scale(const ink_layers& il, double scale_factor);

    ink_layers split_into_layers_adaptive(
        const dimensions<double>& sz,
        const std::vector<gray_polygon>& polys,
        std::span<const ch::brush_expr_ptr> brushes,
        std::span<const uchar> gray_levels);

    ink_layers split_into_layers_simple(
        const dimensions<double>& sz,
        const std::vector<gray_polygon>& polys,
        std::span<const ch::brush_expr_ptr> brushes,
        std::span<const uchar> gray_levels);

    std::vector<ch::gray_polygon> to_polygons(const ink_layer& il);

    void debug_layers(cv::Mat img);
}