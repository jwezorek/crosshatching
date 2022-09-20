#pragma once

#include "geometry.hpp"
#include "raster_to_vector.hpp"
#include <vector>
#include <span>

namespace ch {

    struct ink_layer_item {
        int id;
        uchar value;
        ch::polygon poly;
        ink_layer_item* parent;
    };

    using ink_layer = std::vector<ink_layer_item>;
    using ink_layers = std::vector<ink_layer>;

    ink_layers split_into_layers(const std::vector<gray_polygon>& polys,
        std::span<const uchar> gray_levels);

    void debug_layers();
}