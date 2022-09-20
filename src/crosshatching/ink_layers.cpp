#include "ink_layers.hpp"
#include "point_set.hpp"
#include <range/v3/all.hpp>
#include "qdebug.h"

/*------------------------------------------------------------------------------------------------------*/

namespace r = ranges;
namespace rv = ranges::views;

namespace {
    std::vector<std::pair<uchar, uchar>> gray_ranges(std::span<const uchar> gray_levels) {
        std::vector<uchar> end_of_ranges = gray_levels | r::to_vector;
        if (end_of_ranges.back() != 255) {
            end_of_ranges.push_back(255);
        }
        int n = static_cast<int>(end_of_ranges.size());
        auto start_of_ranges = rv::concat(
            rv::single(static_cast<uchar>(0)),
            end_of_ranges |
            rv::take(n - 1) |
            rv::transform(
                [](uchar v)->uchar {return v - 1; }
            )
        );
        return rv::zip(start_of_ranges, end_of_ranges) | r::to_vector;
    }

    using graph = std::unordered_map<int, std::vector<int>>;

    struct edge_record {
        std::array<int, 2> poly_ids;

        edge_record(int id) : poly_ids{id,-1}
        {}

        int count() const {
            return r::accumulate(
                poly_ids |  rv::transform( [](int i) { return (i >= 0) ? 1 : 0;} ),
                0
            );
        }

        void insert(int id) {
            auto iter = r::find_if(poly_ids, [](int id) {return id >= 0; });
            if (iter == poly_ids.end()) {
                throw std::runtime_error("error constructing polygon graph");
            }
            *iter = id;
        }
    };

    auto enumerated_polygons(const std::vector<ch::gray_polygon>& polys) {
        return rv::enumerate(polys |
            rv::transform([](const auto& gp) {
                return std::get<1>(gp); }
            )
        );
    }

    graph polygons_to_adjacency_graph(const std::vector<ch::gray_polygon>& polys) {

        ch::edge_table<edge_record> edges;
        for (const auto& [index, poly] : enumerated_polygons(polys)) {
        //    for (con)
        }
        return {};
    }

}

ch::ink_layers ch::split_into_layers(const std::vector<gray_polygon>& polys,
        std::span<const uchar> gray_levels) {
    edge_table<int> test;
    test[std::tuple{point{ 40.0,40.0 }, point{ 60.0,120.0 }}] = 42;

    auto hello = test[std::tuple{ point{  60.0,120.0 }, point{40.0,40.0 } }];
    return {};
}

void ch::debug_layers() {
    edge_table<int> test;
    auto poly = ch::make_polygon({
            { 100,300 }, { 200,500 }, { 300,300 }, { 400,500 },
            { 600,500 }, { 700,300 }, { 600,100 }, { 500,300 },
            { 400,100 }, { 200,100 }, { 100,300 }
        }, {
            {{ 400, 200 }, { 450,300 }, { 400,400 }, { 350,300 }, {400,200}},
            {{250,250}, {250,350}, {150,350}, {150,250}, {250,250}}
        }
    );

    for (const auto& [u, v] : all_edges(poly)) {
        qDebug() << "{" << to_string(u).c_str() << " , " << to_string(v).c_str() << "}";
    }
}