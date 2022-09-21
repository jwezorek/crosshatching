#include "ink_layers.hpp"
#include "point_set.hpp"
#include <range/v3/all.hpp>
#include "qdebug.h"


#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
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

        edge_record(int id = -1) : poly_ids{id,-1}
        {}

        int count() const {
            return r::accumulate(
                poly_ids |  rv::transform( [](int i) { return (i >= 0) ? 1 : 0;} ),
                0
            );
        }

        void insert(int id) {
            auto iter = r::find_if(poly_ids, [](int id) {return id < 0; });
            if (iter == poly_ids.end()) {
                throw std::runtime_error("error constructing polygon graph");
            }
            *iter = id;
            if (poly_ids.back() >= 0) {
                r::sort(poly_ids);
            }
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

        ch::edge_table<edge_record> edge_table;
        graph g;
        ch::int_edge_set already_added;

        for (const auto& [index, poly] : enumerated_polygons(polys)) {
            for (const auto& [u, v] : ch::all_edges(poly)) {
                edge_table[{ u,v }].insert(index);
            }
        }

        auto edges = edge_table |
            rv::transform(
                [](const auto& e)->std::tuple<int,int> {
                    const auto& ary = e.second.poly_ids;
                    return { ary[0],ary[1] };
                }
            );

        for (const auto& [u, v] : edges) {
            if (v == -1 || already_added.find({ u,v }) != already_added.end()) {
                continue;
            }
            already_added.insert({ u,v });
            g[u].push_back(v);
            g[v].push_back(u);
        }

        return g;
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
    auto mat = cv::imread("C:\\test\\test_img.png");
    mat = ch::convert_to_1channel_gray(mat);
    auto polys = ch::raster_to_vector_grayscale(mat, 1.5);
    debug_polygons("C:\\test\\test_img_labeled.png", polys);
    auto graph = polygons_to_adjacency_graph(polys);
    for (const auto& [u, adj_list] : graph) {
        std::stringstream ss;
        ss << u << " : ";
        for (int v : adj_list) {
            ss << v << " ";
        }
        qDebug() << ss.str().c_str();
    }
}