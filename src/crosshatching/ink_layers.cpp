#include "ink_layers.hpp"
#include "point_set.hpp"
#include <range/v3/all.hpp>
#include "qdebug.h"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <map>
#include <memory>
#include <unordered_set>
#include <sstream>
#include <functional>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>

/*------------------------------------------------------------------------------------------------------*/

namespace r = ranges;
namespace rv = ranges::views;

namespace {
    auto gray_ranges(std::span<const uchar> gray_levels) {
        std::vector<uchar> end_of_ranges = gray_levels | r::to_vector;
        if (end_of_ranges.back() != 255) {
            end_of_ranges.push_back(255);
        }
        int n = static_cast<int>(end_of_ranges.size());
        auto start_of_ranges = rv::concat(
            rv::single(static_cast<uchar>(1)),
            end_of_ranges |
            rv::take(n - 1) |
            rv::transform(
                [](uchar v)->uchar {return v + 1; }
            )
        );
        auto value_ranges = rv::zip(start_of_ranges, end_of_ranges) |
                rv::transform(
                    [](const auto& pair)->std::tuple<double, double> {
                        auto [from, to] = pair;
                        return { from,to };
                    }
            ) | r::to_vector;
        return value_ranges | rv::reverse | r::to_vector;
    }

    class adjacency_graph {
        std::unordered_map<int, std::set<int>> graph_;

        void replace_neighbor(int u, int v, int new_v) {
            auto iter = graph_.find(u);
            if (iter == graph_.end()) {
                throw std::runtime_error("replace_neighbor");
            }
            auto& neighbors_of_u = iter->second;
            neighbors_of_u.erase(v);
            neighbors_of_u.insert(new_v);
        }

    public:
        adjacency_graph(size_t sz) {
            graph_.reserve(sz);
        }

        void insert_edge(int u, int v) {
            graph_[u].insert(v);
            graph_[v].insert(u);
        }

        void remove_edge(int u, int v) {
            graph_[u].erase(v);
            graph_[v].erase(u);
        }

        bool contains_edge(int u, int v) const {
            auto iter = graph_.find(u);
            if (iter == graph_.end()) {
                return false;
            }
            const auto& neighbors_of_u = iter->second;
            return neighbors_of_u.find(v) != neighbors_of_u.end();
        }
        
        auto adjacent_verts(int v) const {
            static std::set<int> dummy;
            const std::set<int>* adj_list = nullptr;
            if (graph_.find(v) == graph_.end()) {
                adj_list = &dummy;
            } else {
                adj_list = &graph_.at(v);
            }
            return rv::all(*adj_list);
        }

        int merge_vertices(int v1, int v2) {

            remove_edge(v1, v2);

            // order v1 and v2 such that v1 is the lower index.
            if (v1 > v2) {
                std::swap(v1, v2);
            }

            // make all of v2's neighbors point to v1 instead of v2
            for (int v2_neighbor : adjacent_verts(v2)) {
                replace_neighbor(v2_neighbor, v2, v1);
            }

            // merge v2's neighbors into v1's neighbor list
            auto& v1s_neighbors = graph_[v1];
            const auto& v2_neighbors = graph_[v2];
            v1s_neighbors.insert(v2_neighbors.begin(), v2_neighbors.end());

            // delete v2
            graph_.erase(v2);

            return v1;
        }

        void join(r::any_view<int> ids) {
            int id = -1;
            for (int joinee : ids) {
                if (id < 0) {
                    id = joinee;
                    continue;
                }
                id = merge_vertices(id, joinee);
            }
        }
    };

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

    adjacency_graph polygons_to_adjacency_graph(const std::vector<ch::gray_polygon>& polys) {

        ch::edge_table<edge_record> edge_table;
        adjacency_graph g(2*polys.size());

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
            if (v == -1 || g.contains_edge(u,v)) {
                continue;
            } else {
                g.insert_edge(u, v);
            }
        }

        return g;
    }

    class polygon_layer_factory {
        using cluster_info = std::tuple<uchar, std::vector<int>>;

        class cluster {
            int id_;
            uchar color_;
            double area_;
            std::shared_ptr<std::vector<int>> component_ids_;
        public:
            cluster(int id = -1, double area = 0.0, uchar color = 0)  :
                    id_(id), area_(area), color_(color)  {
                if (id_ >= 0) {
                    component_ids_ = std::make_shared<std::vector<int>>();
                    component_ids_->push_back(id_);
                }
            }

            cluster( uchar color, r::any_view<cluster> children) :
                    color_(color),
                    component_ids_(
                        std::make_shared<std::vector<int>>(
                            children | 
                                rv::transform([](const auto& c) {
                                    return rv::all(c.components()); }
                                ) | 
                                rv::join |
                                r::to_vector
                        )
                    ),
                    area_(
                        r::accumulate(
                            children | rv::transform([](const auto& c) {return c.area(); }),
                            0.0
                        )
                    ) {
                id_ = r::min(*component_ids_);
            }

            cluster_info to_cluster_info() const {
                return { color_, *component_ids_ };
            }

            int id() const {
                return id_;
            }

            uchar color() const {
                return color_;
            }

            void set_color(uchar c) {
                color_ = c;
            }

            double area() const {
                return area_;
            }

            const std::vector<int>& components() const {
                return *component_ids_;
            }
        };
        
        class cluster_set {
            std::unordered_map<int, cluster> set_;
            adjacency_graph* graph_;

        public:
            using value_type = std::unordered_map<int, cluster>::value_type;

            cluster_set() : graph_(nullptr)
            {}

            cluster_set(adjacency_graph& graph, r::any_view<value_type> vals = {}) :
                graph_(&graph)
            {
                set_ = vals | r::to_< std::unordered_map<int, cluster>>();
            }

            size_t size() const {
                return set_.size();
            }

            void insert(const cluster& c) {
                set_[c.id()] = c;
            }

            auto contents() const {
                return set_ |
                    rv::transform(
                        [](const value_type& v) {
                            return v.second;
                        }
                );
            }

            bool empty() const {
                return set_.empty();
            }

            std::vector<cluster> clusters_by_area() const
            {
                auto clusters = contents() | r::to_vector;
                r::sort(clusters,
                    [](const auto& lhs, const auto& rhs) {
                        return lhs.area() > rhs.area();
                    }
                );
                return clusters;
            }

            auto adjacent_cluster_ids(int c_id) {
                return graph_->adjacent_verts(c_id) |
                    rv::remove_if(
                        [this](int v) {
                            return set_.find(v) == set_.end();
                        }
                );
            }

            auto adjacent_clusters(int c_id) {
                return adjacent_cluster_ids(c_id) |
                    rv::transform([&](int id) {return this->at(id); });
            }

            cluster at(int c_id) {
                return set_.at(c_id);
            }

            std::vector<cluster> select_merge_targets(const std::vector<cluster>& candidates) {
                std::unordered_map<uchar, double> areas;
                for (const auto& c : candidates) {
                    areas[c.color()] += c.area();
                }
                auto iter = r::max_element(areas,
                    [](const auto& lhs, const auto& rhs) {
                        return lhs.second > rhs.second;
                    }
                );
                uchar selected_color = iter->first;
                return candidates |
                    rv::remove_if(
                        [selected_color](const auto& c) {
                            return c.color() != selected_color;
                        }
                ) | r::to_vector;
            }

            bool contains(const cluster& c) {
                return set_.find(c.id()) != set_.end();
            }

            cluster remove(int id) {
                cluster c = set_.at(id);
                set_.erase(id);
                return c;
            }

            void merge_clusters(r::any_view<cluster> mergees) {
                uchar color = 0;
                for (const auto& c : mergees) {
                    if (c.color() > 0) {
                        color = c.color();
                    }
                    if (contains(c)) {
                        remove(c.id());
                    }
                }

                graph_->join(
                    mergees |  rv::transform([](const auto& c) {return c.id(); })
                );

                auto new_cluster = cluster{ color, mergees };
                insert(new_cluster);
            }

            std::vector<cluster> merge(const cluster_set& src) {
                bool merged_one = false;
                std::vector<cluster> slop;
                slop.reserve(src.size());

                do {
                    bool merged_one = false;
                    for (auto ex_clust : src.clusters_by_area()) {

                        auto neighbors = adjacent_clusters(ex_clust.id()) |  r::to_vector;
                        if (! neighbors.empty()) {
                            auto targets = select_merge_targets(neighbors);
                            auto mergees = rv::concat(rv::single(ex_clust), targets) | r::to_vector;
                            merge_clusters(mergees);
                            merged_one = true;
                        } else {
                            slop.push_back(ex_clust);
                        }
                    }
                } while (merged_one);

                slop.shrink_to_fit();
                return slop;
            }

            void insert(const std::vector<cluster>& clusts) {
                for (const auto& c : clusts) {
                    insert(c);
                }
            }


            void flatten() {
                for (auto& [id, val] : set_) {
                    val.set_color(0);
                }

                using id_set = std::set<int>;
                std::vector<std::vector<cluster>> connected_components;
                std::unordered_set<int> visited;

                auto dfs = [&](int start_id)->id_set {
                    std::vector<int> stack = { start_id };
                    
                    id_set connected_component;
                    while (!stack.empty()) {
                        auto c_id = stack.back();
                        stack.pop_back();

                        if (visited.find(c_id) != visited.end()) {
                            continue;
                        }
                        visited.insert(c_id);
                        connected_component.insert(c_id);

                        for (const auto& neighbor_id : adjacent_cluster_ids(c_id)) {
                            stack.push_back(neighbor_id);
                        }
                    }

                    return connected_component;
                };

                for (const auto& [index, c] : set_) {
                    if (visited.find(c.id()) != visited.end()) {
                        continue;
                    }
                    auto ids = dfs(c.id());
                    connected_components.push_back(
                        ids |
                        rv::transform([this](int id) {return set_.at(id); }) |
                        r::to_vector
                    );
                }

                for (const auto& cc : connected_components) {
                    merge_clusters(cc);
                }
            }
        };

        const std::vector<ch::gray_polygon>* polys_;
        adjacency_graph graph_;
        cluster_set prev_layer_;

        std::vector<int> polys_in_gray_range(uchar from, uchar to) {
            return rv::enumerate(*polys_) |
                rv::transform(
                    [from,to](const auto& id_val_pair)->int {
                        const auto& [id, val] = id_val_pair;
                        const auto& [gray, poly] = val;
                        if (gray >= from && gray <= to) {
                            return id;
                        }
                        return -1;
                    }
                ) |
                rv::remove_if(
                    [](int id) {return id < 0; }
                ) |
                r::to_vector;
        }

        cluster_set seed_clusters(uchar from_gray, uchar to_gray) {
            auto seeds = polys_in_gray_range(from_gray, to_gray);
            auto vals = seeds |
                rv::transform(
                    [this](int id)->cluster_set::value_type {
                        const auto& [color, poly] = polys_->at(id);
                        return { 
                            id, 
                            cluster(
                                id, 
                                boost::geometry::area(poly),
                                color
                            ) 
                        };
                    }
            );
            return cluster_set(graph_, vals);
        }

    public:
        polygon_layer_factory(const std::vector<ch::gray_polygon>& polys) :
                polys_(&polys),
                graph_(polygons_to_adjacency_graph(polys))
        {
        }

        const adjacency_graph& graph() const {
            return graph_;
        }

        std::vector<cluster_info> next_layer(uchar from_gray, uchar to_gray) {
            auto clusters = seed_clusters(from_gray, to_gray);

            std::vector<cluster> slop;
            if (!prev_layer_.empty()) {
                slop = clusters.merge(prev_layer_);
            }

            auto output = clusters.contents() |
                rv::transform(
                    [](const cluster& c)->cluster_info {
                        return c.to_cluster_info();
                    }
                ) | r::to_vector;

            clusters.insert(slop);
            prev_layer_ = clusters;
            prev_layer_.flatten();

            return output;
        }
    };

    ch::polygon join_polygons(const ch::polygon& poly1, const ch::polygon& poly2) {
        ch::polygons output;
        boost::geometry::union_(poly1, poly2, output);
        if (output.empty()) {
            std::vector<std::tuple<ch::color, ch::polygon>> debug = {
                {ch::rgb(255,0,0), poly1},
                {ch::rgb(0,0,255), poly2}
            };
            ch::polygons_to_svg<ch::color>("C:\\test\\debug.svg", debug);
            qDebug() << "boost::geometry::union_ failure";
            return {};
        }
        return output[0];
    }

    ch::polygon join_polygons(const std::vector<ch::polygon>& polys, const adjacency_graph& graph,
            const std::vector<int> ids) {

        std::unordered_set<int> id_set = ids | r::to<std::unordered_set<int>>();
        std::unordered_set<int> visited;

        ch::polygon joined;
        std::function<void(int)> dfs;
        dfs = [&](int id) {
            if (visited.find(id) != visited.end()) {
                return;
            }
            visited.insert(id);
            if (joined.outer().empty()) {
                joined = polys[id];
            } else {
                joined = join_polygons(joined, polys[id]);
            }

            auto neighbors_in_cluster = graph.adjacent_verts(id) |
                rv::remove_if(
                    [&](int v) {
                        return id_set.find(v) == id_set.end(); 
                    }
                );

            for (auto neighbor : neighbors_in_cluster) {
                dfs(neighbor);
            }
        };
        dfs(ids.front());

        return joined;
    }
}

ch::ink_layers ch::split_into_layers(const std::vector<gray_polygon>& cpolys,
        std::span<const ch::brush_expr_ptr> brush_expr_ptrs,
        std::span<const uchar> gray_levels) {

    polygon_layer_factory layer_factory(cpolys);
    auto adj_graph = layer_factory.graph();
    auto polygons = cpolys | rv::transform([](const auto& p) {return std::get<1>(p); }) | r::to_vector;
    int item_id = 0;
    auto ranges = gray_ranges(gray_levels);
    std::vector<ch::brush_expr_ptr> brushes = (!brush_expr_ptrs.empty()) ?
        brush_expr_ptrs | r::to_vector :
        std::vector<ch::brush_expr_ptr>(ranges.size(), nullptr);
    auto layer_content = ranges |
        rv::transform(
            [&](const auto& rng_tup)->std::vector<ink_layer_item> {
                auto [from_gray, to_gray] = rng_tup;
                auto clusters = layer_factory.next_layer(from_gray, to_gray);
                return clusters | 
                    rv::transform(
                        [&](const auto& ci)->ink_layer_item {
                            const auto& [value, poly_ids] = ci;
                            return {
                                item_id++,
                                value,
                                join_polygons(polygons, adj_graph, poly_ids),
                                nullptr
                            };
                        }
                    ) | r::to_vector;
            }
        );

    return rv::zip(brushes, layer_content) |
        rv::transform(
            [](auto&& pair)->ink_layer {
                auto&& [brush, content] = pair;
                return { brush, content };
            }
    ) | r::to_vector;
}

std::vector<ch::gray_polygon> ch::to_polygons(const ink_layer& il) {
    return il.content |
        rv::transform(
            [](const ink_layer_item& ili)->ch::gray_polygon {
                return { ili.value, ili.poly };
            }
    ) | r::to_vector;
}

void ch::debug_layers(cv::Mat img) {
    auto mat = img.clone();
    mat = ch::convert_to_1channel_gray(mat, true);
    auto all_polys = ch::raster_to_vector_grayscale(mat,1.4);
    auto polys = all_polys |
            rv::remove_if(
                [](const gray_polygon& gp) {
                    return std::get<0>(gp) == 0 || std::get<1>(gp).outer().size() == 2;
                }
        ) | r::to_vector;
    /*
    auto layers = ch::split_into_layers(polys, {}, { {60, 128, 255} });
    for (const auto& [index, layer] : rv::enumerate(layers)) {
        std::string fname = "C:\\test\\layers\\test-" + std::to_string(index) + ".svg";
        auto gpolys = layer |
            rv::transform(
                [](const ch::ink_layer_item& ili)->ch::gray_polygon {
                    return {
                        255-ili.value,
                        ili.poly
                    };
                }
            ) | r::to_vector;
        if (!gpolys.empty()) {
            polygons_to_svg<uchar>(fname, gpolys);
        }
    }
    */
}