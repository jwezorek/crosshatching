#include "ink_layers.hpp"
#include "point_set.hpp"
#include "qdebug.h"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <map>
#include <memory>
#include <unordered_set>
#include <sstream>
#include <functional>
#include <boost/geometry.hpp>
#include <chrono>
#include <array>
#include <boost/functional/hash.hpp>

/*------------------------------------------------------------------------------------------------*/

namespace r = ranges;
namespace rv = ranges::views;

namespace {

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

    ch::polygon join_polygons(const ch::polygon& poly1, const ch::polygon& poly2) {
        ch::polygons output;
        boost::geometry::union_(poly1, poly2, output);
        if (output.empty()) {
            qDebug() << "boost::geometry::union_ failure";
            return poly1.outer().size() > poly2.outer().size() ? poly1 : poly2;
        }
        return output[0];
    }
    
    ch::polygon join_polygons(const std::vector<std::tuple<int,ch::polygon>>& id_poly_pairs, const adjacency_graph& graph) {

        std::unordered_set<int> id_set = id_poly_pairs | rv::transform(
            [](const auto& tup) {
                const auto& [id, p] = tup;
                return id;
            }
        ) | r::to<std::unordered_set<int>>();

        std::unordered_map<int, ch::polygon> polys = id_poly_pairs | 
            rv::transform(
                [](const auto& tup)->std::unordered_map<int, ch::polygon>::value_type {
                    const auto& [id, p] = tup;
                    return { id,p };
                }
            ) | r::to<std::unordered_map<int, ch::polygon>>();

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
        dfs(*id_set.begin());

        return joined;
    }
    
    /*
    ch::polygon join_polygons(const std::vector<std::tuple<int, ch::polygon>>& id_poly_pairs, const adjacency_graph& graph) {

        std::unordered_set<int> id_set = id_poly_pairs | rv::transform(
            [](const auto& tup) {
                const auto& [id, p] = tup;
                return id;
            }
        ) | r::to<std::unordered_set<int>>();

        std::unordered_map<int, ch::polygon> polys = id_poly_pairs |
            rv::transform(
                [](const auto& tup)->std::unordered_map<int, ch::polygon>::value_type {
                    const auto& [id, p] = tup;
                    return { id,p };
                }
        ) | r::to<std::unordered_map<int, ch::polygon>>();

        std::unordered_set<int> visited;
        ch::polygon joined;
        
        std::vector<int> stack;
        stack.reserve(id_poly_pairs.size() * 20);
        stack.push_back(*id_set.begin());

        while (!stack.empty()) {
            int id = stack.back();
            stack.pop_back();

            if (visited.find(id) != visited.end()) {
                continue;
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
                stack.push_back(neighbor);
            }
        }

        return joined;
    }
    */

    auto gray_ranges(std::span<const uchar> gray_levels, bool reverse) {
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
                    [](const auto& pair)->std::tuple<uchar, uchar> {
                        auto [from, to] = pair;
                        return { from,to };
                    }
            ) | r::to_vector;

        if (reverse)
            return value_ranges | rv::reverse | r::to_vector;
        else
            return value_ranges | r::to_vector;
    }

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
                return;
                //throw std::runtime_error("error constructing polygon graph");
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
        using cluster_info = std::tuple<uchar, std::vector<int>, ch::polygon>;

        class cluster {
            int id_;
            uchar color_;
            double area_;
            std::shared_ptr<std::vector<int>> component_ids_;
            std::shared_ptr<ch::polygon> poly_ptr_;
        public:
            cluster(const ch::polygon& poly = {}, int id = -1, uchar color = 0) :
                    id_(id), 
                    poly_ptr_( (!poly.outer().empty()) ? std::make_shared<ch::polygon>(poly) : nullptr ),
                    color_(color)  {
                if (id_ >= 0) {
                    component_ids_ = std::make_shared<std::vector<int>>();
                    component_ids_->push_back(id_);
                }
                area_ = (poly_ptr_) ? boost::geometry::area(poly) : 0.0;
            }

            cluster(const ch::polygon& poly, uchar color, r::any_view<cluster> children) :
                    poly_ptr_(std::make_shared<ch::polygon>(poly)),
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
                return { color_, *component_ids_, *poly_ptr_ };
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

            std::shared_ptr<ch::polygon> poly_ptr() const {
                return poly_ptr_;
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

            std::vector<cluster> clusters_by_color() const {
                auto clusters = contents() | r::to_vector;
                r::sort(clusters,
                    [](const auto& lhs, const auto& rhs) {
                        return lhs.color() > rhs.color();
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
                uchar selected_color = r::max(
                    candidates | 
                        rv::transform(
                            [](auto&& c)->uchar {
                                return c.color();
                            }
                        )
                );
                return candidates |
                    rv::remove_if(
                        [selected_color](const auto& c) {
                            return c.color() != selected_color;
                        }
                ) | r::to_vector;
            }

            /*
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
            */

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

                auto ids = mergees | rv::transform([](const auto& c) {return c.id(); });
                std::vector<std::tuple<int,ch::polygon>> polys = mergees | 
                    rv::transform(
                        [](const auto& c)->std::tuple<int, ch::polygon> {
                            return {
                                c.id(),
                                *c.poly_ptr()
                            };
                        }
                    ) | r::to_vector;

                auto poly = join_polygons(polys, *graph_);
                graph_->join(ids);
                auto new_cluster = cluster{ poly, color, mergees };
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
                                poly,
                                id, 
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

        std::vector<cluster_info> next_layer(uchar from_gray, uchar to_gray, size_t n) {
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

            if (n > 1) {
                clusters.insert(slop);
                prev_layer_ = clusters;
                prev_layer_.flatten();
            }
            return output;
        }
    };

    void populate_parents(const std::unordered_map<int, std::vector<int>>& layer_item_to_components,
            ch::ink_layers& layers) {
        std::unordered_map<int, ch::ink_layer_item*> component_id_to_layer_item;
        for (auto& layer : layers.content) {
            for (auto& item : layer) {
                const auto& components = layer_item_to_components.at(item.id);
                for (int comp_id : components) {
                    auto iter = component_id_to_layer_item.find(comp_id);
                    if (iter == component_id_to_layer_item.end()) {
                        component_id_to_layer_item[comp_id] = &item;
                    } else {
                        auto& component_layer_item = iter->second;
                        if (!component_layer_item->parent) {
                            component_layer_item->parent = &item;
                        }
                    }
                }
            }
        }
    }
}

/*------------------------------------------------------------------------------------------------*/

ch::brush_token ch::ink_layer_item::token() const {
    std::array<uchar, 8> values;
    values.fill(0);

    const auto* ili_ptr = this;
    while (ili_ptr) {
        values[ili_ptr->layer_id] = ili_ptr->value;
        ili_ptr = ili_ptr->parent;
    }

    ch::brush_token tok = 0;
    for (int i = 0; i < values.size(); ++i) {
        tok |= values[i] << 8 * i;
    }

    ili_ptr = this;
    while (ili_ptr) {
        if (ili_ptr->brush) {
            boost::hash_combine(tok, ili_ptr->brush.get());
        } else {
            boost::hash_combine(tok, 0xffff);
        }
        ili_ptr = ili_ptr->parent;
    }

    return tok;
}

ch::brush_token ch::ink_layer_item::parent_token() const {
    if (parent == nullptr) {
        return 0;
    } else {
        return parent->token();
    }
}

ch::brush_token ch::ink_layer_item::brush_token() const {
    ch::brush_token seed = parent_token();
    ch::brush_token br_tok = (brush) ? reinterpret_cast<ch::brush_token>(brush.get()) : 0xffff;
    boost::hash_combine(seed, br_tok);
    return seed;
}

ch::json ch::ink_layer_item::to_json(
        const std::unordered_map<brush_expr*, std::string>& brush_name_tbl) const {
    return {
        {"id", id},
        {"layer_id", layer_id},
        {"value", value},
        {"poly", polygon_to_json(poly)},
        {"brush", brush_name_tbl.at(brush.get())},
        {"parent", (parent) ? parent->id : -1},
        {"flow_dir", flow_dir}
    };
}

void ch::ink_layer_item::from_json(const json& js,
        const std::unordered_map<std::string, brush_expr_ptr>& brush_name_tbl,
        const std::unordered_map<int, ink_layer_item*>& parent_tbl) {
    id = js["id"].get<int>();
    layer_id = js["layer_id"].get<int>();
    value = js["value"].get<uint8_t>();
    poly = json_to_polygon(js["poly"]);
    brush = brush_name_tbl.at(js["brush"].get<std::string>());
    auto parent_id = js["parent"].get<int>();
    parent = (parent_id >= 0) ? parent_tbl.at(parent_id) : nullptr;
    flow_dir = js["flow_dir"].get<double>();
}

/*------------------------------------------------------------------------------------------------*/

int ch::ink_layer_index(const ink_layer& il) {
    return il.front().layer_id;
}

ch::ink_layers ch::ink_layers::clone() const {
    std::unordered_map<int, ink_layer_item*> id_to_item;
    auto cloned_ink_layers = *this;
    for (auto& layer : cloned_ink_layers.content) {
        for (auto& item : layer) {
            id_to_item[item.id] = &item;
        }
    }
    for (auto& layer : cloned_ink_layers.content) {
        for (auto& item : layer) {
            if (item.parent) {
                auto parent_id = item.parent->id;
                item.parent = id_to_item[parent_id];
            }
        }
    }
    return cloned_ink_layers;
}

int ch::ink_layers::count() const {
    return static_cast<int>(content.size());
}

bool ch::ink_layers::empty() const {
    return content.empty();
}

void ch::ink_layers::clear() {
    return content.clear();
}

ch::json ch::ink_layers::to_json(
        const std::unordered_map<brush_expr*, std::string>& brush_name_tbl) const {

    json sz_json(2, {});
    sz_json[0] = sz.wd;
    sz_json[1] = sz.hgt;

    json ctnt_json(content.size(), {});
    for (const auto& [i, layer] : rv::enumerate(content)) {
        json layer_json(layer.size(), {});
        for (const auto& [j, item] : rv::enumerate(layer)) {
            layer_json[j] = item.to_json(brush_name_tbl);
        }
        ctnt_json[i] = std::move(layer_json);
    }

    return {
        {"sz" , sz_json},
        {"content", ctnt_json}
    };
}

void ch::ink_layers::from_json(const json& js,
        const std::unordered_map<std::string, brush_expr_ptr>& brush_name_tbl) {
    std::unordered_map<int, ink_layer_item*> id_to_item;
    sz = { js["sz"].at(0).get<double>(), js["sz"].at(1).get<double>() };
    const auto& ink_layers_json = js["content"];
    content.clear();
    content.resize(ink_layers_json.size());
    for (const auto& [i, layer_json] : rv::enumerate(ink_layers_json)) {
        content[i].resize(layer_json.size());
        for (const auto& [j, item_json] : rv::enumerate(layer_json)) {
            auto& ink_layer_item = content[i][j];
            ink_layer_item.from_json(item_json, brush_name_tbl, id_to_item);
            id_to_item[ink_layer_item.id] = &ink_layer_item;
        }
    }
}

/*------------------------------------------------------------------------------------------------*/

ch::ink_layers ch::scale(const ink_layers& il, double scale_factor) {
    auto clone = il.clone();
    clone.sz = scale_factor * clone.sz;
    for (auto& layer : clone.content) { 
        for (auto& item : layer) {
            item.poly = scale(item.poly, scale_factor);
        }
    }
    return clone;
}

ch::ink_layers ch::split_into_layers_adaptive(
        const dimensions<double>& sz,
        const std::vector<gray_polygon>& cpolys,
        std::span<const ch::brush_expr_ptr> brush_expr_ptrs,
        std::span<const uchar> gray_levels) {

    polygon_layer_factory layer_factory(cpolys);
    auto adj_graph = layer_factory.graph();
    auto polygons = cpolys | rv::transform([](const auto& p) {return std::get<1>(p); }) | r::to_vector;
    
    int item_id = 0;
    int layer_id = 0;
    auto ranges = gray_ranges(gray_levels, true);
    auto n = ranges.size();
    std::vector<ch::brush_expr_ptr> brushes = (!brush_expr_ptrs.empty()) ?
        brush_expr_ptrs | rv::reverse | r::to_vector :
        std::vector<ch::brush_expr_ptr>(ranges.size(), nullptr);
    std::unordered_map<int, std::vector<int>> layer_item_to_components;
    auto layer_content = ranges |
        rv::transform(
            [&](const auto& rng_tup)->std::vector<ink_layer_item> {
                auto current_layer_id = layer_id++;
                auto [from_gray, to_gray] = rng_tup;
                auto clusters = layer_factory.next_layer(from_gray, to_gray, n);
                return clusters |
                    rv::transform(
                        [&](const auto& ci)->ink_layer_item {
                            item_id++;
                            auto&& [value, poly_ids, poly] = ci;
                            layer_item_to_components[item_id] = std::move(poly_ids);
                            return {
                                item_id,
                                current_layer_id,
                                value,
                                poly,
                                nullptr
                            };
                        }
                ) | r::to_vector;
            }
        ) | r::to_vector;

    ch::ink_layers layers = {
        sz,
        rv::zip(brushes, layer_content) |
            rv::transform(
                [](auto&& pair)->ink_layer {
                    auto&& [br, content] = pair;
                    for (auto& ili : content) {
                        ili.brush = br;
                    }
                    return  content ;
                }
            ) | r::to_vector
    };

    populate_parents(layer_item_to_components, layers);

    return {
        layers.sz,
        layers.content | rv::reverse | r::to_vector
    };
}

std::vector<ch::gray_polygon> merge_gray_polygons(const std::vector<ch::gray_polygon>& polys) {
    auto adj_graph = polygons_to_adjacency_graph(polys);
    std::vector<ch::gray_polygon> output;
    output.reserve(2 * polys.size());
    std::unordered_set<int> visited;

    auto n = static_cast<int>(polys.size());
    for (int i = 0; i < n; ++i) {
        auto color = std::get<0>(polys[i]);
        std::stack<int> stack;
        stack.push(i);
        ch::polygon cluster_poly;
        while (!stack.empty()) {
            auto poly_id = stack.top();
            stack.pop();

            if (visited.contains(poly_id)) {
                continue;
            }
            visited.insert(poly_id);

            const auto& curr_poly = std::get<1>(polys[poly_id]);
            if (cluster_poly.outer().empty()) {
                cluster_poly = curr_poly;
            } else {
                cluster_poly = join_polygons(cluster_poly, curr_poly);
            }

            for (auto neighbor : adj_graph.adjacent_verts(poly_id)) {
                if (std::get<0>(polys[neighbor]) != color) {
                    continue;
                }
                stack.push(neighbor);
            }
        }
        if (cluster_poly.outer().empty()) {
            continue;
        }
        output.push_back({ color, cluster_poly });
    }

    return output;
}

ch::ink_layer_item to_ink_layer_item(const ch::gray_polygon& gp, const ch::brush_expr_ptr& brush_ptr) {
    const auto& [value, poly] = gp;
    return {
        -1,
        -1,
        value,
        poly,
        brush_ptr,
        nullptr,
        0.0
    };
}

ch::ink_layer to_simple_ink_layer(const std::vector<ch::gray_polygon>& polys,
        const std::tuple<uchar, uchar>& rng, ch::brush_expr_ptr br_expr) {

    auto [low, high] = rng;
    auto layer_polys = polys |
        rv::remove_if(
            [low](const auto& gray_poly) {
                auto val = std::get<0>(gray_poly);
                return val < low;
            }
        ) | rv::transform(
            [high](const auto& gray_poly)->ch::gray_polygon {
                const auto [val, p] = gray_poly;
                if (val <= high) {
                    return gray_poly;
                }
                return { high, p };
            }
        ) | r::to_vector;

    layer_polys = merge_gray_polygons(layer_polys);

    return ch::ink_layer{
        layer_polys |
            rv::transform(
                [br_expr](auto& gp)->ch::ink_layer_item {
                    return to_ink_layer_item(gp, br_expr);
                }
            ) | r::to_vector
    };
}

struct ili_to_poly {
    const ch::polygon& operator()(ch::ink_layer_item* p) const {
        return p->poly;
    }
};

using ink_layer_item_tree = ch::polygon_tree<ch::ink_layer_item*, ili_to_poly>;

void populate_ids(std::vector<ch::ink_layer>& layers) {
    int blob_id = 0;
    for (int layer_id = 0; layer_id < static_cast<int>(layers.size()); ++layer_id) {
        for (auto& ili : layers[layer_id]) {
            ili.id = blob_id++;
            ili.layer_id = layer_id;
        }
    }
}

void populate_parents(std::vector<ch::ink_layer>& layers) {
    if (layers.size() <= 1) {
        return;
    }
    for (int i = 1; i < layers.size(); ++i) {
        auto& curr = layers[i];
        auto& prev = layers[i - 1];
        ink_layer_item_tree rtree;
        for (auto& ili : prev) {
            rtree.insert(&ili);
        }
        for (auto& ili : curr) {
            auto polys = rtree.query(ch::representative_point(ili.poly));
            if (polys.size() < 1) {
                // TODO: this should not happen but apparently does...
                continue;
            }
            ili.parent = polys.front();
        }
    }
}

ch::ink_layers ch::split_into_layers_simple(
        const ch::dimensions<double>& sz,
        const std::vector<ch::gray_polygon>& polys,
        std::span<const ch::brush_expr_ptr> brushes,
        std::span<const uchar> gray_levels) {

    auto layers = ink_layers{ sz, {} };
    auto ranges = gray_ranges( gray_levels, false );
    layers.content = rv::zip(ranges, brushes) |
        rv::transform(
            [&](auto&& tup)->ink_layer {
                const auto& [rng, brush] = tup;
                return to_simple_ink_layer(polys, rng, brush);
            }
    ) | r::to_vector;
    populate_ids(layers.content);
    populate_parents(layers.content);
    return layers;

}

std::vector<ch::gray_polygon> ch::to_polygons(const ink_layer& il) {
    return il |
        rv::transform(
            [](const ink_layer_item& ili)->ch::gray_polygon {
                return { ili.value, ili.poly };
            }
    ) | r::to_vector;
}

void ch::debug_layers(cv::Mat img) {
    auto mat = cv::imread("C:\\test\\squares.png");
    auto vec = ch::raster_to_vector(mat, 0.05);
    auto bw = ch::to_monochrome(vec, true);
    auto start = std::chrono::high_resolution_clock::now();
    auto layers = ch::split_into_layers_adaptive(mat_dimensions(mat), bw, {{nullptr,nullptr,nullptr}}, 
        {{ 130, 242, 255}} );

    int n = 0;
    for (const auto& layer : layers.content) {
        qDebug() << "layer: " << n++;
        for (const auto& item : layer) {
            std::string parent = (item.parent) ? 
                std::to_string(item.parent->id) : 
                std::string("nil");
            qDebug() << "  " << item.id << " " << boost::geometry::area(item.poly) / 1600.0
                << " " << parent.c_str();
        }
        qDebug() << "";
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    //qDebug() << duration.count() << " | " << layers[0].content.size();
}