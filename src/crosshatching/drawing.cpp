#include "drawing.hpp"
#include "geometry.hpp"
#include "point_set.hpp"
#include "brush_lang.hpp"
#include "segmentation.hpp"
#include "qdebug.h"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <range/v3/all.hpp>
#include <stdexcept>
#include <unordered_set>
#include <array>
#include <chrono>
#include <memory>
#include <fstream>

#include "brush_lang.hpp"
#include "raster_to_vector.hpp"

/*------------------------------------------------------------------------------------------------*/

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    constexpr int k_typical_number_of_strokes = 10000;

    class progress {
    private:
        std::string job_name_;
        std::function<void(double)> prog_fn_;
        std::function<void(const std::string&)> stat_fn_;
        std::function<void(const std::string&)> log_fn_;
        size_t total_blobs;
        size_t curr_blob;
        size_t layer_blob;
        size_t total_layer_blobs;

        static bool completed_new_quantile(size_t running_count, size_t total_count, int quantile) {
            if (running_count == 0) {
                return true;
            }
            auto curr_quantile = (running_count * quantile) / total_count;
            auto prev_quantile = ((running_count - 1) * quantile) / total_count;

            return curr_quantile > prev_quantile;
        }

        static int pcnt_complete(size_t running_count, size_t total_count) {
            return static_cast<int>((100.0 * running_count) / total_count);
        }

    public:
        progress(const std::string& job_name, const ch::callbacks& cbs) :
            job_name_(job_name),
            prog_fn_(cbs.update_progress_cb),
            stat_fn_(cbs.update_status_cb),
            log_fn_(cbs.log_message_cb),
            total_blobs(0),
            curr_blob(0)
        {};

        void set_polygon_count(size_t n) {
            total_blobs = n;
        }

        void start_new_layer(size_t n) {
            layer_blob = 0;
            total_layer_blobs = n;
        }

        void tick() {
            if (total_blobs == 0) {
                return;
            }
            ++curr_blob;
            ++layer_blob;
            if (prog_fn_) {
                prog_fn_(static_cast<double>(curr_blob) / static_cast<double>(total_blobs));
            }
            if (log_fn_) {
                if (completed_new_quantile(layer_blob, total_layer_blobs, 10)) {
                    log(std::string("        ") +
                        std::to_string(pcnt_complete(layer_blob, total_layer_blobs)) +
                        "% complete...");
                }
            }
        }

        void log(const std::string& msg) {
            if (log_fn_) {
                log_fn_(msg);
            }
        }

        void status(const std::string& msg) {
            if (stat_fn_) {
                stat_fn_(msg);
                log(std::string("----") + msg + std::string("----"));
            }
        }
    };

    std::tuple<size_t, size_t> drawing_job_stats(const ch::ink_layers& layers) {
        size_t poly_count = r::accumulate(
            layers.content | rv::transform([](auto&& l) {return l.content.size(); }),
            size_t{ 0 }
        );
        size_t vert_count = r::accumulate(
            layers.content | rv::transform(
                [](auto&& l) {
                    return r::accumulate(
                        l.content | rv::transform(
                            [](auto&& ili) {return ch::vert_count(ili.poly); }
                        ),
                        size_t{ 0 }
                    );
                }
            ),
            size_t{ 0 }
        );
        return { poly_count, vert_count };
    }

    using swatch_table = std::unordered_map<ch::brush_token, ch::bkgd_swatches>;
    swatch_table create_swatch_table() {
        swatch_table tbl;
        tbl[0] = {};
        return tbl;
    }

    std::tuple<std::vector<ch::drawn_stroke_cluster>, swatch_table> draw_ink_layer(
            const ch::ink_layer& layer, swatch_table& tbl,
            const ch::parameters& params, progress& prog) {

        size_t n = layer.content.size();
        prog.start_new_layer(n);

        std::vector<ch::drawn_stroke_cluster> output;
        output.reserve(k_typical_number_of_strokes);
        std::unordered_map<ch::brush_token, ch::brush_ptr> brush_table;
        swatch_table output_table = create_swatch_table();

        size_t count = 0;
        for (const auto& blob: layer.content) {
            auto parent_tok = blob.parent_token();
            auto iter = brush_table.find(parent_tok);
            ch::brush_ptr current_brush;
            if (iter != brush_table.end()) {
                current_brush = iter->second;
            } else {
                current_brush = std::make_shared<ch::brush>(
                    layer.brush,
                    params.epsilon,
                    ch::dimensions(static_cast<double>(params.swatch_sz)),
                    tbl.at(parent_tok)
                );
                current_brush->build_n(20);
                brush_table[parent_tok] = current_brush;
            }

            auto value = static_cast<double>(blob.value) / 255.0;
            auto strokes = current_brush->draw_strokes(blob.poly, value);
            std::copy(strokes.begin(), strokes.end(), std::back_inserter(output));

            auto tok = blob.token();
            if (output_table.find(tok) == output_table.end()) {
                output_table[tok] = current_brush->render_swatches(value);
            }
            prog.tick();
        }
        output.shrink_to_fit();

        return { std::move(output), std::move(output_table) };
    }

    ch::drawing draw(const ch::ink_layers& inp_layers, const ch::parameters& params, progress& prog) {

        auto [poly_count, vert_count] = drawing_job_stats(inp_layers);
        prog.log("  (job contains " + std::to_string(poly_count) + " polygons with " +
            std::to_string(vert_count) + " total vertices");
        prog.set_polygon_count(poly_count);

        auto layers = scale(inp_layers, params.scale);
        
        std::vector<std::vector<ch::drawn_stroke_cluster>> layer_strokes(layers.content.size());
        swatch_table tok_to_bkgd = create_swatch_table();
        for (auto [index, layer] : rv::enumerate(layers.content)) {
            prog.log(std::string("  - layer ") + std::to_string(index));
            std::tie(layer_strokes[index], tok_to_bkgd) = draw_ink_layer(
                layer, tok_to_bkgd, params, prog);
        }

        return ch::drawing{
            layer_strokes | rv::join | r::to_vector,
            layers.sz
        };
    }
}

ch::drawing ch::generate_crosshatched_drawing(const ch::crosshatching_job& job, const callbacks& cbs) {
    progress prog(job.title, cbs);

    auto start_time = std::chrono::high_resolution_clock().now();

    prog.status(std::string("drawing ") + job.title);
    auto result = draw(job.layers, job.params, prog);
    prog.status(std::string("complete.")); // TODO: or error

    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock().now() - start_time;
    prog.log(std::string("( ") + std::to_string(elapsed.count()) + " seconds)");

    return result;
}

void ch::write_to_svg(const std::string& filename, const drawing& d,
        std::function<void(double)> update_progress_cb) {

}

ch::drawing ch::scale(const drawing& d, double scale) {
    return {};
}

cv::Mat ch::paint_drawing(const ch::drawing& d, std::function<void(double)> update_progress_cb) {
    return {};
}

void ch::debug_drawing() {

}


    /*
    ch::polylines clip_crosshatching_to_bbox(
            ch::crosshatching_range swatch, const ch::rectangle& bbox) {

        auto input = swatch | r::to_vector;
        size_t n = 0;
        for (const auto& poly : input) {
            n += poly.size() - 1;
        }
        ch::polylines output;
        output.reserve(n);
        for (const auto& poly : input) {
            for (auto rng : poly | rv::sliding(2)) {
                auto p1 = rng[0];
                auto p2 = rng[1];
                auto clipped = ch::linesegment_rectangle_intersection({ p1,p2 }, bbox);
                if (clipped) {
                    auto [clipped_p1, clipped_p2] = *clipped;
                    output.push_back({ clipped_p1, clipped_p2 });
                }
            }
        }
        return output;
    }

    ch::polylines clip_swatch_to_poly(
            const ch::crosshatching_swatch& swatch, const ch::polygon& poly) {

        auto bbox = ch::bounding_rectangle(poly.outer());
        auto cross_hatching_segments = clip_crosshatching_to_bbox(swatch.content, bbox);
        return ch::clip_lines_to_poly(cross_hatching_segments, poly);
    }

    ch::polyline int_polyline_to_polyline(const int_polyline& ip) {
        return ip |
            rv::transform(
                [](const cv::Point p)->cv::Point2d {
                    return {
                        static_cast<double>(p.x),
                        static_cast<double>(p.y)
                    };
                }
            ) |
            r::to<ch::polyline>();
    }

    std::vector<std::tuple<uchar, uchar>> layer_ranges(std::span<const double> thresholds) {
        auto start_of_ranges = thresholds |
            rv::transform(
                [](double val)->uchar {
                    uchar v = 255 - static_cast<uchar>(val * 255.0);
                    return v;
                }
            ) |
            r::to_vector;
                auto n = thresholds.size();
                auto end_of_ranges = rv::concat(rv::single(254),
                    start_of_ranges |
                    rv::take(n - 1) |
                    rv::transform(
                        [](uchar start_of_next)->uchar {
                            return start_of_next - 1;
                        }
                    )
                ) |
                r::to_vector;
                    return rv::zip(start_of_ranges, end_of_ranges) |
                        rv::transform(
                            [](std::pair<uchar, uchar> p)->std::tuple<uchar, uchar> {
                                return { p.first, p.second };
                            }
                        ) |
                r::to_vector;
    }

    std::unordered_map<int, std::vector<int>> build_local_id_to_gobal_id_tbl(
            const ch::segmentation& global_seg, const ch::segmentation& local_seg) {
        std::unordered_map<int, std::vector<int>> local_id_to_global_id;
        for (int global_index = 0; global_index < global_seg.count(); ++global_index) {
            auto pt_in_global_seg = global_seg.component_to_point(global_index);
            auto local_index = local_seg.point_to_component(pt_in_global_seg);
            if (local_index < 0) {
                continue;
            }
            auto iter = local_id_to_global_id.find(local_index);
            if (iter != local_id_to_global_id.end()) {
                iter->second.push_back(global_index);
            } else {
                local_id_to_global_id[local_index] = { global_index };
            }
        }
        return local_id_to_global_id;
    }

    std::vector<int> get_blob_neighbors_global_indices(
            int local_blob, uchar from, uchar to, 
            const std::unordered_map<int, std::vector<int>>& local_id_to_global_id, 
            const ch::segmentation& global_seg) {

        const auto& global_blob_components = local_id_to_global_id.at(local_blob);
        std::unordered_set<int> neighbors;
        for (auto gbc : global_blob_components) {
            const auto& neighbor_indices = global_seg.neighbor_indices(gbc);
            for (auto neigh_index : neighbor_indices) {
                auto neigh_color = global_seg.value(neigh_index);
                if (neigh_color >= from && neigh_color <= to) {
                    neighbors.insert(neigh_index);
                }
            }
        }
        return neighbors | r::to_vector;
    }

    uchar color_of_greatest_area_comp(std::span<ch::segmentation::cc_properties> props) {
        uchar color_of_max;
        int max_area = -1;
        for (const auto& p : props) {
            if (p.area > max_area) {
                max_area = p.area;
                color_of_max = p.color;
            }
        }
        return color_of_max;
    }

    uchar choose_color_for_inherited_blob(
            const ch::segmentation& global_seg, std::span<int> neighbors) {

        //if there is only one neighbor, use its color...
        if (neighbors.size() == 1) {
            return global_seg.value(neighbors.front());
        }

        //strore all the properties in a vector...
        auto neighbor_properties = neighbors |
            rv::transform(
                [&global_seg](int index) {
                    return global_seg.properties(index);
                }
        ) | r::to_vector;

        // make a table mapping shades of gray to total area of neighbors of that color...
        std::map<int, int> color_to_area;
        for (const auto& props : neighbor_properties) {
            color_to_area[props.color] += props.area;
        }

        // choose the shade of gray with greatest area...
        uchar color_with_max_area = 0;
        int max_area = -1;
        for (auto [color, area] : color_to_area) {
            if (area > max_area) {
                max_area = area;
                color_with_max_area = color;
            }
        }

        return color_with_max_area;
    }

    std::tuple<cv::Mat, cv::Mat> blobs_in_range_and_mask(
            const ch::segmentation& seg, uchar from, uchar to) {
        cv::Mat ink_layer_img(seg.dimensions(), CV_8UC1, cv::Scalar(255));
        cv::Mat mask(seg.dimensions(), CV_8UC1, cv::Scalar(0));
        for (int i = 0; i < seg.count(); ++i) {
            auto val = seg.value(i);
            if (val < from || val > to) {
                continue;
            }
            seg.paint_component(i, ink_layer_img);
            seg.paint_component(i, mask, 255);
        }
        return { ink_layer_img, mask };
    }

    std::tuple<cv::Mat, cv::Mat> ink_layer_image_and_mask(
            const ch::segmentation& seg, cv::Mat prev_level_mask, uchar from, uchar to) {

        auto [ink_layer_img, mask] = blobs_in_range_and_mask(seg, from, to);
        if (prev_level_mask.empty()) {
            return { ink_layer_img, mask };
        }

        ch::segmentation prev_lev_seg(prev_level_mask);

        auto prev_lev_id_to_gobal_id = build_local_id_to_gobal_id_tbl(seg, prev_lev_seg);

        for (int i = 0; i < prev_lev_seg.count(); ++i) {
            prev_lev_seg.paint_component(i, mask, 255);
            auto neighbors = get_blob_neighbors_global_indices(i, from, to, prev_lev_id_to_gobal_id, seg);
            if (neighbors.empty()) {
                continue;
            }
            uchar value = choose_color_for_inherited_blob(seg, neighbors);
            prev_lev_seg.paint_component(i, ink_layer_img, value);
        }

        return { ink_layer_img, mask };
    }

    std::vector<cv::Mat> generate_ink_layer_images(
            const ch::segmentation& global_segmentation, 
            std::span<std::tuple<uchar, uchar>> ranges_of_gray) {
        std::vector<cv::Mat> images;
        images.reserve(ranges_of_gray.size());
        cv::Mat prev_mask;
        for (const auto& [from_val, to_val] : ranges_of_gray | rv::reverse) {
            cv::Mat ink_layer_img;
            std::tie(ink_layer_img, prev_mask) = ink_layer_image_and_mask(global_segmentation, prev_mask, from_val, to_val);
            images.push_back(ink_layer_img);
        }
        return images | rv::reverse | r::to_vector;
    }

    template<typename T>
    using blobs_per_gray_t = std::tuple<T, std::vector<ch::polygon>>;

    template<typename T>
    std::vector<std::tuple<T, ch::polygon>> blobs_per_gray_to_layer(
            const std::vector< blobs_per_gray_t<T>>& blobs_per_gray) {
        using blob = std::tuple< T, ch::polygon>;
        return rv::join(
            blobs_per_gray |
            rv::transform(
                [](const auto& bpg) {
                    const auto& [value, polygons] = bpg;
                    return polygons |
                        rv::transform(
                            [value](const auto& poly)->blob {
                                return {
                                    value,
                                    poly
                                };
                            }
                    );
                }
            )
          ) | r::to_vector;
    }

    template<typename T>
    std::vector<std::tuple<T, ch::polygon>> ink_layer_img_to_blobs(const cv::Mat& img) {

        auto contours = perform_find_contours<T>(img);
        std::vector<blobs_per_gray_t<T>> blobs_per_gray = contours |
            rv::transform(
                [](const auto& tup)->blobs_per_gray_t<T> {
                    const auto& [value, contour_info] = tup;
                    auto polygons = contour_info_to_polygons(contour_info);
                    return { value, std::move(polygons) };
                }
            ) |
            r::to_vector;

         return blobs_per_gray_to_layer(blobs_per_gray);
    }

    std::vector<std::vector<std::tuple<uchar, ch::polygon>>> ink_layer_images_to_blobs(
            std::span<const cv::Mat> ink_layer_imgs) {
        return ink_layer_imgs | rv::transform(ink_layer_img_to_blobs<uchar>) | r::to_vector;
    }

    std::vector<cv::Mat> ink_layer_images(const cv::Mat& gray_scale_img, const cv::Mat& label_img, const std::vector<double>& thresholds) {

        auto seg = ch::segmentation(gray_scale_img, label_img);
        auto ranges = layer_ranges(thresholds);
        return generate_ink_layer_images(seg, ranges);
    }

    class ink_layer_stack {

    public:  
        using brush_token = uint64_t;

    private:

        struct blob_properties {
            uchar value;
            ch::polygon poly;
            int layer_id;
            int blob_id;
            blob_properties* parent;

            size_t vert_count() const {
                size_t n = poly.outer().size();
                for (const auto& hole : poly.inners()) {
                    n += hole.size();
                }
                return n;
            }

            cv::Point2i point_in_blob() const {
                return poly.outer().front();
            }

            brush_token tokenize() const {
                std::vector<uchar> values(layer_id + 2, 255);
                const auto* blob_ptr = this;
                while (blob_ptr) {
                    values[blob_ptr->layer_id] = blob_ptr->value;
                    blob_ptr = blob_ptr->parent;
                }
                return make_brush_token(values);
            }

            brush_token parent_token() const {
                return (parent) ?
                    parent->tokenize() :
                    255;
            }

            bool is_white() const {
                return value == 255;
            }

        };

        struct ink_layer {
            ch::brush_fn brush;
            std::vector<blob_properties> blobs;

            void scale(double scale_factor) {
                for (auto& poly : blobs | 
                        rv::transform([](auto& b)->ch::polygon& {return b.poly;})){
                    poly = ch::scale(poly, scale_factor);
                }
            }

            void simplify(double param) {
                auto polygons = blobs | 
                    rv::transform( 
                        []( auto&& blob) {
                            return std::move(blob.poly);
                        } 
                    ) | r::to_vector;
                polygons = ch::simplify_polygons(polygons, param);
                for (size_t i = 0; i < polygons.size(); ++i) {
                    blobs[i].poly = std::move(polygons[i]);
                }
            }
        };

        class blob_map {
        private:

            cv::Mat labels_;
            std::unordered_map<int, int> label_to_blob_id_;

        public:
            blob_map(cv::Mat ink_layer_img, const std::vector<blob_properties>& blobs) {
                labels_ = ch::grayscale_to_label_image(ink_layer_img);
                for (const auto& [id, blob] : rv::enumerate(blobs)) {
                    auto label = labels_.at<int>(blob.point_in_blob());
                    if (label < 0) {
                        throw std::runtime_error("bad ink layer map access.");
                    }
                    label_to_blob_id_[label] = id;
                }
            }

            std::optional<int> id_from_point(const cv::Point2i& pt) const {
                auto label = labels_.at<int>(pt);
                auto iter = label_to_blob_id_.find(label);
                if (iter != label_to_blob_id_.end()) {
                    return iter->second;
                } else {
                    return {};
                }
            }
        };

        std::vector<ink_layer> layers_;

        static ink_layer to_ink_layer(
                const std::vector<std::tuple<uchar, ch::polygon>>& value_poly_pairs,
                int layer_index, ch::brush_fn brush) {
            ink_layer il{
                brush,
                rv::enumerate(value_poly_pairs) |
                    rv::transform([&](const auto& tup)->blob_properties {
                        const auto& [index, value_poly_pair] = tup;
                        int blob_index = static_cast<int>(index);
                        return {
                            std::get<0>(value_poly_pair),
                            std::get<1>(value_poly_pair),
                            layer_index,
                            blob_index,
                            nullptr
                        };
                    }) | r::to_vector
            };
            return il;
        }

        std::tuple<int, int> find_blob_parent(const blob_properties& blob, const std::vector<blob_map>& blob_maps) {
            for (int parent_layer_id = blob.layer_id - 1; parent_layer_id >= 0; --parent_layer_id) {
                auto parent_blob_id = blob_maps.at(parent_layer_id).id_from_point(blob.point_in_blob());
                if (parent_blob_id) {
                    return { parent_layer_id, parent_blob_id.value() };
                }
            }
            return { -1,-1 };
        }

        void populate_blob_parents( std::span<const cv::Mat> layer_images) {
            auto blob_maps = rv::zip(layer_images, layers_) |
                rv::transform(
                    [](const auto& tup) {
                        const auto& [mat, layer] = tup;
                        return blob_map(mat, layer.blobs);
                    }
            ) | r::to_vector;

            for (int layer_index = 0; layer_index < static_cast<int>(layers_.size()); ++layer_index) {
                auto& layer = layers_[layer_index];
                for (auto& blob : layer.blobs) {
                    auto [parent_layer, parent_blob] = find_blob_parent(blob, blob_maps);
                    blob.parent = (parent_blob != -1) ?
                        &(layers_[parent_layer].blobs[parent_blob]) :
                        nullptr;
                }
            }
        }

        static brush_token make_brush_token(const std::vector<uchar>& values) {
            if (values.size() > 8) {
                throw std::runtime_error("too many ink layers");
            }
            brush_token tok = 0;
            for (int i = 0; i < values.size(); ++i) {
                tok |= values[i] << 8 * i;
            }

            return tok;
        }

    public:

        ink_layer_stack(
                std::span<const cv::Mat> layers, std::vector<ch::brush_fn> brushes, 
                double scale_factor, double simplify_param) {
            auto blobs = ink_layer_images_to_blobs(layers);
            layers_ = rv::enumerate(blobs) |
                rv::transform(
                    [&brushes](const auto& tup) {
                        const auto& [layer_index, blobs_per_layer] = tup;
                        return to_ink_layer(blobs_per_layer, layer_index, brushes.at(layer_index));
                    }
                ) | r::to_vector;
            populate_blob_parents(layers);
            qDebug() << simplify_param;
            simplify_polygons(simplify_param);
            scale(scale_factor);
        }

        using brush_blob = std::tuple<ch::brush_fn, blob_properties>;
        using blob_range = r::any_view<brush_blob>;

        std::tuple<size_t,size_t> blob_stats() const {
            size_t blob_count = 0;
            size_t vert_count = 0;
            for (const auto& layer : layers_) {
                blob_count += layer.blobs.size();
                for (const auto& blob : layer.blobs) {
                    vert_count += blob.vert_count();
                }
            }
            return { blob_count, vert_count };
        }

        r::any_view<blob_range> blobs() const {
            return  layers_ |
                rv::transform(
                    [](const auto& layer) {
                        const auto& brush = layer.brush;
                        return layer.blobs |
                            rv::transform(
                                [&brush](const blob_properties& blob)->brush_blob {
                                    return { brush, blob };
                                }
                        );
                    }
                );
        }

        void scale(double scale_factor) {
            for (auto& layer : layers_) {
                layer.scale(scale_factor);
            }
        }

        void simplify_polygons(double param) {
            if (param <= 1.0) {
                return;
            }
            for (auto& layer : layers_) {
                layer.simplify(param);
            }
        }
       
    };

    std::tuple<std::vector<ch::brush_fn>, std::vector<double>> split_brush_thresholds(
            const std::vector<std::tuple<ch::brush_fn, double>>& brush_thresholds) {
        return {
            brush_thresholds | 
                rv::transform([](const auto& tup) {return std::get<0>(tup); }) | r::to_vector,
            brush_thresholds | 
                rv::transform([](const auto& tup) {return std::get<1>(tup); }) | r::to_vector
        };
    }

    std::vector<ch::polyline> crosshatch_polygon(
            const ch::polygon& input, double color, ch::brush_ptr brush) {
        auto [x1, y1, x2, y2] = ch::bounding_rectangle(input.outer());
        cv::Point2d center = { (x1 + x2) / 2.0,(y1 + y2) / 2.0 };
        auto poly = ch::transform(input, ch::translation_matrix(-center));
        auto swatch = brush->get_hatching(color, { x2 - x1,y2 - y1 });
        auto strokes = clip_swatch_to_poly(swatch, poly);
        return ch::transform(strokes, ch::translation_matrix(center));
    }

    using swatch_table = std::unordered_map<ink_layer_stack::brush_token, ch::bkgd_swatches>;

    swatch_table create_swatch_table() {
        swatch_table tbl;
        tbl[255] = {};
        return tbl;
    }

    double gray_byte_to_value(uchar gray) {
        return 1.0 - static_cast<double>(gray) / 255.0;
    }

    std::tuple<std::vector<ch::polyline>, swatch_table> paint_ink_layer(
            ink_layer_stack::blob_range layer, const swatch_table& tbl, 
            const ch::parameters& params, progress_logger& prog) {
        size_t n = r::distance(layer);
        prog.start_new_layer(n);

        std::vector<ch::polyline> output;
        output.reserve(k_typical_number_of_strokes);
        std::unordered_map<ink_layer_stack::brush_token, ch::brush_ptr> brush_table;
        swatch_table output_table = create_swatch_table();

        size_t count = 0;
        for (const auto& [brush_func, blob] : layer) {
            auto parent_tok = blob.parent_token();
            auto iter = brush_table.find(parent_tok);
            ch::brush_ptr current_brush;
            if (iter != brush_table.end()) {
                current_brush = iter->second;
            } else {
                current_brush = std::make_shared<ch::old_brush>(
                    brush_func, 
                    params.stroke_width, 
                    params.epsilon,
                    ch::dimensions(static_cast<double>(params.swatch_sz)), 
                    tbl.at(parent_tok)
                );
                current_brush->build_n(20);
                brush_table[parent_tok] = current_brush;
            }

            auto value = gray_byte_to_value(blob.value);
            auto strokes = crosshatch_polygon( blob.poly, value, current_brush);
            std::copy(strokes.begin(), strokes.end(), std::back_inserter(output));

            auto tok = blob.tokenize();
            if (output_table.find(tok) == output_table.end()) {
                output_table[tok] = current_brush->render_swatches(value);
            }
            prog.tick();
        }
        output.shrink_to_fit();

        return { std::move(output), std::move(output_table) };
    }

    void sanitize_input_images(cv::Mat* img, cv::Mat*label_img) {
        if (img->channels() != 1) {
            *img = ch::convert_to_1channel_gray(*img);
        }

        if (label_img->empty()) {
            *label_img = ch::grayscale_to_label_image(*img);
        }

        if (label_img->type() != CV_32SC1) {
            throw std::runtime_error("bad label image");
        }
    }

    std::tuple<std::vector<ch::brush_fn>, std::vector<cv::Mat>> split_layers_and_brushes(
            std::span<const std::tuple<ch::brush_fn, cv::Mat>> brushes_and_imgs) {
        return {
            brushes_and_imgs | 
                rv::transform([](const auto& p) {return std::get<0>(p); }) | r::to_vector,
            brushes_and_imgs | 
                rv::transform([](const auto& p) {return std::get<1>(p); }) | r::to_vector
        };
    }

    ch::drawing draw_layers( std::span<const std::tuple<ch::brush_fn, cv::Mat>> layers_and_brushes, 
            const ch::parameters& params, progress_logger& ps) {

        auto [brushes, layers] = split_layers_and_brushes(layers_and_brushes);

        ps.log("  constructing ink layer stack...");
        ink_layer_stack stack(layers, brushes, params.scale, params.polygon_simplification_param);

        auto [blob_count, vert_count]  = stack.blob_stats();
        ps.log("  (job contains " + std::to_string(blob_count) + " polygons with " + 
            std::to_string(vert_count) + " total vertices");
        ps.set_total_blobs(blob_count);

        std::vector<std::vector<ch::polyline>> layer_strokes(layers.size());
        swatch_table tok_to_bkgd = create_swatch_table();

        for (auto [index, layer] : rv::enumerate(stack.blobs())) {
            ps.log(std::string("  - layer ") + std::to_string(index));
            std::tie(layer_strokes[index], tok_to_bkgd) = paint_ink_layer(
                    layer,  tok_to_bkgd, params, ps);
        }

        return ch::drawing{
            layer_strokes | rv::join | r::to_vector,
             params.scale * ch::dimensions<double>(ch::mat_dimensions(layers.front())),
            static_cast<double>(params.stroke_width)
        };
    }
}

std::vector<std::tuple<ch::brush_fn, cv::Mat>> ch::generate_ink_layer_images(cv::Mat img, 
        cv::Mat label_img,  
        const std::vector<std::tuple<ch::brush_fn, double>>& brush_intervals) {
    sanitize_input_images(&img, &label_img);
    auto [brushes, thresholds] = split_brush_thresholds(brush_intervals);
    auto layer_images = ink_layer_images(img, label_img, thresholds);
    return rv::zip(brushes, layer_images) | 
        r::to<std::vector<std::tuple<ch::brush_fn,cv::Mat>>>();
}

ch::drawing ch::generate_crosshatched_drawing(const ch::crosshatching_job& job, 
    const ch::callbacks& cbs) {
    
    progress_logger ps(job.title, cbs);

    auto start_time = std::chrono::high_resolution_clock().now();

    ps.status(std::string("drawing ") + job.title);
    auto result = draw_layers(job.layers, job.params, ps);
    ps.status(std::string("complete.")); // TODO: or error

    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock().now() - start_time;
    ps.log(std::string("( ") + std::to_string(elapsed.count()) + " seconds)");

    return result;
}

void ch::write_to_svg(const std::string& filename, const drawing& d, 
        std::function<void(double)> update_progress) {
    std::ofstream outfile(filename);

    outfile << svg_header(static_cast<int>(d.sz.wd), static_cast<int>(d.sz.hgt));
    
    auto n = static_cast<double>(d.strokes.size());
    for (const auto& [index, poly] : rv::enumerate(d.strokes)) {
        outfile << polyline_to_svg(poly, d.stroke_wd) << std::endl;
        if (update_progress) {
            update_progress(index / n);
        }
    }
    outfile << "</svg>" << std::endl;
    outfile.close();
}

ch::drawing ch::scale(const ch::drawing& d, double val) {
    return {
        d.strokes |
            rv::transform([val](const ch::polyline& poly) { return ch::scale(poly, val); }) |
            r::to_vector,
        {val * d.sz.wd, val * d.sz.hgt},
        val * d.stroke_wd
    };
}

cv::Mat ch::paint_drawing(const drawing& d, std::function<void(double)> update_progress_cb) {
    cv::Mat mat(static_cast<int>(d.sz.hgt), static_cast<int>(d.sz.wd), CV_8U, 255);
    auto n = static_cast<double>(d.strokes.size());
    for (const auto& [index, ls] : rv::enumerate(d.strokes)) {
        ch::paint_polyline_aa(mat, ls, d.stroke_wd, 0, point{ 0, 0 });
        if (update_progress_cb) {
            update_progress_cb(index / n);
        }
    }
    return mat;
}

std::vector<std::tuple<uchar,ch::polygon>> ch::detail::to_blobs_from_1channel_image(
    const cv::Mat& input) {
    return ink_layer_img_to_blobs<uchar>(input);
}

std::vector<std::tuple<ch::color, ch::polygon>> ch::detail::to_blobs_from_3channel_image(
    const cv::Mat& input) {
    return ink_layer_img_to_blobs<ch::color>(input);
}

void ch::raster_to_vector(const std::string& fname, const cv::Mat& mat, double scale, double param) {
    auto blobs = to_blobs<color>(mat);
    blobs = ch::simplify_colored_polygons<ch::color>(blobs, param);
    polygons_to_svg<ch::color>(fname, blobs, scale);
}


void ch::debug_drawing(const cv::Mat& mat) {
    debug_brushes();
}

*/