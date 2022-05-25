#include "drawing.hpp"
#include "geometry.hpp"
#include "clipper.hpp"
#include "brush_language.hpp"
#include "segmentation.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <range/v3/all.hpp>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <array>

namespace r = ranges;
namespace rv = ranges::views;
namespace cl = ClipperLib;

namespace {

    struct polygon_with_holes {
        ch::polyline border;
        std::vector<ch::polyline> holes;
    };

    struct blob {
        uchar value;
        polygon_with_holes poly;

        blob(uchar v, const polygon_with_holes& p) :
            value(v), poly(p)
        {}

        cv::Point2i point_in_blob() const {
            return poly.border.front();
        }
    };


    class blob_map {
    private:

        cv::Mat labels_;
        std::unordered_map<int, int> label_to_blob_id_;

    public:
        blob_map(cv::Mat ink_layer_img, const std::vector<blob>& blobs) {
            labels_ = ch::grayscale_to_label_image(ink_layer_img);
            for (const auto& [id, blob] : rv::enumerate(blobs)) {
                auto label = labels_.at<int>( blob.point_in_blob() );
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

    std::vector<std::tuple<uchar, cv::Mat>> gray_masks(const cv::Mat& input) {
        auto levels = ch::unique_gray_values(input);
        return
            levels |
            rv::transform(
                [&input](auto val)->std::map<uchar, cv::Mat>::value_type {
                    cv::Mat output;
                    cv::inRange(input, val, val, output);
                    return { val, output };
                }
            ) |
            r::to<std::vector<std::tuple<uchar, cv::Mat>>>();
    }

    using int_polyline = std::vector<cv::Point>;

    struct find_contour_output {
        std::vector<int_polyline> contours;
        std::vector<cv::Vec4i> hierarchy;
    };

    std::vector<std::tuple<uchar, find_contour_output>> perform_find_contours(const cv::Mat& input) {
        auto planes = gray_masks(input);
        return planes |
            rv::transform(
                [](const auto& pair)->std::tuple<uchar, find_contour_output> {
                    const auto& [gray, mat] = pair;
                    find_contour_output output;
                    cv::findContours(mat, output.contours, output.hierarchy, cv::RETR_CCOMP, cv::CHAIN_APPROX_NONE);

                    return { gray, std::move(output) };
                }
            ) |
            r::to<std::vector<std::tuple<uchar, find_contour_output>>>();
    }

    enum class direction {
        N, NE, E, SE, S, SW, W, NW
    };

    struct pixel_crawler {
        cv::Point loc;
        direction head;
    };

    pixel_crawler roll(const pixel_crawler& pc) {
        static std::unordered_map<direction, direction> dir_to_next_dir = {
            {direction::NW, direction::SW},
            {direction::SW, direction::SE},
            {direction::SE, direction::NE},
            {direction::NE, direction::NW}
        };
        return { pc.loc, dir_to_next_dir[pc.head] };
    }

    direction direction_to(const cv::Point& from_pt, const cv::Point& to_pt) {
        static std::unordered_map<cv::Point, direction, ch::point_hasher> offset_to_direction = {
            {{0,-1}, direction::N },
            {{1,-1}, direction::NE},
            {{1,0},  direction::E },
            {{1,1},  direction::SE},
            {{0,1},  direction::S },
            {{-1,1}, direction::SW},
            {{-1,0}, direction::W },
            {{-1,-1},direction::NW}
        };

        return offset_to_direction[ch::normalize_offset(to_pt - from_pt)];
    }

    cv::Point2d get_vertex(const pixel_crawler& pc) {
        static std::unordered_map<direction, cv::Point2d> dir_to_vert_offset = {
            {direction::NW, {0,0} },
            {direction::NE, {1,0}},
            {direction::SE, {1,1}},
            {direction::SW, {0,1}}
        };
        return cv::Point2d(pc.loc) + dir_to_vert_offset[pc.head];
    }

    bool is_cardinal_dir(direction dir) {
        static std::unordered_set<direction> nesw = {
            direction::N,
            direction::E,
            direction::S,
            direction::W
        };
        return nesw.find(dir) != nesw.end();
    }

    bool is_ordinal_dir(direction dir) {
        return !is_cardinal_dir(dir);
    }

    pixel_crawler flip(const pixel_crawler& pc, const cv::Point& to_pt) {
        auto dir = direction_to(pc.loc, to_pt);
        static std::unordered_map<direction, direction> dir_to_flip_dir = {
            {direction::N , direction::S },
            {direction::NE, direction::SW},
            {direction::E , direction::W },
            {direction::SE, direction::NW},
            {direction::S , direction::N },
            {direction::SW, direction::NE},
            {direction::W , direction::E },
            {direction::NW, direction::SE}
        };
        auto flipped_crawler = pixel_crawler{ to_pt, dir_to_flip_dir[pc.head] };
        if (is_ordinal_dir(dir)) {
            return flipped_crawler;
        } else {
            return  roll(flipped_crawler);
        }
    }

    std::tuple<direction, direction> get_shared_vert_directions(const cv::Point& from_pt, const cv::Point& to_pt) {
        static std::unordered_map<direction, std::tuple<direction, direction>> dir_to_shared_verts = {
            {direction::N , {direction::NW, direction::NE}},
            {direction::NW, {direction::NW, direction::NW}},
            {direction::W , {direction::NW, direction::SW} },
            {direction::SW, {direction::SW, direction::SW}},
            {direction::S , {direction::SW, direction::SE} },
            {direction::SE, {direction::SE, direction::SE}},
            {direction::E , {direction::NE, direction::SE} },
            {direction::NE, {direction::NE, direction::NE}}
        };
        return dir_to_shared_verts[direction_to(from_pt, to_pt)];
    }

    bool can_flip(const pixel_crawler& pc, const cv::Point& to_pt) {
        if (pc.loc == to_pt) {
            return false;
        }
        auto [shared_vert_1, shared_vert_2] = get_shared_vert_directions(pc.loc, to_pt);
        return pc.head == shared_vert_1 || pc.head == shared_vert_2;
    }

    int find_northwest_index(const int_polyline& ip) {
        if (ip.size() == 1) {
            return 0;
        }
        int north_west_index = 0;
        for (int i = 1; i < static_cast<int>(ip.size()); ++i) {
            const auto& current_min = ip[north_west_index];
            const auto& pt = ip[i];
            if (pt.y < current_min.y) {
                north_west_index = i;
            } else if (pt.y == current_min.y && pt.x < current_min.x) {
                north_west_index = i;
            }
        }
        return north_west_index;
    }

    std::tuple<pixel_crawler, int> initialize_contour_crawl(const int_polyline& ip, bool counter_clockwise) {
        int northwest_index = find_northwest_index(ip);
        return { {ip[northwest_index], counter_clockwise ? direction::NW : direction::SW}, northwest_index };
    }

    void push_if_new(ch::polyline& poly, const cv::Point2d& pt) {
        if (!poly.empty() && poly.back() == pt) {
            return;
        }
        poly.push_back(pt);
    }

    auto canonicalized_cyclic_contour_view(const int_polyline& ip) {
        int start_index = find_northwest_index(ip);
        auto start_point = ip.at(start_index);
        return ch::rotate_view(rv::all(ip), start_index) | rv::cycle;
    }

    template<typename T>
    auto neighbor_view(T rng) {
        return rng | rv::sliding(2) |
            rv::transform(
                [](auto pair)->std::tuple<cv::Point, cv::Point> {
                    return { pair[0],pair[1] };
                }
        );
    }

    ch::polyline contour_to_polygon(const int_polyline& ip, bool counter_clockwise = true) {
        ch::polyline poly;
        int n = static_cast<int>(ip.size());
        poly.reserve(n);

        auto contour = canonicalized_cyclic_contour_view(ip);
        auto crawler = pixel_crawler{ contour[0], counter_clockwise ? direction::NW : direction::SW };
        auto neighbors = neighbor_view(contour);
        auto iter = neighbors.begin();

        do {
            const auto& [current, next] = *iter;
            auto current_vert = get_vertex(crawler);
            push_if_new(poly, current_vert);
            if (can_flip(crawler, next)) {
                crawler = flip(crawler, next);
                ++iter;
            }
            crawler = roll(crawler);
        } while (get_vertex(crawler) != poly.front());

        poly.shrink_to_fit();
        return ch::simplify_rectilinear_polygon(poly);
    }

    std::vector<polygon_with_holes> contour_info_to_polygons(const find_contour_output& contours) {
        std::unordered_map<int, int> contour_index_to_poly_index;
        std::vector<polygon_with_holes> blobs;
        for (int i = 0; i < contours.hierarchy.size(); ++i) {
            const cv::Vec4i& h = contours.hierarchy[i];
            int parent = h[3];
            const auto& contour = contours.contours[i];
            if (parent == -1) { // contour is an outer polygon
                int poly_index = static_cast<int>(blobs.size());
                blobs.emplace_back(
                    polygon_with_holes{ contour_to_polygon(contour, true), {} }
                );
                contour_index_to_poly_index[i] = poly_index;
            } else { // contour is a hole, so find its parent...
                auto iter = contour_index_to_poly_index.find(parent);
                if (iter == contour_index_to_poly_index.end()) {
                    throw std::runtime_error("contour hierearchy had a child contour before its parent");
                }
                blobs[iter->second].holes.push_back(contour_to_polygon(contour, false));
            }
        }
        return blobs;
    }

    std::string loop_to_path_commands(const ch::polyline& poly, double scale) {
        std::stringstream ss;
        ss << "M " << scale * poly[0].x << "," << scale * poly[0].y << " L";
        for (const auto& pt : rv::tail(poly)) {
            ss << " " << scale * pt.x << "," << scale * pt.y;
        }
        ss << " Z";
        return ss.str();
    }

    std::string svg_path_commands(const polygon_with_holes& poly, double scale) {
        std::stringstream ss;
        ss << loop_to_path_commands(poly.border, scale);
        for (const auto& hole : poly.holes) {
            ss << " " << loop_to_path_commands(hole, scale);
        }
        return ss.str();
    }

    std::string poly_with_holes_to_svg(uchar gray, const polygon_with_holes& poly, double scale) {
        std::stringstream ss;
        ss << "<path fill-rule=\"evenodd\" stroke=\"none\" fill=\"";
        ss << ch::gray_to_svg_color(gray) << "\" d=\"";
        ss << svg_path_commands(poly, scale);
        ss << "\" />";
        return ss.str();
    }

    /*
    std::string ink_plane_to_svg(const ink_plane& ip, double scale) {
        std::stringstream ss;
        for (const auto& blob : ip.blobs) {
            ss << poly_with_holes_to_svg(static_cast<uchar>((1.0 - blob.value) * 255.0), blob.poly, scale) << "\n";
        }
        return ss.str();
    }
    */

    polygon_with_holes scale(const polygon_with_holes& poly, double scale) {
        return {
            ch::scale(poly.border, scale),
            poly.holes | rv::transform([scale](const auto& p) {return ch::scale(p, scale); }) | r::to_vector
        };
    }

    blob scale(const blob& b, double s) {
        return {
            b.value,
            scale(b.poly, s)
        };
    }

    std::vector<blob> scale_blobs(const std::vector<blob>& blobs, double val) {
        return blobs | rv::transform([val](const auto& blob) { return scale(blob, val); }) | r::to_vector;
    }

    polygon_with_holes transform(const polygon_with_holes& p, const ch::matrix& mat) {
        return {
            ch::transform(p.border, mat),
            p.holes | rv::transform([&mat](const auto& hole) {return ch::transform(hole, mat); }) | r::to_vector
        };
    }

    cl::cInt to_cint(double val, double scale) {
        return static_cast<cl::cInt>(std::round(val * scale));
    }

    double from_cint(cl::cInt val, double scale) {
        return val / scale;
    }

    cl::Path polyline_to_clipper_path(const ch::polyline& poly, double scale) {
        return poly |
            rv::transform(
                [scale](const ch::point& pt)->cl::IntPoint {
                    return {
                        to_cint(pt.x, scale),
                        to_cint(pt.y, scale)
                    };
                }
            ) |
            r::to_vector;
    }

    cl::Paths polylines_to_clipper_paths(const std::vector<ch::polyline>& polys, double scale) {
        return polys |
            rv::transform(
                [scale](const auto& poly)->cl::Path {
                    return polyline_to_clipper_path(poly, scale);
                }
            ) |
            r::to_vector;
    }

    cl::Paths poly_with_holes_to_clipper_paths(const polygon_with_holes& poly, double scale) {
        return polylines_to_clipper_paths(
            rv::concat(rv::single(poly.border), poly.holes) | r::to_vector,
            scale
        );
    }

    ch::polyline clipper_path_to_polyline(const cl::Path& path, double scale) {
        return path |
            rv::transform(
                [scale](const cl::IntPoint& pt)->ch::point {
                    return {
                        from_cint(pt.X, scale),
                        from_cint(pt.Y, scale)
                    };
                }
            ) |
            r::to_vector;
    }

    std::vector<ch::polyline> clipper_paths_to_polylines(const cl::Paths& paths, double scale) {
        return paths |
            rv::transform(
                [scale](const auto& path)->ch::polyline {
                    return clipper_path_to_polyline(path, scale);
                }
            ) |
            r::to_vector;
    }

    std::vector<ch::polyline> clip_crosshatching_to_bbox(ch::crosshatching_range swatch, const ch::rectangle& bbox) {
        auto input = swatch | r::to_vector;
        size_t n = 0;
        for (const auto& poly : input) {
            n += poly.size() - 1;
        }
        std::vector<ch::polyline> output;
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

    std::vector<ch::polyline> clip_swatch_to_poly(ch::crosshatching_swatch swatch, polygon_with_holes& poly, double scale = 100.0) {
        cl::PolyTree polytree;
        cl::Clipper clipper;
        auto bbox = ch::bounding_rectangle(poly.border);
        auto cross_hatching_segments = clip_crosshatching_to_bbox(swatch.content, bbox);
        auto subject = polylines_to_clipper_paths(cross_hatching_segments, scale);
        auto clip = poly_with_holes_to_clipper_paths(poly, scale);

        clipper.AddPaths(subject, cl::ptSubject, false);
        clipper.AddPaths(clip, cl::ptClip, true);
        if (!clipper.Execute(cl::ctIntersection, polytree, cl::pftEvenOdd, cl::pftEvenOdd)) {
            throw std::runtime_error("clipper failed.");
        }

        cl::Paths clipped;
        OpenPathsFromPolyTree(polytree, clipped);

        return clipper_paths_to_polylines(clipped, scale);
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

    std::vector<std::tuple<uchar, uchar>> layer_ranges(const std::vector<double>& thresholds) {
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

    std::unordered_map<int, std::vector<int>> build_local_id_to_gobal_id_tbl(const ch::segmentation& global_seg, const ch::segmentation& local_seg) {
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

    std::vector<int> get_blob_neighbors_global_indices(int local_blob, uchar from, uchar to, const std::unordered_map<int, std::vector<int>>& local_id_to_global_id, const ch::segmentation& global_seg) {
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

    uchar color_of_greatest_area_comp(const std::vector<ch::segmentation::cc_properties>& props) {
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

    uchar choose_color_for_inherited_blob(const ch::segmentation& global_seg, const std::vector<int>& neighbors) {

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

    std::tuple<cv::Mat, cv::Mat> blobs_in_range_and_mask(const ch::segmentation& seg, uchar from, uchar to) {
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

    std::tuple<cv::Mat, cv::Mat> ink_layer_image_and_mask(const ch::segmentation& seg, cv::Mat prev_level_mask, uchar from, uchar to) {

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

    std::vector<cv::Mat> generate_ink_layer_images(const ch::segmentation& global_segmentation, const std::vector<std::tuple<uchar, uchar>>& ranges_of_gray) {
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

    using blobs_per_gray_t = std::tuple<uchar, std::vector<polygon_with_holes>>;
    std::vector<blob> blobs_per_gray_to_layer(const std::vector<blobs_per_gray_t>& blobs_per_gray) {
        return rv::join(
            blobs_per_gray |
            rv::remove_if(
                [](const blobs_per_gray_t& bpg)->bool {
                    return std::get<0>(bpg) == 255;
                }
            ) |
            rv::transform(
                [](const blobs_per_gray_t& bpg) {
                    const auto& [gray, polygons] = bpg;
                    return polygons |
                        rv::transform(
                            [gray](const polygon_with_holes& poly)->blob {
                                return {
                                    gray,
                                    poly
                                };
                            }
                    );
                }
            )
          ) | r::to_vector;
    }

    std::vector<blob> ink_layer_img_to_blobs(const cv::Mat& img) {

        auto contours = perform_find_contours(img);
        std::vector<blobs_per_gray_t> blobs_per_gray = contours |
            rv::transform(
                [](const auto& tup)->blobs_per_gray_t {
                    const auto& [gray, contour_info] = tup;
                    auto polygons = contour_info_to_polygons(contour_info);
                    return { gray, std::move(polygons) };
                }
            ) |
            r::to_vector;

         return blobs_per_gray_to_layer(blobs_per_gray);
    }

    std::vector<std::vector<blob>> ink_layer_images_to_blobs(const std::vector<cv::Mat>& ink_layer_imgs) {
        return ink_layer_imgs | rv::transform(ink_layer_img_to_blobs) | r::to_vector;
    }

    std::vector<std::vector<blob>> layers_of_blobs(const cv::Mat& gray_scale_img, const cv::Mat& label_img, const std::vector<double>& thresholds) {

        auto seg = ch::segmentation(gray_scale_img, label_img);
        auto ranges = layer_ranges(thresholds);
        auto ink_layer_images = generate_ink_layer_images(seg, ranges);

        return ink_layer_images_to_blobs(ink_layer_images);
    }

    std::vector<std::vector<blob>> scale(const std::vector<std::vector<blob>>& layers, double scale) {
        return layers |
            rv::transform(
                [scale](const std::vector<blob>& layer)->std::vector<blob> {
                    return scale_blobs(layer, scale);
                }
            ) |
            r::to_vector;
    }
    /*
    void write_to_svg(const std::string& filename, const std::vector<ink_plane>& planes, int wd, int hgt, double scale) {
        std::ofstream outfile(filename);

        outfile << ch::svg_header(static_cast<int>(scale * wd), static_cast<int>(scale * hgt));

        for (const auto& plane : planes)
            outfile << ink_plane_to_svg(plane, scale);

        outfile << "</svg>" << std::endl;
        outfile.close();
    }
    */
}

std::tuple<std::vector<ch::brush>, std::vector<double>> split_brush_thresholds(const std::vector<std::tuple<ch::brush, double>>& brush_thresholds) {
    return {
        brush_thresholds | rv::transform([](const auto& tup) {return std::get<0>(tup); }) | r::to_vector,
        brush_thresholds | rv::transform([](const auto& tup) {return std::get<1>(tup); }) | r::to_vector
    };
}

ch::drawing ch::generate_crosshatched_drawing(cv::Mat img, cv::Mat label_img, double scale, const std::vector<std::tuple<ch::brush, double>>& brush_thresholds) {

    auto [brushes, thresholds] = split_brush_thresholds(brush_thresholds);
    auto blob_layers = layers_of_blobs(img, label_img, thresholds);
    //TODO
    return {};
}

ch::drawing ch::generate_crosshatched_drawing(cv::Mat img, double scale, const std::vector<std::tuple<ch::brush, double>>& brushes) {
    //TODO
    return {};
}

std::vector<ch::polyline> crosshatched_poly_with_holes(const polygon_with_holes& input, double color, ch::brush& brush) {
    auto [x1, y1, x2, y2] = ch::bounding_rectangle(input.border);
    cv::Point2d center = { (x1 + x2) / 2.0,(y1 + y2) / 2.0 };
    auto poly = ::transform(input, ch::translation_matrix(-center));
    auto swatch = brush.get_hatching(color, { x2 - x1,y2 - y1 });
    auto strokes = clip_swatch_to_poly(swatch, poly);
    return ch::transform(strokes, ch::translation_matrix(center));
}

void ch::to_svg(const std::string& filename, const drawing& d) {
    std::ofstream outfile(filename);

    outfile << svg_header(static_cast<int>(d.sz.wd), static_cast<int>(d.sz.hgt));

    for (const auto& poly : d.strokes)
        outfile << polyline_to_svg(poly, d.stroke_wd) << std::endl;

    outfile << "</svg>" << std::endl;
    outfile.close();
}

void ch::debug(cv::Mat im, cv::Mat la) {

    auto img = cv::imread("C:\\test\\ink_plane_test.png");
    img = ch::convert_to_1channel_gray(img);

    auto labels = grayscale_to_label_image(img);
    ch::write_label_map_visualization(labels, "C:\\test\\test_labels.png");

    std::string script_1 = R"(
        (pipe 
            (lin_brush
                (norm_rnd (lerp 0 800) (lerp 50 100))
                (norm_rnd (lerp 200 0) (lerp 20 0.05))
                (norm_rnd (lerp 7 0.5) (lerp 0.5 0.05))
            )
        	(rot 45)
            (jiggle (norm_rnd 0.0 0.02))
            (dis (ramp 0.20 false true))
        )
	)";
    auto result_1 = ch::parse_brush_language(script_1);
    ch::brush brush_1(std::get<ch::brush_fn>(result_1));
    brush_1.build_n(10);

    std::string script_2 = R"(
        (pipe 
            (lin_brush
                (norm_rnd (lerp 0 800) (lerp 50 100))
                (norm_rnd (lerp 200 0) (lerp 20 0.05))
                (norm_rnd (lerp 7 0.5) (lerp 0.5 0.05))
            )
        	(rot 315)
            (jiggle (norm_rnd 0.0 0.02))
            (dis (ramp 0.20 false true))
        )
	)";
    auto result_2 = ch::parse_brush_language(script_2);
    ch::brush brush_2(std::get<ch::brush_fn>(result_2));
    brush_2.build_n(10);

    auto ch_drawing = generate_crosshatched_drawing(img, labels, 4.0, { {brush_1, 0.5},{brush_2,1.0} });
    ch::to_svg("C:\\test\\test_drawing.svg", ch_drawing);
}

/*
void ch::debug(cv::Mat im, cv::Mat la) {

    auto img = cv::imread("C:\\test\\ink_plane_test.png");
    img = ch::convert_to_gray(img);

    auto labels = grayscale_to_label_image(img);
    ch::write_label_map_visualization(labels, "C:\\test\\test_labels.png");

    segmentation seg(img, labels);
    //std::vector<std::tuple<uchar, uchar>> levels = { {0,0}, {1,156}, {157, 254} };
    std::vector<std::tuple<uchar, uchar>> levels = { {0,156}, {157, 254} };
    auto ink_planes = generate_ink_plane_images(seg, levels);
    int i = 0;
    for (const auto& ip : ink_planes) {
        std::string filename = "C:\\test\\ink_plane_" + std::to_string(++i) + ".png";
        cv::imwrite(filename, ip);
    }

}
*/




