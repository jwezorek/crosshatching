#include "drawing.hpp"
#include "geometry.hpp"
#include "clipper.hpp"
#include "brush_language.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <range/v3/all.hpp>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <array>
#include <boost/functional/hash.hpp>

namespace r = ranges;
namespace rv = ranges::views;
namespace cl = ClipperLib;

namespace {

    struct polygon_with_holes {
        ch::polyline border;
        std::vector<ch::polyline> holes;
    };

    struct ink_plane_blob {
        uchar value;
        uchar bkgd_value;
        polygon_with_holes poly;
    };

    struct ink_plane {
        ch::brush brush;
        std::vector<ink_plane_blob> blobs;
    };

    cv::Rect union_rect_and_pt(const cv::Rect& r, cv::Point2i pt) {
        int x1 = r.x;
        int y1 = r.y;
        int x2 = r.x + r.width - 1;
        int y2 = r.y + r.height - 1;

        x1 = std::min(x1, pt.x);
        y1 = std::min(y1, pt.y);
        x2 = std::max(x2, pt.x);
        y2 = std::max(y2, pt.y);

        return {
            x1, 
            y1,
            x2 - x1 + 1,
            y2 - y1 + 1
        };
    }
    
    int max_val_in_mat(cv::Mat mat) {
        double min_val;
        double max_val;
        cv::Point minLoc;
        cv::Point maxLoc;

        cv::minMaxLoc(mat, &min_val, &max_val, &minLoc, &maxLoc);
        return static_cast<int>(max_val);
    }

    class edge_set {
    private:

        struct hasher {
            size_t operator()(const std::tuple<int, int>& edge) const {
                const auto& [i, j] = edge;
                std::size_t seed = 0;
                boost::hash_combine(seed, i);
                boost::hash_combine(seed, j);
                return seed;
            }
        };

        std::unordered_set<std::tuple<int, int>, hasher> set_;

    public:
        void insert(int i, int j) {
            set_.insert(std::tuple<int, int>(i, j));
        }
        bool contains(int i, int j) {
            return set_.find(std::tuple<int, int>(i, j)) != set_.end();
        }
    };

    class segmentation {
    private:

        struct connected_component {
            int index;
            uchar value;
            cv::Rect bounds;
            cv::Point2i point;
            int area;
            std::vector<connected_component*> neighbors;

            connected_component(int i = -1, cv::Rect b = {}, int a = 0, uchar v = 255, ch::point p = { -1,-1 }) :
                index(i), point(p), area(a), value(v)
            {}
        };

        cv::Mat label_mat_;
        std::vector<connected_component> connected_components_;

        void get_neighbors(int x, int y, std::array<int, 4>& ary) {
            int label = label_mat_.at<int>(y, x);
            std::fill(ary.begin(), ary.end(), -1);
            int i = 0;
            int neigh;
            if ((y > 0) && ((neigh = label_mat_.at<int>(y - 1, x)) != label)) {
                ary[i++] = neigh;
            }
            if ((y < label_mat_.rows - 1) && ((neigh = label_mat_.at<int>(y + 1, x)) != label)) {
                ary[i++] = neigh;
            }
            if ((x > 0) && ((neigh = label_mat_.at<int>(y, x - 1)) != label)) {
                ary[i++] = neigh;
            }
            if ((x < label_mat_.cols - 1) && ((neigh = label_mat_.at<int>(y, x + 1)) != label)) {
                ary[i++] = neigh;
            }
        }

        void insert_edges(int label, int x, int y, edge_set& edges) {
            std::array<int, 4> neighbors;
            get_neighbors(x, y, neighbors);
            for (auto neigh : neighbors) {
                if (neigh == -1) {
                    break;
                }
                if (edges.contains(label, neigh)) {
                    continue;
                }
                edges.insert(label, neigh);
                connected_components_[label].neighbors.push_back(&connected_components_[neigh]);
            }
        }

        void init_new_connected_component(int label, uchar val, int x, int y, edge_set& edges) {
            auto& c = connected_components_[label];
            c.index = label;
            c.value = val;
            c.bounds = cv::Rect(x, y, 1, 1);
            c.area = 1;
            c.point = { x,y };
            insert_edges(label, x, y, edges);
        }

        void update_connected_component(int label, int x, int y, edge_set& edges) {
            auto& c = connected_components_[label];
            c.index = label;
            c.area++;
            c.bounds = union_rect_and_pt(c.bounds, { x,y });
            insert_edges(label, x, y, edges);
        }

        bool is_uninitialized(int index) {
            return connected_components_[index].index == -1;
        }

        

    public:
        segmentation(cv::Mat img, cv::Mat label_image) : 
                label_mat_(label_image) {
            edge_set edges;
            int n = max_val_in_mat(label_image) + 1;
            connected_components_.resize(n);
            for (int y = 0; y < img.rows; y++) {
                for (int x = 0; x < img.cols; x++) {
                    int label = label_image.at<int>(y, x);
                    if (is_uninitialized(label)) {
                        uchar value = img.at<uchar>(y, x);
                        init_new_connected_component(label, value, x, y, edges);
                    } else {
                        update_connected_component(label, x, y, edges);
                    }
                }
            }
        }

        segmentation(cv::Mat mask) {
            cv::Mat stats;
            cv::Mat centroids;
            cv::connectedComponentsWithStats(mask, label_mat_, stats, centroids);

            int n = stats.rows - 1;
            connected_components_.resize(n);
            for (int row = 1; row < n+1; ++row)  {
                int index = row - 1;
                auto& cc = connected_components_[index];
                cc.index = index;
                cc.bounds = cv::Rect(
                    stats.at<int>(cv::Point(0, row)),
                    stats.at<int>(cv::Point(1, row)),
                    stats.at<int>(cv::Point(2, row)),
                    stats.at<int>(cv::Point(3, row))
                );
                cc.area = stats.at<int>(cv::Point(4, row));
            }

            label_mat_.forEach<int>(
                [this,&mask](int& lbl, const int* position) -> void {
                    lbl--;
                    if (lbl >= 0) {
                        int x = position[1];
                        int y = position[0];
                        auto& cc = this->connected_components_[lbl];
                        if (cc.point.x < 0) {
                            cc.point = { x,y };
                            cc.value = mask.at<uchar>(y, x);
                        }
                    }
                }
            );
        }

        uchar value(int index) const {
            return connected_components_[index].value;
        }

        std::vector<int> neighbor_indices(int index) const {
            const auto& neighbors = connected_components_[index].neighbors;
            return neighbors |
                rv::transform(
                    [](const auto* cc)->int {
                        return cc->index;
                    }
                ) |
                r::to_vector;
        }

        cv::Size dimensions() const {
            return { label_mat_.cols, label_mat_.rows };
        }

        int count() const {
            return static_cast<int>(connected_components_.size());
        }

        template<typename T>
        void paint_component(int index, cv::Mat dst, T v) const {
            const auto& c = connected_components_.at(index);

            cv::Mat roi = cv::Mat(label_mat_, c.bounds);

            cv::Mat mask;
            cv::inRange(roi, index, index, mask);

            cv::Mat dst_roi(dst, c.bounds);
            dst_roi.setTo(v, mask);
        }

        void paint_component(int index, cv::Mat dst) const {
            const auto& c = connected_components_.at(index);
            paint_component(index, dst, c.value);
        }

        cv::Point2i component_to_point(int index) const {
            return connected_components_.at(index).point;
        }

        int point_to_component(const cv::Point2i& pt) const {
            return label_mat_.at<int>(pt.y, pt.x);
        }

        struct cc_properties {
            int index;
            uchar color;
            int area;
        };

        cc_properties properties(int index) const {
            const auto& cc = connected_components_.at(index);
            return {
                cc.index,
                cc.value,
                cc.area
            };
        }

        cv::Mat paint_components() const {
            cv::Mat img( dimensions(), CV_8UC1, cv::Scalar(255));
            for (int i = 0; i < count(); ++i) {
                paint_component(i, img);
            }
            return img;
        }
    };

    enum class direction {
        N, NE, E, SE, S, SW, W, NW
    };

    std::vector<std::tuple<uchar, cv::Mat>> gray_planes(const cv::Mat& input) {
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
        auto planes = gray_planes(input);
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
    
    std::vector<polygon_with_holes>  contour_info_to_polygons(const find_contour_output& contours) {
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

    std::string ink_plane_to_svg(const ink_plane& ip, double scale) {
        std::stringstream ss;
        for (const auto& blob : ip.blobs) {
            ss << poly_with_holes_to_svg( static_cast<uchar>((1.0 - blob.value)*255.0), blob.poly, scale) << "\n";
        }
        return ss.str();
    }

    polygon_with_holes scale(const polygon_with_holes& poly, double scale) {
        return {
            ch::scale(poly.border, scale),
            poly.holes | rv::transform([scale](const auto& p) {return ch::scale(p, scale); }) | r::to_vector
        };
    }

    ink_plane_blob scale(const ink_plane_blob& blob, double s) {
        return {
            blob.value,
            blob.bkgd_value,
            scale(blob.poly, s)
        };
    }

    ink_plane scale_ink_plane(const ink_plane& plane, double val) {
        return {
            plane.brush,
            plane.blobs | rv::transform([val](const auto& blob) { return scale(blob, val); }) | r::to_vector
        };
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
                auto clipped = ch::linesegment_rectangle_intersection(  { p1,p2 }, bbox );
                if (clipped) {
                    auto [clipped_p1, clipped_p2] = *clipped;
                    output.push_back({ clipped_p1, clipped_p2 });
                }
            }
        }
        return output;
    }
    //TODO
    std::vector<ch::polyline> clip_swatch_to_poly( ch::crosshatching_swatch swatch, polygon_with_holes& poly, double scale = 100.0) {
        cl::PolyTree polytree;
        cl::Clipper clipper;
        auto bbox = ch::bounding_rectangle(poly.border);
        auto cross_hatching_segments = clip_crosshatching_to_bbox(swatch.content, bbox);
        auto subject = polylines_to_clipper_paths( cross_hatching_segments, scale );
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

    std::vector<std::tuple<uchar, uchar>> ink_plane_ranges(const std::vector<double>& thresholds) {
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

    std::unordered_map<int, std::vector<int>> build_local_id_to_gobal_id_tbl(const segmentation& global_seg, const segmentation& local_seg) {
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

    std::vector<int> get_blob_neighbors_global_indices(int local_blob, uchar from, uchar to, const std::unordered_map<int, std::vector<int>>& local_id_to_global_id, const segmentation& global_seg) {
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

    uchar color_of_greatest_area_comp(const std::vector<segmentation::cc_properties>& props) {
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

    uchar choose_color_for_inherited_blob(const segmentation& global_seg, const std::vector<int>& neighbors) {

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

    std::tuple<cv::Mat, cv::Mat> blobs_in_range_and_mask(const segmentation& seg, uchar from, uchar to) {
        cv::Mat ink_plane_img(seg.dimensions(), CV_8UC1, cv::Scalar(255));
        cv::Mat mask(seg.dimensions(), CV_8UC1, cv::Scalar(0));
        for (int i = 0; i < seg.count(); ++i) {
            auto val = seg.value(i);
            if (val < from || val > to) {
                continue;
            }
            seg.paint_component(i, ink_plane_img);
            seg.paint_component(i, mask, 255);
        }
        return { ink_plane_img, mask };
    }

    std::tuple<cv::Mat, cv::Mat> ink_plane_image_and_mask(const segmentation& seg, cv::Mat prev_level_mask, uchar from, uchar to) {

        auto [ink_plane_img, mask] = blobs_in_range_and_mask(seg, from, to);
        if (prev_level_mask.empty()) {
            return { ink_plane_img, mask };
        }

        if (to == 255) {
            cv::imwrite("C:\\test\\prev_level_mask.png", prev_level_mask);
        }

        segmentation prev_lev_seg(prev_level_mask);

        auto prev_lev_id_to_gobal_id = build_local_id_to_gobal_id_tbl(seg, prev_lev_seg);

        for (int i = 0; i < prev_lev_seg.count(); ++i) {
            prev_lev_seg.paint_component(i, mask, 255);
            auto neighbors = get_blob_neighbors_global_indices(i, from, to, prev_lev_id_to_gobal_id, seg);
            if (neighbors.empty()) {
                continue;
            }
            uchar value = choose_color_for_inherited_blob(seg, neighbors);
            prev_lev_seg.paint_component(i, ink_plane_img, value);
        }

        return { ink_plane_img, mask };
    }

    std::vector<cv::Mat> generate_ink_plane_images(const segmentation& global_segmentation, const std::vector<std::tuple<uchar, uchar>>& ranges_of_gray) {
        std::vector<cv::Mat> images;
        images.reserve(ranges_of_gray.size());
        cv::Mat prev_mask;
        for (const auto& [from_val, to_val] : ranges_of_gray | rv::reverse) {
            cv::Mat ink_plane_img;
            std::tie(ink_plane_img, prev_mask) = ink_plane_image_and_mask(global_segmentation, prev_mask, from_val, to_val);
            images.push_back(ink_plane_img);
        }
        return images | rv::reverse | r::to_vector;
    }

    using blobs_per_gray_t = std::tuple<uchar, std::vector<polygon_with_holes>>;
    ink_plane blobs_per_gray_to_inkplane(const std::vector<blobs_per_gray_t>& blobs_per_gray, ch::brush br, cv::Mat prev_image) {
        std::vector<ink_plane_blob> ink_plane_blobs = rv::join(
                blobs_per_gray |
                rv::remove_if(
                    [](const blobs_per_gray_t& bpg)->bool {
                        return std::get<0>(bpg) == 255;
                    }
                ) |
                rv::transform(
                    [prev_image](const blobs_per_gray_t& bpg) {
                        const auto& [gray, polygons] = bpg;
                        return polygons |
                            rv::transform(
                                [gray, prev_image](const polygon_with_holes& poly)->ink_plane_blob {
                                    uchar bkgd_color = (!prev_image.empty()) ? prev_image.at<uchar>(poly.border.front()) : 255;
                                    return {
                                        gray,
                                        bkgd_color,
                                        poly
                                    };
                                }
                            );
                    }
                )
            ) | r::to_vector;
        return {
            br,
            ink_plane_blobs
        };
    }

    ink_plane ink_plane_img_to_ink_plane(const std::tuple<cv::Mat , ch::brush, cv::Mat>& img_brush_prev) {
        const auto& [img, brush, prev_img] = img_brush_prev;

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

            return blobs_per_gray_to_inkplane(blobs_per_gray, brush, prev_img);
    }

    std::vector<ink_plane> ink_plane_images_to_ink_planes(const std::vector<cv::Mat>& ink_plane_imgs, const std::vector<ch::brush>& brushes) {
        int n = static_cast<int>(ink_plane_imgs.size());
        auto prev_images = rv::concat(rv::single(cv::Mat()), ink_plane_imgs | rv::take(n - 1));
        return rv::zip(ink_plane_imgs, brushes, prev_images) |
            rv::transform(ink_plane_img_to_ink_plane) |
            r::to_vector;
    }

    std::vector<ink_plane> extract_ink_planes_brush_thresh(const cv::Mat& gray_scale_img, const cv::Mat& label_img, const std::vector<ch::brush>& brushes, const std::vector<double>& thresholds) {

        auto seg = segmentation(gray_scale_img, label_img);
        auto ranges = ink_plane_ranges(thresholds);

        auto ink_plane_images = generate_ink_plane_images(seg, ranges);

        for (size_t i = 0; i < ink_plane_images.size(); ++i) {
            std::string name = "C:\\test\\test-ip-" + std::to_string(i + 1) + ".png";
            cv::imwrite(name, ink_plane_images[i]);
        }

        return ink_plane_images_to_ink_planes(ink_plane_images, brushes);
    }
}

std::vector<ink_plane> extract_ink_planes(const cv::Mat& gray_scale_img, const cv::Mat& label_img, const std::vector<std::tuple<ch::brush, double>>& brushes_and_thresholds) {
    auto [brushes, thresholds] = std::tuple{
            brushes_and_thresholds | rv::transform([](const auto& tup) {return std::get<0>(tup); }) | r::to_vector,
            brushes_and_thresholds | rv::transform([](const auto& tup) {return std::get<1>(tup); }) | r::to_vector,
    };
    return extract_ink_planes_brush_thresh(gray_scale_img, label_img, brushes, thresholds);
}

std::vector<ink_plane> scale(const std::vector<ink_plane>& planes, double scale) {
    return planes |
        rv::transform(
            [scale](const ink_plane& plane)->ink_plane {
                return scale_ink_plane(plane, scale);
            }
        ) |
        r::to_vector;
}

void write_to_svg(const std::string& filename, const std::vector<ink_plane>& planes, int wd, int hgt, double scale) {
    std::ofstream outfile(filename);

    outfile << ch::svg_header(static_cast<int>(scale * wd), static_cast<int>(scale * hgt));

    for (const auto& plane : planes)
        outfile << ink_plane_to_svg(plane, scale);

    outfile << "</svg>" << std::endl;
    outfile.close();
}

ch::drawing ch::generate_crosshatched_drawing(cv::Mat img, cv::Mat label_img, double scale, const std::vector<std::tuple<ch::brush, double>>& brushes) {
    auto ink_planes = extract_ink_planes(img, label_img, brushes);
    //TODO
    return {};
}

ch::drawing ch::generate_crosshatched_drawing(cv::Mat img, double scale, const std::vector<std::tuple<ch::brush, double>>& brushes) {
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

cv::Mat grayscale_to_label_image(cv::Mat input) {
    auto values = ch::unique_gray_values(input);
    cv::Mat output(input.size(), CV_32SC1, cv::Scalar(-1));
    int label = 0;
    for (auto val : values) {
        cv::Mat mask;
        cv::inRange(input, val, val, mask);
        segmentation seg(mask);
        for (int i = 0; i < seg.count(); ++i) {
            seg.paint_component(i, output, label++);
        }
    }
    return output;
}

void ch::debug(cv::Mat im, cv::Mat la) {

    auto img = cv::imread("C:\\test\\ink_plane_test.png");
    img = ch::convert_to_gray(img);

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




