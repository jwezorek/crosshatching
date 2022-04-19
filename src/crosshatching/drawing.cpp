#include "drawing.hpp"
#include "geometry.hpp"
#include "clipper.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <range/v3/all.hpp>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <unordered_set>

namespace r = ranges;
namespace rv = ranges::views;
namespace cl = ClipperLib;

namespace {

    enum class direction {
        N, NE, E, SE, S, SW, W, NW
    };

    std::vector<uchar> get_gray_values(const cv::Mat& input) {
        if (input.channels() != 1) {
            throw std::runtime_error("called get_gray_levels on color image");
        }
        std::array<bool, 256> grays = {};
        input.forEach<uchar>([&grays](uchar gray, const int* pos) { grays[gray] = true;  });
        return rv::iota(0) |
            rv::take(256) |
            rv::filter([&grays](int g) {return grays[g]; }) |
            r::to<std::vector<uchar>>();
    }

    std::vector<std::tuple<uchar, cv::Mat>> gray_planes(const cv::Mat& input, bool full_range = false) {
        auto levels = get_gray_values(input);
        return
            levels |
            rv::transform(
                [&input,full_range](auto val)->std::map<uchar, cv::Mat>::value_type {
                    cv::Mat output;
                    cv::inRange(input, !full_range ? val : 0, val, output);
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
        return poly;
    }

    ch::gray_level_plane contour_info_to_gray_level_plane(uchar gray, const find_contour_output& contours) {
        ch::gray_level_plane glp = { gray,{} };
        std::unordered_map<int, int> contour_index_to_poly_index;

        for (int i = 0; i < contours.hierarchy.size(); ++i) {
            const cv::Vec4i& h = contours.hierarchy[i];
            int parent = h[3];
            const auto& contour = contours.contours[i];
            if (parent == -1) {
                int poly_index = static_cast<int>(glp.blobs.size());
                glp.blobs.emplace_back(
                    ch::polygon_with_holes{ contour_to_polygon(contour, true), {} }
                );
                contour_index_to_poly_index[i] = poly_index;
            } else {
                auto iter = contour_index_to_poly_index.find(parent);
                if (iter == contour_index_to_poly_index.end()) {
                    throw std::runtime_error("contour hierearchy had a child contour before its parent");
                }
                glp.blobs[iter->second].holes.push_back(contour_to_polygon(contour, false));
            }
        }
        return glp;
    }

    ch::polyline simplify(const ch::polyline& poly, double param) {
        std::vector<cv::Point2d> output;
        cv::approxPolyDP(poly | rv::transform([](auto p)->cv::Point2f {return p; }) | r::to_vector, output, param, true);
        return output;
    }

    ch::polygon_with_holes simplify(const ch::polygon_with_holes& poly, double param) {
        return {
            simplify(poly.border, param),
            poly.holes | rv::transform([param](const auto& p) {return simplify(p, param); }) | r::to_vector
        };
    }

    ch::gray_level_plane simplify(const ch::gray_level_plane& p, double param) {
        return {
            p.gray,
            p.blobs | rv::transform([param](const auto& p) {return simplify(p, param); }) | r::to_vector
        };
    }

    std::vector<ch::gray_level_plane> simplify(const std::vector<ch::gray_level_plane>& planes, double param) {
        return planes | rv::transform([param](const auto& p) {return simplify(p, param); }) | r::to_vector;
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

    std::string svg_path_commands(const ch::polygon_with_holes& poly, double scale) {
        std::stringstream ss;
        ss << loop_to_path_commands(poly.border, scale);
        for (const auto& hole : poly.holes) {
            ss << " " << loop_to_path_commands(hole, scale);
        }
        return ss.str();
    }

    std::string poly_with_holes_to_svg(uchar gray, const ch::polygon_with_holes& poly, double scale) {
        std::stringstream ss;
        ss << "<path fill-rule=\"evenodd\" stroke=\"none\" fill=\"";
        ss << ch::gray_to_svg_color(gray) << "\" d=\"";
        ss << svg_path_commands(poly, scale);
        ss << "\" />";
        return ss.str();
    }

    std::string gray_level_plane_to_svg(const ch::gray_level_plane& lvl, double scale) {
        std::stringstream ss;
        for (const auto& poly : lvl.blobs) {
            ss << poly_with_holes_to_svg(lvl.gray, poly, scale) << "\n";
        }
        return ss.str();
    }

    ch::polygon_with_holes scale(const ch::polygon_with_holes& poly, double scale) {
        return {
            ch::scale(poly.border, scale),
            poly.holes | rv::transform([scale](const auto& p) {return ch::scale(p, scale); }) | r::to_vector
        };
    }

    ch::gray_level_plane scale(const ch::gray_level_plane& glp, double scale) {
        return {
            glp.gray,
            glp.blobs | rv::transform([scale](const auto& poly) { return ::scale(poly,scale); }) | r::to_vector
        };
    }

    std::tuple<double, double, double, double> bounds(const ch::polyline& poly) {
        auto floats = poly | rv::transform([](const cv::Point2d& pt) {return cv::Point2f(pt); }) | r::to_vector;
        cv::Rect2d rect = cv::boundingRect(floats);
        return {
            static_cast<double>(rect.x),
            static_cast<double>(rect.y),
            static_cast<double>(rect.x + rect.width),
            static_cast<double>(rect.y + rect.height)
        };
    }

    ch::polygon_with_holes transform(const ch::polygon_with_holes& p, const ch::matrix& mat) {
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

    cl::Paths poly_with_holes_to_clipper_paths(const ch::polygon_with_holes& poly, double scale) {
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

    std::vector<ch::polyline> clip_swatch_to_poly( ch::crosshatching_swatch swatch, ch::polygon_with_holes& poly, double scale = 100.0) {
        cl::PolyTree polytree;
        cl::Clipper clipper;

        auto subject = polylines_to_clipper_paths(swatch.content | r::to<std::vector<ch::polyline>>(), scale);
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

    std::vector<ch::polyline> gray_levels_to_strokes(const std::vector<ch::gray_level_plane>& glps, ch::brush& brush) {
        std::vector<ch::polyline> output;
        output.reserve(10000);
        int count = 0;
        for (const auto& glp : glps) {
            std::cout << ++count << " of " << glps.size() << "\n";
            double gray = 1.0 - static_cast<double>(glp.gray) / 255.0;
            if (gray == 0.0) {
                continue;
            }
            int j = 0;
            for (const auto& poly : glp.blobs) {
                std::cout << "   " << ++j << " of " << glp.blobs.size() << "\n";
                auto strokes = ch::crosshatched_poly_with_holes(poly, gray, brush);
                std::copy(strokes.begin(), strokes.end(), std::back_inserter(output));
            }
        }
        return output;
    }
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

std::vector<ch::gray_level_plane> ch::extract_gray_level_planes(const cv::Mat& gray_scale_img)
{
    auto contours = perform_find_contours(gray_scale_img);
    return contours |
        rv::transform(
            [](const auto& tup)->ch::gray_level_plane {
                const auto& [gray, contour_info] = tup;
                return contour_info_to_gray_level_plane(gray, contour_info);
            }
        ) |
        r::to<std::vector<ch::gray_level_plane>>();
}

std::vector<ch::gray_levels> ch::extract_gray_levels(const cv::Mat& gray_scale_img) {
    return {};
}

std::vector<ch::gray_level_plane> ch::scale(const std::vector<gray_level_plane>& planes, double scale)
{
    return planes |
        rv::transform(
            [scale](const ch::gray_level_plane& plane)->ch::gray_level_plane {
                return ::scale(plane, scale);
            }
        ) |
        r::to_vector;
}

void ch::write_to_svg(const std::string& filename, const std::vector<gray_level_plane>& levels, int wd, int hgt, double scale)
{
    std::ofstream outfile(filename);

    outfile << svg_header(static_cast<int>(scale * wd), static_cast<int>(scale * hgt));

    for (const auto& lvl : levels)
        outfile << gray_level_plane_to_svg(lvl, scale);

    outfile << "</svg>" << std::endl;
    outfile.close();
}

void ch::debug() {
    cv::Mat img = cv::imread("C:\\test\\dunjon2.png");
    cv::Mat mat = ch::do_segmentation(img, 8, 3.0f, 12);
    auto gray_levels = extract_gray_level_planes(mat);
    //gray_levels = simplify(gray_levels, 1.0);
    write_to_svg("C:\\test\\dunjon2.svg", gray_levels, mat.cols, mat.rows, 4.0);
}

ch::drawing ch::generate_crosshatched_drawing(const std::string& image_file, segmentation_params params, double scale, brush& br)
{   
    
    cv::Mat img = cv::imread(image_file);
    cv::Mat mat = ch::do_segmentation(img, params.sigmaS, params.sigmaR, params.minSize);
    mat = ch::convert_to_gray(mat);
    auto gray_levels = extract_gray_level_planes(mat);
    gray_levels = ch::scale(gray_levels, scale);

    return ch::drawing{
        gray_levels_to_strokes(gray_levels, br),
        { img.cols * scale, img.rows * scale},
        br.stroke_width()
    };
}

std::vector<ch::polyline> ch::crosshatched_poly_with_holes(const ch::polygon_with_holes& input, double color, ch::brush& brush)
{
    auto [x1, y1, x2, y2] = bounds(input.border);
    cv::Point2d center = { (x1 + x2) / 2.0,(y1 + y2) / 2.0 };
    auto poly = ::transform(input, translation_matrix(-center));
    auto swatch = brush.get_hatching(color, { x2 - x1,y2 - y1 });
    auto strokes = clip_swatch_to_poly(swatch, poly);
    return transform(strokes, translation_matrix(center));
}

void ch::to_svg(const std::string& filename, const drawing& d)
{
    std::ofstream outfile(filename);

    outfile << svg_header(static_cast<int>(d.sz.wd), static_cast<int>(d.sz.hgt));

    for (const auto& poly : d.strokes)
        outfile << polyline_to_svg(poly, d.stroke_wd) << std::endl;

    outfile << "</svg>" << std::endl;
    outfile.close();
}
