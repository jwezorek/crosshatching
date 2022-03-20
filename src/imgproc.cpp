#include "imgproc.hpp"
#include "meanshift.hpp"
#include "util.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <range/v3/all.hpp>
#include <set>
#include <map>
#include <stdexcept>
#include <array>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <fstream>

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    enum class direction {
        N,NE,E,SE,S,SW,W,NW
    };

    std::vector<uchar> get_gray_levels(const cv::Mat& input) {
        if (input.channels() != 1) {
            throw std::runtime_error("called get_gray_levels on color image");
        }
        std::array<bool,256> grays = {};
        input.forEach<uchar>( [&grays](uchar gray, const int* pos) { grays[gray] = true;  } );
        return rv::iota(0) | 
            rv::take(256) | 
            rv::filter([&grays](int g) {return grays[g]; }) | 
            r::to<std::vector<uchar>>();
    }

    std::vector<std::tuple<uchar, cv::Mat>> gray_planes(const cv::Mat& input) {
        auto levels = get_gray_levels(input);
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

    cv::Mat paint_contours(int wd, int hgt, const std::vector<std::tuple<uchar, find_contour_output>>& cont) {
        auto img = cv::Mat(hgt, wd, CV_8UC1, 255);
        for (const auto& [gray, c] : cont) {
            int n = static_cast<int>(c.contours.size());
            for ( int i = 0; i < n; ++i) {
                cv::drawContours(img, c.contours, i, gray, cv::FILLED, 8, c.hierarchy);
            }
        }
        return img;
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
        return offset_to_direction[to_pt - from_pt];
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

    pixel_crawler flip(const pixel_crawler& pc, const cv::Point& to_pt) {
        auto dir = direction_to(pc.loc, to_pt);
        auto old_v = get_vertex(pc);
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
        if (get_vertex(flipped_crawler) != old_v) {
            flipped_crawler = roll(flipped_crawler);
        }

        return flipped_crawler;
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

    bool is_left_of(direction d1, direction d2) {
        static std::unordered_map<direction, direction> right_to_left = {
            {direction::N , direction::NW },
            {direction::NW, direction::W},
            {direction::W , direction::SW },
            {direction::SW, direction::S},
            {direction::S , direction::SE },
            {direction::SE, direction::E},
            {direction::E , direction::NE },
            {direction::NE, direction::N}
        };
        return right_to_left[d1] == d2;
    }

    bool can_flip(const pixel_crawler& pc, const cv::Point& to_pt) {
        auto dir = direction_to(pc.loc, to_pt);
        if (pc.head == dir) {
            return true;
        }
        if (is_ordinal_dir(dir)) {
            return false;
        }
        return is_left_of(pc.head, dir);
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

}

cv::Mat ch::do_segmentation(const cv::Mat& input, int sigmaS, float sigmaR, int minSize) {
    cv::Mat output;
    cv::Mat labels;
    auto mss = createMeanShiftSegmentation(sigmaS, sigmaR, minSize, true);
    mss->processImage(input, output, labels);

    cv::Mat gray;
    cv::cvtColor(output, gray, cv::COLOR_BGR2GRAY);

    return gray;
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

void ch::write_to_svg(const std::string& filename, const std::vector<gray_level_plane>& levels, int wd, int hgt, double scale)
{
    std::ofstream outfile(filename);

    outfile << svg_header( static_cast<int>(scale*wd), static_cast<int>(scale*hgt));

    for (const auto& lvl : levels)
        outfile << gray_level_plane_to_svg(lvl, scale);

    outfile << "</svg>" << std::endl;
    outfile.close();
}

void ch::debug() {
    cv::Mat img = cv::imread("C:\\test\\dunjon.png");
    cv::Mat mat = ch::do_segmentation(img, 8, 3.0f, 12);
    auto gray_levels = extract_gray_level_planes(mat);
    //gray_levels = simplify(gray_levels, 1.0);
    write_to_svg("C:\\test\\dunjon.svg", gray_levels, mat.cols, mat.rows, 4.0);
}

/*
void ch::debug() {
    cv::Mat mat = cv::imread("C:\\test\\big_blob.png");
    cv::Mat gray;
    cv::cvtColor(mat, gray, cv::COLOR_BGR2GRAY);
    auto gray_levels = extract_gray_level_planes(gray);
    write_to_svg("C:\\test\\big_blob.svg", gray_levels, mat.cols, mat.rows, 10.0);
    //cv::imshow("gray", gray);
    //cv::imshow("contours", output);
}
*/