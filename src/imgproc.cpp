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

    cv::Mat dilate_blob(const cv::Mat& blob) {
        static cv::Mat strel = (
            cv::Mat_<uchar>(3, 3) <<
                0, 0, 0,
                0, 1, 1,
                0, 1, 1
        );
        cv::Mat output;
        cv::dilate(blob, output, strel, { 1,1 });
        return output;
    }

    std::vector<std::vector<cv::Point>> apply_douglas_peucker(const std::vector<std::vector<cv::Point>>& polys, float param) {
        return polys |
            rv::transform(
                [param](const std::vector<cv::Point>& poly)->std::vector<cv::Point> {
                    std::vector<cv::Point> output;
                    cv::approxPolyDP(poly, output, param, true);
                    return output;
                }
            ) |
            r::to<std::vector<std::vector<cv::Point>>>();
    }


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

    
    void apply_douglas_peucker(std::map<uchar, find_contour_output>& ct, float param) {
        for (auto& [k, v] : ct) {
            v.contours = apply_douglas_peucker(v.contours, param);
        }
    }

    /*
    void scale(std::map<uchar, contours_trees>& cont_trees, double s) {
        for (auto& [k, cont_tree] : cont_trees) {
            for (auto& polyline : cont_tree.contours) {
                for (auto& pt : polyline) {
                    pt = { 
                        static_cast<int>(std::round(s * pt.x)), 
                        static_cast<int>(std::round(s * pt.y))
                    };
                }
            }
        }
    }
    */

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
    return { to_pt, dir_to_flip_dir[pc.head] };
}

bool is_cardinal_dir( direction dir) {
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

std::tuple<pixel_crawler, int> initialize_contour_crawl(const int_polyline& ip) {
    int northwest_index = find_northwest_index(ip);
    return { {ip[northwest_index], direction::NW}, northwest_index };
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

void push_if_new(ch::polyline& poly, const cv::Point2d& pt) {
    if (!poly.empty() && poly.back() == pt) {
        return;
    }
    poly.push_back(pt);
}

ch::polyline contour_to_polygon(const int_polyline& ip) {
    ch::polyline poly;
    int n = static_cast<int>(ip.size());
    poly.reserve(n);
    
    auto [crawler,i] = initialize_contour_crawl(ip);

    while (poly.empty() || (get_vertex(crawler) != poly.front())) {
        auto current_vert = get_vertex(crawler);
        push_if_new(poly, current_vert);
        int next_i = (i + 1) % n;
        if (can_flip(crawler, ip[next_i])) {
            crawler = roll(flip(crawler, ip[next_i]));
            i = next_i;
        } else {
            crawler = roll(crawler);
        }
    }

    poly.shrink_to_fit();
    return poly;
}

ch::polyline contour_to_polygon_clockwise(const int_polyline& ip) {
    //TODO
    auto contour_cc = rv::reverse(ip) | r::to<int_polyline>();
    auto poly_cc = contour_to_polygon(contour_cc);
    return rv::reverse(poly_cc) | r::to< ch::polyline>();
}

ch::gray_level_plane contour_info_to_gray_level_plane(uchar gray, const find_contour_output& contours) {
    ch::gray_level_plane glp = { gray,{} };
    std::unordered_map<int, int> contour_index_to_poly_index;

    if (gray == 255) {
        int aaa;
        aaa = 5;
    }

    for (int i = 0; i < contours.hierarchy.size(); ++i) {
        const cv::Vec4i& h = contours.hierarchy[i];
        int parent = h[3];
        const auto& contour = contours.contours[i]; 
        if (parent == -1) {
            int poly_index = static_cast<int>(glp.blobs.size());
            glp.blobs.emplace_back(
                ch::polygon_with_holes{ contour_to_polygon(contour), {}  }
            );
            contour_index_to_poly_index[i] = poly_index;
        } else {
            auto iter = contour_index_to_poly_index.find(parent);
            if (iter == contour_index_to_poly_index.end()) {
                throw std::runtime_error("contour hierearchy had a child contour before its parent");
            }
            glp.blobs[iter->second].holes.push_back( contour_to_polygon_clockwise(contour) );
        }
    }
    return glp;
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

std::string loop_to_path_commands(const ch::polyline& poly, double scale) {
    std::stringstream ss;
    ss << "M " << scale*poly[0].x << "," << scale * poly[0].y << " L";
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

void ch::write_to_svg(const std::string& filename, const std::vector<gray_level_plane>& levels, int wd, int hgt, double scale)
{
    std::ofstream outfile(filename);

    outfile << svg_header( static_cast<int>(scale*wd), static_cast<int>(scale*hgt));

    for (const auto& lvl : levels)
        outfile << gray_level_plane_to_svg(lvl, scale);

    outfile << "</svg>" << std::endl;
    outfile.close();
}

std::string ch::to_string(const ch::gray_level_plane& glp) {
    std::stringstream ss;

    ss << "{\n";
    ss << "   " << static_cast<int>(glp.gray) << ",\n";
    for (const auto& pwh : glp.blobs) {
        ss << "   " << to_string(pwh) << "\n";
    }
    ss << "}";

    return ss.str();
}

std::string ch::to_string(const ch::polygon_with_holes& poly) {
    std::stringstream ss;
    if (poly.holes.empty()) {
        ss << poly_to_string(poly.border);
    } else {
        ss << "{ " << poly_to_string(poly.border);
        for (const auto& hole : poly.holes) {
            ss << " , " << poly_to_string(hole);
        }
        ss << "}";
    }
    return ss.str();
}

std::string vec4_to_string(const cv::Vec4i& v) {
    std::stringstream ss;
    ss << "[ ";
    for (int i = 0; i < 4; ++i) {
        ss << v[i] << " ";
    }
    ss << "]\n";
    return ss.str();
}

void debug_contours(const find_contour_output& ct) {
    for (const auto& poly : ct.contours) {
        std::cout << ch::poly_to_string(poly) << "\n";
    }
    for (const auto& h : ct.hierarchy) {
        std::cout << vec4_to_string(h) << " ";
    }
    std::cout << "\n";
}

void debug_contours(uchar gray, const find_contour_output& ct) {
    std::cout << " { " << static_cast<int>(gray) << "\n";
    debug_contours(ct);
    std::cout << "}\n";
}

void debug_contours(const std::vector<std::tuple<uchar, find_contour_output>>& contours) {
    for (const auto [gray, trees] : contours) {
        debug_contours(gray, trees);
    }
}

void ch::debug() {
    cv::Mat mat = cv::imread("C:\\test\\big_blob.png");
    cv::Mat gray;
    cv::cvtColor(mat, gray, cv::COLOR_BGR2GRAY);
    auto gray_levels = extract_gray_level_planes(gray);
    write_to_svg("C:\\test\\big_blob.svg", gray_levels, mat.cols, mat.rows, 10.0);
    //cv::imshow("gray", gray);
    //cv::imshow("contours", output);
}