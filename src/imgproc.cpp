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
#include <fstream>

namespace r = ranges;
namespace rv = ranges::views;

namespace {

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
                    cv::findContours(mat, output.contours, output.hierarchy, cv::RETR_CCOMP, cv::CHAIN_APPROX_SIMPLE);

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

ch::gray_level_plane contour_info_to_gray_level_plane(uchar gray, const find_contour_output& contours) {
    ch::gray_level_plane glp = { gray,{} };
    std::unordered_map<int, int> contour_index_to_poly_index;
    for (int i = 0; i < contours.hierarchy.size(); ++i) {
        const cv::Vec4i& h = contours.hierarchy[i];
        int parent = h[3];
        auto polygon = int_polyline_to_polyline(contours.contours[i]);
        if (parent == -1) {
            int poly_index = static_cast<int>(glp.blobs.size());
            glp.blobs.emplace_back(
                ch::polygon_with_holes{ std::move(polygon), {}  }
            );
            contour_index_to_poly_index[i] = poly_index;
        } else {
            auto iter = contour_index_to_poly_index.find(parent);
            if (iter == contour_index_to_poly_index.end()) {
                throw std::runtime_error("contour hierearchy had a child contour before its parent");
            }
            glp.blobs[iter->second].holes.push_back(std::move(polygon));
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

std::string gray_level_plane_to_svg(const ch::gray_level_plane& lvl) {
    return {};
}

void ch::write_to_svg(const std::string& filename, const std::vector<gray_level_plane>& levels, int wd, int hgt, double scale)
{
    std::ofstream outfile(filename);

    outfile << svg_header(wd, hgt);

    for (const auto& lvl : levels)
        outfile << gray_level_plane_to_svg(lvl);

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
    cv::Mat mat = cv::imread("C:\\test\\holes.png");
    cv::Mat gray;
    cv::cvtColor(mat, gray, cv::COLOR_BGR2GRAY);
    auto gray_levels = extract_gray_level_planes(gray);
    for (const auto& glp : gray_levels) {
        std::cout << to_string(glp) << "\n";
    }

    //cv::imshow("gray", gray);
    //cv::imshow("contours", output);
}