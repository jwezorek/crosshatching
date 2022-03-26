#include "drawing.hpp"
#include "geometry.hpp"
#include "clipper.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <stdexcept>
#include <fstream>
#include <iostream>

namespace r = ranges;
namespace rv = ranges::views;
namespace cl = ClipperLib;

namespace {

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

ch::drawing ch::generate_crosshatched_drawing(const std::string& image_file, segmentation_params params, double scale, brush& br)
{
    cv::Mat img = cv::imread(image_file);
    cv::Mat mat = ch::do_segmentation(img, params.sigmaS, params.sigmaR, params.minSize);
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
