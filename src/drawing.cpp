#include "drawing.hpp"
#include "geometry.hpp"
#include "clipper.hpp"
#include <stdexcept>

namespace r = ranges;
namespace rv = ranges::views;
namespace cl = ClipperLib;

namespace {

    std::tuple<double, double, double, double> bounds(const ch::polyline& poly) {
        cv::Rect2d rect = cv::boundingRect(poly);
        return {
            rect.x,
            rect.y,
            rect.x + rect.width,
            rect.y + rect.height
        };
    }

    ch::polygon_with_holes transform(const ch::polygon_with_holes& p, const ch::matrix& mat) {
        return {
            ch::transform(p.border, mat),
            p.holes | rv::transform([&mat](const auto& hole) {return ch::transform(hole, mat); }) | r::to_vector
        };
    }

    cl::Path polyline_to_clipper_path(const ch::polyline& poly, double scale) {

    }

    cl::Paths polylines_to_clipper_paths(const std::vector<ch::polyline>& polys, double scale) {

    }

    cl::Paths poly_with_holes_to_clipper_paths(const ch::polygon_with_holes& poly, double scale) {

    }

    std::vector<ch::polyline> clipper_paths_to_polyline(const cl::Paths paths, double scale) {

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

        return clipper_paths_to_polyline(clipped, 1.0 / scale);
    }
}

std::vector<ch::polyline> ch::crosshatched_poly_with_holes(const ch::polygon_with_holes& input, uchar color, ch::brush& brush)
{
    auto [x1, y1, x2, y2] = bounds(input.border);
    cv::Point2d center = { (x1 + x2) / 2.0,(y1 + y2) / 2.0 };
    auto poly = ::transform(input, translation_matrix(-center));
    auto swatch = brush.get_hatching(color, { x2 - x1,y2 - y1 });
    auto strokes = clip_swatch_to_poly(swatch, poly);
    return transform(strokes, translation_matrix(-center));
}

void ch::to_svg(const std::string& filename, const drawing& d)
{
}
