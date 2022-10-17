#include "strokes.hpp"
#include "util.hpp"

namespace r = ranges;
namespace rv = ranges::views;

/*------------------------------------------------------------------------------------------------*/

ch::stroke_groups ch::transform(const stroke_groups& s, const ch::matrix& mat) {
    return s |
        rv::transform(
            [&mat](const auto& cluster)->stroke_group {
                return {
                    cluster.strokes | rv::transform(
                        [&mat](const auto& poly)->ch::polyline {
                            return ch::transform(poly, mat);
                        }
                    ) | to_polylines,
                    cluster.thickness
                };
            }
        ) | r::to_vector;
}

void ch::paint_strokes(QPainter& g, const ch::stroke_groups& strks) {

    for (const auto& cluster : strks) {
        g.setPen(create_pen(0, cluster.thickness));
        for (auto strk : cluster.strokes) {
            QList<QPointF> points = strk |
                rv::transform(
                    [](const ch::point& pt)->QPointF {
                        return {
                            static_cast<float>(pt.x),
                            static_cast<float>(pt.y)
                        };
                    }
            ) | r::to<QList>();
            g.drawLines(points);
        }
    }
}

void ch::paint_strokes(QPainter& g, strokes_ptr str_ptr) {
    return paint_strokes(g, *str_ptr);
}

void ch::append(strokes_ptr dst, strokes_ptr src) {
    dst->insert(
        dst->end(),
        std::make_move_iterator(src->begin()),
        std::make_move_iterator(src->end())
    );
}

std::string ch::to_svg(const ch::stroke_groups& s) {
    return {}; //TODO
}

ch::strokes_ptr ch::clip_strokes(const polygon& poly, strokes_ptr strokes) {
    return to_strokes(
        *strokes |
        rv::transform(
            [&poly](const stroke_group& sg)->stroke_group {
                return {
                    clip_polylines_to_poly(sg.strokes, poly),
                    sg.thickness
                };
            }
        )
    );
}

ch::strokes_ptr ch::transform(ch::strokes_ptr strokes, const ch::matrix& mat) {
    return std::make_shared<stroke_groups>(
        transform(*strokes, mat)
    );
}

cv::Mat ch::strokes_to_mat(const stroke_groups& strokes, cv::Mat mat) {
    auto qimg = ch::mat_to_qimage(mat, false);
    QPainter g(&qimg);
    g.setRenderHint(QPainter::Antialiasing, true);
    ch::paint_strokes(g, strokes);
    return mat;
}

cv::Mat ch::strokes_to_mat(strokes_ptr strokes, cv::Mat mat) {
    return strokes_to_mat(*strokes, mat);
}

