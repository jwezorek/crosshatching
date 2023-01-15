#include "strokes.hpp"
#include "util.hpp"
#include "qdebug.h"

namespace r = ranges;
namespace rv = ranges::views;

/*------------------------------------------------------------------------------------------------*/

void paint_strokes(QPainter& g, const ch::drawing_component& sg) {
    if (!sg.is_stippling) {
        g.setPen(ch::create_pen(0, sg.thickness));
        g.setBrush(Qt::BrushStyle::NoBrush);
        for (auto strk : sg.strokes) {
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
    } else {
        g.setPen(Qt::PenStyle::NoPen);
        g.setBrush(QBrush(QColor(0, 0, 0)));
        for (auto pt : sg.strokes.front()) {
            auto r = sg.thickness / 2.0;
            g.drawEllipse(
                pt.x - r, pt.y - r,
                sg.thickness, sg.thickness
            );
        }
    }
}

ch::drawing_component::drawing_component() : 
    thickness(0), 
    is_stippling(false)
{
}

ch::drawing_component::drawing_component(ch::polylines&& strks, double thk, bool isstip) :
    strokes(strks),
    thickness(thk),
    is_stippling(isstip)
{}


ch::drawing_comps ch::transform(const drawing_comps& s, const ch::matrix& mat) {
    return s |
        rv::transform(
            [&mat](const auto& sg)->drawing_component {
                return {
                    sg.strokes | rv::transform(
                        [&mat](const auto& poly)->ch::polyline {
                            return ch::transform(poly, mat);
                        }
                    ) | to_polylines,
                    sg.thickness,
                    sg.is_stippling
                };
            }
        ) | r::to_vector;
}

void ch::paint_strokes(QPainter& g, const ch::drawing_comps& strks) {
    for (const auto& strk_group : strks) {
        paint_strokes(g, strk_group);
    }
}

void ch::paint_strokes(QPainter& g, drawing_comps_ptr str_ptr) {
    return paint_strokes(g, *str_ptr);
}

void ch::append(drawing_comps_ptr dst, drawing_comps_ptr src) {
    dst->insert(
        dst->end(),
        std::make_move_iterator(src->begin()),
        std::make_move_iterator(src->end())
    );
}

std::string ch::to_svg(const ch::drawing_comps& s) {
    return {}; //TODO
}

ch::drawing_comps_ptr ch::clip_strokes(const polygon& poly, drawing_comps_ptr strokes) {
    return
        *strokes |
        rv::transform(
            [&poly](const drawing_component& sg)->drawing_component {
                if (!sg.is_stippling) {
                    return {
                        clip_polylines_to_poly(sg.strokes, poly),
                        sg.thickness,
                        false
                    };
                } else {
                    return sg;
                }
            }
        ) | to_strokes;
}

ch::drawing_comps_ptr ch::transform(ch::drawing_comps_ptr strokes, const ch::matrix& mat) {
    return std::make_shared<drawing_comps>(
        transform(*strokes, mat)
    );
}

cv::Mat ch::strokes_to_mat(const drawing_comps& strokes, cv::Mat mat) {
    auto qimg = ch::mat_to_qimage(mat, false);
    QPainter g(&qimg);
    g.setRenderHint(QPainter::Antialiasing, true);
    ch::paint_strokes(g, strokes);
    return mat;
}

cv::Mat ch::strokes_to_mat(drawing_comps_ptr strokes, cv::Mat mat) {
    return strokes_to_mat(*strokes, mat);
}

