#include "strokes.hpp"
#include "util.hpp"
#include "qdebug.h"

namespace r = ranges;
namespace rv = ranges::views;

/*------------------------------------------------------------------------------------------------*/

void paint_strokes(QPainter& g, const ch::stroke_group& sg) {
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

ch::stroke_group::stroke_group() : 
    thickness(0), 
    is_stippling(false)
{
}

ch::stroke_group::stroke_group(ch::polylines&& strks, double thk, bool isstip) :
    strokes(strks),
    thickness(thk),
    is_stippling(isstip)
{}


ch::stroke_groups ch::transform(const stroke_groups& s, const ch::matrix& mat) {
    return s |
        rv::transform(
            [&mat](const auto& sg)->stroke_group {
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

void ch::paint_strokes(QPainter& g, const ch::stroke_groups& strks) {
    for (const auto& strk_group : strks) {
        paint_strokes(g, strk_group);
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
    return
        *strokes |
        rv::transform(
            [&poly](const stroke_group& sg)->stroke_group {
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

