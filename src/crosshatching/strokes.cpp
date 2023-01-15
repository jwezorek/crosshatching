#include "strokes.hpp"
#include "util.hpp"
#include "qdebug.h"
#include <variant>

namespace r = ranges;
namespace rv = ranges::views;

/*------------------------------------------------------------------------------------------------*/

namespace {
    void paint_strokes(QPainter& g, const ch::drawing_component& dc) {
        std::visit(
            overload{
                [&](const ch::stroke_group& sg) {
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
                },
                [&](const ch::stippling_group& sg) {
                   g.setPen(Qt::PenStyle::NoPen);
                   g.setBrush(QBrush(QColor(0, 0, 0)));
                   for (auto pt : sg.points) {
                       auto r = sg.thickness / 2.0;
                       g.drawEllipse(
                           pt.x - r, pt.y - r,
                           sg.thickness, sg.thickness
                       );
                   }
                },
                [&](const ch::shaded_polygon& sp) {
                   ch::paint_polygon(g, sp.poly, ch::ink_shade_to_color(sp.shade));
                }
            },
            dc
        );
    }

    ch::drawing_component transform(const ch::drawing_component& dc, const ch::matrix& mat) {
        return std::visit(
            overload{
                [&mat](const ch::stroke_group& sg)->ch::drawing_component {
                    return ch::stroke_group {
                        sg.strokes | rv::transform(
                            [&mat](const auto& poly)->ch::polyline {
                                return ch::transform(poly, mat);
                            }
                        ) | to_polylines,
                        sg.thickness
                    };
                },
                [&mat](const ch::stippling_group& sg)->ch::drawing_component {
                    return ch::stippling_group {
                        ch::transform(sg.points, mat),
                        sg.thickness
                    };
                },
                [&mat](const ch::shaded_polygon& sp)->ch::drawing_component {
                    return ch::shaded_polygon {
                        ch::transform(sp.poly, mat),
                        sp.shade
                    };
                }
            },
            dc
        );
    }
}

bool ch::is_empty(const ch::drawing_component& dc) {
    return std::visit(
        overload{
            [](const stroke_group& sg)->bool {
                return sg.strokes.empty();
            },
            [](const stippling_group& sg)->bool {
                return sg.points.empty();
            },
            [](const shaded_polygon& sp)->bool {
                return sp.poly.outer().empty();
            }
        },
        dc
    );
}

ch::drawing_comps ch::transform(const drawing_comps& dcs, const ch::matrix& mat) {
    return dcs |
        rv::transform(
            [&mat](const auto& dc)->drawing_component {
                return ::transform(dc, mat);
            }
        ) | r::to_vector;
}

void ch::paint_strokes(QPainter& g, const ch::drawing_comps& strks) {
    for (const auto& strk_group : strks) {
        ::paint_strokes(g, strk_group);
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
            [&poly](const drawing_component& dc)->drawing_component {
                return std::visit(
                    overload{
                        [&](const stroke_group& sg)->drawing_component {
                            return stroke_group {
                                clip_polylines_to_poly(sg.strokes, poly),
                                sg.thickness
                            };
                        },
                        [](const stippling_group& sg)->drawing_component {
                            return sg;
                        },
                        [](const shaded_polygon& sp)->drawing_component {
                            return sp;
                        }
                    },
                    dc
                );
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

