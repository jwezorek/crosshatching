#include "strokes.hpp"
#include "util.hpp"

namespace r = ranges;
namespace rv = ranges::views;

/*------------------------------------------------------------------------------------------------*/

namespace {
    ch::drawn_stroke_cluster to_drawn_stroke_cluster(ch::stroke_cluster cluster) {
        return {
            ch::to_polylines(cluster.strokes),
            cluster.thickness
        };
    }

    ch::stroke_range transform(ch::stroke_range sr, const ch::matrix& mat) {
        return sr | rv::transform(
            [mat](const auto& pt) {
                return ch::transform(pt, mat);
            }
        );
    }

    ch::stroke_ranges transform_stroke_ranges(ch::stroke_ranges ranges, const ch::matrix& mat) {
        return ranges |
            rv::transform(
                [mat](auto rng)->r::any_view<ch::point> {
                    return rng | rv::transform(
                        [mat](auto pt)->ch::point {
                            return ch::transform(pt, mat);
                        }
                    );
                }
        );
    }
}

ch::stroke_cluster ch::transform(stroke_cluster sc, const ch::matrix& mat) {
    return {
        sc.strokes | rv::transform(
            [mat](auto sr) {
                return ::transform(sr, mat);
            }
        ),
        sc.thickness
    };
}

ch::strokes ch::transform(ch::strokes strokes, const ch::matrix& mat) {
    return strokes |
        rv::transform(
            [mat](auto stroke) {
                return ch::transform(stroke, mat);
            }
    );
}

ch::polylines ch::to_polylines(ch::stroke_ranges ranges) {
    auto polys = ranges |
        rv::transform(
            [](ch::stroke_range sr)->ch::polyline {
                /*
                std::stringstream ss;
                for (auto pt : sr) {
                    ss << ch::to_string(pt) << " ";
                }
                qDebug() << ss.str().c_str();
                */
                return sr | r::to<ch::polyline>();
            }
    );

    ch::polylines output;
    output.resize(r::distance(polys));
    for (auto&& [i, poly] : rv::enumerate(polys)) {
        output[i] = std::move(poly);
    }

    return output;
}

ch::drawn_strokes ch::to_drawn_strokes(ch::strokes strks) {
    return strks | rv::transform(to_drawn_stroke_cluster) | r::to_vector;
}

void ch::paint_strokes(QPainter& g, const ch::drawn_strokes& strks) {

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