#include "crosshatch.hpp"
#include "geometry.hpp"
#include "util.hpp"
#include <opencv2/highgui.hpp>
#include <iostream>
#include <cmath>
#include <numbers>
#include <sstream>
#include <fstream>

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    struct run_and_space {
        double run;
        double space;
        double x;
    };

    auto row_of_crosshatching(double wd, double hgt, double y, const ch::rnd_fn& run_length, const ch::rnd_fn& space_length, const ch::unit_of_hatching_fn& h_fn) {
        auto half_wd = wd / 2.0;
        return rv::join(
            rv::generate([=]() { return run_and_space{ run_length(), space_length(), -half_wd }; }) |
            rv::partial_sum(
                [](run_and_space a, run_and_space b)->run_and_space {
                    return { b.run, b.space, a.x + a.run + a.space };
                }
            ) |
            rv::transform(
                [y, hgt, h_fn](run_and_space rs) {
                    auto x1 = rs.x;
                    auto x2 = rs.x + rs.run - 1;
                    return h_fn(x1, x2, y, hgt);
                }
            ) |
            rv::take_while([half_wd](auto rng) {
                    return (*rng.begin()).front().x < half_wd;
                }
            )
       );
    }

    auto fragment_length(double len, const ch::rnd_fn& run_length) {
        return rv::concat(rv::single(0.0), rv::generate(run_length)) |
            rv::partial_sum |
            rv::take_while([len](double y) { return y < len; });
    }

    auto fragment_line_segment(const ch::point& pt1, const ch::point& pt2, const ch::rnd_fn& run_length) {
        auto len = ch::euclidean_distance(pt1, pt2);
        auto theta = std::atan2(pt2.y - pt1.y, pt2.x - pt1.x);
        auto x_delta = std::cos(theta);
        auto y_delta = std::sin(theta);

        return fragment_length(len, run_length) |
            rv::transform(
                [x_delta, y_delta, pt1](double d)->ch::point {
                    return pt1 + ch::point{ d * x_delta, d * y_delta };
                }
            );
    }

    auto jitter_line_segment(const ch::point& pt1, const ch::point& pt2, const ch::rnd_fn& jit) {
        auto theta = std::atan2(pt2.y - pt1.y, pt2.x - pt1.x);
        ch::point j = { std::cos(theta + std::numbers::pi / 2.0) , std::sin(theta + std::numbers::pi / 2.0) };
        return rv::single(pt1) |
            rv::transform(
                [j,jit](const ch::point& pt)->ch::point {
                    return jit() * j + pt;
                }
            );
    }

    ch::polyline jitter(const ch::polyline& poly, const ch::rnd_fn& jit) {
        return r::to_vector(
            rv::concat(
                rv::join(
                    poly |
                    rv::sliding(2) |
                    rv::transform(
                        [jit](auto rng) {
                            return jitter_line_segment(rng[0], rng[1], jit);
                        }
                    )
                ),
                rv::single( poly.back() )
            )
        );
    }

    ch::polyline fragment(const ch::polyline& poly, const ch::rnd_fn& frag) {
        return r::to_vector(
            rv::concat(
                rv::join(
                    poly |
                    rv::sliding(2) |
                    rv::transform(
                        [frag](auto rng) {
                            return fragment_line_segment(rng[0], rng[1], frag);
                        }
                    )
                ),  
                rv::single(poly.back())
            )
        );
    }
    
    ch::polyline jiggle(const ch::polyline& poly, const ch::rnd_fn& jig) {
        auto p = ch::mean_point(poly);
        auto theta = jig();
        ch::matrix rotate = ch::translation_matrix(p.x, p.y) * ch::rotation_matrix(theta) * ch::translation_matrix(-p.x, -p.y);
        return ch::transform(poly, rotate);
    }

    struct running_sum_item {
        double next_item;
        double sum;

        running_sum_item(double v = 0, double s = 0) : 
            next_item(v), 
            sum(s)
        {}

        double next_sum() const {
            return next_item + sum;
        }
    };  

    auto running_sum(const ch::rnd_fn& gen_fn, double init_sum) {
        return 
            rv::generate(gen_fn) |
            rv::transform(
                [init_sum](double v)->running_sum_item {
                    return { v, init_sum };
                }
            ) |
            rv::partial_sum(
                [](const running_sum_item& item1, const running_sum_item& item2) ->running_sum_item {
                    return { item2.next_item, item1.next_sum() };
                }
            );
    }

    ch::crosshatching_range brick_stroke(double x1, double x2, double y, double hgt) {
        return rv::concat(
            rv::single(ch::polyline{ {x1,y} , {x2,y} }),
            rv::single(ch::polyline{ {(x1 + x2) / 2.0,y} , {(x1 + x2) / 2.0,y + hgt} })
        );
    }

}

ch::unit_of_hatching_fn ch::make_brick_stroke()
{
    return [](double x1, double x2, double y, double hgt) {return brick_stroke(x1, x2, y, hgt); };
}

double shading_stroke_length(double stroke_len, double sz_pcnt, double stddev) {
    auto len = ch::normal_rnd(sz_pcnt * stroke_len, stddev);
    if (len <= 1.0) {
        return 1.0;
    } else if (len > stroke_len) {
        return stroke_len;
    } 
    return len;
}

ch::unit_of_hatching_fn ch::make_shading_stroke(double sz_pcnt, double stddev, bool centered) {
    return [=](double x1, double x2, double y, double hgt)->crosshatching_range {
        auto n = static_cast<int>(hgt)+1;
        auto total_len = x2 - x1;
        return
            rv::iota(0, n) |
            rv::transform(
                [=](int i) {
                    auto len = (i == 0) ? total_len : total_len * std::pow(sz_pcnt, i + 1) + ch::normal_rnd(0, stddev);
                    if (len < 1.0) {
                        len = 1.0;
                    } else if (len > total_len) {
                        len = total_len;
                    }
                    auto x = centered ? (total_len - len) / 2.0 + x1 : x2 - len;
                    auto yy = (n-i) + y;

                    return polyline{
                        {x, yy},{x + len,yy}
                    };
                }
            ); 
    };
}

ch::dimensions operator*(const ch::dimensions& sz, double k) {
    return {
        sz.wd * k,
        sz.hgt * k
    };
}

ch::crosshatching_range ch::one_horz_stroke(double x1, double x2, double y, double hgt) {
    return rv::single(ch::polyline{ {x1,y} , {x2,y} });
}

ch::crosshatching_swatch ch::linear_crosshatching( rnd_fn run_length, rnd_fn space_length, rnd_fn vert_space, 
        unit_of_hatching_fn h_fn, ch::dimensions viewable_swatch_sz) {
    auto dim = (std::sqrt(2.0) * (viewable_swatch_sz.wd + viewable_swatch_sz.hgt)) / 2.0;
    auto swatch_sz = dimensions(dim);
    return {
        rv::join(
            running_sum(vert_space, -swatch_sz.hgt / 2.0) |
            rv::take_while([swatch_sz](const running_sum_item& rsi) { return rsi.sum < swatch_sz.hgt / 2.0; }) |
            rv::transform([=](const running_sum_item& rsi) {return row_of_crosshatching(swatch_sz.wd, rsi.next_item, rsi.sum, run_length, space_length, h_fn); })
        ),
        viewable_swatch_sz
    };
}

ch::crosshatching_swatch ch::fragment(ch::crosshatching_swatch swatch, ch::rnd_fn frag) {
    return  {
        swatch.content |
        rv::transform(
            [frag](const auto& poly) {
                return ::fragment(poly, frag);
            }
         ),
        swatch.sz
    };
}

ch::crosshatching_swatch ch::jitter(ch::crosshatching_swatch swatch, ch::rnd_fn jit)
{
    return {
        swatch.content |
        rv::transform(
            [jit](const auto& poly) {
                return ::jitter(poly, jit);
            }
        ),
        swatch.sz
    };
}

ch::crosshatching_swatch ch::jiggle(ch::crosshatching_swatch swatch,  ch::rnd_fn jig) {
    return {
        swatch.content |
            rv::transform(
                [jig](const auto& poly) {
                    return ::jiggle(poly, jig);
                }
            ),
        swatch.sz
    };
}

ch::crosshatching_swatch ch::rotate(ch::crosshatching_swatch swatch, double theta) {
    matrix rotation = rotation_matrix(theta);
    return {
        ch::transform(swatch.content, rotation),
        swatch.sz
    };
}

ch::crosshatching_swatch ch::disintegrate(ch::crosshatching_swatch swatch, double amount) {

    if (amount >= 1.0) {
        return swatch;
    } else if (amount == 0) {
        return { {}, swatch.sz };
    }

    return { swatch.content |
        rv::filter(
            [amount](const auto& p) {
                return ch::uniform_rnd(0.0, 1.0) < amount;
            }
        ),
        swatch.sz 
    };
}

cv::Mat ch::paint_cross_hatching(int thickness, ch::crosshatching_swatch swatch) {
    cv::Mat mat(static_cast<int>(swatch.sz.hgt), static_cast<int>(swatch.sz.wd), CV_8U, 255);
    for (const auto& ls : swatch.content) {
        ch::paint_polyline(mat, ls, thickness, 0, point{ swatch.sz.wd / 2.0, swatch.sz.hgt / 2.0 });
    }
    return mat;
}

double ch::gray_level(int thickness, crosshatching_swatch swatch)
{
    auto mat = paint_cross_hatching(thickness, swatch);
    auto n = swatch.sz.wd  * swatch.sz.hgt;
    auto white_pixels = cv::countNonZero(mat);
    return static_cast<double>(n - white_pixels) / static_cast<double>(n);
}

//TODO: there needs to be an offset in here because of the new way rotation is handled.
std::string polyline_to_svg(const ch::polyline& poly, int thickness) {
    std::stringstream ss;
    ss << "<polyline points=\"";
    for (const auto& pt : poly) {
       
        ss << " " << pt.x << "," << pt.y;
    }
    ss << "\" style=\"fill:none;stroke:black;stroke-width:" << thickness << "\" />";
    return ss.str();
}

void ch::to_svg(const std::string& filename, int thickness, crosshatching_swatch swatch)
{
    std::ofstream outfile(filename);

    outfile << svg_header(static_cast<int>(swatch.sz.wd), static_cast<int>(swatch.sz.hgt));

    for (const auto& poly : swatch.content)
        outfile << polyline_to_svg(poly, thickness) << std::endl;

    outfile << "</svg>" << std::endl;
    outfile.close();
}

ch::dimensions::dimensions(double d) : wd(d), hgt(d)
{
}

ch::dimensions::dimensions(double w, double h) : wd(w),hgt(h)
{
}
