#include "crosshatch.hpp"
#include "geometry.hpp"
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

    auto row_of_crosshatching(int wd, double hgt, double y, const ch::rnd_fn& run_length, const ch::rnd_fn& space_length, const ch::unit_of_hatching_fn& h_fn) {
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

    auto divide_length(double len, const ch::rnd_fn& run_length) {
        return rv::concat(rv::single(0.0), rv::generate(run_length)) |
            rv::partial_sum |
            rv::take_while([len](double y) { return y < len; });
    }

    auto divide_line_segment(const ch::point& pt1, const ch::point& pt2, const ch::rnd_fn& run_length) {
        auto len = ch::euclidean_distance(pt1, pt2);
        auto theta = std::atan2(pt2.y - pt1.y, pt2.x - pt1.x);
        auto x_delta = std::cos(theta);
        auto y_delta = std::sin(theta);

        return divide_length(len, run_length) |
            rv::transform(
                [x_delta, y_delta, pt1](double d)->ch::point {
                    return pt1 + ch::point{ d * x_delta, d * y_delta };
                }
            );
    }

    auto jitter_line_segment(const ch::point& pt1, const ch::point& pt2, const ch::rnd_fn& run_length, const ch::rnd_fn& jitter) {
        auto theta = std::atan2(pt2.y - pt1.y, pt2.x - pt1.x);
        ch::point j = { std::cos(theta + std::numbers::pi / 2.0) , std::sin(theta + std::numbers::pi / 2.0) };
        return divide_line_segment(pt1, pt2, run_length) |
            rv::transform(
                [j, jitter](const ch::point& pt)->ch::point {
                    return jitter() * j + pt;
                }
            );
    }

    ch::polyline apply_jitter(const ch::polyline& poly, const ch::rnd_fn& run_length, const ch::rnd_fn& jitter) {
        return r::to_vector(
            rv::concat(
                rv::join(
                    poly |
                    rv::sliding(2) |
                    rv::transform(
                        [run_length, jitter](auto rng) {
                            return jitter_line_segment(rng[0], rng[1], run_length, jitter);
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
                            return divide_line_segment(rng[0], rng[1], frag);
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
        auto rotate = ch::translation_matrix(p.x, p.y) * ch::rotation_matrix(theta) * ch::translation_matrix(-p.x, -p.y);
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

    ch::hatching_range brick_stroke(double x1, double x2, double y, double hgt) {
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
    return [=](double x1, double x2, double y, double hgt)->hatching_range {
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

ch::hatching_range ch::one_horz_stroke(double x1, double x2, double y, double hgt) {
    return rv::single(ch::polyline{ {x1,y} , {x2,y} });
}

ch::hatching_range ch::linear_crosshatching(const ch::rnd_fn& run_length, const ch::rnd_fn& space_length,
    const ch::rnd_fn& vert_space, const ch::unit_of_hatching_fn& h_fn, int viewable_swatch_sz)
{
    auto swatch_sz = viewable_swatch_sz * 1.5;
    auto half_swatch_sz = swatch_sz / 2.0;
    return rv::join(
        running_sum(vert_space, -half_swatch_sz) |
        rv::take_while([half_swatch_sz](const running_sum_item& rsi) { return rsi.sum < half_swatch_sz; }) |
        rv::transform([=](const running_sum_item& rsi) {return row_of_crosshatching(swatch_sz, rsi.next_item, rsi.sum, run_length, space_length, h_fn); })
    );
}

ch::hatching_range ch::fragment(ch::hatching_range rng, const ch::rnd_fn& frag) {
    return rng |
        rv::transform(
            [frag](const auto& poly) {
                return ::fragment(poly, frag);
            }
    );
}

ch::hatching_range ch::jitter(ch::hatching_range rng, const ch::rnd_fn& run_length, const ch::rnd_fn& jitter)
{
    return rng |
        rv::transform(
            [run_length, jitter](const auto& poly) {
                return ::apply_jitter(poly, run_length, jitter);
            }
        );
}

ch::hatching_range ch::jiggle(ch::hatching_range rng, const ch::rnd_fn& jig) {
    return rng |
        rv::transform(
            [jig](const auto& poly) {
                return ::jiggle(poly, jig);
            }
    );
}

ch::hatching_range ch::rotate(ch::hatching_range input, double theta) {
    matrix rotation = rotation_matrix(theta);
    return ch::transform(input, rotation);
}

ch::hatching_range ch::disintegrate(ch::hatching_range rng, double amount) {

    if (amount >= 1.0) {
        return rng;
    } else if (amount == 0) {
        return {};
    }

    return rng |
        rv::filter(
            [amount](const auto& p) {
                return ch::uniform_rnd(0.0, 1.0) < amount;
            }
    );
}

cv::Mat ch::paint_cross_hatching(int thickness, ch::hatching_range rng, int swatch_sz) {
    cv::Mat mat(swatch_sz, swatch_sz, CV_8U, 255);
    for (const auto& ls : rng) {
        ch::paint_polyline(mat, ls, thickness, 0, point{ swatch_sz / 2.0,swatch_sz / 2.0 });
    }
    return mat;
}

double ch::gray_level(int thickness, hatching_range rng)
{
    auto mat = paint_cross_hatching(thickness, rng, k_swatch_sz);        
    auto n = k_swatch_sz * k_swatch_sz;
    auto white_pixels = cv::countNonZero(mat);
    return static_cast<double>(n - white_pixels) / static_cast<double>(n);
}

std::string svg_header(int wd, int hgt)
{
    std::stringstream ss;

    ss << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
    ss << "<svg width=\"" + std::to_string(wd) + "px\" height=\"" + std::to_string(hgt) + "px\"  xmlns = \"http://www.w3.org/2000/svg\" version = \"1.1\">\n";
    ss << "<rect width=\"" + std::to_string(wd) + "\" height=\"" + std::to_string(hgt) + "\"  fill=\"white\"/>\n";

    return ss.str();
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

void ch::to_svg(const std::string& filename, int thickness, hatching_range rng, int swatch_sz)
{
    std::ofstream outfile(filename);

    outfile << svg_header(swatch_sz, swatch_sz);

    for (const auto& poly : rng)
        outfile << polyline_to_svg(poly, thickness) << std::endl;

    outfile << "</svg>" << std::endl;
    outfile.close();
}
