#include "brush.hpp"
#include "geometry.hpp"
#include "util.hpp"
#include <opencv2/highgui.hpp>
#include <iostream>
#include <future>
#include <numeric>
#include <cmath>
#include <numbers>
#include <sstream>
#include <fstream>
#include <queue>

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
                [j, jit](const ch::point& pt)->ch::point {
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
                rv::single(poly.back())
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

    constexpr auto k_num_samples = 4;

    class pipeline_visitor {
    private:
        double stroke_wd_;
        ch::dimensions sz_;
        double t_;
        ch::crosshatching_swatch hatching_;
    public:
        pipeline_visitor(double thickness, ch::dimensions sz, double t) : stroke_wd_(thickness), sz_(sz), t_(t) {
        }

        void operator()(const ch::param_adapter_fn& param_adapter) {
            t_ = param_adapter(t_);
        }

        void operator()(const ch::brush_adapter_fn& brush_adapter) {
            hatching_ = brush_adapter(hatching_, t_);
        }

        void operator()(const ch::brush_fn& brush) {
            hatching_ = brush(stroke_wd_,sz_,t_);
        }

        ch::crosshatching_swatch output() const {
            return hatching_;
        }
    };

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
        auto n = static_cast<int>(hgt) + 1;
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
                    auto yy = (n - i) + y;

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

ch::crosshatching_swatch ch::linear_crosshatching(rnd_fn run_length, rnd_fn space_length, rnd_fn vert_space,
    unit_of_hatching_fn h_fn, ch::dimensions viewable_swatch_sz, double stroke_wd) {
    auto dim = (std::sqrt(2.0) * (viewable_swatch_sz.wd + viewable_swatch_sz.hgt)) / 2.0;
    auto swatch_sz = dimensions(dim);
    return {
        rv::join(
            running_sum(vert_space, -swatch_sz.hgt / 2.0) |
            rv::take_while([swatch_sz](const running_sum_item& rsi) { return rsi.sum < swatch_sz.hgt / 2.0; }) |
            rv::transform([=](const running_sum_item& rsi) {return row_of_crosshatching(swatch_sz.wd, rsi.next_item, rsi.sum, run_length, space_length, h_fn); })
        ),
        viewable_swatch_sz,
        stroke_wd
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
        swatch.sz,
        swatch.stroke_wd
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
        swatch.sz,
        swatch.stroke_wd
    };
}

ch::crosshatching_swatch ch::jiggle(ch::crosshatching_swatch swatch, ch::rnd_fn jig) {
    return {
        swatch.content |
            rv::transform(
                [jig](const auto& poly) {
                    return ::jiggle(poly, jig);
                }
            ),
        swatch.sz,
        swatch.stroke_wd
    };
}

ch::crosshatching_swatch ch::rotate(ch::crosshatching_swatch swatch, double theta) {
    matrix rotation = rotation_matrix(theta);
    return {
        ch::transform(swatch.content, rotation),
        swatch.sz,
        swatch.stroke_wd
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
        swatch.sz ,
        swatch.stroke_wd
    };
}

cv::Mat ch::paint_cross_hatching(ch::crosshatching_swatch swatch) {
    cv::Mat mat(static_cast<int>(swatch.sz.hgt), static_cast<int>(swatch.sz.wd), CV_8U, 255);
    for (const auto& ls : swatch.content) {
        ch::paint_polyline(mat, ls, swatch.stroke_wd, 0, point{ swatch.sz.wd / 2.0, swatch.sz.hgt / 2.0 });
    }
    return mat;
}

double ch::gray_level(crosshatching_swatch swatch)
{
    auto mat = paint_cross_hatching(swatch);
    auto n = swatch.sz.wd * swatch.sz.hgt;
    auto white_pixels = cv::countNonZero(mat);
    return static_cast<double>(n - white_pixels) / static_cast<double>(n);
}

void ch::to_svg(const std::string& filename, crosshatching_swatch swatch)
{
    std::ofstream outfile(filename);

    outfile << svg_header(static_cast<int>(swatch.sz.wd), static_cast<int>(swatch.sz.hgt));

    for (const auto& poly : swatch.content)
        outfile << polyline_to_svg(poly, swatch.stroke_wd) << std::endl;

    outfile << "</svg>" << std::endl;
    outfile.close();
}

ch::dimensions::dimensions(double d) : wd(d), hgt(d)
{
}

ch::dimensions::dimensions(double w, double h) : wd(w), hgt(h)
{
}


ch::crosshatching_swatch ch::run_brush_pipeline(const brush_pipeline& pipeline, double thickness, dimensions sz, double t)
{
    pipeline_visitor visitor(thickness,sz,t);
    for (const auto& item : pipeline) {
        std::visit(visitor, item);
    }
    return visitor.output();
}

ch::param_adapter_fn ch::make_constant_fn(double k) {
    return [k](double t) {return k; };
}

ch::param_adapter_fn ch::make_lerp_fn(double v1, double v2)
{
    return [v1, v2](double t) { return std::lerp(v1, v2, t);  };
}

ch::param_adapter_fn ch::make_lerped_normal_dist_fn(double mu1, double sigma1, double mu2, double sigma2)
{
    return [mu1, sigma1, mu2, sigma2](double t) {
        auto mu = std::lerp(mu1, mu2, t);
        auto sigma = std::lerp(sigma1, sigma2, t);
        return ch::normal_rnd(mu, sigma);
    };
}

ch::param_unit_of_hatching_fn ch::make_default_hatching_unit() {
    return [](double t, double x1, double x2, double y, double hgt) {
        return ch::one_horz_stroke(x1, x2, y, hgt);
    };
}

ch::brush_adapter_fn ch::make_one_param_brush_adaptor(brush_adapter_fn fn, param_adapter_fn param) {
    return [fn, param](crosshatching_swatch input, double t)->crosshatching_swatch {
        return fn(input, param(t));
    };
}

ch::brush_adapter_fn ch::make_random_brush_adaptor(random_brush_adaptor_fn fn, rnd_fn rnd) {
    return [fn, rnd](crosshatching_swatch input, double t)->crosshatching_swatch {
        return fn(input, rnd);
    };
}

/*
ch::param_adapter_fn ch::make_one_parameter_param_adapter(ch::param_adapter_fn fn, ch::param_adapter_fn arg) {
    return [fn, arg]( double t)->double {
        return fn(arg(t));
    };
}
*/

ch::brush_fn ch::make_linear_hatching_brush_fn(const param_adapter_fn& run_length, const param_adapter_fn& space_length,
    const param_adapter_fn& vert_space, const param_unit_of_hatching_fn& h_fn)
{
    return [=](double stroke, ch::dimensions sz, double t)->crosshatching_swatch {
        return linear_crosshatching(
            [run_length, t]()->double {return run_length(t); },
            [space_length, t]()->double {return space_length(t); },
            [vert_space, t]()->double {return vert_space(t); },
            [t, h_fn](double a, double b, double c, double d)->crosshatching_range { return h_fn(t, a, b, c, d); },
            sz,
            stroke
        );
    };
}

ch::brush_fn ch::make_run_pipeline_fn(const ch::brush_pipeline& pipeline) {
    return [pipeline](double stroke_wd, ch::dimensions sz, double t)->crosshatching_swatch {
        return run_brush_pipeline(pipeline, stroke_wd, sz, t);
    };
}

ch::brush_fn ch::make_merge_fn(const std::vector<ch::brush_fn>& brushes) {
    return [brushes](double stroke_wd, dimensions sz, double t)->crosshatching_swatch {
        return {
            rv::join(
                brushes |
                rv::transform([stroke_wd,sz,t](auto fn)->crosshatching_range { return fn(stroke_wd, sz,t).content; })
            ),
            sz,
            stroke_wd
        };
    };
}

ch::param_adapter_fn ch::make_ramp_fn(double k, bool right, bool up)
{
    return [k, up, right](double t) {return ch::ramp(t, k, right, up); };
}

//std::map<double, double> gray_to_param_;
//brush_fn brush_fn_;
//int line_thickness_;
ch::brush::brush(brush_fn fn, int line_thickness, double epsilon, dimensions sz) :
    brush_fn_(fn),
    line_thickness_(line_thickness),
    swatch_sz_{ sz },
    eps_(epsilon)
{}

struct range_queue_item {
    double left;
    double right;
    double gray_level_left;
    double gray_level_right;
};

double sample(ch::dimensions sz, ch::brush_fn fn, double t, int n, int thickness) {
    std::vector<std::future<double>> samples(n);

    auto compute_gray_level = [sz](int th, ch::brush_fn f, double t)->double {
        return ch::gray_level(f(th, sz, t));
    };
    std::generate(samples.begin(), samples.end(),
        [=]() {return std::async(std::launch::async, compute_gray_level, thickness, fn, t); }
    );
    double sum = std::accumulate(samples.begin(), samples.end(), 0.0,
        [](double s, std::future<double>& fut) {
            return s + fut.get();
        }
    );
    return sum / n;
}

bool compare_items(const range_queue_item& v1, const range_queue_item& v2) {
    return v1.gray_level_right - v1.gray_level_left < v2.gray_level_right - v2.gray_level_left;
}

using range_queue = std::priority_queue<
    range_queue_item, 
    std::vector<range_queue_item>,
    std::function<bool(const range_queue_item&, const range_queue_item&)>
>;

void debug(const range_queue& input) {
    range_queue q = input;
    while (!q.empty()) {
        auto item = q.top();
        q.pop();
        std::cout << item.gray_level_right - item.gray_level_left << " ";
    }
    std::cout << "\n";
}

bool ch::brush::is_uinitiailized() const
{
    return gray_to_param_.empty();
}

double ch::brush::get_or_sample_param(double param) {
    auto iter = param_to_gray_.find(param);
    if (iter != param_to_gray_.end()) {
        return iter->second;
    }
    auto gray = sample(swatch_sz_, brush_fn_, param, k_num_samples, line_thickness_);
    gray_to_param_[gray] = param;
    param_to_gray_[param] = gray;

    return gray;
}

double ch::brush::build_between(double gray, gray_map_iter left, gray_map_iter right)
{
    double mid_param = (left->second + right->second) / 2.0;
    auto mid_gray_value = get_or_sample_param(mid_param);
    if (std::abs(mid_gray_value - gray) < eps_) {
        return mid_param;
    }
    auto mid_iter = gray_to_param_.find(mid_gray_value);
    if (gray < mid_gray_value) {
        return build_between(gray, left, mid_iter);
    } else {
        return build_between(gray, mid_iter, right);
    }
}

void ch::brush::build() {
    if (!is_uinitiailized()) {
        throw std::runtime_error("brush is already built");
    }
    range_queue queue(compare_items);

    auto make_item = [this](double left, double right)->range_queue_item {
        auto left_gray = get_or_sample_param(left);
        auto right_gray = get_or_sample_param(right);
        return {
            left,
            right,
            left_gray,
            right_gray
        };
    };
    
    queue.push(make_item(0, 1));
    while (!queue.empty()) {
        //debug(queue);
        auto item = queue.top();
        queue.pop();
        if (item.gray_level_right - item.gray_level_left > eps_) {
            double mid = (item.left + item.right) / 2.0;
            queue.push(make_item(item.left, mid));
            queue.push(make_item(mid, item.right));
        }
    }
}

void ch::brush::build_n(int n)
{
    if (!is_uinitiailized()) {
        throw std::runtime_error("brush is already built");
    }
    double delta = 1.0 / static_cast<double>(n);
    for (int i = 0; i <= n; ++i) {
        get_or_sample_param(delta * i);
    }
}

double ch::brush::stroke_width() const
{
    return line_thickness_;
}

double ch::brush::build_to_gray_level(double gray_level) {
    auto right = gray_to_param_.lower_bound(gray_level);
    auto left = (right != gray_to_param_.begin()) ? std::prev(right) : right;
    if (std::abs(gray_level - left->first) < eps_) {
        return left->second;
    }
    if (std::abs(gray_level - right->first) < eps_) {
        return right->second;
    }
    return build_between(gray_level, left, right);
}
    
ch::crosshatching_swatch ch::brush::get_hatching(double gray_level, dimensions sz) {
    if (is_uinitiailized()) {
        throw std::runtime_error("brush is not built");
    }
    if (gray_level < min_gray_level()) {
        return get_hatching(ch::brush::min_gray_level(), sz);
    } else if (gray_level > max_gray_level()) {
        return get_hatching(ch::brush::max_gray_level(), sz);
    }
    auto param = build_to_gray_level(gray_level);
    return brush_fn_(line_thickness_, sz, param);
}

double ch::brush::min_gray_level() const {
    if (is_uinitiailized()) {
        throw std::runtime_error("brush is not built");
    }
    return gray_to_param_.begin()->first;
}

double ch::brush::max_gray_level() const {
    if (is_uinitiailized()) {
        throw std::runtime_error("brush is not built");
    }
    return std::prev(gray_to_param_.end())->first;
}