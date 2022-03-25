#include "brush.hpp"
#include "geometry.hpp"
#include "crosshatch.hpp"
#include <queue>
#include <iostream>
#include <future>
#include <numeric>

namespace r = ranges;
namespace rv = ranges::views;

namespace {

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