#include "brush.hpp"
#include "geometry.hpp"
#include "crosshatch.hpp"
#include <queue>
#include <iostream>

namespace r = ranges;
namespace rv = ranges::views;

namespace {

    constexpr auto k_num_samples = 4;

    class pipeline_visitor {
    private:
        double t_;
        ch::hatching_range hatching_;
    public:
        pipeline_visitor(double t) : t_(t) {
        }

        void operator()(const ch::param_adapter_fn& param_adapter) {
            t_ = param_adapter(t_);
        }

        void operator()(const ch::brush_adapter_fn& brush_adapter) {
            hatching_ = brush_adapter(hatching_, t_);
        }

        void operator()(const ch::brush_fn& brush) {
            hatching_ = brush(t_);
        }

        ch::hatching_range output() const {
            return hatching_;
        }
    };

}

ch::hatching_range ch::run_brush_pipeline(const brush_pipeline& pipeline, double t)
{
    pipeline_visitor visitor(t);
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

ch::brush_adapter_fn ch::make_one_param_brush_adaptor_fn(std::function< ch::hatching_range(ch::hatching_range, double)> fn, param_adapter_fn param) {
    return [fn, param](hatching_range input, double t)->hatching_range {
        return fn(input, param(t));
    };
}

ch::brush_fn ch::make_linear_hatching_brush_fn(const param_adapter_fn& run_length, const param_adapter_fn& space_length,
    const param_adapter_fn& vert_space, const param_unit_of_hatching_fn& h_fn)
{
    return [=](double t)->hatching_range {
        return linear_crosshatching(
            [run_length, t]()->double {return run_length(t); },
            [space_length, t]()->double {return space_length(t); },
            [vert_space, t]()->double {return vert_space(t); },
            [t, h_fn](double a, double b, double c, double d)->hatching_range { return h_fn(t, a, b, c, d); }
        );
    };
}

ch::brush_fn ch::make_run_pipeline_fn(const ch::brush_pipeline& pipeline) {
    return [pipeline](double t)->hatching_range {
        return run_brush_pipeline(pipeline, t);
    };
}

ch::brush_fn ch::make_merge_fn(const std::vector<ch::brush_fn>& brushes) {
    return [brushes](double t)->hatching_range {
        return rv::join(
            brushes |
            rv::transform([t](auto fn)->hatching_range {return fn(t); })
        );
    };
}

//std::map<double, double> gray_to_param_;
//brush_fn brush_fn_;
//int line_thickness_;
ch::brush::brush(brush_fn fn, int line_thickness) :
    brush_fn_(fn),
    line_thickness_(line_thickness)
{}

struct range_queue_item {
    double left;
    double right;
    double gray_level_left;
    double gray_level_right;
};

double sample(ch::brush_fn fn, double t, int n, int thickness) {
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += ch::gray_level(thickness, fn(t));
    }
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

double ch::brush::get_or_sample(double param) {
    auto iter = param_to_gray_.find(param);
    if (iter != param_to_gray_.end()) {
        return iter->second;
    }
    auto gray = sample(brush_fn_, param, k_num_samples, line_thickness_);
    gray_to_param_[gray] = param;
    param_to_gray_[param] = gray;

    return gray;
}

void ch::brush::build(double epsilon) {
    range_queue queue(compare_items);

    auto make_item = [this](double left, double right)->range_queue_item {
        auto left_gray = get_or_sample(left);
        auto right_gray = get_or_sample(right);
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
        if (item.gray_level_right - item.gray_level_left > epsilon) {
            double mid = (item.left + item.right) / 2.0;
            queue.push(make_item(item.left, mid));
            queue.push(make_item(mid, item.right));
        }
    }
}
    
ch::hatching_range ch::brush::get_hatching(double gray_level, int wd) const {
    auto left = gray_to_param_.lower_bound(gray_level);
    auto right = std::next(left);
    return {};
}