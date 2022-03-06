#include "brush.hpp"
#include "geometry.hpp"
#include "crosshatch.hpp"

namespace r = ranges;
namespace rv = ranges::views;

namespace {

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

