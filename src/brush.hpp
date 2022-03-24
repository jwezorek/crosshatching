#pragma once

#include "crosshatch.hpp"
#include <functional>
#include <variant>
#include <vector>
#include <map>

namespace ch {

    constexpr int k_swatch_sz = 512;

    using brush_fn = std::function<ch::crosshatching_swatch(ch::dimensions sz, double t)>;
    using param_adapter_fn = std::function<double (double t)>;
    using brush_adapter_fn = std::function<ch::crosshatching_swatch(ch::crosshatching_swatch, double t)>;
    using brush_pipeline_item = std::variant<brush_fn, param_adapter_fn, brush_adapter_fn>;
    using brush_pipeline = std::vector<brush_pipeline_item>;
    using param_unit_of_hatching_fn = std::function<crosshatching_range(double, double, double, double, double)>;
    using random_brush_adaptor_fn = std::function<crosshatching_swatch(crosshatching_swatch rng, ch::rnd_fn)>;
    ch::crosshatching_swatch run_brush_pipeline(const brush_pipeline& pipeline, dimensions sz, double t);

    param_adapter_fn make_constant_fn(double k);
    param_adapter_fn make_lerp_fn(double v1, double v2);
    param_adapter_fn make_lerped_normal_dist_fn(double mu1, double sigma1, double mu2, double sigma2);

    param_unit_of_hatching_fn make_default_hatching_unit();
    brush_adapter_fn make_one_param_brush_adaptor(brush_adapter_fn fn, param_adapter_fn param);
    brush_adapter_fn make_random_brush_adaptor(random_brush_adaptor_fn fn, rnd_fn rnd);

    brush_fn make_linear_hatching_brush_fn(const param_adapter_fn& run_length, const param_adapter_fn& space_length, const param_adapter_fn& vert_space,
        const param_unit_of_hatching_fn& h_fn);

    brush_fn make_run_pipeline_fn(const brush_pipeline& pipeline);
    brush_fn make_merge_fn(const std::vector<brush_fn>& brushes);

    param_adapter_fn make_ramp_fn(double k, bool right, bool up);
    //param_adapter_fn make_one_parameter_param_adapter(param_adapter_fn f, param_adapter_fn arg);

    class brush {
    private:
        using gray_map_iter = std::map<double, double>::iterator;
        std::map<double, double> gray_to_param_;
        std::unordered_map<double, double> param_to_gray_;
        brush_fn brush_fn_;
        int line_thickness_;
        dimensions swatch_sz_;

        bool is_uinitiailized() const;
        double get_or_sample_param(double param);
        double build_between(double v, gray_map_iter left, gray_map_iter right, double epsilon);
        double build_to_gray_level(double gray_level, double epsilon);

    public:
        brush(brush_fn fn, int line_thickness = 1, dimensions swatch_sz = { k_swatch_sz });
        void build(double epsilon);
        void build_n(int n);
        crosshatching_swatch get_hatching(double gray_level, dimensions sz, double epsilon);
        double min_gray_level() const;
        double max_gray_level() const;
    };

}