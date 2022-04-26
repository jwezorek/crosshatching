#pragma once
#include <map>
#include "geometry.hpp"
#include "util.hpp"
#include <opencv2/core.hpp>
#include <range/v3/all.hpp>
#include <tuple>
#include <functional>
#include <variant>
#include <vector>

namespace ch {
    struct dimensions {
        double wd;
        double hgt;

        dimensions(double d = 0.0);
        dimensions(double w, double h);
    };

    using crosshatching_range = ranges::any_view<polyline>;
    struct crosshatching_swatch {
        crosshatching_range content;
        dimensions sz;
        double stroke_wd;
    };

    using unit_of_hatching_fn = std::function<crosshatching_range(double, double, double, double)>;

    unit_of_hatching_fn make_shading_stroke(double sz_pcnt, double variance, bool centered);
    unit_of_hatching_fn make_brick_stroke();
    crosshatching_range one_horz_stroke(double x1, double x2, double y, double hgt);

    crosshatching_swatch linear_crosshatching(rnd_fn run_length, rnd_fn space_length, rnd_fn vert_space,
        unit_of_hatching_fn h_fn = one_horz_stroke, dimensions sz = { 512,512 }, double stroke_wd = 1);
    crosshatching_swatch fragment(crosshatching_swatch rng, ch::rnd_fn frag);
    crosshatching_swatch jitter(crosshatching_swatch rng, ch::rnd_fn jitter);
    crosshatching_swatch jiggle(crosshatching_swatch rng, ch::rnd_fn jiggle);
    crosshatching_swatch rotate(crosshatching_swatch rng, double theta);
    crosshatching_swatch disintegrate(crosshatching_swatch rng, double amount);

    cv::Mat paint_cross_hatching(crosshatching_swatch swatch, cv::Mat bkgd);
    double measure_gray_level(crosshatching_swatch rng, cv::Mat bkgd);
    void to_svg(const std::string& filename, crosshatching_swatch swatch);

    constexpr int k_swatch_sz = 512;
    constexpr double k_epsilon = 0.001;

    using brush_fn = std::function<ch::crosshatching_swatch(double, ch::dimensions, double)>;
    using param_adapter_fn = std::function<double (double t)>;
    using brush_adapter_fn = std::function<ch::crosshatching_swatch(ch::crosshatching_swatch, double t)>;
    using brush_pipeline_item = std::variant<brush_fn, param_adapter_fn, brush_adapter_fn>;
    using brush_pipeline = std::vector<brush_pipeline_item>;
    using param_unit_of_hatching_fn = std::function<crosshatching_range(double, double, double, double, double)>;
    using random_brush_adaptor_fn = std::function<crosshatching_swatch(crosshatching_swatch rng, ch::rnd_fn)>;
    ch::crosshatching_swatch run_brush_pipeline(const brush_pipeline& pipeline, double thickness, dimensions sz, double t);

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
        double eps_;
        std::vector<cv::Mat> bkgds_;

        bool is_uinitiailized() const;
        double get_or_sample_param(double param);
        double build_between(double v, gray_map_iter left, gray_map_iter right);
        double build_to_gray_level(double gray_level);

    public:
        brush(brush_fn fn, int line_thickness = 1, double epsilon = k_epsilon, 
            dimensions swatch_sz = { k_swatch_sz }, const std::vector<cv::Mat>& bkgds = {});

        void build();
        void build_n(int n);
        double stroke_width() const;
        double gray_value_to_param(double gray_val);
        crosshatching_swatch get_hatching(double gray_val, dimensions sz);
        std::vector<cv::Mat> render_swatches(double gray_level);
        double min_gray_level() const;
        double max_gray_level() const;
        static int num_samples();
    };

    using hierarchical_brush_component = std::function< brush_fn(double) >;

    class hierarchical_brush {
    private:
        std::map<double, brush_fn> gray_val_to_brush_; 
        std::map<double, hierarchical_brush_component> gray_val_to_brush_gen_;
        double epsilon_;
        double line_thickness_;
        dimensions swatch_sz_;
        std::map<double, std::function<ch::crosshatching_swatch(dimensions)>> brushes_;

        brush_fn gray_val_to_brush(double gray_val);

    public:
        hierarchical_brush(const std::vector<brush_fn>& brush_fns, const std::vector<double>& gray_intervals, 
            int line_thickness = 1, double epsilon = k_epsilon, dimensions swatch_sz = { k_swatch_sz });
        hierarchical_brush(const std::vector<hierarchical_brush_component>& brush_fns, const std::vector<double>& gray_intervals,
            int line_thickness = 1, double epsilon = k_epsilon, dimensions swatch_sz = { k_swatch_sz });
        void build(const std::vector<double>& gray_values);
        crosshatching_swatch get_hatching(double gray_level, dimensions sz);
        double stroke_width() const;
    };

}