#include "geometry.hpp"
#include "crosshatch.hpp"
#include "brush.hpp"
#include <opencv2/highgui.hpp>
#include <iostream>
#include <numbers>
#include <chrono>

int main() {
    auto test = ch::jiggle(
        ch::linear_crosshatching(ch::normal_rnd_fn(100, 30), ch::normal_rnd_fn(40, 20), ch::const_rnd_fn(6), ch::one_horz_stroke, 1024), 
        ch::normal_rnd_fn(0, 0.005)
    );
    auto mat = ch::paint_cross_hatching(1, test, 1024);
    cv::imshow("crosshatch", mat);
    int k = cv::waitKey(0);
}

/*
ch::brush_fn make_pipeline_fn(double theta) {
    return ch::make_run_pipeline_fn(
        ch::brush_pipeline{
            ch::make_linear_hatching_brush_fn(
                ch::make_lerped_normal_dist_fn(10, 25, 100, 20),
                ch::make_lerped_normal_dist_fn(100, 20, 0, 0.05),
                ch::make_lerped_normal_dist_fn(7, 0.5, 0.75, 0.05),
                ch::make_default_hatching_unit()
            ),
            [](ch::hatching_range rng, double t)->ch::hatching_range {
                return ch::jitter(rng, ch::normal_rnd_fn(3.0*t+1.0, 0.3+1.5*t), ch::normal_rnd_fn(0.0, 0.0001+3*t));
            } ,
            ch::make_one_param_brush_adaptor_fn(ch::rotate, ch::make_constant_fn(theta))
        }
    );
}

int main()
{
    auto brush_fn = ch::make_merge_fn({ make_pipeline_fn(std::numbers::pi / 4.0), make_pipeline_fn(-std::numbers::pi / 4.0) });

    ch::brush br(brush_fn);
    
    std::cout << "building...\n";
    auto start = std::chrono::high_resolution_clock::now();
    br.build(0.1);
    auto finish = std::chrono::high_resolution_clock::now();;
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count() << "\n";

    auto hatching = br.get_hatching(0.5, 512);

    return 0;
}
*/