#include "geometry.hpp"
#include "brush.hpp"
#include "drawing.hpp"
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <numbers>
#include <chrono>
#include <iostream>


ch::brush_fn make_pipeline_fn(double theta, double start) {
    return ch::make_run_pipeline_fn(
        ch::brush_pipeline{    
            ch::make_ramp_fn(start, true,true),
            ch::make_linear_hatching_brush_fn(
                ch::make_lerped_normal_dist_fn(0, 50, 800, 100),
                ch::make_lerped_normal_dist_fn(200, 20, 0, 0.05),
                ch::make_lerped_normal_dist_fn(7, 0.5, 0.5, 0.05),
                ch::make_default_hatching_unit()
            ),
            ch::make_one_param_brush_adaptor(ch::rotate, ch::make_constant_fn(theta)),
            ch::make_ramp_fn(0.20, false, true),
            ch::disintegrate
        }
    );
}

int main()
{
    ch::brush br(
        ch::make_run_pipeline_fn(
            ch::brush_pipeline{
                ch::make_merge_fn({
                    make_pipeline_fn(std::numbers::pi / 4.0, 0),
                    make_pipeline_fn(-std::numbers::pi / 4.0, 0.65)
                }),
                ch::make_random_brush_adaptor(ch::jiggle, ch::normal_rnd_fn(0.0, 0.02))
            }
        ),
        1
    );
    br.build_n(10);

    auto drawing = ch::generate_crosshatched_drawing("C:\\test\\paris.png", { 8, 3.0f, 12 }, 4.0, br);
    ch::to_svg("C:\\test\\paris.svg", drawing);

    return 0;
}
