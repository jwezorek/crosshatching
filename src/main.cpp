#include "geometry.hpp"
#include "crosshatch.hpp"
#include "brush.hpp"
#include "drawing.hpp"
#include "imgproc.hpp"
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

    auto drawing = ch::generate_crosshatched_drawing("C:\\test\\grace2.png", { 8, 3.0f, 12 }, 4.0, br);
    ch::to_svg("C:\\test\\grace2-test.svg", drawing);
    /*
    ch::polygon_with_holes poly = { {
            {10,10},
            {10,60},
            {20,60},
            {20,40},
            {30,40},
            {30,60},
            {40,60},
            {40,10},
            {30,10},
            {30,30},
            {20,30},
            {20,10}
        },
        {}
    };

    ch::drawing drawing = {
        ch::crosshatched_poly_with_holes(poly, 0.75, br),
        {70,70},
        1
    };
    */

    

    return 0;
}

/*
int main() {
    auto test = ch::jitter(
        ch::fragment(
            ch::linear_crosshatching(ch::normal_rnd_fn(100, 30), ch::const_rnd_fn(1), ch::const_rnd_fn(6), ch::one_horz_stroke, 512),
            ch::normal_rnd_fn(4, 1)
        ),
        ch::normal_rnd_fn(0, 0.5)
    );
    auto mat = ch::paint_cross_hatching(1, test, 512);
    cv::imshow("crosshatch", mat);
    int k = cv::waitKey(0);
}
*/

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