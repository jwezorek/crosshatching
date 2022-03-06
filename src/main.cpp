#include "geometry.hpp"
#include "crosshatch.hpp"
#include "brush.hpp"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <range/v3/all.hpp>
#include <iostream>
#include <functional>
#include <tuple>
#include <vector>
#include <random>
#include <numbers>

int main()
{
    auto pipeline = ch::brush_pipeline{
        ch::make_linear_hatching_brush_fn(
            ch::make_lerped_normal_dist_fn(10, 25, 100, 20),
            ch::make_lerped_normal_dist_fn(100, 20, 0, 0.05),
            ch::make_lerped_normal_dist_fn(7, 0.5, 0.75, 0.05),
            ch::make_default_hatching_unit()
        ),
        ch::make_one_param_brush_adaptor_fn( ch::rotate, ch::make_constant_fn(std::numbers::pi/4.0))
    };

    int n = 100;
    for (int i = 0; i <= n; i++) {
        auto crosshatching = ch::run_brush_pipeline(pipeline, i * (1.0 / n));
        auto mat = ch::paint_cross_hatching(2, crosshatching);
        cv::imshow("crosshatch", mat);
        int k = cv::waitKey(0);
    }

    return 0;
}
