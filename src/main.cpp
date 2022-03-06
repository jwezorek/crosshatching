#include "geometry.hpp"
#include "crosshatch.hpp"
#include "brush.hpp"
#include <opencv2/highgui.hpp>
#include <iostream>
#include <numbers>

ch::brush_fn make_pipeline_fn(double theta) {
    return ch::make_run_pipeline_fn(
        ch::brush_pipeline{
            ch::make_linear_hatching_brush_fn(
                ch::make_lerped_normal_dist_fn(10, 25, 100, 20),
                ch::make_lerped_normal_dist_fn(100, 20, 0, 0.05),
                ch::make_lerped_normal_dist_fn(7, 0.5, 0.75, 0.05),
                ch::make_default_hatching_unit()
            ),
            ch::make_one_param_brush_adaptor_fn(ch::rotate, ch::make_constant_fn(theta))
        }
    );
}

int main()
{
    auto brush = ch::make_merge_fn({ make_pipeline_fn(std::numbers::pi / 4.0), make_pipeline_fn(-std::numbers::pi / 4.0) });
    int n = 100;
    
    for (int i = 0; i <= n; i++) {
        auto crosshatching = brush(i * (1.0 / n));
        auto mat = ch::paint_cross_hatching(2, crosshatching);
        cv::imshow("crosshatch", mat);
        int k = cv::waitKey(50);
    }

    return 0;
}
