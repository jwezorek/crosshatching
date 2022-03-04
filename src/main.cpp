#include "geometry.hpp"
#include "crosshatch.hpp"
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
    
    
    //auto crosshatching = ch::apply_jitter(
    //    ch::horizontal_crosshatching(256, ch::normal_rnd_fn(100, 30.0), ch::const_rnd_fn(0), ch::const_rnd_fn(5), ch::make_shading_stroke(1, 0.02, true)),
    //    ch::normal_rnd_fn(2,0.1), ch::normal_rnd_fn(0, 0.05)
    //);
    auto crosshatching = ch::linear_crosshatching( 512, ch::normal_rnd_fn(80,40.0), ch::normal_rnd_fn(15,5.0), ch::normal_rnd_fn(2, 0.2));

    auto vec =  ranges::to_vector(crosshatching);
    auto mat = ch::paint_cross_hatching(1, 512, vec);
    cv::imshow("crosshatch", mat);
    std::cout << "crosshatch: " << ch::gray_level(1, crosshatching) << "\n";

    auto vec2 = ranges::to_vector(ch::disintegrate(vec, 0.5 /*, ch::fragment_param{4,ch::normal_rnd_fn(2, 0.2)}*/));
    cv::imshow("disintegrate", ch::paint_cross_hatching(1, 512, vec2));
    std::cout << "disintegrate: " << ch::gray_level(1, vec2) << "\n";

    ch::to_svg("C:\\test\\ch.svg", 512, 1, vec);

    int k = cv::waitKey(0); // Wait for a keystroke in the window
    

    return 0;
}

/*
{
    {hz run:[norm(30,15)...norm(100,5)] 
        space:[norm(15,5)...norm(0,1)] 
        vert:[norm(7,0.2)...norm(0.5,0.1)]
    },
    floor_left(
        0.85,
        {hz 
            rot: 90
            run:[norm(30,15)...norm(100,5)]
            space:[norm(15,5)...norm(0,1)]
            vert:[norm(7,0.2)...norm(0.5,0.1)]
        }
    )
}
*/

/*
int main()
{
    auto crosshatching = ch::apply_jitter(
        ch::horz_crosshatching(1000, ch::normal_rnd_fn(100, 30.0), ch::normal_rnd_fn(15, 5.0), ch::normal_rnd_fn(6, 0.2)),
        ch::normal_rnd_fn(7.0, 1.0), ch::normal_rnd_fn(0.0, 0.6)
    );
    auto mat = ch::paint_cross_hatching(1, 1000, crosshatching);
    //auto crosshatching = ch::horz_crosshatching_mat(2, 1000, ch::normal_rnd_fn(40,10.0), ch::normal_rnd_fn(15,5.0), ch::normal_rnd_fn(7, 0.2));
    cv::imshow("crosshatch", mat);
    int k = cv::waitKey(0); // Wait for a keystroke in the window

    return 0;
}
*/

/*
int main()
{
    auto pi_over_4 = std::numbers::pi / 3.0;
    auto crosshatching = ranges::views::concat(
        ch::rotated_crosshatching( pi_over_4, 1000, ch::normal_rnd_fn(100, 30.0), ch::normal_rnd_fn(15, 5.0), ch::normal_rnd_fn(6, 0.2)),
        ch::rotated_crosshatching( -pi_over_4, 1000, ch::normal_rnd_fn(100, 30.0), ch::normal_rnd_fn(15, 5.0), ch::normal_rnd_fn(6, 0.2)),
        ch::horz_crosshatching(1000, ch::normal_rnd_fn(100, 30.0), ch::normal_rnd_fn(15, 5.0), ch::normal_rnd_fn(6, 0.2))
    );
    ch::test();
    auto mat = ch::paint_cross_hatching(1, 1000, crosshatching);
    //auto crosshatching = ch::horz_crosshatching_mat(2, 1000, ch::normal_rnd_fn(40,10.0), ch::normal_rnd_fn(15,5.0), ch::normal_rnd_fn(7, 0.2));
    cv::imshow("crosshatch", mat);
    int k = cv::waitKey(0); // Wait for a keystroke in the window

    return 0;
}
*/