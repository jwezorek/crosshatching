#pragma once

#include "image_box.h"
#include "settingctrls.hpp"
#include <QWidget>

namespace ui {

    class flow_direction_panel : public image_box {

        Q_OBJECT

    public:
        flow_direction_panel();
    };

    class rgn_properties_panel : public QWidget {

        Q_OBJECT

    public:
        rgn_properties_panel();
    };


    class rgn_map_ctrl : public image_box {

        Q_OBJECT

    private:
        ch::dimensions<int> base_sz_;
        double scale_;
        std::vector<ch::colored_polygon> scaled_regions_;

        void display();

    public:
        rgn_map_ctrl();
        void set_regions(vector_graphics_ptr gfx);
        void set_scale(double sc);
    };

}