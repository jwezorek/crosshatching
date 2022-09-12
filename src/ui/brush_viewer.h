#pragma once

#include <QDialog>
#include <QLabel>
#include <string>
#include "float_value_slider.h"
#include "../crosshatching/old_brush.hpp"

namespace ui {

    class brush_viewer : public QDialog
    {
        Q_OBJECT

    public:
        brush_viewer(const std::string& brush_name, ch::brush_fn brush,  QWidget* parent = nullptr);

    private:

        void update_swatch();

        QLabel* swatch_box_;
        float_value_slider* slider_;
        ch::old_brush brush_;
    };
}