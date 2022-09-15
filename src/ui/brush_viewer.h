#pragma once

#include <QDialog>
#include <QLabel>
#include <string>
#include "float_value_slider.h"
#include "../crosshatching/brush.hpp"

namespace ui {

    class brush_viewer : public QDialog
    {
        Q_OBJECT

    public:
        brush_viewer(const std::string& brush_name, ch::brush_expr_ptr brush,  QWidget* parent = nullptr);

    private:

        void update_swatch();

        QLabel* swatch_box_;
        float_value_slider* slider_;
        ch::brush brush_;
    };
}