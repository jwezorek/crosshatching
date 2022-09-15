#include "brush_viewer.h"
#include <QtWidgets>

ui::brush_viewer::brush_viewer(const std::string& brush_name, ch::brush_expr_ptr b_expr, QWidget *parent)
    : QDialog(parent)
{
    setWindowTitle(brush_name.c_str());
    auto* main_layout = new QVBoxLayout();
    auto* row = new QHBoxLayout();
    this->setLayout(main_layout);
    row->addWidget(swatch_box_ = new QLabel());
    swatch_box_->setFixedSize(QSize(ch::k_swatch_sz, ch::k_swatch_sz));
    row->addWidget(slider_ = new ui::float_value_slider("val", 0, 1.0, 0, 100));
    connect(slider_, &float_value_slider::slider_released,
        [this]() {
            this->update_swatch();
        }
    );
    main_layout->addLayout(row);

    brush_ = ch::brush(b_expr);
    brush_.build_n(10);
    update_swatch();
}

void ui::brush_viewer::update_swatch() {
    auto value = slider_->value();
    auto swatch = brush_.swatch(value);

    swatch = ch::convert_to_3channel_grayscale(swatch);
    int wd = swatch.cols;
    int hgt = swatch.rows;
    int stride = swatch.step;

    swatch_box_->setPixmap(QPixmap::fromImage(QImage(swatch.data, wd, hgt, stride, QImage::Format_BGR888)));
}
