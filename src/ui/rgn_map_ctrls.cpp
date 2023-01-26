#include "rgn_map_ctrls.hpp"
#include <QtWidgets>
#include "../crosshatching/util.hpp"

ui::flow_direction_panel::flow_direction_panel() {
    set_image(ch::blank_monochrome_bitmap(100));
}

ui::rgn_properties_panel::rgn_properties_panel() {
    QVBoxLayout* layout = new QVBoxLayout(this);
    layout->addWidget(new QLabel("brush"));
    layout->addWidget(new QComboBox());
    layout->addWidget(new QLabel("flow"));
    layout->addWidget(new flow_direction_panel());
    layout->addStretch();
}

ui::rgn_map_ctrl::rgn_map_ctrl() : 
        base_sz_(-1), scale_(1.0) {
}

void ui::rgn_map_ctrl::set_regions(vector_graphics_ptr gfx) {
    scaled_regions_ = (scale_ != 1.0) ?
        ch::scale(gfx->polygons, scale_) :
        gfx->polygons;
    base_sz_ = gfx->sz;
    auto sz = scale_ * base_sz_;
    setFixedSize({ sz.wd, sz.hgt });
    display();
}

void ui::rgn_map_ctrl::display() {
    auto mat = ch::paint_polygons(scaled_regions_, scale_ * base_sz_);
    set_image(mat);
}

void ui::rgn_map_ctrl::set_scale(double sc) {
    scale_ = sc;
    display();
}