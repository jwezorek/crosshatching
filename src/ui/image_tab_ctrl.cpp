#include "image_tab_ctrl.h"
#include <QScrollArea>

namespace {
    QScrollArea* create_img_scroller(const cv::Mat& mat) {
        auto scroller = new QScrollArea();
        auto img_box = new ui::image_box(mat);
        scroller->setWidget(img_box);
        return scroller;
    }
}

image_tab_ctrl::image_tab_ctrl()
{}

void image_tab_ctrl::set_content(const std::vector<std::tuple<std::string, cv::Mat>>& content)
{
    clear_and_delete();
    for (const auto& [name, img] : content) {
        this->addTab(create_img_scroller(img), name.c_str());
    }
}

bool image_tab_ctrl::has_images() const {
    return this->count() > 0;
}

void image_tab_ctrl::clear_and_delete() {
    int n = this->count();
    std::vector<QWidget*> tabs(n);
    for (int i = 0; i < n; ++i) {
        tabs[i] = this->widget(i);
    }
    this->clear();
    for (auto tab : tabs) {
        delete tab;
    }
}