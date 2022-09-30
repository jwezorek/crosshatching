#pragma once

#include <QTabWidget>
#include <opencv2/core.hpp>
#include <vector>
#include <tuple>
#include <string>
#include "image_box.h"

class image_tab_ctrl  : public QTabWidget
{
    Q_OBJECT

public:
    image_tab_ctrl();
    void set_content(const std::vector<std::tuple<std::string, cv::Mat>>& content);
    bool has_images() const;
private:
    void clear_and_delete();
};
