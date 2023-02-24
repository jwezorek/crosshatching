#pragma once

#include "treepanel.h"
#include "../crosshatching/util.hpp"

namespace ui {

    class layer_panel : public ui::list_panel {

        Q_OBJECT

    public:
        layer_panel();
        void set_brush_names(const std::vector<std::string>& brush_names);
        std::vector<std::tuple<std::string, double>> layers() const;
        bool has_content() const;
        ch::json to_json() const;
        void from_json(const ch::json& json);

    signals:
        void layers_changed();

    private:

        std::vector<std::string> brush_names_;
        std::map<double, std::string> layers_;

        std::tuple<std::string, double> row(int n) const;
        void add_layer();
        void cell_double_clicked(int r, int col);
        void delete_layer();
        void insert_layer(const std::string& brush, double end_of_range);
        void setRowText(int row, const std::string& brush, const std::string& from, const std::string& to);
        void sync_layers_to_ui();
    };

}